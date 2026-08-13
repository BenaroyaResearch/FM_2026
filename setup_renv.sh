#!/usr/bin/env bash
#
# One-time setup of the R reference library that munge_gwas needs.
#
# munge_sumstats.R loads four SNPlocs/BSgenome reference data packages. They come to
# roughly 13 GB, too large to ship inside the container image, so they live in this
# project's renv library and get mounted into the job pods instead. The image carries
# MungeSumstats; this library carries the reference data. Both are required.
#
#   ./setup_renv.sh          restore the library, then verify it
#   ./setup_renv.sh check    verify only, change nothing
#
# Environment overrides:
#   RENV_PATHS_CACHE   shared renv package cache (default: /renv)
#
# Safe to re-run: renv::restore() is idempotent, and 'check' is read-only.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

SNAKEFILE="pipeline_munge/Snakefile"
SETTINGS="renv/settings.json"
LOG_DIR="$PROJECT_ROOT/logs"

# Keep in sync with the pre-flight check in munge_gwas.
REF_PACKAGES=(
    SNPlocs.Hsapiens.dbSNP155.GRCh38
    SNPlocs.Hsapiens.dbSNP155.GRCh37
    BSgenome.Hsapiens.NCBI.GRCh38
    BSgenome.Hsapiens.1000genomes.hs37d5
)

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
info() { printf '%s\n' "$*"; }
ok() { printf '  ok    %s\n' "$*"; }
bad() { printf '  FAIL  %s\n' "$*"; }

# --- version expectations --------------------------------------------------

# Read the required R version out of the Snakefile rather than hardcoding it here, so this
# script cannot drift from what the Snakefile actually checks the library against.
required_r() {
    local v
    v=$(sed -n 's/^CONTAINER_R_VERSION[[:space:]]*=[[:space:]]*"R-\([0-9][0-9.]*\)".*/\1/p' \
        "$SNAKEFILE" | head -1)
    [[ -n "$v" ]] || die "could not read CONTAINER_R_VERSION from $SNAKEFILE"
    printf '%s' "$v"
}

bioc_version() {
    sed -n 's/.*"bioconductor.version"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' \
        "$SETTINGS" 2>/dev/null | head -1
}

# major.minor only — R refuses to load packages built by a newer R *minor* version, so that
# is the granularity that actually matters.
installed_r() {
    command -v R >/dev/null 2>&1 || die "R not found on PATH."
    R --version 2>/dev/null | sed -n '1s/^R version \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p'
}

check_r_version() {
    local want have
    want=$(required_r)
    have=$(installed_r)
    [[ -n "$have" ]] || die "could not parse the output of 'R --version'."

    if [[ "$have" != "$want" ]]; then
        die "R $have is on PATH but the container image ships R $want.
  R refuses to load packages built under a different R minor version, so a library built
  here would fail inside the job pods — hours into a munge, not at launch.
  Switch to R $want before running this (module load, conda, rig, ...), then re-run."
    fi
    ok "R $have matches the image (R $want), Bioconductor $(bioc_version)"
}

# --- library verification -------------------------------------------------

# Mirrors resolve_ref_lib() in the Snakefile, including its refusal to guess between
# multiple libraries, so a pass here means the Snakefile will resolve the same path.
find_libraries() {
    local d
    for d in renv/library/*/*/*; do
        [[ -d "$d" ]] && printf '%s\n' "$d"
    done
}

verify() {
    local failed=0 libs lib count pkg r_component

    info "Verifying reference library:"

    mapfile -t libs < <(find_libraries)
    count=${#libs[@]}

    if (( count == 0 )); then
        bad "no renv library found at renv/library/<platform>/<R-version>/<arch>"
        info ""
        info "Run './setup_renv.sh' (without 'check') to build it."
        return 1
    fi

    if (( count > 1 )); then
        bad "$count renv libraries found — the Snakefile refuses to guess between them:"
        printf '        %s\n' "${libs[@]}"
        info ""
        info "Delete the stale one, or name the right one with 'ref_lib:' in"
        info "config/config_munge.yaml."
        return 1
    fi

    lib="${libs[0]}"
    ok "single library: $(cd "$lib" && pwd)"

    # The Snakefile only warns about this; treat it as fatal here, since the whole point of
    # this script is to catch it before a munge burns hours to discover it.
    r_component=$(printf '%s' "$lib" | tr '/' '\n' | grep -E '^R-[0-9]+\.[0-9]+$' | head -1)
    if [[ -n "$r_component" && "$r_component" != "R-$(required_r)" ]]; then
        bad "library is keyed to $r_component but the image ships R-$(required_r)"
        failed=1
    fi

    # renv library entries are symlinks into the shared cache, so the dangling case has to be
    # tested before -d: a dangling symlink fails -d (which follows the link) and would
    # otherwise be reported as simply missing. The distinction matters — an unmounted cache
    # and a never-installed package need entirely different fixes.
    for pkg in "${REF_PACKAGES[@]}"; do
        if [[ -L "$lib/$pkg" && ! -e "$lib/$pkg" ]]; then
            bad "$pkg is a dangling symlink -> $(readlink "$lib/$pkg")"
            info "        the shared cache is not mounted at ${RENV_PATHS_CACHE:-/renv}"
            failed=1
        elif [[ -d "$lib/$pkg" ]]; then
            ok "$pkg"
        elif [[ -e "$lib/$pkg" ]]; then
            bad "$pkg exists but is not a directory"
            failed=1
        else
            bad "$pkg missing"
            failed=1
        fi
    done

    if (( failed )); then
        info ""
        info "Re-run './setup_renv.sh' to repair, or see the Dockerfile header for the"
        info "standalone BiocManager install if you would rather not use renv."
        return 1
    fi

    info ""
    info "Library is ready. Confirm the Snakefile discovers it with:"
    # Not 'head': it closes the pipe early and snakemake turns that into a BrokenPipeError
    # traceback instead of exiting quietly. grep reads the stream to the end.
    info "    ./run_munge.sh dryrun 2>&1 | grep 'Reference package library:'"
    return 0
}

# --- restore ---------------------------------------------------------------

do_restore() {
    local cache log stamp
    cache="${RENV_PATHS_CACHE:-/renv}"
    export RENV_PATHS_CACHE="$cache"

    # The container image sets this FALSE so the project's renv .Rprofile cannot hijack
    # .libPaths() inside job pods and hide the image's MungeSumstats. That is correct there
    # and wrong here: restore needs the autoloader to bootstrap renv itself.
    export RENV_CONFIG_AUTOLOADER_ENABLED=TRUE

    if [[ -d "$cache" ]]; then
        ok "package cache: $cache"
    else
        info "  warn  cache $cache not present — packages will be compiled from source,"
        info "        which can take hours rather than minutes for ~13 GB of reference data."
    fi

    mkdir -p "$LOG_DIR"
    stamp=$(date +%Y%m%dT%H%M%S)
    log="$LOG_DIR/renv_setup_$stamp.log"

    info ""
    info "Running renv::restore() — logging to $log"
    info ""

    # tee so a failed restore leaves a record to read afterwards, while still showing
    # progress live.
    if ! R --no-save --no-restore -q -e 'renv::restore(prompt = FALSE)' 2>&1 | tee "$log"; then
        die "renv::restore() failed. Full output: $log"
    fi

    info ""
}

usage() {
    sed -n '3,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

case "${1:-}" in
    "")
        check_r_version
        do_restore
        verify
        ;;
    check)
        check_r_version
        verify
        ;;
    -h|--help|help) usage ;;
    *) die "unknown argument '$1'. Run '$0 --help'." ;;
esac
