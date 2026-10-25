#!/usr/bin/env bash
#
# Launch and control the pipeline_finemap Snakemake run.
#
# A full-studies run outlives any terminal session, and unlike pipeline_munge this pipeline
# has no Kubernetes/container profile: every job runs locally, sharing this one workspace's
# memory. get_associations loads two large Bioconductor annotation packages
# (SNPlocs.Hsapiens.dbSNP155.GRCh37/GRCh38) per study to build-check lead SNP positions, and
# with no concurrency cap, enough of those jobs running at once can exceed this workspace's
# memory cgroup and take the whole container down (tmux server included) rather than just
# failing a job. This wrapper starts the run detached under nohup (not tmux — the tmux server
# is itself a single point of failure, and we watched it die independently of the pipeline),
# records the driver PID, and caps that rule's concurrency by default.
#
#   ./run_finemap.sh start [extra snakemake args...]   launch detached
#   ./run_finemap.sh dryrun [extra snakemake args...]  -n in the foreground
#   ./run_finemap.sh status                            is it alive, what's pending
#   ./run_finemap.sh log                               follow the live log
#   ./run_finemap.sh cancel                            graceful stop (SIGINT)
#   ./run_finemap.sh unlock                            release a stale .snakemake lock
#
# Environment overrides:
#   SNAKEMAKE      path to the snakemake binary   (default: first on PATH)
#   CONDA_ENV      conda env to activate first    (default: bri-snakemake, activated only
#                                                  when snakemake is not already on PATH.
#                                                  Takes a name or a full prefix path.)
#
# Concurrency caps:
#   JOBS           jobs in flight overall         (default: 8)
#   MAX_JOBS       hard ceiling on JOBS            (default: 16)
#   MAX_SNPLOCS    simultaneous get_associations jobs loading SNPlocs.Hsapiens.dbSNP155
#                                                  (default: 2)
#   MAX_SUBSET     simultaneous subset_gwas_locus / rescue_missing_ld_neighbors jobs, each of
#                                                  which reads a whole study's munged sumstats
#                                                  (or cut-SNP) file into memory (default: 1)
#
# MAX_SNPLOCS and MAX_SUBSET exist because --jobs cannot express them. Every get_associations job
# that needs it loads a good deal of Bioconductor machinery (GenomicRanges/BSgenome/Biostrings/
# DelayedArray/SummarizedExperiment); every subset_gwas_locus job re-reads the entire munged
# sumstats file fresh (observed at ~9.9 GB RSS for the largest study) with no reuse across a
# study's loci. With no per-job pod isolation, running many of either at once is what took the
# workspace down — subset_gwas_locus specifically caused a real OOM pod restart (not just a close
# call), which is why MAX_SUBSET defaults to 1 rather than 2. The effective limit for a rule is
# the smaller of its own cap and JOBS.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

# pipeline_finemap/Snakefile declares `workdir: "../"`, which expects to be launched from
# inside pipeline_finemap/ (unlike pipeline_munge/Snakefile, which derives an absolute workdir
# from workflow.basedir). Every function below that actually invokes snakemake cds there first.
FINEMAP_DIR="$PROJECT_ROOT/pipeline_finemap"
SNAKEFILE="Snakefile"
CONFIG="../config/config_finemap.yaml"
LOG_DIR="$PROJECT_ROOT/logs"
PID_FILE="$LOG_DIR/finemap.pid"
LATEST_LOG="$LOG_DIR/finemap_latest.log"

# Activated only when snakemake is not already available — see activate_conda_env.
DEFAULT_CONDA_ENV="${DEFAULT_CONDA_ENV:-bri-snakemake}"

# Concurrency caps.
#
# JOBS is snakemake's --jobs: the ceiling on jobs in flight overall. MAX_JOBS is a hard
# ceiling on JOBS itself, so a typo or an optimistic setting cannot flood the machine.
#
# MAX_SNPLOCS and MAX_SUBSET are what usually matter more, since --jobs alone cannot express
# them. Both are enforced by the driver's scheduler via --resources.
JOBS="${JOBS:-8}"
MAX_JOBS="${MAX_JOBS:-16}"
MAX_SNPLOCS="${MAX_SNPLOCS:-2}"   # simultaneous get_associations SNPlocs loads
MAX_SUBSET="${MAX_SUBSET:-1}"     # simultaneous subset_gwas_locus / rescue_missing_ld_neighbors jobs

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
info() { printf '%s\n' "$*"; }

# --- concurrency caps ------------------------------------------------------

require_count() {
    local name=$1 val=$2
    [[ "$val" =~ ^[1-9][0-9]*$ ]] \
        || die "$name must be a positive integer, got '$val'"
}

# Validated and clamped before any jobs run, so an over-large JOBS is corrected rather than
# silently honoured.
apply_caps() {
    require_count MAX_JOBS "$MAX_JOBS"
    require_count JOBS "$JOBS"
    require_count MAX_SNPLOCS "$MAX_SNPLOCS"
    require_count MAX_SUBSET "$MAX_SUBSET"

    if (( JOBS > MAX_JOBS )); then
        info "note: JOBS=$JOBS exceeds the MAX_JOBS=$MAX_JOBS cap; using $MAX_JOBS."
        info "      Raise the cap with MAX_JOBS=<n> if that is really what you want."
        JOBS=$MAX_JOBS
    fi

    # Nothing breaks if a per-rule cap exceeds JOBS — snakemake takes the minimum — but saying
    # so avoids the impression that a higher number is in effect.
    if (( MAX_SNPLOCS > JOBS )); then
        info "note: MAX_SNPLOCS=$MAX_SNPLOCS is above JOBS=$JOBS, so JOBS is the real limit."
    fi
    if (( MAX_SUBSET > JOBS )); then
        info "note: MAX_SUBSET=$MAX_SUBSET is above JOBS=$JOBS, so JOBS is the real limit."
    fi
}

# --- discovery -------------------------------------------------------------

# Report a failed activation. Fatal when the caller named the env, advisory when we were
# only falling back to the default — in that case find_snakemake reports the missing binary
# in general terms, which is more useful than a conda error to someone who never mentioned
# conda.
activation_failed() {
    local explicit=$1; shift
    (( explicit )) && die "$*"
    info "note: $*"
    return 1
}

# Activation, so a plain `./run_finemap.sh start` works in a fresh shell with no setup.
#
# Order of preference:
#   1. SNAKEMAKE — an explicit binary path, nothing to activate.
#   2. CONDA_ENV — activate exactly what the caller named (a name or a prefix path).
#   3. snakemake already on PATH — you activated something yourself, so stay out of the way.
#   4. DEFAULT_CONDA_ENV — the project's env, activated on your behalf.
#
# 3 before 4 is what keeps the default from being intrusive: it only ever fires when the run
# would otherwise fail outright, so it cannot silently swap out an interpreter you chose.
#
# `conda activate` is a shell function that conda's hook defines, and non-interactive bash
# never sources ~/.bashrc, so the hook has to be sourced explicitly here.
#
# This must run in the current shell rather than a command substitution: the PATH and
# CONDA_PREFIX it exports have to be inherited by the driver cmd_start launches.
activate_conda_env() {
    if [[ -n "${SNAKEMAKE:-}" ]]; then
        [[ -n "${CONDA_ENV:-}" ]] \
            && info "note: SNAKEMAKE is set, so CONDA_ENV=$CONDA_ENV is ignored."
        return 0
    fi

    local env explicit=1
    env="${CONDA_ENV:-}"
    if [[ -z "$env" ]]; then
        command -v snakemake >/dev/null 2>&1 && return 0
        env="$DEFAULT_CONDA_ENV"
        explicit=0
    fi

    # Already active — re-activating would only stack a second layer of the same env.
    [[ "${CONDA_DEFAULT_ENV:-}" == "$env" ]] && return 0

    local conda_exe base hook
    conda_exe="${CONDA_EXE:-$(command -v conda || true)}"
    [[ -n "$conda_exe" && -x "$conda_exe" ]] \
        || activation_failed $explicit "no conda binary was found, so '$env' cannot be activated.
  Set CONDA_EXE=/path/to/bin/conda, or skip conda entirely:
  SNAKEMAKE=/path/to/envs/$env/bin/snakemake $0 <subcommand>" || return 0

    base=$("$conda_exe" info --base) \
        || activation_failed $explicit "'$conda_exe info --base' failed, so conda's shell hook cannot be located" \
        || return 0
    hook="$base/etc/profile.d/conda.sh"
    [[ -r "$hook" ]] \
        || activation_failed $explicit "conda's shell hook is not readable at $hook" || return 0

    # conda's hook and activation scripts are not written for strict mode: they dereference
    # unset variables, which aborts under the set -u at the top of this file. Both flags go
    # off for the duration and back on immediately after.
    set +eu
    # shellcheck source=/dev/null
    source "$hook" || {
        set -eu
        activation_failed $explicit "failed to source conda's shell hook at $hook" || return 0
    }
    # A name only resolves if the env sits in one of this installation's envs_dirs, so an env
    # belonging to a second conda install shows up in `conda info --envs` as a bare path with
    # no name and cannot be activated by name. CONDA_ENV takes a prefix path too, which is
    # the way out of that: conda activate accepts either form.
    conda activate "$env" || {
        set -eu
        activation_failed $explicit "'conda activate $env' failed (conda base: $base).
  Check the name against '$conda_exe info --envs'. Envs listed there without a name live
  outside this install's envs_dirs and only activate by full path — pass that instead:
  CONDA_ENV=/full/path/to/envs/<env> $0 <subcommand>" || return 0
    }
    set -eu

    info "conda env:  $env"
}

# Snakemake builds a run header that calls getpass.getuser(), which reads LOGNAME/USER/LNAME/
# USERNAME first and only falls back to /etc/passwd. Workspace pods here run as an arbitrary
# uid with no passwd entry, so when those variables are also unset the driver aborts with
# "OSError: No username set in the environment" before it even builds the DAG. The value is
# only ever printed, so any stable non-empty string will do.
ensure_username() {
    [[ -n "${USER:-}" || -n "${LOGNAME:-}" ]] && return 0

    local name
    name=$(id -un 2>/dev/null) || name=""          # fails for the same reason getuser() does
    [[ -n "$name" ]] || name=$(basename "${HOME:-}" 2>/dev/null)
    [[ -n "$name" && "$name" != "/" ]] || name="uid-$(id -u)"

    export USER="$name" LOGNAME="$name"
    info "note: USER and LOGNAME were unset, which snakemake treats as fatal; using USER=$name."
}

find_snakemake() {
    if [[ -n "${SNAKEMAKE:-}" ]]; then
        [[ -x "$SNAKEMAKE" ]] || die "SNAKEMAKE=$SNAKEMAKE is not executable"
        printf '%s' "$SNAKEMAKE"
        return
    fi
    if command -v snakemake >/dev/null 2>&1; then
        command -v snakemake
        return
    fi
    die "snakemake not found on PATH, and activating '$DEFAULT_CONDA_ENV' did not provide it
  (any reason why is noted above). Check that the env exists and has snakemake installed:
    conda info --envs
    ls \$CONDA_PREFIX/bin/snakemake
  Then name a different env (CONDA_ENV=<name|path> $0 <subcommand>) or point at the
  binary directly (SNAKEMAKE=/path/to/bin/snakemake $0 <subcommand>)."
}

# --- pid handling ----------------------------------------------------------

# A recorded PID is only trusted if it is still alive AND still looks like our snakemake
# run. PIDs get recycled, and this script's whole job is to send signals to it.
running_pid() {
    [[ -f "$PID_FILE" ]] || return 1
    local pid
    pid=$(<"$PID_FILE")
    [[ "$pid" =~ ^[0-9]+$ ]] || return 1
    kill -0 "$pid" 2>/dev/null || return 1
    tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null | grep -q 'snakemake' || return 1
    printf '%s' "$pid"
}

# --- shared snakemake arguments -------------------------------------------

common_args() {
    printf '%s\n' \
        --snakefile "$SNAKEFILE" \
        --configfile "$CONFIG" \
        --jobs "$JOBS" \
        --keep-going \
        --rerun-incomplete \
        --resources "snplocs_slots=$MAX_SNPLOCS" "subset_slots=$MAX_SUBSET"
}

# --- subcommands -----------------------------------------------------------

cmd_start() {
    local pid snakemake_bin stamp log
    if pid=$(running_pid); then
        die "a run is already active (pid $pid). Use '$0 cancel' first, or '$0 log' to watch it."
    fi

    apply_caps
    activate_conda_env
    ensure_username
    snakemake_bin=$(find_snakemake)
    mkdir -p "$LOG_DIR"

    stamp=$(date +%Y%m%dT%H%M%S)
    log="$LOG_DIR/finemap_$stamp.log"

    local -a args
    mapfile -t args < <(common_args)

    info "snakemake:  $snakemake_bin"
    info "concurrency: $JOBS jobs (cap $MAX_JOBS), $MAX_SNPLOCS get_associations SNPlocs loads, $MAX_SUBSET subset_gwas_locus loads"
    info "log:        $log"

    # The subshell cds into FINEMAP_DIR (see the note above SNAKEFILE) and then execs into
    # nohup, which itself execs into snakemake — no fork in that chain, so $! (the backgrounded
    # subshell's PID) ends up being exactly the snakemake process's PID, same as run_munge.sh's
    # plain `nohup ... &` achieves when no cd is needed first.
    (cd "$FINEMAP_DIR" && exec nohup "$snakemake_bin" "${args[@]}" "$@" > "$log" 2>&1) &
    local driver=$!
    printf '%s\n' "$driver" > "$PID_FILE"
    ln -sfn "$log" "$LATEST_LOG"
    disown "$driver" 2>/dev/null || true

    info ""
    info "started, pid $driver. Follow with '$0 log', stop with '$0 cancel'."
}

cmd_dryrun() {
    local snakemake_bin
    apply_caps
    activate_conda_env
    ensure_username
    snakemake_bin=$(find_snakemake)
    local -a args
    mapfile -t args < <(common_args)
    ( cd "$FINEMAP_DIR" && "$snakemake_bin" "${args[@]}" -n "$@" )
}

cmd_status() {
    local pid
    if pid=$(running_pid); then
        info "RUNNING  pid $pid  (started $(ps -o lstart= -p "$pid" | sed 's/^ *//'))"
    else
        info "NOT RUNNING"
        # A pid file with no live snakemake behind it is left over from a finished or killed
        # run. Prune it so status stays quiet, rather than reporting it after every run.
        if [[ -f "$PID_FILE" ]]; then
            info "  (previous run has ended; clearing $PID_FILE)"
            rm -f "$PID_FILE"
        fi
    fi

    if [[ -e "$LATEST_LOG" ]]; then
        info ""
        info "last lines of $(readlink -f "$LATEST_LOG"):"
        tail -5 "$LATEST_LOG" | sed 's/^/  /'
    fi
}

cmd_log() {
    [[ -e "$LATEST_LOG" ]] || die "no log yet at $LATEST_LOG"
    tail -f "$LATEST_LOG"
}

# SIGINT rather than SIGTERM/SIGKILL: snakemake traps it and releases the .snakemake lock on
# the way out. Killing it outright leaves the lock behind, and the next start then fails on a
# stale lock.
cmd_cancel() {
    local pid
    if pid=$(running_pid); then
        info "sending SIGINT to pid $pid (snakemake will finish its running jobs and exit)..."
        kill -INT "$pid" 2>/dev/null || true

        local waited=0
        while (( waited < 120 )) && kill -0 "$pid" 2>/dev/null; do
            sleep 2
            waited=$((waited + 2))
        done

        if kill -0 "$pid" 2>/dev/null; then
            info "still alive after ${waited}s, escalating to SIGTERM..."
            kill -TERM "$pid" 2>/dev/null || true
            sleep 5
            kill -0 "$pid" 2>/dev/null && kill -KILL "$pid" 2>/dev/null || true
            info "killed."
        else
            info "driver exited cleanly after ${waited}s."
        fi
    else
        info "no active run recorded."
    fi
    rm -f "$PID_FILE"
}

cmd_unlock() {
    local snakemake_bin
    if running_pid >/dev/null; then
        die "a run is still active — cancel it before unlocking."
    fi
    activate_conda_env
    ensure_username
    snakemake_bin=$(find_snakemake)
    ( cd "$FINEMAP_DIR" && "$snakemake_bin" --snakefile "$SNAKEFILE" --configfile "$CONFIG" --unlock )
}

usage() {
    sed -n '3,44p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

case "${1:-}" in
    start)  shift; cmd_start "$@" ;;
    dryrun) shift; cmd_dryrun "$@" ;;
    status) shift; cmd_status ;;
    log)    shift; cmd_log ;;
    cancel) shift; cmd_cancel ;;
    unlock) shift; cmd_unlock ;;
    ""|-h|--help|help) usage ;;
    *) die "unknown subcommand '$1'. Run '$0 --help'." ;;
esac
