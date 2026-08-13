#!/usr/bin/env bash
#
# Launch and control the pipeline_munge Snakemake run.
#
# The munge pipeline downloads and munges whole GWASs, so a run outlives any terminal
# session: expect tens of minutes per download and up to munge_timeout_hours per munge.
# This wrapper starts it detached under nohup, records the driver PID, and gives you a
# cancel path that stops the Kubernetes jobs too instead of orphaning them.
#
#   ./run_munge.sh start [extra snakemake args...]   launch detached
#   ./run_munge.sh dryrun [extra snakemake args...]  -n in the foreground, no pods
#   ./run_munge.sh status                            is it alive, what is pending
#   ./run_munge.sh log                               follow the live log
#   ./run_munge.sh cancel                            graceful stop (SIGINT)
#   ./run_munge.sh cancel --force                    also delete leftover k8s jobs
#   ./run_munge.sh unlock                            release a stale .snakemake lock
#
# Environment overrides:
#   SNAKEMAKE      path to the snakemake binary   (default: first on PATH)
#   CONDA_ENV      conda env to activate first    (default: bri-snakemake, activated only
#                                                  when snakemake is not already on PATH.
#                                                  Takes a name or a full prefix path.)
#   PROFILE        snakemake profile name         (default: coder)
#   IMAGE          container image                (default: the container: key in
#                                                  config/config_munge.yaml, so the flag the
#                                                  k8s executor needs cannot drift from the
#                                                  directive the Snakefile declares)
#
# Concurrency caps:
#   JOBS           jobs in flight overall         (default: 4)
#   MAX_JOBS       hard ceiling on JOBS           (default: 8)
#   MAX_DOWNLOADS  simultaneous GWAS downloads    (default: 2)
#   MAX_MUNGES     simultaneous munges            (default: 2)
#
# The per-rule caps exist because --jobs cannot express them: every munge asks for 64 GB,
# and every download is a multi-GB pull from the same EBI FTP server. The effective limit
# for a rule is the smaller of its own cap and JOBS.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

CONFIG="config/config_munge.yaml"
SNAKEFILE="pipeline_munge/Snakefile"
LOG_DIR="$PROJECT_ROOT/logs"
TMP_DIR="$PROJECT_ROOT/.tmp"
PID_FILE="$LOG_DIR/munge.pid"
LATEST_LOG="$LOG_DIR/munge_latest.log"

PROFILE="${PROFILE:-coder}"

# Activated only when snakemake is not already available — see activate_conda_env.
DEFAULT_CONDA_ENV="${DEFAULT_CONDA_ENV:-bri-snakemake}"

# Concurrency caps.
#
# JOBS is snakemake's --jobs: the ceiling on jobs in flight overall. MAX_JOBS is a hard
# ceiling on JOBS itself, so a typo or an optimistic setting cannot flood the cluster with
# pods — every job here is a multi-GB download or a 64 GB munge, not a cheap task.
#
# The two per-rule caps are what usually matter more, since --jobs alone cannot express
# them. They are enforced by the driver's scheduler via --resources, so they apply on any
# executor. Snakemake takes the effective limit as the minimum of these and JOBS.
JOBS="${JOBS:-4}"
MAX_JOBS="${MAX_JOBS:-8}"
MAX_DOWNLOADS="${MAX_DOWNLOADS:-2}"   # simultaneous EBI FTP pulls
MAX_MUNGES="${MAX_MUNGES:-2}"         # simultaneous munges, i.e. 2 x munge_gwas mem_mb

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
info() { printf '%s\n' "$*"; }

# --- concurrency caps ------------------------------------------------------

require_count() {
    local name=$1 val=$2
    [[ "$val" =~ ^[1-9][0-9]*$ ]] \
        || die "$name must be a positive integer, got '$val'"
}

# Validated and clamped before any pods are created, so an over-large JOBS is corrected
# rather than silently honoured.
apply_caps() {
    require_count MAX_JOBS "$MAX_JOBS"
    require_count JOBS "$JOBS"
    require_count MAX_DOWNLOADS "$MAX_DOWNLOADS"
    require_count MAX_MUNGES "$MAX_MUNGES"

    if (( JOBS > MAX_JOBS )); then
        info "note: JOBS=$JOBS exceeds the MAX_JOBS=$MAX_JOBS cap; using $MAX_JOBS."
        info "      Raise the cap with MAX_JOBS=<n> if that is really what you want."
        JOBS=$MAX_JOBS
    fi

    # Nothing breaks if a per-rule cap exceeds JOBS — snakemake takes the minimum — but
    # saying so avoids the impression that a higher number is in effect.
    if (( MAX_DOWNLOADS > JOBS )); then
        info "note: MAX_DOWNLOADS=$MAX_DOWNLOADS is above JOBS=$JOBS, so JOBS is the real limit."
    fi
    if (( MAX_MUNGES > JOBS )); then
        info "note: MAX_MUNGES=$MAX_MUNGES is above JOBS=$JOBS, so JOBS is the real limit."
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

# Activation, so a plain `./run_munge.sh start` works in a fresh shell with no setup.
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

# The k8s executor plugin ignores the Snakefile's container: directive and needs the image
# named on the command line instead. Reading it back out of the config keeps the two in
# agreement automatically.
find_image() {
    if [[ -n "${IMAGE:-}" ]]; then
        printf '%s' "$IMAGE"
        return
    fi
    local img
    img=$(sed -n 's/^[[:space:]]*container:[[:space:]]*["'"'"']\{0,1\}\([^"'"'"'#[:space:]]\{1,\}\).*/\1/p' \
          "$CONFIG" | head -1)
    [[ -n "$img" ]] || die "no 'container:' key found in $CONFIG (set IMAGE=... to override)"
    printf '%s' "$img"
}

kube_namespace() {
    local pcfg="$HOME/.config/snakemake/$PROFILE/config.yaml"
    if [[ -r "$pcfg" ]]; then
        sed -n 's/^[[:space:]]*kubernetes-namespace:[[:space:]]*\([^[:space:]#]\{1,\}\).*/\1/p' \
            "$pcfg" | head -1
    fi
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

# tmpdir is deliberately on the shared filesystem: check_gwas_coverage decompresses the
# entire GWAS through mktemp, which would otherwise fill the pod's ephemeral storage.
# Snakemake exports resources.tmpdir as $TMPDIR and mktemp honors it. The nested quoting
# is required — snakemake evaluates default-resources values as Python expressions.
#
# munge_gwas's own mem_mb comes from the escalating schedule in the Snakefile, so it is
# deliberately not set here; mem_mb below is only the floor for the small rules.
common_args() {
    printf '%s\n' \
        --snakefile "$SNAKEFILE" \
        --profile "$PROFILE" \
        --jobs "$JOBS" \
        --keep-going \
        --resources "downloads=$MAX_DOWNLOADS" "munge_slots=$MAX_MUNGES" \
        --default-resources "tmpdir='$TMP_DIR'" mem_mb=4000 disk_mb=20000
}

# --- subcommands -----------------------------------------------------------

cmd_start() {
    local pid image snakemake_bin stamp log
    if pid=$(running_pid); then
        die "a run is already active (pid $pid). Use '$0 cancel' first, or '$0 log' to watch it."
    fi

    apply_caps
    activate_conda_env
    ensure_username
    snakemake_bin=$(find_snakemake)
    image=$(find_image)
    mkdir -p "$LOG_DIR" "$TMP_DIR"

    stamp=$(date +%Y%m%dT%H%M%S)
    log="$LOG_DIR/munge_$stamp.log"

    local -a args
    mapfile -t args < <(common_args)
    args+=(--container-image "$image")

    info "snakemake:  $snakemake_bin"
    info "image:      $image"
    info "profile:    $PROFILE"
    info "concurrency: $JOBS jobs (cap $MAX_JOBS), $MAX_DOWNLOADS downloads, $MAX_MUNGES munges"
    info "tmpdir:     $TMP_DIR"
    info "log:        $log"

    nohup "$snakemake_bin" "${args[@]}" "$@" > "$log" 2>&1 &
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
    mkdir -p "$TMP_DIR"
    local -a args
    mapfile -t args < <(common_args)
    # No --container-image: a dry run creates no pods, and requiring the profile's
    # executor here would mean needing cluster access just to check the DAG.
    "$snakemake_bin" "${args[@]}" -n "$@"
}

cmd_status() {
    local pid ns
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

    ns=$(kube_namespace)
    info ""
    if [[ -z "$ns" ]]; then
        info "profile $PROFILE declares no kubernetes-namespace; skipping job listing."
    elif ! command -v kubectl >/dev/null 2>&1; then
        info "kubectl not on PATH — cannot list jobs in namespace $ns."
    else
        info "kubernetes jobs in namespace $ns:"
        kubectl get jobs -n "$ns" 2>/dev/null | grep -E 'NAME|snakejob-' | sed 's/^/  /' \
            || info "  none"
    fi
}

cmd_log() {
    [[ -e "$LATEST_LOG" ]] || die "no log yet at $LATEST_LOG"
    tail -f "$LATEST_LOG"
}

# SIGINT rather than SIGTERM/SIGKILL: snakemake traps it, cancels the Kubernetes jobs it
# submitted, and releases the .snakemake lock on the way out. Killing it outright leaves
# both behind, and the next start then fails on a stale lock.
cmd_cancel() {
    local force=0
    [[ "${1:-}" == "--force" ]] && force=1

    local pid
    if pid=$(running_pid); then
        info "sending SIGINT to pid $pid (snakemake will cancel its running jobs)..."
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
            info "killed. Jobs may have been orphaned — check '$0 status'."
        else
            info "driver exited cleanly after ${waited}s."
        fi
    else
        info "no active run recorded."
    fi
    rm -f "$PID_FILE"

    # Orphan sweep is opt-in: the name pattern cannot distinguish this run's jobs from a
    # concurrent one's, so deleting them by default could take out someone else's work.
    local ns leftover
    ns=$(kube_namespace)
    if [[ -n "$ns" ]] && ! command -v kubectl >/dev/null 2>&1; then
        info ""
        info "Note: kubectl is not on PATH, so leftover jobs in namespace $ns could not be"
        info "checked. A clean SIGINT exit means snakemake cancelled them itself; after an"
        info "escalated kill, verify with kubectl from a shell that has it."
    elif [[ -n "$ns" ]]; then
        leftover=$(kubectl get jobs -n "$ns" -o name 2>/dev/null | grep 'snakejob-' || true)
        if [[ -n "$leftover" ]]; then
            if (( force )); then
                info "deleting leftover jobs in $ns:"
                printf '%s\n' "$leftover" | sed 's/^/  /'
                printf '%s\n' "$leftover" | xargs -r kubectl delete -n "$ns"
            else
                info ""
                info "leftover jobs still in namespace $ns:"
                printf '%s\n' "$leftover" | sed 's/^/  /'
                info "These may belong to another run. Remove them with '$0 cancel --force'."
            fi
        fi
    fi
}

cmd_unlock() {
    local snakemake_bin
    if running_pid >/dev/null; then
        die "a run is still active — cancel it before unlocking."
    fi
    activate_conda_env
    ensure_username
    snakemake_bin=$(find_snakemake)
    "$snakemake_bin" --snakefile "$SNAKEFILE" --unlock
}

usage() {
    sed -n '3,37p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

case "${1:-}" in
    start)  shift; cmd_start "$@" ;;
    dryrun) shift; cmd_dryrun "$@" ;;
    status) shift; cmd_status ;;
    log)    shift; cmd_log ;;
    cancel) shift; cmd_cancel "${1:-}" ;;
    unlock) shift; cmd_unlock ;;
    ""|-h|--help|help) usage ;;
    *) die "unknown subcommand '$1'. Run '$0 --help'." ;;
esac
