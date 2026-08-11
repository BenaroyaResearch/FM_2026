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
#   SNAKEMAKE   path to the snakemake binary   (default: first on PATH)
#   PROFILE     snakemake profile name         (default: coder)
#   JOBS        concurrent jobs                (default: 4)
#   IMAGE       container image                (default: the container: key in
#                                               config/config_munge.yaml, so the flag the
#                                               k8s executor needs cannot drift from the
#                                               directive the Snakefile declares)

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
JOBS="${JOBS:-4}"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
info() { printf '%s\n' "$*"; }

# --- discovery -------------------------------------------------------------

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
    die "snakemake not found on PATH.
  Activate the environment first (conda activate snakemake),
  or point at it directly: SNAKEMAKE=/path/to/bin/snakemake $0 $*"
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
        --default-resources "tmpdir='$TMP_DIR'" mem_mb=4000 disk_mb=20000
}

# --- subcommands -----------------------------------------------------------

cmd_start() {
    local pid image snakemake_bin stamp log
    if pid=$(running_pid); then
        die "a run is already active (pid $pid). Use '$0 cancel' first, or '$0 log' to watch it."
    fi

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
    info "profile:    $PROFILE (jobs=$JOBS)"
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
    snakemake_bin=$(find_snakemake)
    if running_pid >/dev/null; then
        die "a run is still active — cancel it before unlocking."
    fi
    "$snakemake_bin" --snakefile "$SNAKEFILE" --unlock
}

usage() {
    sed -n '3,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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
