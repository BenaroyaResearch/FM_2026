#!/usr/bin/env python3
"""Read-only audit of every study's pipeline_finemap status.

Classifies each study in config/study_table.tsv as not-run / partial / complete, and for
complete studies, tallies each eligible locus's QC_STATUS as a proxy for confidence (the old
composite FINAL SCORE / HIGH-MEDIUM-LOW-UNRELIABLE system this used to rely on is gone).

Purely reads results/ and config/ -- never writes to either, never triggers or suggests a
rerun. Safe to run anytime, as often as wanted. Run from the repo root:

    python3 pipeline_finemap/scripts/audit_study_status.py
"""

import argparse
import glob
import os

import pandas as pd
import yaml

# The old scoring system that used to define "confident" fine-mapping is gone (see
# pipeline_finemap/scripts/locus_summary.R and build_final_report.R, rewritten 2026-08-21).
# This is a placeholder proxy -- the narrowest, highest-precision definition available from
# what run_susie.R still computes (QC_STATUS, driven by config's lambda_warn_threshold /
# lambda_fail_threshold) -- pending a definition the user actually wants. Adjust here if a
# broader definition (e.g. also counting QC_WARN, or factoring in R2_VS_OPENTARGETS /
# N_SNPS_PIP_APPROX_1 from each locus's *_summary.txt) is wanted later.
NOT_CONFIDENT_QC_STATUSES = {"QC_FAIL"}

STUDY_TABLE_PATH = "config/study_table.tsv"
CONFIG_PATH = "config/config_finemap.yaml"
RESULTS_DIR = "results"


def load_config_studies(config_path):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    return set(cfg.get("studies", []))


def munge_done(gcst):
    return os.path.isfile(f"{RESULTS_DIR}/{gcst}/munged/{gcst}_fm_decision.txt")


def finemap_status(gcst):
    assoc_dir = f"{RESULTS_DIR}/{gcst}/associations"
    if not os.path.isdir(assoc_dir):
        return "NOT_RUN"

    final_report = f"{RESULTS_DIR}/{gcst}/susieR/{gcst}_final_report.txt"
    finalized_marker = f"{RESULTS_DIR}/{gcst}/susieR/study_finalized.txt"

    report_exists = os.path.isfile(final_report)
    finalized_exists = os.path.isfile(finalized_marker)

    report_clean = True
    if report_exists:
        with open(final_report, encoding="utf-8", errors="replace") as f:
            # The exact marker build_final_report.R writes when a locus's summary file goes
            # missing unexpectedly -- a report containing it is self-flagged as inconsistent.
            report_clean = "unexpected — investigate" not in f.read()

    if report_exists and finalized_exists and report_clean:
        return "COMPLETE"
    return "PARTIAL"


def confidence_breakdown(gcst):
    counts = {"QC_PASS": 0, "QC_WARN": 0, "QC_FAIL": 0, "NOT_RUN": 0, "OTHER": 0}
    pattern = f"{RESULTS_DIR}/{gcst}/*/susieR/*_susie_summary.txt"
    for path in glob.glob(pattern):
        status = None
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("QC_STATUS:"):
                    status = line.split(":", 1)[1].strip()
                    break
        counts[status if status in counts else "OTHER"] += 1

    n_loci = sum(counts.values())
    has_not_confident = any(
        counts.get(status, 0) > 0 for status in NOT_CONFIDENT_QC_STATUSES
    )
    return counts, n_loci, has_not_confident


def audit(study_table_path, config_path):
    study_table = pd.read_csv(study_table_path, sep="\t")
    active_studies = load_config_studies(config_path)

    rows = []
    for _, study_row in study_table.iterrows():
        gcst = study_row["gcst_id"]
        row = {
            "gcst_id": gcst,
            "trait": study_row.get("Trait"),
            "in_config": gcst in active_studies,
            "munge_done": munge_done(gcst),
        }
        row["finemap_status"] = finemap_status(gcst)

        if row["finemap_status"] == "COMPLETE":
            counts, n_loci, has_not_confident = confidence_breakdown(gcst)
            row["n_eligible_loci"] = n_loci
            row["n_qc_pass"] = counts["QC_PASS"]
            row["n_qc_warn"] = counts["QC_WARN"]
            row["n_qc_fail"] = counts["QC_FAIL"]
            row["n_qc_not_run"] = counts["NOT_RUN"]
            row["has_not_confident_loci"] = has_not_confident
        else:
            row["n_eligible_loci"] = None
            row["n_qc_pass"] = None
            row["n_qc_warn"] = None
            row["n_qc_fail"] = None
            row["n_qc_not_run"] = None
            row["has_not_confident_loci"] = None

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-table", default=STUDY_TABLE_PATH)
    parser.add_argument("--config", default=CONFIG_PATH)
    parser.add_argument(
        "--output",
        default=f"{RESULTS_DIR}/checking/study_pipeline_audit.csv",
        help="Where to write the audit CSV (default: results/checking/study_pipeline_audit.csv)",
    )
    args = parser.parse_args()

    result = audit(args.study_table, args.config)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, index=False)

    n = len(result)
    by_status = result["finemap_status"].value_counts().to_dict()
    n_not_confident = int((result["has_not_confident_loci"] == True).sum())  # noqa: E712 (NaN-safe: != True for both False and NaN)
    n_not_in_config = int((~result["in_config"]).sum())
    n_munge_pending = int((~result["munge_done"]).sum())

    print(f"Audited {n} studies -> {args.output}")
    print(f"  FINEMAP_STATUS: {by_status}")
    print(f"  not in config_finemap.yaml's studies list: {n_not_in_config}")
    print(f"  pipeline_munge not yet complete: {n_munge_pending}")
    print(f"  COMPLETE studies with >=1 QC_FAIL locus: {n_not_confident}")


if __name__ == "__main__":
    main()
