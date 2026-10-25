# FM_2026: automated GWAS fine-mapping with SuSiE

## Purpose

A genome-wide association study (GWAS) finds regions of the genome linked to a disease, but each
region usually holds dozens to thousands of SNPs that are inherited together (they are in
"linkage disequilibrium", or LD). The GWAS can tell you a region matters; it cannot tell you
which SNP in that region is actually responsible. **Fine-mapping** answers that question
statistically. For every SNP in a region it estimates the probability that the SNP is causal
(its **PIP**), and it groups the most likely candidates into **credible sets** — small lists of
SNPs that together very probably contain a causal variant. To do this it needs two things: the
GWAS results for the region, and a picture of how the SNPs in that region are correlated, which
we take from a reference population (the 1000 Genomes Project).

Doing this by hand for one study is tedious; doing it for hundreds is not realistic. Every
study in the GWAS Catalog is formatted a little differently, some are on an older genome build,
some are missing columns fine-mapping needs, the files are several gigabytes each, and every
study has tens of regions that each need their own reference data, alignment and checks. This
repository automates all of it. You give it a list of study IDs; it downloads each study,
checks whether it is usable, cleans it into one standard format, fine-maps every lead SNP,
runs quality checks, compares the answers to the published fine-mapping in Open Targets, and
writes a plain-text report saying what was fine-mapped and, for anything that was not, why.

The work is split into two pipelines that run one after the other, each started with a single
command:

```
study list  ──►  munge pipeline (run_munge.sh)  ──►  fine-mapping pipeline (run_finemap.sh)  ──►  results/{study}/
                 download, check, clean             per-region LD, SuSiE, QC, Open Targets
```

- **Munge pipeline** (`pipeline_munge/`): downloads each study's summary statistics, checks
  they contain the data required for fine-mapping, converts them to a standard format on genome build
  GRCh38, and consults parameters to decide if the study can be fine-mapped.
- **Fine-mapping pipeline** (`pipeline_finemap/`): for every lead SNP, subsets a window of
  the genome, builds the matching LD matrix from 1000 Genomes, fine-maps it with SuSiE, checks the
  result, and compares it with Open Targets (if possible).
- **What you get per study:** a SNP-level results table, a human-readable final report, QC
  summaries for every region, Open Targets comparison tables and plots, and a log of every
  step.

### Contents

1. [Setup: R reference library (renv)](#setup-r-reference-library-renv)
2. [Setup: other software and reference files](#setup-other-software-and-reference-files)
3. [Running: run_munge.sh](#running-run_mungesh)
4. [Running: run_finemap.sh](#running-run_finemapsh)
5. [Running a batch from start to finish](#running-a-batch-from-start-to-finish)
6. [Config files](#config-files)
7. [Choosing studies: the study table](#choosing-studies-the-study-table)
8. [The munge pipeline, step by step](#the-munge-pipeline-step-by-step)
9. [The fine-mapping pipeline, step by step](#the-fine-mapping-pipeline-step-by-step)
10. [File tree map](#file-tree-map)
11. [Final files, explained](#final-files-explained)
12. [Known issues and ideas](#known-issues-and-ideas)

## Setup: R reference library (renv)

`munge_sumstats.R` needs four SNPlocs/BSgenome reference data packages. They come to
roughly 13 GB, which is too large to ship inside the container image, so they live in the
project's renv library and are mounted into the job pods instead. The container carries
MungeSumstats; the renv library carries the reference data. Both are needed.

The library **must be built against R 4.4 / Bioconductor 3.20**, matching the image
(`bioconductor/bioconductor_docker:RELEASE_3_20`, R 4.4.2). R refuses to load packages
built by a newer R minor version, so a mismatch is a guaranteed failure — and it surfaces
hours in, at the end of a long munge, not at launch.

Run `setup_renv.sh` once per checkout, from the project root:

```bash
./setup_renv.sh            # restore the library, then verify it
./setup_renv.sh check      # verify an existing setup, changes nothing
```

It checks that R matches the image before doing any work, restores the library, then
verifies the result and exits non-zero if anything is wrong. Both modes are safe to re-run.
`renv::restore()` output is logged to `logs/renv_setup_<timestamp>.log`.

A healthy `check` looks like:

```
  ok    R 4.4 matches the image (R 4.4), Bioconductor 3.20
Verifying reference library:
  ok    single library: /home/<you>/FM_2026/renv/library/linux-ubuntu-noble/R-4.4/x86_64-pc-linux-gnu
  ok    SNPlocs.Hsapiens.dbSNP155.GRCh38
  ok    SNPlocs.Hsapiens.dbSNP155.GRCh37
  ok    BSgenome.Hsapiens.NCBI.GRCh38
  ok    BSgenome.Hsapiens.1000genomes.hs37d5
```

Then confirm the Snakefile discovers the same path:

```bash
./run_munge.sh dryrun 2>&1 | grep -E 'Running [0-9]+ studies:|Reference package library:'
# Running 1 studies: ['GCST90132226']
# Reference package library: /home/<you>/FM_2026/renv/library/linux-ubuntu-noble/R-4.4/x86_64-pc-linux-gnu
```

Do not pipe a dryrun through `head`. It exits after its line count and closes the pipe, and
because Python ignores `SIGPIPE`, snakemake raises `BrokenPipeError` from inside its logger
and buries the output you wanted in two tracebacks. `grep` and `sed -n '1,20p'` read to the
end of the stream, so they are safe. The `2>&1` matters too: snakemake logs to stderr, so an
unredirected pipe sees only part of the output.

#### What it does, and what it catches

The manual equivalent is three steps — confirm R is 4.4.x, point `RENV_PATHS_CACHE` at the
shared cache (without it, restore compiles ~13 GB from source instead of symlinking out of
the cache in minutes), then `R -e 'renv::restore(prompt = FALSE)'` from the project root,
where `.Rprofile` bootstraps renv automatically.

The value is in the checks around that, each of which otherwise surfaces hours into a munge
rather than at setup time:

- **Wrong R.** Refuses to restore at all if R's major.minor does not match the image, since
  R will not load packages built under a different R minor version. The expected version is
  read out of `CONTAINER_R_VERSION` in the Snakefile, so the script cannot drift from what
  the Snakefile checks against.
- **No library** — `renv::restore()` never ran, or ran outside the project root.
- **Multiple libraries** matching `renv/library/*/*/*`, usually left over from an R upgrade.
  `resolve_ref_lib()` in the Snakefile refuses to guess between them, so this is fatal;
  delete the stale one or name the right one via `ref_lib:` in `config/config_munge.yaml`.
- **A library keyed to the wrong R**, e.g. `R-4.3` when the image ships `R-4.4`. The
  Snakefile only warns about this; the script treats it as a failure.
- **Dangling symlinks.** Library entries point into the shared `/renv` cache, so job pods
  must mount both the project directory and `/renv` — the `coder` profile does both. These
  are reported as dangling rather than merely missing, because an unmounted cache and a
  never-installed package need completely different fixes.

To use a standalone BiocManager library instead of renv, see the install command in the
Dockerfile header and set `ref_lib:` in `config/config_munge.yaml` (or pass
`--config ref_lib=/path`).

Note: keep credentials such as `GITHUB_PAT` in `~/.Renviron`, never in a project-level
`.Renviron` — the latter sits inside the repo and gets committed.

## Setup: other software and reference files

Do the renv setup above first. The rest of this section covers everything else the two
pipelines need. Most of it only has to be done once.

### Software

| What | Used by | Notes |
|---|---|---|
| **Snakemake 9.25.1** in a conda env | both pipelines | The launchers activate `bri-snakemake` for you if Snakemake isn't already available (see "Which snakemake it uses" below). Keep this version matched to the one in the `Dockerfile`. |
| **Kubernetes access** with the `coder` Snakemake profile (`~/.config/snakemake/coder/config.yaml`) | munge pipeline | The munge steps run as Kubernetes jobs inside the container image `ghcr.io/benaroyaresearch/fm_2026` (built from `Dockerfile` by `.github/workflows/build-image.yaml`). Nothing else needs installing for munging. |
| **PLINK 1.9** on your `PATH` | fine-mapping pipeline | Used to cut reference data and compute LD. Check with `plink --version`. |
| **R** with `susieR`, `data.table`, `ggplot2`, `ggrepel`, `stringr`, `readr`, `jsonlite`, `rtracklayer` | fine-mapping pipeline | The fine-mapping steps run on your own machine, not in the container, using the `R` on your `PATH`. The dbSNP reference packages are loaded from the renv library above. |
| **Python 3 with pandas** | study-table helper scripts | The Snakemake conda env already has it. |

Install the R packages once from an R session:

```r
install.packages(c("susieR", "data.table", "ggplot2", "ggrepel", "stringr", "readr", "jsonlite", "BiocManager"))
BiocManager::install("rtracklayer")
```

### Reference files

The fine-mapping pipeline needs four kinds of reference data. Their locations are set in
`config/config_finemap.yaml`, so they can live anywhere you like.

#### 1. The 1000 Genomes LD reference panel

LD can't be calculated from GWAS summary statistics, so we use an out-of-reference LD matrix 
estimated from 1000 Genomes individuals of the same ancestry as the study. We use the 
high-coverage (30x) release on genome build GRCh38, in PLINK format, split by ancestry 
(EUR and EAS) and chromosome.

A ready-made copy already exists on the shared drive at
`/home/gdenonno/fine_mapping/1000_genome_files/1000g_hg38_30x_dedup_plink`, and the config
points at it. To build your own:

1. **Download the VCFs.** Get the phased VCF files (and their `.tbi` index files) for
   chromosomes 1–22 from the 1000 Genomes FTP site:
   `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/`
   (files named `1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`).
   They are large; check each one downloaded completely with `gzip -t <file>`.
2. **Make the ancestry sample lists.** Download the sample information file
   `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel`.
   Write the sample IDs whose `super_pop` is `EUR` to `EUR.1000g.sample.ID.txt` and those whose
   `super_pop` is `EAS` to `EAS.1000g.sample.ID.txt`, one ID per line (503 people in each).
3. **Convert to PLINK, one file set per ancestry and chromosome**, keeping only that ancestry's
   people and only simple SNPs:

   ```bash
   plink --vcf 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
         --keep EUR.keep.txt \
         --snps-only just-acgt --biallelic-only strict --maf 0.01 \
         --keep-allele-order --make-bed \
         --out 1000g_hg38_30x_plink/EUR/1000g_hg38_30x_EUR_chr1
   ```

   (`--keep` takes a two-column file of family ID and sample ID, so give each sample ID twice
   per line.) Repeat for every chromosome and for EAS.
4. **Remove duplicated positions.** The 30x files store some sites as several records at the
   same chromosome position. The pipeline names SNPs by position, so these collide and break the
   alignment step. For each chromosome, find every position that appears more than once in the
   `.bim` file and drop all records at those positions:

   ```bash
   awk '{print $1"\t"$4}' in.bim | sort | uniq -d | awk '{print $1"\t"$2"\t"$2"\tdup"NR}' > dups.txt
   plink --bfile in --exclude range dups.txt --keep-allele-order --make-bed --out out
   ```

The finished panel must be laid out like this, which is what the config expects:

```
{thousand_genomes_plink_dir}/
├── EUR/1000g_hg38_30x_EUR_chr1.bed .bim .fam   ... chr22
└── EAS/1000g_hg38_30x_EAS_chr1.bed .bim .fam   ... chr22
{thousand_genomes_sample_dir}/
├── EUR.1000g.sample.ID.txt
└── EAS.1000g.sample.ID.txt
```

The file-name prefix (`1000g_hg38_30x`) is set by `thousand_genomes_file_prefix`.

There is no chrX in the panel, so studies with a chrX lead SNP can't finalize. See
[Known issues](#known-issues-and-ideas).

#### 2. Genome-build conversion file (liftOver chain)

Some lead SNPs in the study table are recorded on the older genome build (GRCh37). Converting
them to GRCh38 needs the UCSC chain file `hg19ToHg38.over.chain`. Download
`hg19ToHg38.over.chain.gz` from UCSC, **unzip it** (the R function that reads it rejects `.gz`),
and put it at `pipeline_finemap/reference/hg19ToHg38.over.chain` (or point
`grch37_to_grch38_chain` at it).

#### 3. Open Targets credible sets

To compare our answers with published fine-mapping, the pipeline reads the Open Targets
credible-set table (`credible_set_merged_25.03.tsv.gz`, about 2 GB), set by
`opentargets_raw_table`. The first run builds a smaller cache of just our studies
(`results/checking/opentargets_gwas_native_cache.tsv`) and reuses it after that.

#### 4. The lead-variant library

`config/2026_library_leadVariantTable_CSannotated.csv` lists every GWAS association in our
variant library and whether it has already been fine-mapped. It is the starting point for
choosing studies (see [Choosing studies](#choosing-studies-the-study-table)).


## Running: run_munge.sh

`run_munge.sh` at the project root launches the munge pipeline detached under `nohup` and
gives you a clean way to stop it. Runs are long — tens of minutes per GWAS download, up to
`munge_timeout_hours` per munge — so they need to survive the terminal closing.

Before launching: run `./setup_renv.sh` once (above), list the GCST IDs you want in
`studies:` in `config/config_munge.yaml`, and make sure each one has a row in
`config/study_table.tsv`.

```bash
./run_munge.sh dryrun             # check the DAG and library discovery, creates no pods
./run_munge.sh start              # launch detached
./run_munge.sh log                # follow the live log
./run_munge.sh status             # alive? last log lines? which k8s jobs?
./run_munge.sh cancel             # graceful stop
```

Results land in `results/{study}/` inside the project directory, and each study's log at
`results/{study}/{study}.log`. The driver's own log goes to `logs/munge_<timestamp>.log`,
with `logs/munge_latest.log` symlinked to the newest.

Extra arguments are passed straight through to snakemake, so anything not covered by the
wrapper still works:

```bash
./run_munge.sh start --rerun-incomplete
./run_munge.sh start --config munge_mem_mb_schedule="[65536, 131072]"
./run_munge.sh dryrun --config ref_lib=/tmp      # skip the library check entirely
```

Environment overrides: `SNAKEMAKE` (binary path), `CONDA_ENV` (env to activate first, see
below), `PROFILE` (default `coder`), `IMAGE` (default: the `container:` key in
`config/config_munge.yaml`).

#### Which snakemake it uses

No activation step is needed: if `snakemake` is not already available, the wrapper activates
`bri-snakemake` itself. It picks an interpreter in this order, and the first hit wins.

| Order | Source | When |
|---|---|---|
| 1 | `SNAKEMAKE` | set — an explicit binary path, no activation happens |
| 2 | `CONDA_ENV` | set — activates exactly that env (a name or a full prefix path) |
| 3 | `snakemake` on `PATH` | you already activated something; the wrapper leaves it alone |
| 4 | `bri-snakemake` | nothing else resolved, so the project's env is activated for you |

Step 3 sits ahead of step 4 deliberately: the default only fires when the run would fail
anyway, so it can never silently swap out an env you chose yourself. If you have a different
env active, that is what runs. Override the default name with `DEFAULT_CONDA_ENV` if the env
is called something else in your setup.

Activating for step 2 or 4 means sourcing conda's shell hook first, because `conda activate`
is a shell function and non-interactive bash never reads `~/.bashrc`. Set `CONDA_EXE` as well
if `conda` is not on `PATH`. A failure to activate is fatal when you named the env yourself,
and only a warning when the wrapper was falling back to the default — in that case the
"snakemake not found" error reports the real problem.

`CONDA_ENV` accepts a full prefix path as well as a name, which is what you need when more
than one conda is installed. A name only resolves if the env is inside the active install's
`envs_dirs`; `conda info --envs` prints such envs with a name and everything else as a bare
path. If yours is listed without a name, pass the path — or bypass conda with `SNAKEMAKE`:

```bash
CONDA_ENV=/home/<you>/miniforge3/envs/bri-snakemake ./run_munge.sh start
SNAKEMAKE=/home/<you>/miniforge3/envs/bri-snakemake/bin/snakemake ./run_munge.sh start
```

Activation is only needed at launch. `start` resolves the binary to an absolute path before
`nohup`ing it, so the detached driver keeps running after you close the terminal or
`conda deactivate`. `status`, `log`, and `cancel` never need the env at all — they work from
the PID file, so you can stop a run from a fresh shell.

#### Username in the environment

Snakemake's run header calls `getpass.getuser()`, which reads `LOGNAME`/`USER`/`LNAME`/
`USERNAME` and only then falls back to `/etc/passwd`. Workspace pods run as an arbitrary uid
with no passwd entry, so if those variables are unset too, the driver dies with
`OSError: No username set in the environment` before building the DAG — and the traceback
points at `getpass`, not at anything in this pipeline.

The wrapper exports `USER`/`LOGNAME` when they are missing, preferring `id -un`, then the
basename of `$HOME`, then `uid-<n>`. The value is only printed in the header, so any stable
string works. Set `USER` yourself if you want a specific name in the logs.

#### Concurrency caps

With a long `studies:` list the pipeline will happily run many studies at once, which is
rarely what you want here: each munge asks for 64 GB, and each download is a multi-GB pull
from the same EBI FTP server.

| Variable | Default | Caps |
|---|---|---|
| `JOBS` | 4 | jobs in flight overall (snakemake `--jobs`) |
| `MAX_JOBS` | 8 | hard ceiling on `JOBS` itself |
| `MAX_DOWNLOADS` | 2 | simultaneous GWAS downloads |
| `MAX_MUNGES` | 2 | simultaneous munges, i.e. 2 x 64 GB |

```bash
MAX_MUNGES=1 ./run_munge.sh start          # one munge at a time
JOBS=8 MAX_DOWNLOADS=4 ./run_munge.sh start
MAX_JOBS=16 JOBS=16 ./run_munge.sh start   # deliberately raise the ceiling
```

`JOBS` above `MAX_JOBS` is clamped with a note rather than silently honoured, and all four
are validated as positive integers before any pods are created. The two per-rule caps exist
because `--jobs` cannot express them; they are enforced by the driver's scheduler through
`--resources`, so they apply on any executor. The effective limit for a rule is the smaller
of its own cap and `JOBS`, and the wrapper says so if you set one above `JOBS`.

#### What the wrapper sets for you

- `--container-image`, read back out of `config/config_munge.yaml`. The Kubernetes executor
  ignores the Snakefile's `container:` directive and needs the image named on the command
  line, so reading the same config key keeps the two from drifting.
- `--default-resources tmpdir=<project>/.tmp`, keeping `check_gwas_coverage`'s full
  decompression of the GWAS on shared storage instead of filling the pod's ephemeral disk.
- Baseline `mem_mb`/`disk_mb` for the small rules. `munge_gwas` is deliberately left alone
  so it uses the escalating schedule declared in the Snakefile — see below.

#### Cancelling

`cancel` sends **SIGINT**, which snakemake traps: it cancels the Kubernetes jobs it
submitted and releases the `.snakemake` lock on the way out. It escalates to SIGTERM then
SIGKILL only if the driver is still alive after 120s, and in that case jobs may be
orphaned. `cancel --force` additionally deletes any leftover `snakejob-*` jobs — opt-in
because the name pattern cannot tell your jobs apart from a concurrent run's.

If a run is killed rather than interrupted, the next `start` may fail on a stale lock.
Clear it with `./run_munge.sh unlock`.

Listing and sweeping Kubernetes jobs needs `kubectl` on PATH; the wrapper says so
explicitly when it is missing rather than silently reporting nothing to clean.

#### munge_gwas memory

`munge_gwas` requests 64 GB on its first attempt and escalates on retry — 64, then 96, then
128 GB — via a `mem_mb` callable keyed to snakemake's `attempt` counter, with `retries: 2`
declared on the rule. Retrying is cheap because `download_gwas` skips a download that is
already on disk, so a retry re-runs only the munge.

Tune the schedule without editing the Snakefile:

```bash
./run_munge.sh start --config munge_mem_mb_schedule="[65536, 131072, 196608]"
```

Two caveats. On Kubernetes `mem_mb` is both the pod's memory request *and* its limit, so an
entry larger than any node's allocatable memory leaves the pod `Pending` indefinitely
rather than failing fast — keep the top entry under what a node can actually satisfy. And
`retries` applies to every failure of the rule, not just OOM kills: the pre-flight checks
fail in seconds so retrying them is cheap, but a run that hits the `munge_timeout_hours`
limit is retried too. Set `retries: 0` on the rule to disable escalation.


## Running: run_finemap.sh

`run_finemap.sh` at the project root starts the fine-mapping pipeline in the background and
gives you simple commands to watch or stop it. It works the same way as `run_munge.sh`, with
two differences: **every job runs on your own machine** (no Kubernetes), and it uses
`config/config_finemap.yaml`.

Before launching: the munge pipeline must have finished for the studies you want (each needs
`results/{study}/munged/{study}_fm_decision.txt`), the study IDs must be listed under
`studies:` in `config/config_finemap.yaml`, and each must have a row in
`config/study_table.tsv`.

```bash
./run_finemap.sh dryrun           # show what would run, without running anything
./run_finemap.sh start            # launch in the background
./run_finemap.sh log              # follow the live log (Ctrl-C stops watching, not the run)
./run_finemap.sh status           # is it still running? last few log lines
./run_finemap.sh cancel           # stop it cleanly
./run_finemap.sh unlock           # clear a stale lock after a crash
```

The driver's own log goes to `logs/finemap_<timestamp>.log`, with `logs/finemap_latest.log`
always pointing at the newest one. Each study also gets its own log at
`results/{study}/{study}.log`. A study is done when
`results/{study}/checking/susieR/{window}/study_finalized.txt` exists (for example
`.../susieR/350kb_window/study_finalized.txt`).

Anything you add after `start` or `dryrun` is passed straight to Snakemake, for example to
run just one study without editing the config:

```bash
./run_finemap.sh start --config 'studies=["GCST90476277"]'
```

#### Concurrency caps

Because everything shares one machine's memory, the wrapper limits how many heavy jobs run at
once. Set these as environment variables in front of the command.

| Variable | Default | Limits |
|---|---|---|
| `JOBS` | 8 | jobs running at once overall |
| `MAX_JOBS` | 16 | hard ceiling on `JOBS` |
| `MAX_SNPLOCS` | 2 | lead-SNP checking jobs, which each load the large dbSNP reference packages |
| `MAX_SUBSET` | 1 | jobs that read a whole study's sumstats into memory (up to ~10 GB each) |

```bash
JOBS=4 ./run_finemap.sh start
```

`JOBS=4 MAX_SNPLOCS=1 MAX_SUBSET=1` keeps fine-mapping to roughly 20 GB of memory in use, which
leaves room for other work in the same workspace.

**Reading memory use.** The workspace's memory counter (`/sys/fs/cgroup/memory.current`, and the
percentage most monitors show) includes the kernel's file cache. A run reads large sumstats and
panel files, so that number can sit at 70% or more while programs are using only a few GB. The
kernel hands cache back as soon as anything needs the memory. To see what's really in use, check
the `anon` line instead:

```bash
awk '$1=="anon"{printf "%.1f GB in use\n", $2/1e9}' /sys/fs/cgroup/memory.stat
```

`SNAKEMAKE`, `CONDA_ENV` and `DEFAULT_CONDA_ENV` work exactly as described for
`run_munge.sh`.

#### Running a different window or settings without editing the config

To fine-map with a different `locus_window_bp` (or other setting) while leaving
`config/config_finemap.yaml` alone, put the overrides in a small YAML file and pass it after the
main config:

```yaml
# logs/finemap_300kb_120.yaml
locus_window_bp: 300000
studies:
  - GCST90476277
  - ...
```

```bash
JOBS=4 MAX_SNPLOCS=1 MAX_SUBSET=1 ./run_finemap.sh dryrun \
  --configfile ../config/config_finemap.yaml ../logs/finemap_300kb_120.yaml \
  --rerun-triggers mtime
```

Then run the same command with `start` in place of `dryrun`.

- **Paths are relative to `pipeline_finemap/`**, because the wrapper starts Snakemake from
  there; hence the `../`.
- **List both config files.** A second `--configfile` replaces the wrapper's own, so leaving out
  the main config would drop its settings.
- **`--rerun-triggers mtime` protects earlier results.** The window-specific outputs go into
  their own `{N}kb_window` folders, but a few files are shared by every window, such as
  `results/{study}/associations/{study}_lead_snps.tsv`. Without this flag, any edit to a script
  or setting since the last run can make Snakemake rebuild those shared files, and so change what
  earlier windows' results were built from.
- **Check the dry run before starting.** Every planned output path should contain the new
  window's folder, and `get_associations` should be scheduled only for studies that don't yet
  have a `lead_snps.tsv`.

The 300 kb results in `results/` are from such a run (2026-10-24, 120 studies, 30x dedup panel).
Earlier 300 kb test results from September, made on the wrong panels, were moved to
`archive/300kb_window_sep2026/` before that run so they couldn't mix with it.


## Running a batch from start to finish

1. **Pick the studies.** Find their GCST IDs (for example from `config/1008_study_table.tsv`,
   see [Choosing studies](#choosing-studies-the-study-table)).
2. **Make sure each study has a row in `config/study_table.tsv`.** Both pipelines only read
   this file; an ID that isn't in it is silently skipped.
3. **List the IDs** under `studies:` in `config/config_munge.yaml`.
4. **Munge:** `./run_munge.sh dryrun`, then `./run_munge.sh start`. Watch it with
   `./run_munge.sh log`. This is the slow part (downloads and memory-heavy cleaning).
5. **Check what passed.** Each study now has `results/{study}/munged/{study}_fm_decision.txt`
   (or a reason in its log for why not). Studies marked `none` are skipped by fine-mapping.
6. **List the same IDs** under `studies:` in `config/config_finemap.yaml`.
7. **Fine-map:** `./run_finemap.sh dryrun`, then `./run_finemap.sh start`.
8. **Read the results.** Start with each study's final report
   (`results/{study}/checking/susieR/{window}/{study}_final_report.txt`) — see
   [Your final files, explained](#your-final-files-explained).

To see where every study stands, open `susie_fm_status/fm_pipeline_run_status.tsv` (or the
`.csv`). Steps 4 and 7 each refresh it automatically when the run ends; to refresh it at any
other time, run `python susie_fm_status/build_fm_status_table.py`. See
[Tracking progress](#tracking-progress).


## Config files

Each pipeline reads one config file in `config/`. These hold every setting you might want to
change, so you shouldn't need to edit the pipeline code.

### config/config_munge.yaml

| Setting | Current value | What it does |
|---|---|---|
| `studies` | list of GCST IDs | Which studies to download and munge. Each must also be in `config/study_table.tsv`. |
| `container` | `ghcr.io/benaroyaresearch/fm_2026:latest` | The software image the munge jobs run in. |
| `ref_lib` | `renv/library/linux-ubuntu-noble/R-4.4/x86_64-pc-linux-gnu` | Where the reference R packages live (see the renv setup). |
| `finemap_thresholds: snp_mid_low` | 500000 | A study with fewer SNPs than this after cleaning is marked `none` (not fine-mapped). |
| `finemap_thresholds: n_high`, `snp_high` | 100000, 2000000 | A study with more people than `n_high` **and** more SNPs than `snp_high` is marked `both` (see the munge decision below). |
| `finemap_thresholds: n_mid_low` | 30000 | Not currently used. |
| `munge_timeout_hours` | 4 | A munge job that runs longer than this is stopped. |
| `munge_mem_mb_schedule` | not set (default 64, 96, 128 GB) | Memory for the first, second and third munge attempt. Set it in a YAML file (see [Known issues](#known-issues-and-ideas)). |
| `use_functional_annotations` | false | Not currently used. |

### config/config_finemap.yaml

| Setting | Current value | What it does |
|---|---|---|
| `studies` | list of GCST IDs | Which studies to fine-map. Each must be munged and in `config/study_table.tsv`. |
| `finemapping_method` | `susieR` | Fine-mapping method: `susieR` (standard SuSiE), or the variants `susie-inf` / `susie-ash`, which also model many tiny background effects. The method name becomes part of the output folder. |
| `susie_unmappable_estimate_residual_variance` | false | Only used by `susie-inf` / `susie-ash`. |
| `locus_window_bp` | 350000 | How far either side of each lead SNP to look (350,000 bases = a 700 kb region). Results go in a folder named after it, e.g. `350kb_window`, so runs with different windows never overwrite each other. |
| `susie_max_causal` | 10 | The most causal signals SuSiE may find in one region (SuSiE's `L`). |
| `lambda_warn_threshold`, `lambda_fail_threshold` | 0.1, 0.3 | Cut-offs for the **lambda** quality score: above 0.1 a region is flagged `QC_WARN`, above 0.3 `QC_FAIL`. |
| `thousand_genomes_plink_dir` | `.../1000g_hg38_30x_dedup_plink` | The LD reference panel (see setup). |
| `thousand_genomes_file_prefix` | `1000g_hg38_30x` | File-name prefix inside that folder. |
| `thousand_genomes_sample_dir` | `/home/gdenonno/fine_mapping/1000_genome_files` | Where `EUR.1000g.sample.ID.txt` / `EAS.1000g.sample.ID.txt` live. |
| `ld_r2_threshold`, `ld_window_kb` | 0.6, 300 | For the "missing neighbours" check: SNPs within 300 kb of the lead SNP with r² ≥ 0.6 count as its close neighbours. |
| `rescue_missing_ld_neighbors` | false | If true, tries to recover close neighbours that the munge step removed. |
| `ref_lib` | renv library path | Where the dbSNP reference packages live. |
| `grch37_to_grch38_chain` | `pipeline_finemap/reference/hg19ToHg38.over.chain` | The build-conversion file (see setup). |
| `opentargets_raw_table` | `FM_checking/credible_set_merged_25.03.tsv.gz` | Open Targets credible sets used for comparison. |
| `opentargets_cs_size_exclude_threshold` | 50 | Credible sets bigger than this are left out of the set-by-set comparison. |
| `locus_summary_thresholds: pip_approx_1_threshold` | 1.00 | A SNP at or above this PIP is reported as "PIP ≈ 1". |
| `locus_summary_thresholds: opentargets_high_pip_threshold` | 0.9 | What counts as a "high" PIP when comparing with Open Targets. |
| `lead_variant_table` | `config/2026_library_leadVariantTable_CSannotated.csv` | Only used by an optional comparison step that is currently switched off. |
| `locus_output_root` | not set (`results`) | Optional: write the per-region outputs somewhere other than `results/`. |


## Choosing studies: the study table

*The munging and fine-mapping pipelines have been separated into two separate snakemake pipelines. I chose to do this because 1) I see them as two separate processes and thought a munging pipeline could be useful in the future, even unrelated to FM, 2) it makes the debugging process much easier, and 3) many studies don't make it to FM post-munge.

### Building the candidate list

`pipeline_munge/scripts/build_1008_study_table.py` picks the studies worth fine-mapping from
the lead-variant library (`config/2026_library_leadVariantTable_CSannotated.csv`). Run it from
the project root:

```bash
python pipeline_munge/scripts/build_1008_study_table.py
```

It keeps a study only if:

1. its summary statistics are available to download,
2. it is from a single ancestry group (EUR or EAS, not a mix),
3. it is a case/control (binary) trait, and
4. it has at least one association that **hasn't been fine-mapped yet** — meaning the SNP was
   not fine-mapped in an earlier immune study and is not in an existing immune credible set.

It then **ranks** the studies so the most useful come first: it repeatedly picks the study that
adds the most not-yet-fine-mapped SNPs that no higher-ranked study already covers (ties go to
the larger study). Running the studies in this order covers the most new associations soonest.

The result is `config/1008_study_table.tsv` (currently 158 studies covering 1,046 unique
not-yet-fine-mapped SNPs). Every SNP of a chosen study is kept for fine-mapping, not only the
new ones.

### The table the pipelines read: config/study_table.tsv

Both pipelines read `config/study_table.tsv` and nothing else, so **a study has to be in this
file before it can run**. It has the same first 12 columns as `1008_study_table.tsv`; to add
a study, copy its row across (the extra ranking columns at the end of the 1008 table aren't
needed). The columns are:

| Column | Meaning |
|---|---|
| `STUDY_KEY` | Study identifier from the library (PubMed ID, trait, sample description) |
| `Trait` | Disease or trait name |
| `Ancestry` | `EUR` or `EAS` — picks which 1000 Genomes reference people are used |
| `n` | Total sample size |
| `GWAS_dir_url` | GWAS Catalog folder the summary statistics are downloaded from |
| `gcst_id` | GWAS Catalog accession (e.g. `GCST90476277`) — the ID you list in the configs |
| `lead_snps` | Comma-separated lead SNPs (rsIDs, or `chr:position` when there is no rsID) |
| `pubmed_id` | Publication |
| `pct_not_finemapped` | Fraction of the study's associations not yet fine-mapped |
| `lead_snp_positions` | Comma-separated `chr:position` for each lead SNP, in the same order |
| `study_id` | Short label (`study001`, …) |
| `log_file` | Where the study's log is written (`results/{study}/{study}.log`) |

### Tracking progress

Everything that tracks fine-mapping status lives in `susie_fm_status/`; its `README.md` lists
each file and whether it's current. `susie_fm_status/build_fm_status_table.py` is the one script
that produces the status table. It looks at what every study in `config/1008_study_table.tsv`
has produced so far and writes `susie_fm_status/fm_pipeline_run_status.tsv` (plus a `.csv`
copy), with one row per study.

Both pipelines refresh the table automatically at the end of every real run, whether the run
succeeds or fails (`onsuccess`/`onerror` in `pipeline_munge/Snakefile` and
`pipeline_finemap/Snakefile`). Dry runs and `--no-hooks` skip this, and only the driver
refreshes, never a job pod. To refresh by hand, run it from `FM_2026/`:

```bash
python susie_fm_status/build_fm_status_table.py
```

**Every fine-mapping window gets its own columns.** The script finds every
`results/{study}/checking/{method}/{N}kb_window/` folder on disk, plus the window of a finemap
run in progress. For each, it adds a `{method}_{N}kb_status` and `{method}_{N}kb_detail` column
pair, for example `susieR_300kb_status`. A new window size adds its columns the next time the
table is refreshed, with nothing to edit. Per-window statuses:

| Window status | Meaning |
|---|---|
| `finalized` | Finished at this window; the detail gives how many regions were fine-mapped and excluded |
| `running` | The finemap run in progress is on this window and includes this study |
| `incomplete` | Some outputs at this window, but not finalized |
| `not run` | Fine-mappable, but never run at this window |
| (blank) | The study never reached fine-mapping |

`fm_status` is the headline: the furthest stage reached in **any** window. `windows_finalized`
lists the windows where the study is finalized (e.g. `susieR_300kb;susieR_500kb`). Headline
statuses, from furthest along to least:

| Status | Meaning |
|---|---|
| `finalized` | Fine-mapping finished in at least one window (see the per-window detail columns) |
| `fine-mapping running` / `fine-mapping incomplete` | Fine-mapping is running, or started but didn't finish |
| `ready to fine-map` | Munged and approved for fine-mapping |
| `not fine-mappable` | Munged, but too few SNPs (decision `none`) |
| `munge failed` / `munge did not finish` | The cleaning step errored or stopped partway |
| `failed coverage check` | The file is missing something fine-mapping needs, or has too few SNPs |
| `no sumstats file` | No usable summary statistics file was found on the GWAS Catalog |
| `munge pipeline running`, `downloaded`, `download incomplete`, `started`, `not started` | Earlier stages |

### Fixing problem studies

Some studies fail only because the GWAS Catalog's copy of the file is broken (for example, the
cleaned-up "harmonised" file has empty effect-size columns). You can fix these by hand and let
the pipeline carry on. `pipeline_munge/scripts/prepare_GCST90010715.py` is a worked example: it
downloads the study's original file instead, renames columns to standard names, and writes it
where the pipeline expects it.

The important trick is the order of writing: write `results/{study}/raw/{study}_filename.txt`
**first** and `results/{study}/raw/{study}_sumstats.tsv.gz` **second**. Snakemake re-downloads
a file whenever its input is newer, so writing them the other way round would make the pipeline
download over your fixed file. If the source file is on GRCh37, put `GRCh37` in the name you
write to `_filename.txt` so the munge step knows to convert it.


## The munge pipeline, step by step

**Goal:** turn each study's GWAS Catalog file into one clean, standard file on genome build
GRCh38, or record clearly why that wasn't possible.

**Inputs:** `config/config_munge.yaml` (which studies) and `config/study_table.tsv` (where to
download each one).

**Outputs** (all in `results/{study}/`):

| File | What it is |
|---|---|
| `raw/{study}_trait.txt` | The trait name, for reference |
| `raw/{study}_filename.txt` | Which file on the GWAS Catalog was chosen |
| `raw/{study}_sumstats.tsv.gz` | The downloaded summary statistics, untouched |
| `raw/{study}_coverage_check.txt` | `PASS` if the file has what fine-mapping needs |
| `munged/{study}_munged.tsv.gz` | **The cleaned file** the fine-mapping pipeline reads |
| `munged/{study}_fm_decision.txt` | `none`, `regular_only` or `both` (see step e) |
| `munged/*.tsv.gz` (others) | SNPs the cleaning step removed, one file per reason (e.g. `snp_strand_ambiguous.tsv.gz`, `dup_row.tsv.gz`) |
| `{study}.log` | Everything that happened, step by step |

The code is `pipeline_munge/Snakefile` plus `pipeline_munge/scripts/munge_sumstats.R`. Each
step below is one Snakemake rule.

### a. Write trait file

Writes the study's trait name to a small text file so every results folder says what it is.

- Output: `results/{study}/raw/{study}_trait.txt`

### b. Find the summary statistics file

Each GWAS Catalog study folder is laid out a little differently, so this step works out which
file to download and records its name.

**Logic:**

1. If the folder has a `harmonised/` subfolder (the Catalog's cleaned-up version), take the
   file in it whose name contains `.h.`.
2. Otherwise take a file in the main folder ending in `.sumstats.tsv.gz`.
3. Otherwise take the first `.tsv` file whose name contains the study ID, and write a warning
   to the log (it may not be the right file).
4. If nothing is found, log "no sumstats file found" and stop this study.

- Output: `results/{study}/raw/{study}_filename.txt`

### c. Download

Downloads the chosen file to a temporary name, checks it isn't empty, and moves it into place,
compressing it first if the Catalog served it uncompressed. Snakemake only re-runs this if the
file is missing or older than `_filename.txt`.

- Input: `raw/{study}_filename.txt`
- Output: `raw/{study}_sumstats.tsv.gz`

### d. Coverage check

Checks the file has everything fine-mapping needs **before** spending hours cleaning it. The
header line and SNP count are always written to the log first.

**Logic:**

1. Look for these kinds of columns (each kind has a list of accepted names, e.g. `p_value`,
   `PVAL`, `P`, `p_val`). A column only counts if it actually contains values — a column that
   is entirely `NA` or `#NA` is treated as missing.
   - **Position:** chromosome and base-pair position
   - **P-value**
   - **Effect signal**, any one of: a Z-score; **or** an effect size (beta or odds ratio) plus
     its standard error; **or** an odds ratio plus its 95% confidence interval
2. **Fail** if position, p-value or effect signal is missing, or if the file has **fewer than
   500,000 SNPs**. The downloaded file is kept so you can look at it.
3. **Warn** (but pass) if it has 500,000–1,000,000 SNPs — fine-mapping will work but with
   lower coverage.
4. Otherwise write `PASS`.

- Input: `raw/{study}_sumstats.tsv.gz`
- Output: `raw/{study}_coverage_check.txt`

### e. Munge (clean and standardise)

Uses the R package MungeSumstats to put the file into one standard format on GRCh38, then
decides whether the study can be fine-mapped.

**Logic:**

1. Refuse to start unless the coverage check says `PASS`, and check the reference R packages
   are in place.
2. **Work out the genome build.** Look for a hint in the file name (`GRCh38`, `hg38`, `b38`,
   `.h.` mean GRCh38; `GRCh37`, `hg19`, `b37` mean GRCh37). If there's no hint, or the hints
   disagree, MungeSumstats works the build out from the data itself.
3. **Clean** with MungeSumstats: convert to GRCh38 if needed, standardise column names and
   alleles, and remove problem SNPs (strand-ambiguous A/T and C/G SNPs, duplicates, SNPs whose
   alleles don't match the reference genome, and so on). Removed SNPs are saved by reason in
   `munged/`.
4. **Decide** what the study is good for:
   - fewer than 500,000 SNPs left → `none` (not fine-mapped)
   - more than 100,000 people **and** more than 2,000,000 SNPs → `both` (suitable for regular
     and, in future, functional fine-mapping)
   - otherwise → `regular_only`

   Functional fine-mapping isn't built yet, so `both` and `regular_only` are currently
   fine-mapped the same way.

**Memory and time limits:** the first attempt gets 64 GB of memory. If it fails, Snakemake
retries with 96 GB, then 128 GB. Each attempt is stopped after `munge_timeout_hours` (4).
The log says which kind of failure happened:

| Exit code in the log | Meaning |
|---|---|
| 124 | Ran out of time |
| 137 | Probably ran out of memory |
| 2 | MungeSumstats itself errored (the R error is in the log) |
| 3 | The cleaned file came out empty or unreadable |
| 4 | The genome build couldn't be worked out |

- Inputs: `raw/{study}_sumstats.tsv.gz`, `raw/{study}_coverage_check.txt`,
  `raw/{study}_filename.txt`
- Outputs: `munged/{study}_munged.tsv.gz`, `munged/{study}_fm_decision.txt`
- Script: `pipeline_munge/scripts/munge_sumstats.R`


## The fine-mapping pipeline, step by step

**Goal:** fine-map every lead SNP of every study, check the results, and compare them with
Open Targets.

**Inputs:** the munge outputs (`munged/{study}_munged.tsv.gz` and `_fm_decision.txt`),
`config/config_finemap.yaml`, `config/study_table.tsv` (the lead SNPs and ancestry), and the
reference files from setup.

**Outputs:** see [Where everything ends up](#where-everything-ends-up) and
[Your final files, explained](#your-final-files-explained).

**How it works, in plain terms.** Each lead SNP defines a **locus** (region): the lead SNP
plus `locus_window_bp` bases either side. For each locus the pipeline:

1. cuts that window out of the study's GWAS results,
2. cuts the same window out of the 1000 Genomes reference for the study's ancestry,
3. checks whether the lead SNP's close LD partners are present in the GWAS,
4. builds an LD matrix (every SNP's correlation with every other SNP),
5. lines the GWAS results up with the LD matrix SNP by SNP,
6. runs SuSiE, which uses the GWAS signal plus LD to estimate each SNP's PIP and credible sets,
7. runs quality checks and compares with Open Targets.

Then it gathers every locus into study-level tables and a final report.

A note on paths: three settings are part of most output paths. Below, `{window}` is the window
folder (e.g. `350kb_window`, from `locus_window_bp`) and `{method}` is the fine-mapping method
(e.g. `susieR`). Each lead SNP gets its own folder, `results/{study}/{locus}/`, where `{locus}`
is the lead SNP's name. **Output paths below are inside `results/{study}/`** unless they start
with `results/`.

The code is `pipeline_finemap/Snakefile` plus the scripts in `pipeline_finemap/scripts/`.

### Stage 1: Get the lead SNPs ready

#### a. get_associations — check each lead SNP's position

Lead SNPs come from the study table, which mixes genome builds and sometimes has no rsID. This
step makes sure every lead SNP has a correct GRCh38 position, using the dbSNP database
(version 155) for both builds.

**Logic:**

1. If the munge decision is `none`, write an empty table and stop (the study isn't fine-mapped).
2. For each lead SNP **with an rsID**: look it up in dbSNP on GRCh38 and GRCh37.
   - Position matches GRCh38 → keep it.
   - Position matches GRCh37 → replace it with the GRCh38 position.
   - Matches neither → **drop** it ("does not match either GRCh37 or GRCh38 dbSNP position").
3. For each lead SNP **without an rsID** (named like `chr16:77507275`): check which build that
   position is a known SNP on.
   - Known only on GRCh38 → keep it.
   - Known only on GRCh37 → convert it to GRCh38 with the liftOver chain file; drop it if the
     conversion fails or is ambiguous.
   - Known on **both** builds, or neither → **drop** it, because we can't tell which build it
     is on ("genome build could not be confidently resolved").
4. Write the kept lead SNPs and log a summary ("N of M resolved to GRCh38").

Dropped lead SNPs are listed in the study log with their reason.

- Output: `results/{study}/associations/{study}_lead_snps.tsv` (columns SNP, CHR, BP)
- Script: `get_associations.R`

#### b. gather_loci — one job per lead SNP

Creates one empty marker file per lead SNP so Snakemake knows how many loci there are.

- Output: `results/{study}/loci_pending/`

### Stage 2: Build each locus

#### c. subset_gwas_locus — cut the GWAS window

Keeps only the GWAS rows on the lead SNP's chromosome within ±`locus_window_bp` of it.

- Output: `{locus}/{window}/{locus}_gwas_subset.tsv`
- Script: `subset_gwas_locus.R`

#### d. subset_1000g_and_plink — cut the reference window

Takes the same window from the 1000 Genomes panel, keeping only people of the study's
ancestry (EUR or EAS), and renames every reference SNP to `chr:position` so it can be matched
to the GWAS by position.

- Output: `{locus}/{window}/LD_block_check/{locus}.region.chrpos.bed/.bim/.fam`

#### e. compute_ld_r2 and ld_neighbor_check — are the lead SNP's close partners present?

If SNPs strongly correlated with the lead SNP are missing from the GWAS (often removed during
cleaning), fine-mapping can't consider them and may point at the wrong SNP. This check flags
that. It finds every reference SNP within `ld_window_kb` (300 kb) of the lead with
r² ≥ `ld_r2_threshold` (0.6), then sees which are missing from the GWAS window (missing indels
aren't counted).

| Flag | Meaning |
|---|---|
| `PASS` | No close partners missing |
| `MOD_LOW` | One missing partner, r² below 0.8 |
| `MOD_HIGH` | One missing partner, r² 0.8 or higher |
| `HIGH` | More than one missing partner |
| `NOT_IN_REFERENCE` | The lead SNP itself isn't in the reference panel — **locus excluded** |
| `CHECK_FAILED` | The check errored — **locus excluded** |

Only the last two stop a locus; the others are recorded and fine-mapping goes ahead.

- Outputs: `{locus}/{window}/LD_block_check/{locus}_1000g_0.6_ld_snps.ld`,
  `{locus}/{window}/LD_block_check/{locus}_ld_check.txt`
- Script: `ld_neighbor_check.R`

**Optional — rescue_missing_ld_neighbors:** if `rescue_missing_ld_neighbors: true`, missing
partners are looked up in the SNPs the munge step removed and added back when their alleles
match the reference. Off by default.

#### f. build_ld_flag_table and filter_eligible_loci

Collects every locus's flag into one table, then lists the loci that may continue (everything
except `NOT_IN_REFERENCE` and `CHECK_FAILED`).

- Outputs: `associations/{window}/{study}_ld_flags.tsv`, `{window}/loci_eligible/`

#### g. ld_matrix_and_freq — build the LD matrix

Keeps only the reference SNPs that are also in the GWAS window, have a minor allele frequency
of at least 1% and have only two alleles, then calculates the LD matrix and allele frequencies.

- Outputs: `{locus}/{window}/LD_matrix/{locus}.locus.matrix.ld` (the matrix),
  `.locus.matrix.snplist` (its SNP order), `.locus.frq`, `.locus.plink.bim`

#### h. align_gwas_bim — line the GWAS up with the LD matrix

SuSiE needs the GWAS rows in exactly the same order as the LD matrix, with alleles counted the
same way round.

**Logic:**

1. If the GWAS has no Z-score, make one: from beta ÷ standard error, or from the odds ratio
   and its confidence interval. If neither is possible, the locus fails.
2. Drop positions that have more than one allele in the reference.
3. Put the GWAS rows in the LD matrix's order. If the counts still don't match, the locus fails.
4. Where the GWAS alleles are the reverse of the reference, flip the sign of the effect and Z.
5. Final safety check: the row order must exactly match the LD matrix, otherwise the locus
   fails. (This guards against a past bug where rows were silently re-sorted.)

A failed locus gets `ALIGNMENT_STATUS: FAILED` and is excluded.

- Output: `{locus}/{window}/{locus}_gwas_subset_aligned.tsv`
- Script: `align_gwas_bim.R`

#### i. build_exclusions_table — which loci were excluded and why

One row per lead SNP with its LD flag, alignment status, and whether it is excluded:

| Reason | Meaning |
|---|---|
| `none` | Not excluded — fine-mapped |
| `lead SNP not in 1000G reference panel` | Flag `NOT_IN_REFERENCE` |
| `LD-neighbor check errored; ...` | Flag `CHECK_FAILED` |
| `GWAS/LD-matrix dimension mismatch after alignment` | Alignment failed (any cause) |
| `did not reach alignment step (check pipeline for gaps)` | Shouldn't happen — investigate |

- Output: `associations/{window}/{study}_exclusions.tsv`
- Script: `build_exclusions_table.R`

### Stage 3: Fine-map

#### j. run_susie

Runs SuSiE (`susie_rss`, from the R package susieR) on each eligible locus.

**Logic:**

1. Skip if alignment failed. Stop safely (status `FAILED`) if the GWAS and LD matrix don't
   have the same SNPs in the same order, or the LD matrix has missing values.
2. Calculate **lambda**, a score for how well the GWAS Z-scores agree with the reference LD
   (0 = perfect; larger = they disagree, often because the reference people don't match the
   study's people well).
3. Run SuSiE with up to `susie_max_causal` (10) signals. With `susie-inf` / `susie-ash` it also
   models many tiny background effects.
4. Record each SNP's PIP and credible set. SuSiE only reports a credible set if it reaches 95%
   coverage and its SNPs are reasonably correlated with each other (susieR's defaults).
5. Flag individual SNPs whose Z-score disagrees with what LD predicts (a "kriging" check): these
   are the **outlier** ("red") SNPs.
6. Give the locus a QC status: `QC_FAIL` if lambda > 0.3; `QC_WARN` if lambda > 0.1 or SuSiE
   didn't converge; otherwise `QC_PASS`.

- Outputs in `{locus}/{window}/{method}/`: `{locus}_susie_results.tsv`,
  `{locus}_susie_summary.txt`, `{locus}_susie_plot.pdf`, `{locus}_lambda_diagnostic.pdf`
- Script: `run_susie.R`

### Stage 4: Compare with Open Targets

Open Targets publishes fine-mapping for many GWAS Catalog studies. Comparing ours with theirs is
a sanity check.

#### k. build_opentargets_cache and subset_opentargets_native

Pulls every Open Targets credible set for the studies in `study_table.tsv` out of the large Open
Targets table once (a shared cache), then keeps just this study's rows.

- Outputs: `results/checking/opentargets_gwas_native_cache.tsv`,
  `checking/{window}/{study}_opentargets_native.tsv`

#### l. build_opentargets_comparison

Matches our SNPs to theirs **by position** and compares PIPs and credible sets.

- Outputs in `checking/{method}/{window}/`: `{study}_pip_full.tsv`, `{study}_cs_alignment.tsv`,
  `{study}_pip_matched.pdf`, `{study}_pip_all.pdf`

#### m. build_locus_pip_plots and build_opentargets_alignment_summary

Per-locus agreement (R²) and a study-level summary of the comparison, plus one PIP-vs-PIP plot
per locus.

- Outputs: `checking/{method}/{window}/{study}_locus_r2_summary.txt`,
  `{study}_alignment_summary.txt`, and `{locus}/{window}/{method}/{locus}_opentargets_pip_alignment.pdf`

### Stage 5: Summarise

#### n. build_manhattan_plot

One plot per locus: every SNP's p-value along the window, coloured by its LD with the lead SNP.

- Output: `{locus}/{window}/{method}/{locus}_manhattan_plot.pdf` (once every locus has one,
  `collect_manhattan_plots` writes the marker `{window}/manhattan_plots_done.txt`)

#### o. locus_summary

A plain-text summary of each locus that combines the LD check, SuSiE and Open Targets results.

- Output: `{locus}/{window}/{method}/{locus}_summary.txt`
- Script: `locus_summary.R`

#### p. build_snps_table

Stacks every fine-mapped locus's SNPs into one study-wide table.

- Output: `checking/{method}/{window}/{study}_snps_table.tsv`

#### q. build_final_report

The human-readable report: one entry per lead SNP saying whether it was fine-mapped and, if not,
why.

- Output: `checking/{method}/{window}/{study}_final_report.txt`

#### r. finalize_study

Writes a summary banner (total, excluded and fine-mapped loci) and a list of output files into
the study log, then writes the "done" marker.

- Output: `checking/{method}/{window}/study_finalized.txt`

**Not currently used:** `build_reference_comparison` (compares each lead SNP's PIP with the
lead-variant library's earlier fine-mapping) is still in the Snakefile, but nothing asks for its
output, so it never runs. `aggregate_ld_matrices` and `aggregate_susie_results` are debugging
helpers.


## File tree map

Simplified layout for one study, using `300kb_window` and `susieR` as the window and method.
Files marked ★ are the ones most people need.

```
results/
├── checking/
│   └── opentargets_gwas_native_cache.tsv     Open Targets credible sets for all our studies (shared)
└── {study}/                                  e.g. results/GCST90476277/
    ├── {study}.log                           ★ every step's log messages, in order
    ├── raw/                                  ── munge pipeline ──
    │   ├── {study}_trait.txt                 trait name
    │   ├── {study}_filename.txt              which Catalog file was used
    │   ├── {study}_sumstats.tsv.gz           the original download
    │   └── {study}_coverage_check.txt        PASS, or missing if the check failed
    ├── munged/
    │   ├── {study}_munged.tsv.gz             cleaned GRCh38 summary statistics
    │   ├── {study}_fm_decision.txt           none / regular_only / both
    │   └── *.tsv.gz                          SNPs removed during cleaning, by reason
    ├── associations/                         ── fine-mapping pipeline ──
    │   ├── {study}_lead_snps.tsv             lead SNPs with checked GRCh38 positions
    │   └── 350kb_window/
    │       ├── {study}_ld_flags.tsv          missing-neighbour flag per lead SNP
    │       └── {study}_exclusions.tsv        ★ which loci were excluded and why
    ├── loci_pending/                         one marker per lead SNP (internal)
    ├── 300kb_window/
    │   ├── loci_eligible/                    one marker per locus that continues (internal)
    │   └── manhattan_plots_done.txt          internal marker
    ├── {locus}/                              one folder per lead SNP, e.g. rs10486567/
    │   └── 300kb_window/
    │       ├── {locus}_gwas_subset.tsv              GWAS rows in the window
    │       ├── {locus}_gwas_subset_aligned.tsv      same rows, lined up with the LD matrix
    │       ├── LD_block_check/
    │       │   ├── {locus}.region.chrpos.*          1000 Genomes window (PLINK files)
    │       │   └── {locus}_ld_check.txt             missing-neighbour check
    │       ├── LD_matrix/
    │       │   └── {locus}.locus.matrix.ld, ...     the LD matrix and its SNP list
    │       └── susieR/
    │           ├── {locus}_susie_results.tsv        ★ every SNP's PIP and credible set
    │           ├── {locus}_susie_summary.txt        SuSiE run and QC status
    │           ├── {locus}_summary.txt              ★ plain-text summary of the locus
    │           ├── {locus}_susie_plot.pdf           PIP plot
    │           ├── {locus}_lambda_diagnostic.pdf    QC plot (outlier SNPs)
    │           ├── {locus}_manhattan_plot.pdf       p-values coloured by LD
    │           └── {locus}_opentargets_pip_alignment.pdf   our PIPs vs Open Targets
    └── checking/
        ├── 300kb_window/
        │   └── {study}_opentargets_native.tsv       Open Targets credible sets for this study
        └── susieR/300kb_window/
            ├── study_finalized.txt                  "DONE" when the study is finished
            ├── {study}_final_report.txt             ★ the report to read first
            ├── {study}_snps_table.tsv               ★ all fine-mapped SNPs, study-wide
            ├── {study}_pip_full.tsv                 our PIPs and Open Targets', SNP by SNP
            ├── {study}_cs_alignment.tsv             credible-set overlap with Open Targets
            ├── {study}_alignment_summary.txt        ★ Open Targets comparison summary
            ├── {study}_locus_r2_summary.txt         agreement with Open Targets per locus
            └── {study}_pip_matched.pdf, _pip_all.pdf   PIP comparison plots

logs/
├── munge_<timestamp>.log, munge_latest.log          munge driver logs
└── finemap_<timestamp>.log, finemap_latest.log      fine-mapping driver logs
```

Running again with a different `locus_window_bp` or `finemapping_method` writes to new
`{window}` / `{method}` folders and leaves earlier results alone. The study log
(`{study}.log`) and the files above the window folders are shared between runs. The status
table picks up the new folders on its own, adding a column pair for the new method/window (see
[Tracking progress](#tracking-progress)).


## Final files, explained

Start with the **final report**, then use the **SNPs table** for analysis. The QC and Open
Targets files explain *why* a result looks the way it does.

### Results files

#### `{study}_snps_table.tsv` — every fine-mapped SNP in the study

One row per SNP for every fine-mapped locus (excluded loci aren't in it). This is the main
table for downstream analysis.

| Column | Meaning |
|---|---|
| `TRAIT`, `STUDY` | Trait name and GCST ID |
| `LOCUS` | The lead SNP this region was built around |
| `CS_BPS` | The window's edges: "lead position + window, lead position − window" |
| `SNP`, `CHR`, `BP` | SNP ID and GRCh38 position |
| `REF`, `ALT` | The SNP's two alleles, lined up with the reference panel |
| `MAF` | Minor allele frequency in the 1000 Genomes reference people |
| `PVAL`, `Z` | The GWAS p-value and Z-score |
| `PIP` | **Posterior inclusion probability**: the chance this SNP is a causal variant (0–1) |
| `CS` | Which credible set the SNP belongs to (`L1`, `L2`, …), or `NA` if none |
| `CS_SIZE` | How many SNPs are in that credible set (smaller = more precise) |
| `LAMBDA` | The locus's lambda QC score (same for every SNP in a locus) |

#### `{locus}_susie_results.tsv` — one locus, every SNP

The same columns as above (with `A1`/`A2` for the alleles), plus four QC columns:

| Column | Meaning |
|---|---|
| `logLR`, `z_std_diff`, `KRIGING_Z` | Numbers from the outlier check: how far the SNP's Z-score is from what its LD with other SNPs predicts |
| `QC_OUTLIER` | `TRUE` if the SNP's Z-score clearly disagrees with the LD — a "red" SNP. A few are normal; many suggest a mismatch between the study and the reference people. |

### Summary files

#### `{study}_final_report.txt` — read this first

1. **Header:** study ID, date, number of lead associations.
2. **Munging summary:** SNP counts before and after cleaning.
3. **Top 20 SNPs by PIP** across the study (QC-passing loci only).
4. **Loci with QC_WARN or QC_FAIL**, with the reason (e.g. "lambda 0.134 exceeds warn
   threshold 0.100", "did not converge").
5. **One entry per lead SNP:** `FINE-MAPPED: YES` followed by that locus's summary (below), or
   `FINE-MAPPED: NO` with the reason (excluded, SuSiE didn't run, and so on).

#### `{locus}_summary.txt` — one locus in plain text

| Field | Meaning |
|---|---|
| `STATUS` | `OK`, `NOT_RUN` (SuSiE didn't run) or `ERROR` |
| `SUSIE_CONVERGED_UNDER_50_ITER`, `niter` | Whether SuSiE settled on an answer, and how many rounds it took |
| `MISSING_LD_NEIGHBORS_R2_GE_0.8`, `MISSING_LD_NEIGHBORS_R2_0.6_TO_0.8` | How many of the lead SNP's close LD partners are missing from the GWAS. Missing partners with r² ≥ 0.8 matter most: the true causal SNP could be one of them. |
| `SNPS_RESCUED` | Partners added back by the optional rescue step |
| `N_CREDIBLE_SETS`, `CS_SIZES` | How many credible sets SuSiE found and their sizes (e.g. `L1=2`) |
| `N_SNPS_PIP_APPROX_1`, `SNPS_PIP_APPROX_1` | SNPs with PIP ≈ 1. Occasionally real, but several at one locus often signal a problem. |
| `LAMBDA` | Agreement between GWAS and reference LD (near 0 is good; high values mean the reference panel doesn't match the study well) |
| `PCT_SNPS_FAILING_QC`, `N_RED_SNPS`, `RED_SNPS` | Outlier ("red") SNPs from the outlier check |
| `R2_VS_OPENTARGETS` | How well our PIPs agree with Open Targets' at this locus (needs ≥ 3 shared SNPs) |
| Four `..._HIGH_PIP_...` counts | SNPs one side calls high-PIP (> 0.9) but the other calls low or doesn't have at all |
| Final table | SNP, credible set and PIP for every credible-set member |

#### `study_finalized.txt`

Just says `DONE`. Its existence means every step for the study finished.

### QC files

| File | What to look at |
|---|---|
| `{study}_exclusions.tsv` | One row per lead SNP: `LD_FLAG`, `ALIGN_STATUS`, `EXCLUDE` (TRUE/FALSE) and `EXCLUSION_REASON`. Lead SNPs dropped by get_associations aren't in this table — their reasons are in the study log. |
| `{study}_ld_flags.tsv` | Just the missing-neighbour flag per lead SNP (see stage 2e). |
| `{locus}_ld_check.txt` | The flag, the number of missing partners, and a list of them with their r². |
| `{locus}_susie_summary.txt` | `SUSIE_STATUS` (OK / FAILED / SKIPPED), `QC_STATUS` (QC_PASS / QC_WARN / QC_FAIL / NOT_RUN), lambda, convergence, number of credible sets, largest credible set, number and rate of outlier SNPs, and any warnings. |
| `{locus}_lambda_diagnostic.pdf` | The outlier check as a plot: each SNP's observed vs. expected Z-score. Red points are outliers. |
| `{locus}_susie_plot.pdf` | PIP for every SNP in the window, with credible sets highlighted. |
| `{locus}_manhattan_plot.pdf` | GWAS p-values along the window, coloured by LD with the lead SNP. Useful to see whether the signal is one clean peak. |

**"Fine-mapped, but no credible set"** means SuSiE ran fine but couldn't narrow the signal
down to a small, tightly correlated group of SNPs with 95% confidence. Common causes: a weak
signal, a very long stretch of LD, or GWAS results that don't fit the reference LD well (check
lambda and the red SNPs).

### Open Targets comparison files

Open Targets fine-maps many of the same studies with its own methods (mostly SuSiE-inf, and
PICS for older studies). Close agreement is reassuring; big differences are worth a look.

| File | What it contains |
|---|---|
| `{study}_opentargets_native.tsv` | Open Targets' credible sets for this study, as downloaded. |
| `{study}_pip_full.tsv` | Every SNP from either side, matched by position: `my_pip`, `my_cs`, `ot_pip`, `ot_cs_id`, `ot_finemapping_method` and `match_status` — `both`, `mine_only` (only we have it, usually a low-PIP SNP Open Targets didn't report) or `OT_only`. |
| `{study}_cs_alignment.tsv` | Each of our credible sets paired with each overlapping Open Targets set: their sizes, shared SNPs and **Jaccard** overlap (shared ÷ combined; 1 = identical). Sets bigger than 50 SNPs are skipped. |
| `{study}_alignment_summary.txt` | Study-level summary: SNP counts by match status; agreement statistics per Open Targets method (n, R, R², mean absolute error, RMSE); the 10 biggest PIP disagreements; how many loci agree with R² > 0.8; mean and median credible-set Jaccard. |
| `{study}_locus_r2_summary.txt` | Agreement (R²) for each locus with at least 3 shared SNPs. |
| `{study}_pip_matched.pdf`, `{study}_pip_all.pdf`, `{locus}_opentargets_pip_alignment.pdf` | Our PIPs plotted against Open Targets', study-wide and per locus. |

### Logs

- `results/{study}/{study}.log` — every step of both pipelines adds lines here, so it is the
  place to find out why a study or lead SNP stopped. Search it for `WARNING`, `ERROR` or `FAIL`.
- `logs/munge_latest.log`, `logs/finemap_latest.log` — the Snakemake driver output for the
  latest run (which jobs ran, failed or are waiting).


## Known issues and ideas

- **Setting munge memory:** passing `munge_mem_mb_schedule` with `--config` (as in the
  `run_munge.sh` examples above) passes the numbers as text and the run fails with a
  `ResourceConversionError`. Put it in a small YAML file instead and pass both config files:
  `./run_munge.sh start --configfile config/config_munge.yaml my_override.yaml`.
- **Exclusion reasons are coarse:** every alignment failure is reported as "GWAS/LD-matrix
  dimension mismatch", even when the cause was something else (e.g. no usable effect size). The
  study log has the specific reason.
- **`locus_output_root`:** if you set it, also give Snakemake explicit target files — the
  pipeline's default target list still points at `results/`.
- **Hard-coded "0.6":** the LD file name `{locus}_1000g_0.6_ld_snps.ld` always says 0.6, and
  the summary buckets are fixed at 0.8 and 0.6, whatever `ld_r2_threshold` is set to.
- **Unused settings:** `n_mid_low` and `use_functional_annotations` in `config_munge.yaml`
  don't do anything yet. Functional fine-mapping (for `both` studies) is planned.
- **Status table's headline mixes windows:** `fm_status` in
  `susie_fm_status/fm_pipeline_run_status.tsv` is the furthest stage in any window, so
  "finalized" there doesn't say which window. Use the `{method}_{N}kb_status` columns or
  `windows_finalized` for that.
- **No chromosome X in the reference panel:** every panel on disk (`1000g_hg38_plink`,
  `1000g_hg38_30x_plink`, `1000g_hg38_30x_dedup_plink`) has only chr1–22, for both EUR and EAS.
  For a lead SNP on chrX, `subset_1000g_and_plink` fails with "Failed to open
  `..._chrX.bed`". `build_ld_flag_table` needs every locus of a study, so that one locus stops the
  **whole study** from finalizing; `--keep-going` lets other studies carry on. This is why five
  studies never finalize in any window. In `susie_fm_status/fm_pipeline_run_status.tsv` each
  window column shows them as `incomplete` (`running` while a run is on). The five are
  GCST007800 (rs5953283, rs850637), GCST90435202 (rs5912761), GCST90446793 (rs175737),
  GCST90018627 (rs2982584) and GCST90018690 (rs144238538). It happened in the 500 kb run (Oct 19) and again in
  the 300 kb run (Oct 24). Not fixed yet. One option is to drop chrX lead SNPs in
  `get_associations`, which finalizes these studies without those loci. The other is to build
  chrX panel files from the 30x VCFs, which keeps the loci. chrX LD needs care because males
  carry one copy.
- **Ideas:** run the munge pipeline from inside the fine-mapping pipeline for a fully hands-off
  run; switch on the lead-variant
  library comparison (`build_reference_comparison`).
