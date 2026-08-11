# 2026 FM Pipeline Documentation

Sections 1-4 below document how the pipeline works, stage by stage. The two sections
immediately following are the practical "how do I run this" instructions: set up the R
reference library once, then launch with `run_munge.sh`.


## Setup: R reference library (renv)

`munge_sumstats.R` needs four SNPlocs/BSgenome reference data packages. They come to
roughly 13 GB, which is too large to ship inside the container image, so they live in the
project's renv library and are mounted into the job pods instead. The container carries
MungeSumstats; the renv library carries the reference data. Both are needed.

The library **must be built against R 4.4 / Bioconductor 3.20**, matching the image
(`bioconductor/bioconductor_docker:RELEASE_3_20`, R 4.4.2). R refuses to load packages
built by a newer R minor version, so a mismatch is a guaranteed failure — and it surfaces
hours in, at the end of a long munge, not at launch.

Run once per checkout, from the project root:

```bash
# 1. Confirm R is 4.4.x. If it is not, switch R before going further.
R --version                      # expect: R version 4.4.2

# 2. Point renv at the shared package cache. Without this, restore compiles ~13 GB from
#    source; with it, packages are symlinked out of the cache in minutes. The workspace
#    image normally sets this already.
export RENV_PATHS_CACHE=/renv

# 3. Restore. .Rprofile bootstraps renv automatically, so no install step is needed.
R -e 'renv::restore(prompt = FALSE)'
```

Then verify the four data packages actually landed — this is exactly what `munge_gwas`
pre-flight checks before it starts:

```bash
ls -d renv/library/*/*/*/{SNPlocs.Hsapiens.dbSNP155.GRCh38,SNPlocs.Hsapiens.dbSNP155.GRCh37,BSgenome.Hsapiens.NCBI.GRCh38,BSgenome.Hsapiens.1000genomes.hs37d5}
```

And confirm the Snakefile discovers it, which prints the resolved absolute path:

```bash
./run_munge.sh dryrun | head -2
# Running 1 studies: ['GCST90132226']
# Reference package library: /home/<you>/FM_2026/renv/library/linux-ubuntu-noble/R-4.4/x86_64-pc-linux-gnu
```

Things that go wrong here:

- **"no reference package library found"** — `renv::restore()` has not been run, or it was
  run from somewhere other than the project root.
- **"multiple renv libraries found"** — more than one directory matches
  `renv/library/*/*/*`, usually left over from an R upgrade. Delete the stale one, or name
  the right one via `ref_lib:` in `config/config_munge.yaml`.
- **A warning that the library is `R-4.3` (or similar) but the image ships `R-4.4`** — the
  library was built under the wrong R. Rebuild it under R 4.4.
- **Dangling symlinks.** renv library entries point into the shared `/renv` cache, so job
  pods must mount both the project directory and `/renv`. The `coder` profile does both.

To use a standalone BiocManager library instead of renv, see the install command in the
Dockerfile header and set `ref_lib:` in `config/config_munge.yaml` (or pass
`--config ref_lib=/path`).

Note: keep credentials such as `GITHUB_PAT` in `~/.Renviron`, never in a project-level
`.Renviron` — the latter sits inside the repo and gets committed.


## Running: run_munge.sh

`run_munge.sh` at the project root launches the munge pipeline detached under `nohup` and
gives you a clean way to stop it. Runs are long — tens of minutes per GWAS download, up to
`munge_timeout_hours` per munge — so they need to survive the terminal closing.

Before launching: list the GCST IDs you want in `studies:` in `config/config_munge.yaml`,
and make sure each one has a row in `config/study_table.tsv`.

```bash
conda activate snakemake          # or: export SNAKEMAKE=/path/to/bin/snakemake

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

Environment overrides: `SNAKEMAKE` (binary path), `PROFILE` (default `coder`), `JOBS`
(default 4), `IMAGE` (default: the `container:` key in `config/config_munge.yaml`).

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


## 1: Build study table

Build fine-mapping study table from LeAnn's original MPRA variant table. Filters are regularly adjusted to determine which studies should be FM.

(07/31) Currently, they are as follows:
- Sum stats available: TRUE
- Ancestry: EUR or EAS (not EUR/EAS)
- Has FM: less than 50% of associations FM

Associations are grouped by study. Additional columns such as CHR:BP location, percent FM, and log file location are added.

Table is used to collect all GCST IDs eligible for FM. GCST IDs pasted into config files.

#### Relevant files:
	Script: pipeline_munge/scripts/build_study_table.py
	Input: config/07272026_library_leadVariantTable.tsv
	Output: config/study_table.tsv

Potential changes:
Collect studies that have only been FM by PICS


*The munging and fine-mapping pipelines have been separated into two separate snakemake pipelines. I chose to do this because 1) I see them as two separate processes and thought a munging pipeline could be useful in the future, even unrelated to FM, 2) it makes the debugging process much easier, and 3) many studies don't make it to FM post-munge.

Potential changes:
Determine if there is a way to call pipeline_munge within pipeline_fm for complete hands-off fm.



## 2: Adjust config files

Config files contain all the relevant input information for snakemake pipelines.

"studies:" will always contain a list of GCST IDs that require FM

-> config_munge.yaml - for pipeline_munge
use functional annotations: T/F
7/31: This is currently irrelevant. Hope to implement functional fine-mapping if studies pass n and SNP thresholds. Functional fine-mapping is only accurate when GWASs possess a certain level of statistical power. I am using these n and SNP thresholds as a proxy power assessment. The thresholds are currently n >= 100,000 and snp_count >= 2,000,000. Studies will always also undergo regular fine-mapping. *** need literature backing

finemap thresholds:
n high: #
snp hig: #
n mid low: #
snp mid low: #
More on these in munge_sumstats.R

munge timeout hours: #


config_finemap.yaml - for pipeline_finemap


#### Relevant files:
	Input: config/config_munge.yaml
	       config/config_finemap.yaml

Potential changes:
Config files contain many filters. SNP and samples size filters could be adjused within config_munge. Many QC filters could be adjusted in config_finemap.



## 3: pipeline_munge/Snakefile

Snakefile that will automatically:
Download GWAS
Check its coverage/column availability
Munge (if applicable)
Make a decision on how to proceed (FM or not)

a. write trait file

Didn't organize by trait type this year, so recording study trait in file for organizational purposes.
	
Relevant files:
	output: results/{study}/raw/{study}_trait.txt


### b. resolve gwas filename

Determines and records GWAS sum stats file name as listed in GWAS catalog. Checks for harmonised sum stats first, then sum stats under common naming convention (.sumstats.tsv.gz) in main directory. If neither are found, collects the first tsv file that contains the GCST_ID and flags (no obvious sumstats file was found). Filename is written to file. If no filename is found, log and exit pipeline.
	
#### Logic:

	If (directory url contains folder /harmonized/):
		Look inside and collect filename that contains ".h."
	else:
		Check the main directory and collect the filename that contains ".sumstats.tsv.gz"

	If (no filename was found):
		Collect first filename that contains .tsv and the GCST ID.

	If (no filename found):
		Log that no sum stats file was found and exit.

#### Relevant files:
	output: results/{study}/raw/{study}_filename.txt


### c. download gwas

Downloads the resolved GWAS summary statistics file. Skips download if the output already exists. Ensures the final file is gzip-compressed regardless of how the source served it (some GWAS Catalog files are already gzipped, others are plain text with a misleading extension).

Logic:

	If (output file already exists and is non-empty):
    		Skip download
	Else:
   		 Build full download URL from directory + resolved filename
    		Download to a temporary file (never directly to the final output path)
    	If (download failed or file is empty):
        	Log failure and exit
   	If (downloaded file is already gzip-compressed):
        	Move to final output path
   	Else:
        	Compress, then move to final output path

#### Relevant files:
	input:  results/{study}/raw/{study}_filename.txt
	output: results/{study}/raw/{study}_sumstats.tsv.gz

### d. check gwas coverage

Checks the downloaded sum stats for the minimum columns and SNP density required to proceed with fine-mapping. Separately flags (and blocks) studies whose SNP count is large enough to warrant reviewing session memory before munging.

#### Logic:

	Decompress sum stats to a temp file
	Record header and SNP count in log (unconditionally, before any pass/fail decision)

	Search header for column name matches across five categories:
    		position, effect size, standard error, Z-score, p-value

	Determine signal adequacy:
    		Z-score column alone is sufficient
    		OR effect size + standard error together are sufficient

	If (position column missing) OR (signal inadequate) OR (SNP count < 1,000,000):
    		Log failure reason
    		Retain sum stats file for manual review
   	 	Exit

	If (SNP count > 13,000,000):
    		Log memory flag: "Study has over 10,000,000 SNPs. Adjust session memory requirements."
    		Exit

    	(Study is blocked here. To resume: manually verify/adjust session memory,
    	then manually write "PASS" to this study's coverage_check.txt to signal
     	the check has been reviewed, and rerun the pipeline — munge_gwas will
     	then proceed normally.)

	Else:
    		Write "PASS"

#### Relevant files:
	input:  results/{study}/raw/{study}_sumstats.tsv.gz
	output: results/{study}/raw/{study}_coverage_check.txt

Potential changes:
Adjust large SNP count value -> fooling around with what session can handle
work with memory level inputs - unsure on how this works with shared workspaces


### e. munge gwas

Munges the downloaded sum stats into standardized GRCh38 format via MungeSumstats, determines the genome build beforehand (from filename hints, falling back to full inference), and classifies the study into a fine-mapping decision tier (`none` / `regular_only` / `both`) based on final SNP count and sample size. The rule itself is a shell wrapper around `munge\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_sumstats.R`, separating shell-level failures (timeout, out-of-memory) from R-level failures for clearer logging.

#### Logic (shell wrapper):
	
	Re-verify coverage_check file says PASS (defensive check; refuses to munge otherwise)

	Run munge_sumstats.R, wrapped with a configurable timeout

	If (exit code 124): log as timeout (OOT), exit

	If (exit code 137): log as likely out-of-memory (OOM), exit

	If (any other non-zero exit code): log generic failure, point to R-level log detail, exit


#### Logic (munge_sumstats.R):
	
	Determine genome build:
    		Check resolved filename for explicit build tokens
        		(".h.", "GRCh38", "GRCh37", "hg38", "hg19", "b38", "b37")
 
		If (conflicting signals found, e.g. both 37 and 38 tokens present):
        		Fall through to full inference

    		If (no reliable hint found in filename):
        		Run full genome build inference (get_genome_builds)

        	If (inference fails or returns nothing): log and exit

	Munge sum stats via MungeSumstats::format_sumstats():
    		Convert to GRCh38 (liftover applied only if source build differs)
   		Retain indels; do not filter to biallelic-only sites;
        	correct allele frequency orientation for multi-allelic sites
    		Filter out strand-ambiguous SNPs (A/T, C/G)

	If (munging errors): log actual R error message, exit

	Save raw munge result to disk immediately, before any further processing

	Extract the sum stats table from the raw result; verify it is non-empty
 
	If (table missing, malformed, or zero rows): log and exit
	
	Write final munged file
	
	Append MungeSumstats' own internal log messages into the study log

	Classify fine-mapping decision:
    		If (final SNP count < snp_mid_low):
        		decision = none — insufficient coverage to fine-map
    		Else if (sample size > n_high) AND (SNP count > snp_high):
        		decision = both — eligible for regular and functional fine-mapping
    		Else:
        		decision = regular_only

	Write decision to file; log full reasoning

#### Relevant files:
	input:  results/{study}/raw/{study}_sumstats.tsv.gz,
        	results/{study}/raw/{study}_coverage_check.txt,
        	results/{study}/raw/{study}_filename.txt
	output: results/{study}/munged/{study}_munged.tsv.gz,
        	results/{study}/munged/{study}_fm_decision.txt
	script: pipeline_munge/scripts/munge_sumstats.R



## 4: pipeline_finemap/Snakefile

Snakefile that will automatically:
Build the associations table from munged sumstats
Discover loci and subset the GWAS to each
Check LD-neighbor completeness against the 1000G reference
Build LD matrices and align GWAS to them
Run SuSiE fine-mapping and score confidence per locus
Compare against prior fine-mapping where available
Produce final per-study summary tables


### a. get associations

Creates the associations table listing every lead SNP still present in the munged sumstats, matched by rsID first, falling back to CHR:BP position matching when rsID fails (e.g. dbSNP version mismatches between the catalog table and munged data). Refuses to proceed for studies flagged "none" by the munge decision

#### Logic (get_associations.R):
	
	Read fm_decision file

	If (decision == "none"):
    		Write empty associations table (SNP, CHR, BP headers only)
    		Exit cleanly (status 0)

	For each lead SNP in study_table's lead_snps list:
   		Try exact match on SNP == rsID in munged sumstats
 
		If no match AND a valid position exists:
        		Fall back to match on CHR == chr AND BP == bp
        		Log if this fallback match succeeded
   		If still no match:
        		Record as unmatched, log it

	Write matched (SNP, CHR, BP) rows to associations table
	Log count of unmatched SNPs

#### Relevant files:
	input:  results/{study}/munged/{study}_munged.tsv.gz,
        	results/{study}/munged/{study}_fm_decision.txt
	output: results/{study}/associations/{study}_lead_snps.tsv
	script: pipeline_finemap/scripts/get_associations.R


### b. gather loci (checkpoint) + subset gwas locus

The number of loci isn't known until the associations table exists, so a checkpoint creates one marker file per lead SNP; a helper function (get_loci) reads those markers back to give Snakemake the real locus list. Each discovered locus then gets its own GWAS subset — the munged sumstats windowed to ±locus_window_bp around that lead SNP's position.

#### Logic:

	[checkpoint] For each row in associations table:
   		Create empty marker file named after that SNP

	[subset_gwas_locus.R] For this locus:
    		Look up lead SNP's CHR/BP in associations table

    		If not found: log error, exit (real failure — should never happen)

    		Filter munged sumstats to same CHR, BP within ±window_bp

    		Log SNP count retained; warn if subset is empty
    		Write subset

#### Relevant files:
	[checkpoint gather_loci]
		input:  results/{study}/associations/{study}_lead_snps.tsv
		output: directory(results/{study}/loci_pending)

	[subset_gwas_locus]
		input:  results/{study}/munged/{study}_munged.tsv.gz,
        		results/{study}/associations/{study}_lead_snps.tsv
		output: results/{study}/{locus}/{locus}_gwas_subset.tsv
		script: pipeline_finemap/scripts/subset_gwas_locus.R


### c. subset 1000g and plink

Subsets the hg38 1000 Genomes PLINK reference to a window around each locus's lead SNP, restricted to the study's ancestry, and rewrites the resulting .bim's SNP IDs to chr:pos format to avoid rsID mismatches between the GWAS and reference panel.

#### Logic (shell):
	Check ancestry sample-ID file exists and is non-empty
	
	If missing: log error, exit (no reference panel for this ancestry)

	Look up lead SNP's CHR/BP in associations table
	Compute window: [BP - window_bp, BP + window_bp]

	Extract 1000G PLINK region for that window, restricted to ancestry samples
	Rewrite .bim SNP ID column to chr:pos format

#### Relevant files:
	input:  results/{study}/associations/{study}_lead_snps.tsv
	output: results/{study}/{locus}/LD_block_check/{locus}.region.chrpos.bed,
        	results/{study}/{locus}/LD_block_check/{locus}.region.chrpos.bim,
        	results/{study}/{locus}/LD_block_check/{locus}.region.chrpos.fam

Potential changes:
Currently EUR-only reference panel available; EAS (and other ancestries) will fail this rule's guard check until an hg38 reference panel is obtained for them.


### d. compute ld r2

Computes LD (r²) between the lead SNP and every neighboring variant within the configured window, using the ancestry-matched 1000G subset. Handles the case where the lead SNP itself isn't genotyped in the reference panel.

#### Logic (shell):

	Look up lead SNP's CHR/BP
	If lead SNP position not found in region's .bim file:
    		Log warning
    		Write "LEAD_SNP_NOT_IN_REFERENCE" sentinel to output
    		Exit cleanly (status 0)

	Else:
    		Run PLINK --r2 anchored on lead SNP, filtered to r2_threshold and window_kb

#### Relevant files:
	input:  results/{study}/{locus}/LD_block_check/{locus}.region.chrpos.bed/.bim/.fam,
        	results/{study}/associations/{study}_lead_snps.tsv
	output: results/{study}/{locus}/LD_block_check/{locus}_1000g_0.6_ld_snps.ld

Potential changes:
Output filename hardcodes "0.6" regardless of the actual configured threshold — cosmetic only, doesn't affect correctness, but misleading if threshold is changed and files are inspected by eye.


### e. ld neighbor check

Classifies each locus by how well its GWAS data covers the lead SNP's high-LD neighbors. Rather than dropping a locus outright, flags it into a category and lets fine-mapping proceed regardless. Indels are exempted from counting as missing (they were deliberately excluded from GWAS elsewhere, so their absence isn't a real problem).

#### Logic (ld_neighbor_check.R):

	If LD file is the NOT_IN_REFERENCE sentinel:
    		Flag NOT_IN_REFERENCE, exit
	If LD file is empty (no neighbors found at all):
    		Flag PASS, exit

	For each high-LD neighbor PLINK found:
    		If present in GWAS subset: skip, no problem
    		If missing and is an indel: skip, not a real problem
   		Else: record as a genuinely missing neighbor, with its R2

	Classify:
    		0 missing -> PASS
    		1 missing, R2 < 0.8 -> MOD_LOW
    		1 missing, R2 >= 0.8 -> MOD_HIGH
    		>1 missing (any R2) -> HIGH

	Write flag, count, and full list of missing neighbors with R2 values

#### Relevant files:
	input:  results/{study}/{locus}/LD_block_check/{locus}_1000g_0.6_ld_snps.ld,
        	results/{study}/{locus}/{locus}_gwas_subset.tsv,
       		results/{study}/{locus}/LD_block_check/{locus}.region.chrpos.bim
	output: results/{study}/{locus}/LD_block_check/{locus}_ld_check.txt
	script: pipeline_finemap/scripts/ld_neighbor_check.R


### f. build ld flag table

Collapses every locus's individual LD flag into one per-study table, feeding directly into the eligibility checkpoint

#### Logic (build_ld_flag_table.R):

	For each locus's LD check file:
    		Extract FLAG value
	Merge onto full associations list (every original lead SNP present, NA where missing)
	
	Log summary counts per category
	Write combined table

#### Relevant files:
	input:  results/{study}/associations/{study}_lead_snps.tsv,
        	results/{study}/{locus}/LD_block_check/{locus}_ld_check.txt (all discovered loci)
	output: results/{study}/associations/{study}_ld_flags.tsv
	script: pipeline_finemap/scripts/build_ld_flag_table.R


### g. filter eligible loci (checkpoint)

Filters the LD flag table down to every locus except NOT_IN_REFERENCE, creating one marker file per eligible locus. A helper function (get_eligible_loci) reads these back for every downstream rule.

#### Logic (shell):

	For each row in ld_flags table where LD_FLAG != "NOT_IN_REFERENCE":
    		Create empty marker file named after that SNP

#### Relevant files:
	input:  results/{study}/associations/{study}_ld_flags.tsv
	output: directory(results/{study}/loci_eligible)


### h. ld matrix and freq

Takes the wider 1000G window subset and narrows it to the GWAS's exact positions, applying MAF (>0.01) and biallelic filters. Builds the final LD matrix, the SNP list matching its row/column order, and per-SNP allele frequencies.

#### Logic (shell):

	Build extract list from GWAS subset's CHR:BP positions

	Filter reference panel to those positions, MAF > 0.01, biallelic only

	Build square LD matrix + SNP list from filtered panel

	Build allele frequency file from same filtered panel

#### Relevant files:
	input:  results/{study}/{locus}/LD_block_check/{locus}.region.chrpos.bed/.bim/.fam,
        	results/{study}/{locus}/{locus}_gwas_subset.tsv
	output: results/{study}/{locus}/LD_matrix/{locus}.locus.matrix.ld,
        	results/{study}/{locus}/LD_matrix/{locus}.locus.matrix.snplist,
        	results/{study}/{locus}/LD_matrix/{locus}.locus.frq,
        	results/{study}/{locus}/LD_matrix/{locus}.locus.plink.bim


### i. align gwas bim

Reconciles the GWAS subset with the final LD matrix: excludes multi-allelic reference positions, forces row order to exactly match the LD matrix's row/column order, resolves count mismatches with a fallback, and corrects allele orientation (flipping BETA sign where the GWAS's alleles are swapped relative to the reference).

#### Logic (align_gwas_bim.R):

	Split GWAS subset into duplicate/non-duplicate positions
	Resolve duplicates by matching position AND alleles against reference
	Recombine

	Exclude any position that is multi-allelic in the reference panel
    		(prevents many-to-many merge duplication)

	Force row order to exactly match reference .bim's row order
    		(rbindlist does not preserve this on its own)

	If row count still mismatches LD matrix:
    		Attempt fallback: keep only GWAS rows whose position exists in reference
	If still mismatched after fallback:
    		Write ALIGNMENT_STATUS: FAILED, log loudly, exclude locus, exit cleanly

	Check allele orientation against reference:
    		If GWAS alleles exactly swapped relative to reference: flip BETA sign, swap A1/A2, log count
    		If alleles match neither orientation: flag as unresolved, log warning

	Write aligned GWAS subset

#### Relevant files:
	input:  results/{study}/{locus}/{locus}_gwas_subset.tsv,
        	results/{study}/{locus}/LD_matrix/{locus}.locus.plink.bim,
        	results/{study}/{locus}/LD_matrix/{locus}.locus.matrix.ld
	output: results/{study}/{locus}/{locus}_gwas_subset_aligned.tsv
	script: pipeline_finemap/scripts/align_gwas_bim.R

### j. build exclusions table

Combines the LD-neighbor flag and alignment outcome into one authoritative per-locus exclusion record, distinguishing genuine exclusions from an unexpected pipeline gap (a locus that should have reached alignment but didn't).

#### Logic (ld_exclusions_table.R):

	Extract LD flag per locus (all discovered loci)
	Extract alignment status per locus (only eligible loci — FAILED or OK)
	Merge both onto full associations list

	Missing LD_FLAG -> UNKNOWN (should not normally occur)
	Missing ALIGN_STATUS, LD_FLAG == NOT_IN_REFERENCE -> N/A (expected, never reached alignment)
	Missing ALIGN_STATUS, LD_FLAG != NOT_IN_REFERENCE -> NOT_REACHED (unexpected — investigate)

	EXCLUDE = (LD_FLAG == NOT_IN_REFERENCE) OR (ALIGN_STATUS == FAILED)
	Assign specific EXCLUSION_REASON per case

	Log exclusion summary
	If any locus is NOT_REACHED: log loud warning, print directly to console
	
	Write combined table

#### Relevant files:
	input:  results/{study}/associations/{study}_lead_snps.tsv,
       		results/{study}/{locus}/LD_block_check/{locus}_ld_check.txt (all discovered loci),
        	results/{study}/{locus}/{locus}_gwas_subset_aligned.tsv (eligible loci only)
	output: results/{study}/associations/{study}_exclusions.tsv
	script: pipeline_finemap/scripts/build_exclusions_table.R


### k. run susie

Runs SuSiE-RSS fine-mapping on the aligned GWAS/LD data. Computes the lambda consistency diagnostic, per-SNP outlier flags via kriging_rss, credible sets, and PIPs. Skips gracefully if alignment failed upstream.

#### Logic (run_susie.R):

	If aligned GWAS file shows ALIGNMENT_STATUS: FAILED:
    		Write SUSIE_STATUS: SKIPPED, exit cleanly

	Compute Z-scores from BETA/SE
	Compute lambda (estimate_s_rss) — consistency between Z-scores and LD matrix
	If GWAS row count != LD matrix row count: write FAILED, exit cleanly

	Run susie_rss (z, R, n, L)
	Extract PIP, credible set membership, credible set size per SNP
	
	Determine largest credible set size across all sets

	Run kriging_rss diagnostic:
    		Save PIP plot and lambda diagnostic plot
    		Flag SNPs where logLR > 2 AND |z_std_diff| > 2 as outliers
    		Compute count and rate of flagged SNPs

	Classify QC_STATUS:
   		lambda > lambda_fail_threshold -> QC_FAIL
    		lambda > lambda_warn_threshold OR did not converge -> QC_WARN
    		else -> QC_PASS

	Write results table (SNP, CHR, BP, A1, A2, MAF, PVAL, Z, PIP, CS, CS_SIZE, LAMBDA)
	Write summary (status, QC status, n_credible_sets, largest_cs_size, lambda, converged, qc_failed_snps, qc_failed_rate)

#### Relevant files:
	input:  results/{study}/{locus}/{locus}_gwas_subset_aligned.tsv,
        	results/{study}/{locus}/LD_matrix/{locus}.locus.matrix.ld/.snplist/.frq
	output: results/{study}/{locus}/{locus}_susie_results.tsv,
        	results/{study}/{locus}/{locus}_susie_summary.txt,
        	results/{study}/{locus}/{locus}_susie_plot.pdf,
        	results/{study}/{locus}/{locus}_lambda_diagnostic.pdf
	script: pipeline_finemap/scripts/run_susie.R

Potential changes:
z_ld_weight parameter exists in config but is currently commented out/unused — available to test on specific unstable loci (small reference panel relative to locus SNP density) if needed. kriging_rss column names (logLR, z_std_diff) assumed per susieR's documented diagnostic convention but not independently verified against installed version beyond one test locus.


### l. build reference comparison

Compares each locus's lead SNP PIP against a prior fine-mapping PIP recorded for this exact study in the lead variant table (parallel finemappingGWAS/PIP list columns). Informational only — not currently used in scoring.

#### Logic (build_reference_comparison.R):

	For each row in lead variant reference table:
    		Find this study's GCST ID within its finemappingGWAS list
    		Extract PIP at same list position -> reference_pip

	For each locus, read own susie_results.tsv, get lead SNP's own PIP -> my_pip
	Merge onto associations list
	Log count of matched loci and how many agree within 0.2 PIP

#### Relevant files:
	input:  results/{study}/associations/{study}_lead_snps.tsv,
        	results/{study}/{locus}/{locus}_susie_results.tsv (eligible loci)
	output: results/{study}/associations/{study}_reference_comparison.tsv
	script: pipeline_finemap/scripts/build_reference_comparison.R

Potential changes:
OpenTargets credible-set table integration (build-confirmed GRCh38, native GCST004131 match confirmed present) designed but not yet implemented — would follow this same pattern with full per-variant credible set comparison rather than lead-SNP-only. Incorporate into scoring.

### m. locus confidence

Scores each locus's fine-mapping reliability from 0-100 based on LD-neighbor completeness, per-SNP QC outlier rate, lambda, credible set count, and largest credible set size. Two hard disqualifiers bypass scoring entirely.

#### Logic (locus_confidence.R):

	If QC_STATUS == NOT_RUN: score = 0, category UNRELIABLE, exit
	If any SNP has PIP >= 1.00: score = 0, category UNRELIABLE, exit

	score = 100
	Subtract (count of missing R2>=0.8 neighbors) * ld_missing_high_weight
	Subtract (count of missing R2 0.6-0.8 neighbors) * ld_missing_low_weight
	Subtract (qc_failed_rate) * qc_failed_rate_weight
	Subtract accelerating penalty if lambda > lambda_ok_threshold (capped)
	Subtract discrete penalty by credible set count (1/2/3/4+)
	Subtract accelerating penalty if largest_cs_size > cs_size_ok_threshold 	(capped)

	Clamp score to [0, 100]
	Categorize: >=80 HIGH, >=50 MEDIUM, >=25 LOW, else UNRELIABLE

	Record prior fine-mapping comparison line (informational, not scored)
	Log full per-metric penalty breakdown
	Write score, category, breakdown, reasons

#### Relevant files:
	input:  results/{study}/{locus}/LD_block_check/{locus}_ld_check.txt,
        	results/{study}/{locus}/{locus}_susie_summary.txt,
        	results/{study}/{locus}/{locus}_susie_results.tsv,
        	results/{study}/associations/{study}_reference_comparison.tsv
	output: results/{study}/{locus}/{locus}_confidence.txt
	script: pipeline_finemap/scripts/locus_confidence.R

Potential changes:
All weight/threshold/scale constants are first-pass placeholders, not empirically validated — intended to be revisited using real per-metric penalty distributions once run across full batch. See scoring reference table for full formula documentation.


### n. finalize study

Writes a summary banner into the study log once every locus has been scored — total/excluded/fine-mapped counts and confidence category breakdown.

#### Logic (finalize_study.R):

	Read confidence category from every eligible locus's confidence file
	Tabulate counts per category (HIGH/MEDIUM/LOW/UNRELIABLE)
	Log banner with total loci, excluded count, fine-mapped count, category 	breakdown

#### Relevant files:
	input:  results/{study}/associations/{study}_exclusions.tsv,
        	results/{study}/associations/{study}_reference_comparison.tsv,
        	results/{study}/{locus}/{locus}_confidence.txt (eligible loci)
	output: results/{study}/study_finalized.txt
	script: pipeline_finemap/scripts/finalize_study.R


### o. build snps table

Builds the final flat, SNP-level output table for the study — one row per SNP across every fine-mapped locus, in the target reporting format.

#### Logic ():

	For each locus's susie_results.tsv:
    		Skip gracefully if empty/missing (SKIPPED or FAILED locus)
    		Look up lead SNP's BP to compute CS_BPS window bounds
    		Build one row per SNP: TRAIT, STUDY, LOCUS, CS_BPS, SNP, CHR, 			BP, REF, ALT, MAF, PVAL, Z, PIP, CS, CS_SIZE, LAMBDA
		Combine all loci's rows into one table

Relevant files:
---


### n. finalize study

Writes a summary banner into the study log once every locus has been scored — total/excluded/fine-mapped counts and confidence category breakdown.

#### Logic:

Read confidence category from every eligible locus's confidence file
Tabulate counts per category (HIGH/MEDIUM/LOW/UNRELIABLE)
Log banner with total loci, excluded count, fine-mapped count, category breakdown

#### Relevant files:

input:  results/{study}/associations/{study}_exclusions.tsv,
        results/{study}/associations/{study}_reference_comparison.tsv,
        results/{study}/{locus}/{locus}_confidence.txt (eligible loci)
output: results/{study}/study_finalized.txt
script: pipeline_finemap/scripts/finalize_study.R

### o. build snps table

Builds the final flat, SNP-level output table for the study — one row per SNP across every fine-mapped locus, in the target reporting format.

#### Logic:

For each locus's susie_results.tsv:
    Skip gracefully if empty/missing (SKIPPED or FAILED locus)
    Look up lead SNP's BP to compute CS_BPS window bounds
    Build one row per SNP: TRAIT, STUDY, LOCUS, CS_BPS, SNP, CHR, BP, REF, ALT,
        MAF, PVAL, Z, PIP, CS, CS_SIZE, LAMBDA
Combine all loci's rows into one table

#### Relevant files:

input:  results/{study}/associations/{study}_lead_snps.tsv,
        results/{study}/{locus}/{locus}_susie_results.tsv (eligible loci)
output: results/{study}/{study}_snps_table.tsv
script: pipeline_finemap/scripts/build_snps_table.R

### p. build final report

Produces the detailed, human-readable final report for the study — one entry per original lead association, stating whether it was fine-mapped, and if not, why; if so, its full score breakdown including specific missing high-LD neighbors and their R2 values.

#### Logic:

For each original lead SNP:
    If excluded (per exclusions table): report NO, with exclusion reason
    Else, read confidence file:
        If a hard-stop reason or error is recorded: report NO, with that reason
        Else (real scoring occurred):
            Report YES, final score/category
            Re-parse per-metric penalty breakdown from study log
            Pull specific missing-neighbor SNPs/R2 from ld_check file
            Report prior fine-mapping comparison line (informational)
Write full report

#### Relevant files:
input:  results/{study}/associations/{study}_lead_snps.tsv,
        results/{study}/associations/{study}_exclusions.tsv,
        results/{study}/{locus}/{locus}_confidence.txt (eligible loci),
        results/{study}/{locus}/LD_block_check/{locus}_ld_check.txt (all discovered loci)
output: results/{study}/{study}_final_report.txt
script: pipeline_finemap/scripts/build_final_report.R

### q. aggregator rules

Two lightweight touch-only rules exist purely to force checkpoint-dependent expansion to resolve for testing/debugging in isolation — not load-bearing for the final DAG, since run_susie/align_gwas_bim already depend on ld_matrix_and_freq's output directly.

Relevant files:

[aggregate_ld_matrices]
input:  results/{study}/{locus}/LD_matrix/{locus}.locus.matrix.ld (eligible loci)
output: results/{study}/ld_matrix_done.txt

[aggregate_susie_results]
input:  results/{study}/{locus}/{locus}_susie_summary.txt (eligible loci)
output: results/{study}/susie_done.txt
