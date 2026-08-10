# Container image for pipeline_munge/Snakefile.
#
# Published to ghcr.io/benaroyaresearch/fm_2026 by .github/workflows/build-image.yaml.
#
# ---------------------------------------------------------------------------
# Why this base, and not neurogenomicslab/mungesumstats
# ---------------------------------------------------------------------------
# The official MungeSumstats image is ~17.6 GB compressed / ~40 GB on disk, because it
# bakes in the SNPlocs and BSgenome reference data packages. That is too large to build
# on a GitHub-hosted runner, and under the Kubernetes executor every node would pull it
# before the first munge job could start.
#
# This image instead starts from the same R/Bioconductor foundation the official image
# is built on — bioconductor/bioconductor_docker RELEASE_3_20 (R 4.4.2, Bioconductor
# 3.20, Ubuntu 22.04 "jammy") — and installs MungeSumstats from that frozen Bioconductor
# release. So MungeSumstats resolves to the same version the official container ships;
# only the multi-GB reference data is left out. That data lives in a shared R library on
# NFS instead (see "Reference data" below), which keeps this image around 2-3 GB.
#
# ---------------------------------------------------------------------------
# Why Snakemake is installed here
# ---------------------------------------------------------------------------
# The Kubernetes executor does not read the Snakefile's `container:` directive. It puts
# a single image — whatever `--container-image` names — on every job pod, and runs
# Snakemake *itself* inside that pod to execute the rule (the plugin's default for this
# flag is snakemake/snakemake:v<version>). The pod therefore re-parses the Snakefile,
# which imports pandas at module scope, so both Snakemake and pandas must be present.
#
# Snakemake 8+ requires Python >= 3.11 and jammy ships 3.10, so uv provides a managed
# 3.12 rather than pulling in a PPA. SNAKEMAKE_VERSION must match the driver running in
# the Coder workspace (`conda activate snakemake; snakemake --version`) — a mismatched
# pod version is the usual cause of jobs failing immediately with argument errors.
#
# ---------------------------------------------------------------------------
# Reference data (one-time setup, outside this image)
# ---------------------------------------------------------------------------
# munge_sumstats.R loads four data packages that are NOT in this image. Install them
# once into a shared NFS directory, using this same image so they are built against the
# identical R 4.4 / Bioconductor 3.20:
#
#   docker run --rm -v /nfs/<your-path>:/refdata ghcr.io/benaroyaresearch/fm_2026:latest \
#     Rscript -e 'BiocManager::install(c(
#         "SNPlocs.Hsapiens.dbSNP155.GRCh37", "SNPlocs.Hsapiens.dbSNP155.GRCh38",
#         "BSgenome.Hsapiens.NCBI.GRCh38", "BSgenome.Hsapiens.1000genomes.hs37d5"),
#       lib = "/refdata", ask = FALSE, update = FALSE)'
#
# Budget ~13 GB and a long download. Then point `ref_lib` in config/config_munge.yaml at
# that directory; the munge_gwas rule exports it as R_LIBS. Because the packages are
# compiled against a specific R minor version, bumping this image's Bioconductor release
# means reinstalling them.
#
# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
#   docker build --platform linux/amd64 -t ghcr.io/benaroyaresearch/fm_2026:dev .
#
# The base is pinned by digest so a rebuild resolves the same R and Bioconductor.

FROM bioconductor/bioconductor_docker:RELEASE_3_20@sha256:5053f59de1c1a79b5a05d820ba722bae62ded0c1042e547189aa0226f111ea7e

LABEL org.opencontainers.image.title="fm_2026" \
      org.opencontainers.image.description="pipeline_munge: GWAS download, coverage QC, and MungeSumstats munge to GRCh38" \
      org.opencontainers.image.source="https://github.com/BenaroyaResearch/FM_2026"

# ---------------------------------------------------------------------------
# 1. Command-line tools used by the shell rules
# ---------------------------------------------------------------------------
# resolve_gwas_filename : curl, grep -oP (PCRE)
# download_gwas         : curl, mktemp, file, gzip
# check_gwas_coverage   : zcat (gzip), head, tail, wc, mktemp
# munge_gwas            : timeout (coreutils), Rscript
#
# coreutils/grep/gzip are Ubuntu essential packages and already present; they are listed
# so the dependency is recorded rather than assumed. `file` is the one genuinely absent
# from most R images, and download_gwas relies on it to tell whether the GWAS Catalog
# returned a gzip or a plain TSV.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        coreutils \
        curl \
        file \
        grep \
        gzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Fail the build if any tool is missing, or if grep lacks PCRE — resolve_gwas_filename's
# lookbehind patterns silently match nothing on a grep built without it.
RUN for tool in curl file gzip zcat grep head tail wc mktemp timeout date Rscript; do \
        command -v "$tool" >/dev/null || { echo "MISSING TOOL: $tool" >&2; exit 1; }; \
    done \
    && echo 'href="GCST123.h.tsv.gz"' | grep -oP '(?<=href=")[^"]*\.h\.[^"]*' \
    && echo "command-line tool check OK"

# ---------------------------------------------------------------------------
# 2. R packages
# ---------------------------------------------------------------------------
# BiocManager resolves MungeSumstats' Imports (GenomicRanges, rtracklayer,
# VariantAnnotation, BSgenome, ...) automatically. The base sets BIOCONDUCTOR_VERSION and
# BIOCONDUCTOR_USE_CONTAINER_REPOSITORY=TRUE, so these arrive as prebuilt binaries rather
# than source builds.
#
# The four SNPlocs/BSgenome data packages are deliberately excluded — see the header.
RUN Rscript -e ' \
        pkgs <- c("MungeSumstats", "data.table", "R.utils"); \
        have <- pkgs %in% rownames(installed.packages()); \
        cat("already present:", paste(pkgs[have], collapse = " "), "\n"); \
        cat("to install:", paste(pkgs[!have], collapse = " "), "\n"); \
        if (any(!have)) BiocManager::install(pkgs[!have], ask = FALSE, update = FALSE, Ncpus = 4) \
    ' \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

# ---------------------------------------------------------------------------
# 3. Snakemake + pandas
# ---------------------------------------------------------------------------
# Everything lands under /usr/local and /opt rather than /root: the Coder workspaces map
# POSIX UIDs from Azure AD, so job pods do not necessarily run as root and must still be
# able to read these.
ARG SNAKEMAKE_VERSION=9.25.1
ENV UV_TOOL_BIN_DIR=/usr/local/bin \
    UV_TOOL_DIR=/opt/uv-tools \
    UV_PYTHON_INSTALL_DIR=/opt/uv-python

# uv comes from its own published image rather than `curl | sh`. That is what uv's docs
# recommend for Dockerfiles: the version is pinned and visible, and it is a registry pull
# rather than an HTTPS fetch from inside the build, which is one less thing to fail behind
# a TLS-intercepting proxy.
COPY --from=ghcr.io/astral-sh/uv:0.12.3 /uv /usr/local/bin/uv

RUN uv python install 3.12 \
    && uv tool install --python 3.12 "snakemake==${SNAKEMAKE_VERSION}" --with pandas \
    && chmod -R a+rX /opt/uv-tools /opt/uv-python \
    && rm -rf /root/.cache/uv

# ---------------------------------------------------------------------------
# 4. Build-time verification
# ---------------------------------------------------------------------------
# Loads the two libraries munge_sumstats.R needs that are actually in this image, and
# checks the entry points it calls. The SNPlocs/BSgenome packages cannot be checked here
# because they live on NFS — the munge_gwas rule fails loudly at runtime if R_LIBS is
# wrong.
#
# Also confirms the chain file MungeSumstats bundles is present: liftover to GRCh38 first
# tries to fetch a chain from Ensembl and falls back to this copy, which is what lets the
# munge work on a pod with no outbound network.
RUN Rscript -e ' \
        library(MungeSumstats); \
        library(data.table); \
        stopifnot(is.function(MungeSumstats::format_sumstats), \
                  is.function(MungeSumstats::get_genome_builds), \
                  is.function(data.table::fwrite)); \
        chain <- system.file("extdata", "GRCh37_to_GRCh38.chain.gz", package = "MungeSumstats"); \
        stopifnot(nzchar(chain), file.exists(chain)); \
        cat("R dependency check OK\n"); \
        cat("  MungeSumstats", as.character(packageVersion("MungeSumstats")), "\n"); \
        cat("  R", as.character(getRversion()), "\n"); \
        cat("  bundled chain file:", chain, "\n") \
    '

# The pod runs `sh -c "snakemake ..."`, so both must resolve on a bare PATH, and the
# Snakefile's own `import pandas` must succeed under that interpreter.
RUN snakemake --version \
    && /opt/uv-tools/snakemake/bin/python -c \
        "import pandas, snakemake; print('pandas', pandas.__version__, '| snakemake', snakemake.__version__)" \
    && echo "snakemake/pandas check OK"

# The base is an RStudio Server image whose CMD is the s6 init supervisor. Snakemake
# overrides the command anyway, but a shell is the sane default for a batch image.
CMD ["/bin/bash"]
