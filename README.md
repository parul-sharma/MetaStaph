# MetaStaph

MetaStaph is a Nextflow workflow for mixed-strain-aware analysis of *Staphylococcus* sequencing reads.

## Installation

### 1. Install Git LFS, then clone

The database files are large and stored with [Git LFS](https://git-lfs.com/). Install Git LFS before cloning so the actual database files are downloaded instead of small pointer files.

On Ubuntu or Debian:

```bash
git lfs install
git clone https://github.com/parul-sharma/MetaStaph.git
cd MetaStaph
git lfs pull
```

If Git LFS is not available from your system package manager, install it with Conda:

```bash
conda install -c conda-forge git-lfs
git lfs install
```

Then run the clone commands above. If you already cloned the repository without Git LFS installed, install and initialize Git LFS, change into the repository, and run:

```bash
git lfs pull
```

To check the LFS download status, run:

```bash
git lfs ls-files | head
```
In this listing, `*` means the file content is present; `-` means only its pointer is checked out. If database files show `-`, rerun `git lfs pull`. You can check the total database size with `du -sh databases`.

### 2. Create MetaStaph env using the yaml file

```bash
conda env create -f metastaph_env.yml
conda activate metastaph
```

The included [nextflow.config](nextflow.config) enables Conda environments and prefers Mamba. Nextflow will create the tool environments defined under `envs/` on the first run.

## Input data

The workflow currently expects **paired-end FASTQ** reads. For each sample, provide two matching files named with `_R1` and `_R2`, for example:

```text
sample/
  isolateA_R1.fastq.gz
  isolateA_R2.fastq.gz
  isolateB_R1.fastq.gz
  isolateB_R2.fastq.gz
```

## Run

From the repository directory, run the workflow with the default input pattern (`sample/*_R{1,2}.fastq.gz`):

```bash
nextflow run main.nf
```
To specify the read directory/pattern, database directory, and output directory:

```bash
nextflow run main.nf \
  --reads 'reads/*_R{1,2}.fastq.gz' \
  --database_dir "$PWD/databases" \
  --outdir "$PWD/output"
```

Use `nextflow run main.nf --help` to display the workflow's usage summary. Results are written to the output directory, including per-sample reports and combined summary/evidence files.

Relative database paths are resolved from the directory where you launch Nextflow before tasks run in their work directories. With older pipeline versions, pass an absolute path such as `--database_dir "$PWD/databases"` to avoid database lookup failures.

## Reports and tool diagnostics

The summary columns are `Sample`, `Classified_reads`, `Kraken_note`, `Tools_reporting_n`, `Tools_reporting`, `Tool_hit_counts`, `Total_tool_hits`, and `Notes`.

`Classified_reads` is the total classified fragments from the **first** Kraken2 pass (`--kraken_confidence`, default 0.1), across all classified taxa. Paired-end reads count as one fragment per pair. The count sums direct assignments, excluding unclassified reads, so nested taxonomy clade counts are not counted twice. It does not require taxid 1280 to appear in the database. The first-pass report is preserved as `<sample_id>.firstpass.kreport`; `<sample_id>.kreport` is the stringent second-pass report used by LINtax.

StrainScan retains `final_report.txt` alongside `strainscan.stdout.log`, `strainscan.stderr.log`, and `strainscan_status.tsv` in the sample's `strainscan_output/` directory. The status file never replaces the report. A valid report with a zero exit code is successful; a report with a nonzero exit code is marked `report_with_error` and remains available as evidence. Header-only reports are `no_hits`; a zero exit without a report is `no_report`. Failed executions retain their real exit code in the status note instead of assuming low sequencing depth. StrainScan receives absolute input/output paths, reserves eight CPUs for its bundled Jellyfish, and runs at most one sample at a time to limit simultaneous resource use.

StrainGE runs with a fixed internal output prefix, then renames results to retain the complete sample ID, including dots. Its stdout/stderr logs are published with the sample results. Failed or empty-input runs produce header-only strain tables with zero hits and a reason in the stats file. Historical `NA` placeholder rows are also excluded from hit counts.

After updating the pipeline, copy the updated `bin/` scripts as well as `main.nf` and `nextflow.config` to the execution machine. Rerun with `-resume`; changed tool tasks and reports need to run again. The tests can be run with `python -B -m unittest discover -s tests -v` in an environment containing pandas.

If Kraken fails, inspect `<sample_id>.kraken.stderr` in the failing task's work directory for its first-pass error, or `<sample_id>.extract_read.log` for a read-extraction error. After correcting the cause, rerun the same command with `-resume` to reuse cached successful tasks.
