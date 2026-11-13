# MetaStaph

MetaStaph is a Nextflow workflow for mixed-strain-aware analysis of *Staphylococcus* sequencing reads.

## Installation

### 1. Install Git LFS, then clone

The bundled database files are large and stored with [Git LFS](https://git-lfs.com/). Install Git LFS before cloning so the actual database files are downloaded instead of small pointer files. The StrainGE database is hosted separately on Figshare; download it in step 3 below.

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

With Nextflow 26.04 or later, select the legacy parser before running this workflow (also used in the local tests):

```bash
export NXF_SYNTAX_PARSER=v1
```

### 3. Download the StrainGE database

Download `pan-genome-db_9.hdf5` (approximately 1.67 GB) from [Figshare, version 1](https://doi.org/10.6084/m9.figshare.34999555.v1). From the repository directory, run:

```bash
mkdir -p databases
wget -c -O databases/pan-genome-db_9.hdf5.part https://ndownloader.figshare.com/files/71401612
echo 'd330826cc16216248c97aae331af4fa5  databases/pan-genome-db_9.hdf5.part' | md5sum -c - && \
  mv databases/pan-genome-db_9.hdf5.part databases/pan-genome-db_9.hdf5
```

The `.part` file allows interrupted downloads to resume and is renamed only after checksum verification. On macOS, use `md5 databases/pan-genome-db_9.hdf5.part` to compare against `d330826cc16216248c97aae331af4fa5`, then rename the verified file. The pipeline expects `databases/pan-genome-db_9.hdf5`; this file is downloaded separately and is not tracked in Git.

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

`Classified_reads` is the total classified fragments from the **first** Kraken2 pass (`--kraken_confidence`, default 0.1), across all classified taxa. Paired-end reads count as one fragment per pair. The count sums direct assignments, excluding unclassified reads, so nested taxonomy clade counts are not counted twice. It does not require taxid 1280 to appear in the database. The first-pass report is preserved as `<sample_id>.firstpass.kreport`; `<sample_id>.kreport` is the second-pass report used by LINtax (`--lintax_confidence`, default 0.4).

StrainScan retains `final_report.txt` alongside `strainscan.stdout.log`, `strainscan.stderr.log`, and `strainscan_status.tsv` in the sample's `strainscan_output/` directory. The status file never replaces the report. A valid report with a zero exit code is successful; a report with a nonzero exit code is marked `report_with_error` and remains available as evidence. Header-only reports are `no_hits`; a zero exit without a report is `no_report`. Failed executions retain their real exit code in the status note instead of assuming low sequencing depth. StrainScan receives absolute input/output paths, reserves eight CPUs for its bundled Jellyfish, and runs at most one sample at a time to limit simultaneous resource use.

The StrainScan database root must contain `Tree_database/tree_structure.txt`. MetaStaph checks this before launching tasks and recognizes the bundled nested layout `databases/strainscan_db/db/`. To use the database from a standalone StrainScan run, pass `--strainscan_db '/absolute/path/to/database/root'`.

LINtax's `Assigned_reads` includes assignments to the group's descendants. `Direct_Assigned_reads` counts only assignments exactly to the group node. `Unique_Assigned_reads` counts final Kraken assignments to the group or its descendants that belong to **exactly one listed LINgroup**. Shared ancestors outside a group are excluded; assignments covered by overlapping or duplicate groups are excluded from their unique counts to prevent double counting. For disjoint groups, unique and clade counts normally agree. These are group-specific taxonomic assignments, not a measure of sequence-level uniqueness or proof of a particular strain within a multi-strain group.

`Total_reads_length` includes both mates and descendant assignments. The `Total_reads` row reports all classified fragments, their total bases, the sum of exclusive group assignments, and direct group-node assignments. These totals include groups below the reporting cutoff; the raw report retains their individual counts. Percentages use all classified fragments as the denominator, including assignments above the listed groups. LIN prefixes are matched by comma-separated components so `0,1` cannot match `0,10`. Missing group prefixes are flagged in `LINtax_note` in the raw report.

LINtax reports require at least **two assigned fragments per group** by default (`--lintax_min_reads 2`). This removes singleton groups from `<sample_id>.LINreport.txt` and summary/evidence hit counts. The unfiltered counts remain in `<sample_id>.LINreport.raw.csv`. Increase the threshold to require more support, or set it to 1 to include singleton hits. This is a reporting cutoff, not a validated accuracy threshold. Multiple qualifying groups remain separate evidence hits.

StrainGE runs with a fixed internal output prefix, then renames results to retain the complete sample ID, including dots. Its stdout/stderr logs are published with the sample results. Failed or empty-input runs produce header-only strain tables with zero hits and a reason in the stats file. Historical `NA` placeholder rows are also excluded from hit counts.

After updating the pipeline, copy the updated `bin/` scripts as well as `main.nf` and `nextflow.config` to the execution machine. Rerun with `-resume`; changed tool tasks and reports need to run again. The tests can be run with `python -B -m unittest discover -s tests -v` in an environment containing pandas.

If Kraken fails, inspect `<sample_id>.kraken.stderr` in the failing task's work directory for its first-pass error, or `<sample_id>.extract_read.log` for a read-extraction error. After correcting the cause, rerun the same command with `-resume` to reuse cached successful tasks.
