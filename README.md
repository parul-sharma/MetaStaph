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
