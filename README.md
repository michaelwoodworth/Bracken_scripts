# Bracken scripts
Kraken2/Bracken result analytic workflow

## Motivation
[kraken2](https://github.com/DerrickWood/kraken2) is a tool that takes sequence data as input and returns an output table of taxonomic classifications.  [bracken](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual#step2) uses bayesian probabilities to estimate species (or other taxonomic level) relative abundance using each read-length kmer information from kraken2.  

This workflow documents input parameters for upstream steps to produce summary bracken relative abundance tables for reproducibility as well as custom python scripts written to process output tables for downstream analysis and visualization.

## Overall workflow:
1. Prepare reads
2. Classify reads with kraken2
3. Refine relative abundance estimates with bracken

	1. estimate relative abundance of taxa ([01_bracken_summarize_relabundance.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_bracken_summarize_relabundance.py))

4. Analysis in R
5. Produce heatmaps in R

## Workflow detail

### 1. Prepare reads for assembly
Since I work with human-associated metagenomes, in addition to trimming with trimmomatic, I usually attempt to remove reads aligning to reference human genomes.  This can be efficiently done in one step with the Huttenhower lab tool [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/).

### 2. Use kraken2 to classify reads by taxonomy
This step assumes that you have already installed kraken2 & bracken.

*default kraken2 installation instructions suggest:*

```bash
conda activate kraken
kraken2-build --standard --threads 20 --db ${db_path}/${db_name}
```

I have had difficulty using this command on a cluster (e.g. rsync errors, multi-threading errors).  As a workaround, the following files can be downloaded individually: 
- ncbi taxonomy
- kraken databases

```bash
#activate conda environment
conda activate kraken

#download NCBI taxonomy
kraken2-build --download-taxonomy --db kraken2 

#download kraken2 database components
kraken2-build --download-library bacteria --db kraken2 --threads 20
kraken2-build --download-library archaea --db kraken2
kraken2-build --download-library plasmid --db kraken2
kraken2-build --download-library viral --db kraken2
kraken2-build --download-library fungi --db kraken2
kraken2-build --download-library protozoa --db kraken2
```

- Followed by running the build command as follows:

```bash
#build kraken2 database
kraken2-build --build --db kraken2 --threads 20

#build bracken database
bracken-build -d /nv/hmicro1/mwoodworth8/shared3/DB/kraken2 -t 20 -l 150 -t 20

#clean up intermediate files
kraken2-build --clean
```

**Run kraken2**
```bash
kraken2 --db ${db_path}/${db_name} --output ${out_path}/${sample}.kraken2 --report ${out_path}/${sample}.report --paired ${sample}*R1.fastq ${sample}*2.fastq
```

**Run bracken**
```bash
bracken -d ${db_path}/${kraken_db_name} -i ${sample}.kreport -o ${sample}.bracken -r ${read_length} -l ${classification_level} -t ${threshold}
```

### 3. Post-processing

**01_bracken_summarize_relabundance.py** - this python script takes a directory containing ${sample}.bracken files as input and produces a summary tsv file with taxa as rows and samples as columns with relative abundance cell values.

```console
usage: 01_bracken_summarize_relabundance.py [-h] -a  -m  -o  [-v] [-V]

Summarize relative abundance of kraken2/bracken classified 
taxa for plots and analysis.

This script takes the following input:
- directory containing ${sample}.bracken tsv files

The script output is a tsv matrix with unique taxa across all samples as rows, samples as columns, and bracken-estimated relative abundance of reads.

optional arguments:
  -h, --help            show this help message and exit
  -a , --amrfinder_tsv_path 
                        Please specify directory path containing filtered
                        AMRFinder tsv files.
  -m , --coverage_magic_path 
                        Please specify directory path containing coverage
                        magic path.
  -o , --output         Please specify output file path (& optional prefix).
  -v, --verbose         Toggle volume of printed output.
  -V, --validate        Write genes_to_validate.tsv and deduplicated.tsv.
```