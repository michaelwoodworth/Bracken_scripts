# Bracken scripts
Kraken2/Bracken result analytic workflow

## Motivation
[kraken2](https://github.com/DerrickWood/kraken2) is a tool that takes sequence data as input and returns an output table of taxonomic classifications.  [bracken](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual#step2) uses bayesian probabilities to estimate species (or other taxonomic level) relative abundance using each read-length kmer information from kraken2.  

This workflow documents input parameters for upstream steps to produce summary bracken relative abundance tables for reproducibility as well as custom python scripts written to process output tables for downstream analysis and visualization.

## Overall workflow:
1. Prepare reads
2. Classify reads with kraken2
3. Refine relative abundance estimates with bracken

	2. create binary presence/absence matrix ([01_bracken_binary_matrix.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_bracken_binary_matrix.py))
	3. estimate relative abundance of taxa ([02_bracken_validate_and_summarize_relabundance_RPKM.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/02_bracken_validate_and_summarize_relabundance_RPKM.py))

4. Analysis in R
5. Produce heatmaps in R

## Workflow detail

### 1. Prepare reads for assembly
Since I work with human-associated metagenomes, in addition to trimming with trimmomatic, I usually attempt to remove reads aligning to reference human genomes.  This can be efficiently done in one step with the Huttenhower lab tool [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/).

### 2. Use kraken2 to classify reads by taxonomy

**- 2.a: install kraken2 using conda**
**- 2.b: install bracken using conda**
**- 2.c: build kraken2 database(s)**

*default kraken2 installation instructions suggest:*

```bash
conda activate kraken
kraken2-build --standard --threads 20 --db ${db_path}/${db_name}
```

I have had difficulty using this command on a cluster (rsync errors, threading errors, etc).  As a workaround, the following files can be downloaded individually: 
- ncbi taxonomy
- kraken databases
Followed by running the build command as follows:

```bash
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

#build kraken2 database
kraken2-build --build --db kraken2 --threads 20

#build bracken database
bracken-build -d /nv/hmicro1/mwoodworth8/shared3/DB/kraken2 -t 20 -l 150 -t 20

#clean up intermediate files
kraken2-build --clean
```

**- 2.d: run kraken2**
```bash
kraken2 --db ${db_path}/${db_name} --output ${out_path}/${sample}.kraken2 --report ${out_path}/${sample}.report --paired ${sample}*R1.fastq ${sample}*2.fastq
```

**- 2.e: run bracken**
```bash
bracken -d ${db_path}/${kraken_db_name} -i ${sample}.kreport -o ${sample}.bracken -r ${read_length} -l ${classification_level} -t ${threshold}
```

### 3. Post-processing

**1. 00_bracken_filter.py** - this python script was written to optionally filter bracken output tables for minimum relative abundance.

```console
usage: 00_amrfinder_filter.py [-h] -i  -o  [-m] [-j]

Filter AMRFinder Plus results for high confidence matches.

This script filters AMRFinder output tables for matches, with
default criteria focused on high quality & complete matches.
e.g. >90% identity, >90% match length.

Script options also allow filtering for just AMR determinants.

optional arguments:
  -h, --help      show this help message and exit
  -i , --input    Please specify AMRFinder input tsv file name & path.
  -o , --output   Please specify AMRFinder filtered filename & path.
  -m , --method   Please specify filtered AMRFinder output tsv file name &
                  path. Select from: complete -or- add_partial_end
  -j, --just_AMR  Flag to remove non-AMR AMRFinder results (e.g. virulence,
                  stress)

```

**2. 01_amrfinder_binary_matrix.py** - this python script searches for .tsv files in an input directory and produces a binary presence/absence matrix for all genes across all samples coded as 0 for absent and 1 as present.

```console
usage: 01_amrfinder_binary_matrix.py [-h] -i INPUT -o OUTPUT [-v]

Create summary matrix of AMRFinder Plus results for plots & analysis.

This script sumamrizes filtered tables from 00_amrfinder_filter.

optional arguments:
  -h, --help      show this help message and exit
  -i , --input    Please specify input directory path.
  -o , --output   Please specify output filename & path.
  -v, --verbose   Increase output messaging detail, print results.
```

**3. 02_amrfinder_validate_and_summarize_RPKM.py** - this python script performs two main tasks.

- First, a validation step is completed to confirm the presence of all AMRFinder-detected genes in the [in situ gene coverage workflow](https://github.com/rotheconrad/00_in-situ_GeneCoverage/tree/6812ebd32c5127ce8b72ba8e520799b75f45c895) ${unique_id}gene_RPKM.tsv output file.

- Second, [RPKM](https://sites.google.com/site/wiki4metagenomics/pdf/definition/rpkm-calculation) values are tallied for duplicate gene names within a sample (which is printed to STDOUT if the verbose option is selected) and then used to construct a matrix of RPKM values with rows ov unique genes by columns of metagenomes.

```console
usage: 02_amrfinder_validate_and_summarize_RPKM.py [-h] -a  -m  -o  [-v] [-V]

Validate and summarize RPKM of AMRFinder-detected genes 
for plots & analysis.

This script takes the following inputs:
- directory containing filtered AMRFinder tsv files
  (output from step 00)
- directory containing ${uniqueID}_gene_RPKM.tsv files
  (https://github.com/rotheconrad/00_in-situ_GeneCoverage)

With intermediate validation steps (option -V):

- all genes input in AMRFinder tables are tested against all genes 
in the coverage_magic tsv files. If there are any genes that are 
not in the submitted coverage_magic tsv files, these are optionally 
output as: genes_to_validate.tsv

- all duplicated gene RPKM values are summed by sample.  
Input contigs/scaffolds hosting the detected genes are listed with
the summed (deduplicated) RPKM values and gene sequence name and
output as: deduplicated_RPKM.tsv

Genes that have RPKM values are returned with following output:
- specified output file & path containting a single tsv file with 
length / effort normalized AMR gene abundance (RPKM)

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