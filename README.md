# Saccharomycotina_subphylum_analysis
An analysis of mutation bias and natural selection shaping codon usage bias across the 327 Saccharomycotina budding yeast subphylum.


## Obtaining the data

Sequence data, ortholog matrices, and the species were obtained from [Figshare](https://doi.org/10.6084/m9.figshare.5854692.v1). CDS files were cleaned to remove any sequences that do not start with a canonical start codon (ATG), had a length that was not a multiple of 3, or showed high homology to mitochondiral protein-coding sequences. The mitochondrial protein-coding sequences were obtained from [MitoFun](http://mitofun.biol.uoa.gr/). Unfortunately, this resource is no longer available. We have provided FASTA files containing the sequences removed from the files analyzed via ROC-SEMPPR. 

## Running ROC-SEMPPR via `AnaCoDa`

AnaCoDa is currently available on CRAN. However, as developers of the software project, we have implemented new functionality, performance improvements, and bug fixes that have not yet been incorporated into the version on CRAN. We recommend using the version currently available on the `master` branch for [the RibModelFramework GitHub repo](https://github.com/acope3/RibmodelFramework). We also note that a version was implemented on branch `CTG_Ser` that is able to handle the case of yeasts with CTG coding for serine as opposed to leucine. The functionality of this branch is consistent with the `master` branch, aside from some changes to the data structures for mapping codons to amino acids. Specifically, this code changes a static data structure (i.e. initialized at compile time) to a non-static data structure, but this seems to increase run times. Due to this performance issue, this code has not been merged into the `master` branch.

For running ROC-SEMPPR with the `CTG_Ser` branch to account for the CTG-Ser yeasts, specifify the command line argument `--codon_table 12`. See `R_scripts/rocAnalysis.R` for descriptions of the command line arguments or run `
```
INPUT=<path to FASTA file containing CDS>
OUTPUT=<directory for output>
Rscript --vanilla R_scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 20 -s 20000 -a 20 -t 10 -n 8 --est_csp --est_phi --est_hyp  --max_num_runs 2
Rscript --vanilla R_scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 20 -s 20000 -a 20 -t 10 -n 8 --est_csp --est_phi --est_hyp  --max_num_runs 2 --codon_table 12
```

Model fits can be found in `Final_runs/`. `Result_k_1` and `Results_k_2_selectionShared` correspond to the ConstMut and VarMut model fits, respectively.  

## Empirical gene expresison

Empirical gene expression data were taken from RNA-seq data as described in [Cope and Shah Plos Genetics 2022](https://doi.org/10.1371/journal.pgen.1010256). Upon downloading the raw data using your preferred tool (we used `fasterq-dump`), the data were processed as follows. 


For paired-end data, RNA-seq was processed using the following format.

```bash
fastp -i "${SPECIES_ACC}_1.fastq" -I "${SPECIES_ACC}_2.fastq" -o "${SPECIES_ACC}_1_trimmed.fastq" -O "${SPECIES_ACC}_2_trimmed.fastq" -w 8 -j "${SPECIES_ACC}_fastp.json" -h "${SPECIES_ACC}_fastp.html"
kallisto quant -i ${1}/${1}.index -o "${SPECIES_ACC}_tpm_kallisto" --bias -t 8 "${SPECIES_ACC}_1_trimmed.fastq" "${SPECIES_ACC}_2_trimmed.fastq"
```		

For single-end data, RNA-seq was processed using the following format.

```bash
fastp -i "${SPECIES_ACC}.fastq" -o "${SPECIES_ACC}_trimmed.fastq" -w 8 -j "${SPECIES_ACC}_fastp.json" -h "${SPECIES_ACC}_fastp.html"
kallisto quant -i ${1}/${1}.index -o "${SPECIES_ACC}_tpm_kallisto" --bias --single -l 200 -s 25 -t 8 "${SPECIES_ACC}_trimmed.fastq"
```

| Software      	| version     	|
|---------------	|-------------	|
| fasterq-dump     	| 2.9.2         |
| fastp      	    | 0.21.0        |
| kallisto      	| 0.46.2      	|


## tRNA genes

We identified tRNA genes using `tRNAscan-SE 2.0` using the default settings for eukaryotes. The tRNAscan-SE output can be found in `Data/new_tGCNs`.

## Re-creating results from Cope and Shah 2025

Results from our manuscript can be recreated using the R markdown files found in the `R_notebooks/` directory. Running these in order should re-create the the analyses presented ingit  the main text and SI Appendix. There may also be some analyses that did not make it to the final version. To run these files, you will need will need the following R packages (as well as any associated dependencies). Note that `R_notebooks/00_create_csp_data_matrices.R`, `R_notebooks/01_determine_best_model_fit.R`, and `R_notebooks/02_get_elongation_waiting_times.R` are mostly intended for setting up downstream analyses and some basic sanity checks, such as the correlation between amino acid frequency and tGCN. These notebooks can likely be skipped, as the output of these notebooks should be included in `Post_analysis/`.

| Software          | version       |
|---------------    |-------------  |
| AnaCoDa           | 0.1.4.4       |
| tidytext          | 0.4.2         |
| preprocessCore    | 1.66.0        |
| phylolm           | 2.6.5         |
| Biostrings        | 2.72.1        |
| ggnewscale        | 0.5.0         |
| geiger            | 2.0.11        |
| ggmagnify         | 0.4.1         |
| ComplexHeatmap    | 2.20.0        |
| scales            | 1.3.0         |
| patchwork         | 1.3.0         |
| ggtree            | 3.12.0        |
| broom             | 1.0.7         |
| ggrepel           | 0.9.6         |
| ggpubr            | 0.6.0         |
| reshape2          | 1.4.4         |
| cowplot           | 1.1.3         |
| RColorBrewer      | 1.1-3         |
| dendextend        | 1.18.1        |
| factoextra        | 1.0.7         |
| tidyverse         | 2.0.0         |
| ggplot2           | 3.5.1         |
| cluster           | 2.1.6         |
| phylogram         | 2.1.0         |
| phytools          | 2.3-0         |
| ape               | 5.8           |


