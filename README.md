# nORFs, HARs and TEs
Instructions, R code, selected input files and supplementary data for a project studying regulation of nORF transcript expression in brain disorders.

The latest version of R is available at: https://cran.r-project.org/.
It may be helpful to use RStudio, which is available at: https://rstudio.com/

The PLINK 1.90 beta is available at: https://www.cog-genomics.org/plink2
INRICH is available at: https://atgu.mgh.harvard.edu/inrich/

featureCounts can easily be installed on Ubuntu by ```sudo apt-get install featureCounts``` in the Terminal.
For other systems, the instructions for installation are available at: http://bioinf.wehi.edu.au/subread-package/

Data that were not publicly available and were not generated in this project is zipped and password-protected. To access these, send @tonilogbo a message.

The following files are too large to be stored in this repository, which affects several code files. To obtain these, send @tonilogbo a message.
> 2.txt

> GRCh37_GENCODE_rmsk_TE.gtf

> rmsk.saf


# Enrichment of HAR-associated nORFs in disorder-linked loci
Schizophrenia SNPs are available at: https://walters.psycm.cf.ac.uk/

Bipolar disorder SNPs are available at: https://www.med.unc.edu/pgc/download-results/.

Genotype data from 1000 Genomes' EUR population is available at: https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip

To set up files for use in PLINK and INRICH, use pre-process.R

## Plink commands:

Clumping:
``` 
./plink --bfile g1000_eur \
--maf 0.02 \
--clump [input] \
--clump-p1 0.1 \
--clump-p2 1 \
--clump-r2 0.25 \
--clump-kb 500 \
--out [output]
```

LD intervals for INRICH:
```
./plink --bfile g1000_eur \
--hwe 0.001 \
--maf 0.05 \
--show-tags [snp_list] \
--list-all \
--tag-r2 0.5 \
--tag-kb 250 \
--out [output] && \
gawk ' NR>1 { print $2,$5,$6 } ' [output].tags.list > [interval_file].int
```

## INRICH command:
```
./inrich \
-a [interval_file].int \
-g norfs_rgf.gene.map \
-t [set_file].set \
-m [snp_map_file].snp.map \
-p 1 \
-o [output]
```

# Weighted Transcript Co-expression Network Analysis
WTCNA.R contains code used to perform the WTCNA in this project. 

# Correlation of nORF and TE transcript expression
To obtain TE counts from BAM files using featureCounts:
```
featureCounts -M -F SAF -T 1 -s 2 -p -a rmsk.saf -o outfeatureCounts.txt *.bam
```
Pre-process the featureCounts output with counts.R and use correlation.R contains code to identify correlations between TE and nORF expression.

# Other supplementary data for reference
allStrongCorrelations.csv contains all strong TE-nORF correlations with an absolute Spearman/Pearson correlation > 0.8
