# markerQuant
Determine transcriptome expression based on unique markers

## Dependencies

  * Java
  * Snakemake
  * Conda

## Usage

The example provided will run for approximately 3.4 minutes on a single core machine

    git clone https://github.com/thiesgehrmann/markerQuant.git
    cd markerQuant

    snakemake --use-conda --configfile testData/config.json quantifyTargets # Unstranded RNA-Seq
    snakemake --use-conda --configfile testData/config_stranded.json quantifyTargets # Stranded RNA-Seq

You can compare the results from a normal RNA-Seq pipeline (STAR + htseq_count)
  snakemake --use-conda --configfile testData/config.json --snakefile align.Snakefile quantifyTargets
  snakemake --use-conda --configfile testData/config_stranded.json --snakefile align.Snakefile quantifyTargets

### Tasks

There are several tasks that you can run

 * **generateMarkers**: Generate markers for the targets specified in the config file
 * **quantifyTargets**: Quantify the targets specified in the config file
 * **deseqTests**: Perform the differential expression tests specified in the config file.
 * **deseqNorm**: Normalize the expression based on the library size
 * **deseq**: Normalize and do the tests

## Output files

Each of these tasks generates a different set of output files

### generateMarkers

 * *run/markers/markers.fasta*: Aggregated markers to use in the tree
 * *run/markers/gaps.fasta*: Gaps in the aggregated markers that cannot be used for identification

### quantifyTargets

 * *run/quantification/target_counts.`sample`.tsv*: A table describing, in sample `sample`, the number of markers per gene, the number of counts per marker, and the average and standard deviation of these.
 * *run/quantification/quantification.tsv`*: Essentially the table of expression counts that need to go into DESeq. Rows are genes, columns are samples.

### deseqTests

* *run/diffex/tests.tsv*: A table of differential expression tests

### deseqNorm

* *run/diffex/quantification.normalized.tsv*: Normalized transcript abundance
  
