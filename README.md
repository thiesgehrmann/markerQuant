# markerQuant
Determine transcriptome expression based on unique markers.
Useful in deconvolving the expression of several highly homologous transcripts

This pipeline is made up of three steps:
1. Marker identification

   For each given target, identify markers that are unique to that target with respect to a) the other targets, b) the genome and c) the transcriptome (optionally).
2. Marker quantification

   Using those markers, quantify the expression of each target in the dataset
3. Differential expression

   Use DESeq2 to normalize the expression of the targets, and determine differentially expressed targets.

## Dependencies

  * [Snakemake](http://snakemake.readthedocs.io)
  * [Conda](https://conda.io/miniconda.html)
    * With [Bioconda](https://bioconda.github.io/)

## Usage

We provide two pipelines, written using [Snakemake](http://snakemake.readthedocs.io):
 * *Snakefile*: The MarkerQuant pipeline
 * *align.Snakefile*: A [STAR](https://github.com/alexdobin/STAR)/[HTSeq-count](http://www-huber.embl.de/HTSeq/doc/overview.html)/[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) rnaseq pipeline

### Example case

The example provided will run for approximately 3.4 minutes on a single core machine

    git clone https://github.com/thiesgehrmann/markerQuant.git
    cd markerQuant

    ./mq -a quantifyTargets testData/config.json            # Unstranded RNA-Seq
    ./mq -a quantifyTargets testData/config_stranded.json   # Stranded RNA-Seq
    ./mq -a quantifyTargets estData/config_errorfree.json  # Unstranded, error-free RNA-Seq


You can compare the results from a normal RNA-Seq pipeline:

    ./mq -A -a quantifyTargets testData/config.json          # Unstranded RNA-Seq
    ./mq -A -a quantifyTargets testData/config_stranded.json # Stranded RNA-Seq

### Actions

With the `-a` option, you can specify different actions:

 * **generateMarkers**: Generate markers for the targets specified in the config file
 * **quantifyTargets**: Quantify the targets specified in the config file
 * **deseq**: Normalize and do the tests

## Output files

Each of these tasks generates a different set of output files

### generateMarkers

 * *run/markers/markers.fasta*: Aggregated markers to use in the tree
 * *run/markers/gaps.fasta*: Gaps in the aggregated markers that cannot be used for identification

### quantifyTargets

 * *run/quantification/target_counts.`sample`.tsv*: A table describing, in sample `sample`, the number of markers per gene, the number of counts per marker, and the average and standard deviation of these.
 * *run/quantification/quantification.tsv*: Essentially the table of expression counts that need to go into DESeq. Rows are genes, columns are samples.

### deseq

* *run/diffex/tests.tsv*: A table of differential expression tests
* *run/diffex/quantification.normalized.tsv*: Normalized transcript abundance


## Configuration

Configuration is given in a json file, (example in `testData/config.json`).

```python
{
  "outdir" : "testOutput/unstranded",            # Give the output directory you want to use, should be absolute
  "targets" : "testData/tdh_genes.fasta",        # A fasta file of the targets you wish to generate markers for
  "genomes" : "testData/genome.fasta",           #   which are unique relative to this genome
  "transcriptome" : "testData/tdh_genes.fasta",  #   and this transcriptome (can be omitted if same as targets file)
  "genes"   : "testData/tdh_genes.gff",          # Needed for the traditional pipeline in align.Snakefile
  "strandSpecific" : 0,                          # 0 if not strand specific, 1 it yes
  "targetMap": "",                               # If your transcript names are esoteric, you can map them to useful names with this file if necessary
  "minQual" : 25,                                # The minimum PHRED33 quality score to use for read regions that hit a marker
  # Provide here your fastq files, defining for each the sample identifier (e.g. sample_1_r1), a replicate group (e.g. sample_1), and a list of fastq files.
  # Can be single or paired end, just provide an array with one element if single ended.
  "samples" : { "sample_1_r1" : { "replicate_group" : "sample_1",  "fastq" : [ "testData/fastq/unstranded/sample_01_1.fastq", "testData/fastq/unstranded/sample_01_2.fastq" ]} ,
                "sample_1_r2" : { "replicate_group" : "sample_1",  "fastq" : [ "testData/fastq/unstranded/sample_02_1.fastq", "testData/fastq/unstranded/sample_02_2.fastq" ]} ,
                "sample_2_r1" : { "replicate_group" : "sample_2",  "fastq" : [ "testData/fastq/unstranded/sample_03_1.fastq", "testData/fastq/unstranded/sample_03_2.fastq" ]} ,
                "sample_2_r2" : { "replicate_group" : "sample_2",  "fastq" : [ "testData/fastq/unstranded/sample_04_1.fastq", "testData/fastq/unstranded/sample_04_2.fastq" ]} },
  "tests" : [ ["sample_1", "sample_2"] ], # A list of differential expression tests to perform
  "htseq_t" : "gene", # For HTSeq-count in align.Snakefile
  "htseq_i" : "gene"  # For HTSeq-count in align.Snakefile
}
```
  
