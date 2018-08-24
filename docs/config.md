# Configuration

The configuration file that you define takes the form of a JSON file.
It contains several variables, which are outlined here.
For a complete, functioning configuration file, that works for both the markerquant and alignment pipelines, see [the example JSON file](https://github.com/thiesgehrmann/markerQuant/blob/master/testData/config.json).
 

## Data/Output variables

|Variable|Description|
|--------|-----------|
|outdir|Where all output will be located|
|targets|The targets to quantify using the marker quant pipeline.|
|genomes|The genome in which markers should be unique to in the marker quant pipeline, and to which reads are aligned in the alignment pipeline.|
|transcriptome|The transcriptome to which markers should be unique (if same as targets, can be left empty).|
|genes|The GFF3 file describing genes on the genome.|
|samples|A structure of samples, which contain locations for FASTQ files. (Note, they should not be gzipped).|
|tests| A list of conditions between differential expression tests will be conducted.| 

## Configuration variables

To see the default values for these variables, see [the defaults.JSON file](https://github.com/thiesgehrmann/markerQuant/blob/master/pipeline_components/defaults.json).

|Variable|Type|Description|
|--------|----|-----------|
|strandSpecific|Boolean|Are the libraries strand specific, or not?|
|minQual|Integer|Minimum base quality to permit match in marker quant pipeline.|
|k|Integer|Marker length.|
|targetMap|File|A tsv file with two columns which contains a mapping between a target name and another name.|
|knockout|File|A file containing a list of targets that are supposed to be knocked out. Reads mapping to these targets will be output in a seperate file.|
|remove_pcr_duplicates|Boolean|Remove PCR duplicates in alignment pipeline?|
|remove_mismatched_reads|Boolean|Remove reads that have mismatches.|
|remove_mismatched_reads_max|Integer|Together with above parameter, allow for at most how many mismatches before removing?|
|star_index_params|String|Additional parameters for the STAR aligner index generation step.|
|star_params|String|Additional parameters for the STAR alignment step.|
|htseq_params|String|Additional parameters for htseq-count.|
|htseq_t|String|The GFF feature type to count in htseq-count.|
|htseq_i|String|The GFF attribute to aggregate features in htseq-count.|
|fragLengthMeanSubsampleSize|Integer|To calculate FPKM or TPM normalizations, the number of  alignments to subsample from the alignment to estimate the mean fragment length.|
|norm_method|String|tmm/fpkm/tmm. Which normalization method to use.|  

