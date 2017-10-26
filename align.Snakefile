import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__
__JAR__         = "%s/markerQuant/markerQuant.jar" % __INSTALL_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)

###############################################################################

__RUN_DIR__ = os.path.abspath(dconfig["outdir"]) + "/align_run"
__STAR_OUTDIR__ = "%s/star_align" % __RUN_DIR__
__ALNFILTER_OUTDIR__ = "%s/align_filter" % __RUN_DIR__
__QUANT_OUTDIR__ = "%s/quantification" % __RUN_DIR__
__DIFF_OUTDIR__ = "%s/diffex" % __RUN_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)


###############################################################################

rule all:
  input:
    quant = "%s/quantification.tsv" % __QUANT_OUTDIR__,
    norm  = "%s/quantification.normalized.tsv" % __DIFF_OUTDIR__,
    tests = "%s/tests.tsv" % __DIFF_OUTDIR__

###############################################################################

rule star_index:
  input:
    genomes = dconfig["genomes"]
  output:
    index = "%s/star.idx" % __STAR_OUTDIR__
  conda: "%s/env.yaml" % __PC_DIR__
  params:
    star_index_params = dconfig['star_index_params']
  shell: """
    mkdir {output.index}
    STAR --runMode genomeGenerate {params.star_index_params} --genomeDir {output.index} --genomeFastaFiles {input.genomes}
  """

rule star_pass_one_sample:
  input:
    fastq = lambda wildcards: dconfig["samples"][wildcards.sample]["fastq"],
    index = rules.star_index.output.index
  output:
    sj = "%s/aln1.{sample}.SJ.out.tab" % __STAR_OUTDIR__
  conda: "%s/env.yaml" % __PC_DIR__
  threads: 5
  params:
    rule_outdir = __STAR_OUTDIR__,
    star_params = dconfig["star_params"]
  shell: """
    STAR --runMode alignReads --runThreadN {threads} {params.star_params} --genomeDir {input.index} --readFilesIn {input.fastq} --outFileNamePrefix {params.rule_outdir}/aln1.{wildcards.sample}. --outSAMmode None
  """
    

rule star_pass_one:
  input:
    sj = expand("%s/aln1.{sample}.SJ.out.tab" % __STAR_OUTDIR__, sample=dconfig["samples"].keys())
    
    
rule star_pass_two_sample:
  input:
    index = rules.star_index.output.index,
    sj    = rules.star_pass_one.input.sj,
    fastq = lambda wildcards: dconfig["samples"][wildcards.sample]["fastq"]
  output:
    stats = "%s/aln2.{sample}.Log.final.out" % __STAR_OUTDIR__,
    sam   = "%s/aln2.{sample}.Aligned.out.sam" % __STAR_OUTDIR__
  conda: "%s/env.yaml" % __PC_DIR__
  threads: 5
  params:
    rule_outdir = __STAR_OUTDIR__,
    star_params = dconfig["star_params"]
  shell: """
    STAR --runMode alignReads \
         --runThreadN {threads} \
         {params.star_params} \
         --genomeDir {input.index} \
         --readFilesIn {input.fastq} \
         --outFileNamePrefix {params.rule_outdir}/aln2.{wildcards.sample}. \
         --outSAMtype SAM \
         --sjdbFileChrStartEnd {input.sj}
  """

rule all_star:
  input:
    aln = expand("%s/aln2.{sample}.Aligned.sam" % __STAR_OUTDIR__, sample=dconfig["samples"].keys())

###############################################################################

rule remove_multimapped_sort:
  input:
    sam = lambda wildcards: "%s/aln2.%s.Aligned.out.sam" % (__STAR_OUTDIR__, wildcards.sample)
  output:
    bam = "%s/nomultimap.{sample}.bam" % __ALNFILTER_OUTDIR__
  threads: 5
  conda: "%s/env.yaml" % __PC_DIR__
  shell: """
    cat "{input.sam}" \
     | grep -e "\(^@\)\|\(NH:i:1\)" \
     | samtools sort -O BAM -@ {threads} \
     > "{output.bam}"
  """

rule remove_pcr_duplicates:
  input:
    bam = lambda wildcards: "%s/nomultimap.%s.bam" % (__ALNFILTER_OUTDIR__, wildcards.sample)
  output:
    bam = "%s/nodup.{sample}.bam" % __ALNFILTER_OUTDIR__,
    metric = "%s/nodup.{sample}.metric" % __ALNFILTER_OUTDIR__
  conda: "%s/env.yaml" % __PC_DIR__
  shell: """
    picard MarkDuplicates \
      I="{input.bam}" \
      O="{output.bam}" \
      REMOVE_DUPLICATES=true \
      METRICS_FILE="{output.metric}"
  """

###############################################################################

def remove_mismatched_input_toggle(wildcards):
  if dconfig["remove_pcr_duplicates"]:
    return "%s/nodup.%s.bam" % (__ALNFILTER_OUTDIR__, wildcards.sample),
  else:
    return "%s/nomultimap.%s.bam" % (__ALNFILTER_OUTDIR__, wildcards.sample),
  #fi
#edef

rule remove_mismatched:
  input:
    bam = lambda wildcards: remove_mismatched_input_toggle(wildcards)
  output:
    bam = "%s/nomismatches.{sample}.bam" % __ALNFILTER_OUTDIR__
  conda: "%s/env.yaml" % __PC_DIR__
  params:
    max = dconfig["remove_mismatched_reads_max"] + 1
  shell: """
    bamtools filter -tag "nM:<{params.max}" -in "{input.bam}" -out "{output.bam}"
  """
  

###############################################################################

def htseq_count_input_nodup_toggle(wildcards):
  if dconfig["remove_mismatched_reads"]:
    return "%s/nomismatches.%s.bam" % (__ALNFILTER_OUTDIR__, wildcards.sample)
  elif dconfig["remove_pcr_duplicates"]:
    return "%s/nodup.%s.bam" % (__ALNFILTER_OUTDIR__, wildcards.sample),
  else:
    return "%s/nomultimap.%s.bam" % (__ALNFILTER_OUTDIR__, wildcards.sample),
  #fi
#edef

rule htseq_count:
  input:
    bam = lambda wildcards: htseq_count_input_nodup_toggle(wildcards),
    gff = dconfig["genes"]
  output:
    count = "%s/count.{sample}.tsv" % __QUANT_OUTDIR__
  conda: "%s/env.yaml" % __PC_DIR__
  params:
    stranded = "--stranded=yes" if dconfig["strandSpecific"] == 1 else "--stranded=no",
    tparam = "-t %s"% dconfig["htseq_t"],
    iparam = "-i %s"% dconfig["htseq_i"],
    htseq_params = dconfig["htseq_params"]
  shell: """
    htseq-count {params.stranded} {params.htseq_params} {params.tparam} {params.iparam} -f bam -r pos {input.bam} {input.gff} > {output.count}
  """

rule quantifyTargets:
  input:
    count = expand("%s/count.{sample}.tsv" % __QUANT_OUTDIR__, sample=dconfig["samples"].keys())
  output:
    quant = "%s/quantification.tsv" % __QUANT_OUTDIR__
  run:
    import csv
    data = {}
    for sample in dconfig["samples"].keys():
      with open("%s/count.%s.tsv" % (__QUANT_OUTDIR__, sample), "r") as ifd:
        data[sample] = { target : counts for (target, counts) in csv.reader(ifd, delimiter="\t") }
      #ewith
    #efor

    if len(data) > 0:
      samples = sorted(list(dconfig["samples"].keys()), key=lambda x: (dconfig["samples"][x]["replicate_group"], x))
      with open(output.quant, "w") as ofd:
        ofd.write("target\t%s\n" % '\t'.join(samples))
        for target in sorted(list(data[list(data.keys())[0]].keys())):
          ofd.write("%s\t%s\n" % (target, '\t'.join([(data[sample][target]) for sample in samples])))
        #efor
      #ewith
    #fi

###############################################################################

include: "%s/deseq.Snakefile" % __PC_DIR__

rule deseq:
  input:
    tests = rules.deseqTests.output.tests,
    norm  = rules.deseqNorm.output.norm,
    map   = rules.deSeqMapTargetNames.output

###############################################################################


