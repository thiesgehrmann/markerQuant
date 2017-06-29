import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__
__JAR__         = "%s/markerQuant/markerQuant.jar" % __INSTALL_DIR__
__RUN_DIR__ = os.path.abspath(config["outdir"]) + "/align_run"

tconf = {
  "star_params" : ""
}

__STAR_OUTDIR__ = "%s/star_align" % __RUN_DIR__
__HTSEQ_OUTDIR__ = "%s/htseq_count" % __RUN_DIR__

###############################################################################

rule all:
  input:
    aln = expand("%s/aln2.{sample}.Aligned.sortedByCoord.out.bam" % __STAR_OUTDIR__, sample=config["samples"].keys()),
    count = expand("%s/count.{sample}.tsv" % __HTSEQ_OUTDIR__, sample=config["samples"].keys())


###############################################################################

rule star_index:
  input:
    genomes = config["genomes"]
  output:
    index = "%s/star_.idx" % __STAR_OUTDIR__
  conda: "%s/align_env.yaml" % __PC_DIR__
  shell: """
    mkdir {output.index}
    STAR --runMode genomeGenerate --genomeDir {output.index} --genomeFastaFiles {input.genomes}
  """

rule star_pass_one_sample:
  input:
    fastq = lambda wildcards: config["samples"][wildcards.sample]["fastq"],
    index = rules.star_index.output.index
  output:
    sj = "%s/aln1.{sample}.SJ.out.tab" % __STAR_OUTDIR__
  conda: "%s/align_env.yaml" % __PC_DIR__
  threads: 5
  params:
    rule_outdir = __STAR_OUTDIR__,
    star_params = tconf["star_params"]
  shell: """
    STAR --runMode alignReads --runThreadN {threads} {params.star_params} --genomeDir {input.index} --readFilesIn {input.fastq} --outFileNamePrefix {params.rule_outdir}/aln1.{wildcards.sample}. --outSAMmode None
  """
    

rule star_pass_one:
  input:
    sj = expand("%s/aln1.{sample}.SJ.out.tab" % __STAR_OUTDIR__, sample=config["samples"].keys())
    
    
rule star_pass_two_sample:
  input:
    index = rules.star_index.output.index,
    sj    = rules.star_pass_one.input.sj,
    fastq = lambda wildcards: config["samples"][wildcards.sample]["fastq"]
  output:
    stats = "%s/aln2.{sample}.Log.final.out" % __STAR_OUTDIR__,
    bam   = "%s/aln2.{sample}.Aligned.sortedByCoord.out.bam" % __STAR_OUTDIR__
  conda: "%s/align_env.yaml" % __PC_DIR__
  threads: 5
  params:
    rule_outdir = __STAR_OUTDIR__,
    star_params = tconf["star_params"]
  shell: """
    STAR --runMode alignReads \
         --runThreadN {threads} \
         {params.star_params} \
         --genomeDir {input.index} \
         --readFilesIn {input.fastq} \
         --outFileNamePrefix {params.rule_outdir}/aln2.{wildcards.sample}. \
         --outSAMtype BAM SortedByCoordinate \
         --sjdbFileChrStartEnd {input.sj}
  """

rule all_star:
  input:
    aln = expand("%s/aln2.{sample}.Aligned.sortedByCoord.out.bam" % __STAR_OUTDIR__, sample=config["samples"].keys())

###############################################################################

rule htseq_count:
  input:
    bam = lambda wildcards: "%s/aln2.%s.Aligned.sortedByCoord.out.bam" % (__STAR_OUTDIR__, wildcards.sample),
    gff = config["genes"]
  output:
    count = "%s/count.{sample}.tsv" % __HTSEQ_OUTDIR__
  conda: "%s/align_env.yaml" % __PC_DIR__
  params:
    stranded = "--stranded=yes" if config["stranded"] == 1 else "--stranded=no"
  shell: """
    htseq-count {params.stranded} -f bam -r pos -t gene -i gene {input.bam} {input.gff} > {output.count}
  """

rule quantification:
  input:
    count = expand("%s/count.{sample}.tsv" % __HTSEQ_OUTDIR__, sample=config["samples"].keys())
  output:
    quant = "%s/quantification.tsv" % __HTSEQ_OUTDIR__
  run:
    import csv
    data = {}
    for sample in config["samples"].keys():
      with open("%s/count.%s.tsv" % (__HTSEQ_OUTDIR__, sample), "r") as ifd:
        data[sample] = { target : counts for (target, counts) in csv.reader(ifd, delimiter="\t") }
      #ewith
    #efor

    if len(data) > 0:
      samples = sorted(list(config["samples"].keys()), key=lambda x: (config["samples"][x]["replicate_group"], x))
      with open(output.quant, "w") as ofd:
        ofd.write("target\t%s\n" % '\t'.join(samples))
        for target in sorted(list(data[list(data.keys())[0]].keys())):
          ofd.write("%s\t%s\n" % (target, '\t'.join([(data[sample][target]) for sample in samples])))
        #efor
      #ewith
    #fi
