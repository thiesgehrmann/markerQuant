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
__DIFF_OUTDIR__ = "%s/deseq_outdir" % __RUN_DIR__

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
    stranded = "--stranded=yes" if config["strandSpecific"] == 1 else "--stranded=no",
    tparam = "-t %s"% config["htseq_t"] if "htseq_t" in config else "-t gene",
    iparam = "-i %s"% config["htseq_i"] if "htseq_i" in config else "-i gene"
  shell: """
    htseq-count {params.stranded} {params.tparam} {params.iparam} -f bam -r pos {input.bam} {input.gff} > {output.count}
  """

rule quantifyTargets:
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

###############################################################################

rule mapTargetNames:
  input:
    targetMap = config["targetMap"] if "targetMap" in config else "",
    quant = rules.quantifyTargets.output.quant
  output:
    quant = "%s/quantification.map.tsv" % __HTSEQ_OUTDIR__
  run:
    import csv
    mapN = {}
    with open(input.targetMap, "r") as ifd:
      reader = csv.reader(ifd, delimiter="\t")
      for row in reader:
        if len(row) < 2:
          continue
        #fi
        mapN[row[0]] = row[1]
      #efor
    #ewith

    with open(input.quant, "r") as ifd:
      with open(output.quant, "w") as ofd:
        reader = csv.reader(ifd, delimiter="\t")
        for row in reader:
          if row[0] in mapN:
            ofd.write("%s\t%s\n" % (mapN[row[0]], "\t".join(row[1:])))
          else:
            ofd.write("%s\n" % '\t'.join(row))
          #fi
        #efor
      #ewith
    #ewith

###############################################################################

rule deseqSampleInfoTable:
  output:
    table = "%s/sample_info.tsv"% __DIFF_OUTDIR__
  run:
    samples = sorted(list(config["samples"].keys()), key=lambda x: (config["samples"][x]["replicate_group"], x))
    with open(output.table, "w") as ofd:
      ofd.write("sample\tcondition\n")
      for sample in samples:
        ofd.write("%s\t%s\n" % (sample, config["samples"][sample]["replicate_group"]))
      #efor
    #ewith

def deseqTestQuantInput():
  if "targetMap" in config:
    return rules.mapTargetNames.output.quant
  else:
    return rules.quantifyTargets.output.quant
  #fi
#edef

rule deseqTest:
  input:
    quant       = deseqTestQuantInput(),
    sample_info = rules.deseqSampleInfoTable.output.table
  output:
    diff = "%s/test.{test}.output.tsv" %__DIFF_OUTDIR__
  conda : "%s/pipeline_components/env.yaml"% __INSTALL_DIR__
  params:
    sample1 = lambda wildcards: wildcards.test.split('-')[0],
    sample2 = lambda wildcards: wildcards.test.split('-')[1],
    deseq_wrapper = "%s/deseq_wrapper.R" % __PC_DIR__
  shell: """
    Rscript {params.deseq_wrapper} '{input.quant}' {input.sample_info} {params.sample1} {params.sample2} {output.diff}
  """

rule deseqTests:
  input:
    tests = expand("%s/test.{test}.output.tsv" % __DIFF_OUTDIR__, test=[ "%s-%s" % (a,b) for (a,b) in config["tests"] ])
  output:
    tests = "%s/tests.tsv"% __DIFF_OUTDIR__
  shell: """
    cat {input.tests} > {output.tests}
  """

rule deseqNorm:
  input:
    quant       = rules.quantifyTargets.output.quant,
    sample_info = rules.deseqSampleInfoTable.output.table
  output:
    norm = "%s/quantification.normalized.tsv" %__DIFF_OUTDIR__
  conda : "%s/pipeline_components/env.yaml"% __INSTALL_DIR__
  params:
    deseq_wrapper = "%s/deseq_norm.R" % __PC_DIR__
  shell: """
    Rscript {params.deseq_wrapper} '{input.quant}' {input.sample_info} {output.norm}
  """

rule deseq:
  input:
    tests = rules.deseqTests.output.tests,
    norm  = rules.deseqNorm.output.norm

###############################################################################


