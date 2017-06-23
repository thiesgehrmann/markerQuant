#configfile: "config.json"

import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__
__JAR__         = "%s/markerQuant/markerQuant.jar" % __INSTALL_DIR__
__RUN_DIR__ = os.path.abspath(config["outdir"]) + "/run"


__MARKER_OUTDIR__ = "%s/markers" % __RUN_DIR__
rule generateMarkers:
  input:
    targets = config["targets"],
    genomes = config["genomes"]
  output:
    markers = "%s/markers.fasta" % __MARKER_OUTDIR__,
    gaps    = "%s/gaps.fasta"    % __MARKER_OUTDIR__
  params:
    jar = __JAR__,
    strandSpecific = "-s" if config["strandSpecific"] == 1 else ""
  shell: """
    outdir=`dirname {output.markers}`
    mkdir -p "$outdir"
    java -jar {params.jar} markers -g {input.genomes} -t {input.targets} -o "$outdir/" {params.strandSpecific}
  """

__QUANT_OUTDIR__ = "%s/quantification" % __RUN_DIR__
rule quantifyMarkers:
  input:
    markers = rules.generateMarkers.output.markers,
    gaps    = rules.generateMarkers.output.gaps,
    fastq   = lambda wildcards: config["samples"][wildcards.sample]["fastq"]
  output:
    quant = "%s/marker_counts.{sample}.tsv" % __QUANT_OUTDIR__
  params:
    jar = __JAR__,
     fastq = lambda wildcards: ','.join(config["samples"][wildcards.sample]["fastq"]),
    strandSpecific = "-s" if config["strandSpecific"] == 1 else ""
  shell: """
    outdir=`dirname {output.quant}`
    mkdir -p "$outdir"
    java -jar {params.jar} quant -m {input.markers} -g {input.gaps} -f {params.fastq} -o {output.quant} {params.strandSpecific}
  """

rule averageMarkers:
  input:
    quant = lambda wildcards: "%s/marker_counts.%s.tsv" % (__QUANT_OUTDIR__, wildcards.sample)
  output:
    quant = "%s/target_counts.{sample}.tsv" % __QUANT_OUTDIR__
  run:
    import csv
    rawQ = {}
    with open(input.quant, "r") as ifd:
      reader = csv.reader(ifd, delimiter="\t")
      for row in reader:
        if len(row) != 2:
          continue
        #fi
        target   = ':'.join(row[0].split(':')[:-1])
        location = int(row[0].split(':')[-1])
        count    = int(row[1])
        if target not in rawQ:
          rawQ[target] = {}
        #fi
        rawQ[target][location] = count
      #efor
    #ewith

    import statistics
    with open(output.quant, "w") as ofd:
      for target in sorted(rawQ.keys()):
        counts = rawQ[target].values()
        mean  = statistics.mean(counts)
        stdev = statistics.stdev(counts)
        nmarkers = len(counts)
        countString = ','.join([ str(rawQ[target][location]) for location in sorted(rawQ[target].keys())])
        ofd.write("%s\t%f\t%f\t%d\t%s\n" % (target, mean, stdev, nmarkers, countString))
      #efor
    #ewith
     

rule quantifiedMarkers:
  input:
    quant = expand( "%s/target_counts.{sample}.tsv" % (__QUANT_OUTDIR__), sample=config["samples"].keys())
  output:
    quant = "%s/quantification.tsv" % __QUANT_OUTDIR__
  run:
    import csv
    data = {}
    for sample in config["samples"].keys():
      with open("%s/target_counts.%s.tsv" % (__QUANT_OUTDIR__, sample), "r") as ifd:
        data[sample] = { target : (mean, stdev, nmarkers, counts) for (target, mean, stdev, nmarkers, counts) in csv.reader(ifd, delimiter="\t") }
      #ewith
    #efor

    if len(data) > 0:
      samples = sorted(list(config["samples"].keys()), key=lambda x: (config["samples"][x]["replicate_group"], x))
      with open(output.quant, "w") as ofd:
        ofd.write("target\t%s\n" % '\t'.join(samples))
        for target in sorted(list(data[list(data.keys())[0]].keys())):
          ofd.write("%s\t%s\n" % (target, '\t'.join([ data[sample][target][0] for sample in samples])))
        #efor
      #ewith
    #fi

__DIFF_OUTDIR__ = "%s/diffex" % __RUN_DIR__

rule deseq_sample_info_table:
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

rule deseq_test:
  input:
    quant       = rules.quantifiedMarkers.output.quant,
    sample_info = rules.deseq_sample_info_table.output.table
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

rule allDeseq:
  input:
    deseqs = expand("%s/test.{test}.output.tsv" % __DIFF_OUTDIR__, test=[ "%s-%s" % (a,b) for (a,b) in config["tests"] ])


