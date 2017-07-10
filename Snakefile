#configfile: "config.json"

import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__
__JAR__         = "%s/markerQuant/markerQuant.jar" % __INSTALL_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)

###############################################################################

__RUN_DIR__ = os.path.abspath(dconfig["outdir"]) + "/run"
__MARKER_OUTDIR__ = "%s/markers" % __RUN_DIR__
__QUANT_OUTDIR__  = "%s/quantification" % __RUN_DIR__
__DIFF_OUTDIR__   = "%s/diffex" % __RUN_DIR__

###############################################################################

rule all:
  input:
    quant = "%s/quantification.tsv" % __QUANT_OUTDIR__,
    stats = "%s/targetCountStatistics.tsv" % __QUANT_OUTDIR__,
    norm  = "%s/quantification.normalized.tsv"% __DIFF_OUTDIR__,
    tests = "%s/tests.tsv"% __DIFF_OUTDIR__,

###############################################################################

rule generateMarkers:
  input:
    targets = dconfig["targets"],
    genomes = dconfig["genomes"]
  output:
    markers = "%s/markers.fasta" % __MARKER_OUTDIR__,
    gaps    = "%s/gaps.fasta"    % __MARKER_OUTDIR__
  params:
    jar = __JAR__,
    strandSpecific = "-s" if dconfig["strandSpecific"] == 1 else "",
    transcriptome = "-T %s" % dconfig["transcriptome"] if ("transcriptome" in config) else "",
    k = dconfig["k"]
  shell: """
    java -jar {params.jar} markers -g {input.genomes} -t {input.targets} -k {params.k} {params.strandSpecific} {params.transcriptome} -M {output.markers} -G {output.gaps}
  """

###############################################################################

rule quantifyMarkers:
  input:
    markers = rules.generateMarkers.output.markers,
    gaps    = rules.generateMarkers.output.gaps,
    fastq   = lambda wildcards: dconfig["samples"][wildcards.sample]["fastq"]
  output:
    markerQuant = "%s/marker_counts.{sample}.tsv" % __QUANT_OUTDIR__,
    targetQuant = "%s/target_counts.{sample}.tsv" % __QUANT_OUTDIR__,
    logfile     = "%s/quantification.{sample}.log"% __QUANT_OUTDIR__
  params:
    jar = __JAR__,
    fastq = lambda wildcards: ','.join(dconfig["samples"][wildcards.sample]["fastq"]),
    strandSpecific = "-s" if dconfig["strandSpecific"] == 1 else "",
    minQual = dconfig["minQual"],
    k = dconfig["k"],
    knockouts = "-K %s" % dconfig["knockouts"] if dconfig["knockouts"] != "" else ""
  shell: """
    java -Xmx100G -jar {params.jar} quant -m {input.markers} -g {input.gaps} -f {params.fastq} -k {params.k} -M {output.markerQuant} -T {output.targetQuant} {params.knockouts} {params.strandSpecific} -q {params.minQual} 2>&1 | tee {output.logfile}
  """

rule MarkerStatisticsPerTargetPerSample:
  input:
    quant = lambda wildcards: "%s/marker_counts.%s.tsv" % (__QUANT_OUTDIR__, wildcards.sample)
  output:
    quant = "%s/targetCountStatistics.{sample}.tsv" % __QUANT_OUTDIR__
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
        location = int(row[0].split(':')[-1]) if len(row[0].split(':')) > 1 else 0
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
        stdev = statistics.stdev(counts) if len(counts) > 1 else 0
        nmarkers = len(counts)
        countString = ','.join([ str(rawQ[target][location]) for location in sorted(rawQ[target].keys())])
        ofd.write("%s\t%f\t%f\t%d\t%s\n" % (target, mean, stdev, nmarkers, countString))
      #efor
    #ewith
     

rule MarkerStatisticsPerTarget:
  input:
    quant = expand("%s/targetCountStatistics.{sample}.tsv" % (__QUANT_OUTDIR__), sample=dconfig["samples"].keys())
  output:
    quant = "%s/targetCountStatistics.tsv" % __QUANT_OUTDIR__
  run:
    import csv
    data = {}
    for sample in dconfig["samples"].keys():
      with open("%s/targetCountStatistics.%s.tsv" % (__QUANT_OUTDIR__, sample), "r") as ifd:
        data[sample] = { target : {"mean": mean, "stdev": stdev, "nmarkers": nmarkers, "counts":counts} for (target, mean, stdev, nmarkers, counts) in csv.reader(ifd, delimiter="\t") }
      #ewith
    #efor

    if len(data) > 0:
      samples = sorted(list(dconfig["samples"].keys()), key=lambda x: (dconfig["samples"][x]["replicate_group"], x))
      with open(output.quant, "w") as ofd:
        ofd.write("#target\tnmarkers\tmeans\tstdevs\trawcounts\n")
        ofd.write("#Samples: %s\n" % ' '.join(samples))
        for target in sorted(list(data[list(data.keys())[0]].keys())):
          ofd.write("%s\t%s\t%s\t%s\t%s\n" % (target,
                                    data[samples[0]][target]["nmarkers"], 
                                    ','.join([data[s][target]["mean"] for s in samples]),
                                    ','.join([data[s][target]["stdev"] for s in samples]),
                                    ';'.join([data[s][target]["counts"] for s in samples])
                                    ))
        #efor
      #ewith
    #fi

rule quantifyTargets:
  input:
    quant = expand( "%s/target_counts.{sample}.tsv" % (__QUANT_OUTDIR__), sample=dconfig["samples"].keys())
  output:
    quant = "%s/quantification.tsv" % __QUANT_OUTDIR__
  run:
    import csv
    data = {}
    for sample in dconfig["samples"].keys():
      with open("%s/target_counts.%s.tsv" % (__QUANT_OUTDIR__, sample), "r") as ifd:
        data[sample] = { target : counts for (target, counts) in csv.reader(ifd, delimiter="\t") }
      #ewith
    #efor

    if len(data) > 0:
      samples = sorted(list(dconfig["samples"].keys()), key=lambda x: (dconfig["samples"][x]["replicate_group"], x))
      with open(output.quant, "w") as ofd:
        ofd.write("target\t%s\n" % '\t'.join(samples))
        for target in sorted(list(data[list(data.keys())[0]].keys())):
          ofd.write("%s\t%s\n" % (target, '\t'.join([ str(int(float(data[sample][target]))) for sample in samples])))
        #efor
      #ewith
    #fi

###############################################################################

rule deseqSampleInfoTable:
  output:
    table = "%s/sample_info.tsv"% __DIFF_OUTDIR__
  run:
    samples = sorted(list(dconfig["samples"].keys()), key=lambda x: (dconfig["samples"][x]["replicate_group"], x))
    with open(output.table, "w") as ofd:
      ofd.write("sample\tcondition\n")
      for sample in samples:
        ofd.write("%s\t%s\n" % (sample, dconfig["samples"][sample]["replicate_group"]))
      #efor
    #ewith

rule deseqTest:
  input:
    quant       = rules.quantifyTargets.output.quant,
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
    tests = expand("%s/test.{test}.output.tsv" % __DIFF_OUTDIR__, test=[ "%s-%s" % (a,b) for (a,b) in dconfig["tests"] ])
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


rule mapTargetNames:
  input:
    targetMap = dconfig["targetMap"],
    quant = rules.deseqNorm.output.norm
  output:
    quant = "%s/quantification.normalized.map.tsv" % __DIFF_OUTDIR__
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

rule deseq:
  input:
    tests = rules.deseqTests.output.tests,
    norm  = rules.deseqNorm.output.norm

###############################################################################


