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
    strandSpecific = "-s" if dconfig["strandSpecific"] else "",
    transcriptome = "-T %s" % dconfig["transcriptome"] if ("transcriptome" in config) else "",
    k = dconfig["k"]
  conda: "%s/env.yaml" % __PC_DIR__
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
    logfile     = "%s/quantification.{sample}.log"% __QUANT_OUTDIR__,
    kohits      = "%s/ko_hits.{sample}.fastq"% __QUANT_OUTDIR__
  params:
    jar = __JAR__,
    fastq = lambda wildcards: ','.join(dconfig["samples"][wildcards.sample]["fastq"]),
    strandSpecific = "-s" if dconfig["strandSpecific"] else "",
    minQual = dconfig["minQual"],
    k = dconfig["k"],
    knockouts = "-K %s" % dconfig["knockouts"] if dconfig["knockouts"] != "" else ""
  conda: "%s/env.yaml" % __PC_DIR__
  shell: """
    java -Xmx100G -jar {params.jar} quant -m {input.markers} -g {input.gaps} -f {params.fastq} -k {params.k} -M {output.markerQuant} -T {output.targetQuant} {params.knockouts} {params.strandSpecific} -q {params.minQual} 2>&1 \
     | tee {output.logfile} \
     | grep "^#Found Knockout:" -A 4 \
     > "{output.kohits}"
  """

rule quantifyAllMarkers:
  input:
    markerQuants = expand("%s/marker_counts.{sample}.tsv" % __QUANT_OUTDIR__, sample=dconfig["samples"].keys()),
    targetQuants = expand("%s/target_counts.{sample}.tsv" % __QUANT_OUTDIR__, sample=dconfig["samples"].keys()),
    logfiles     = expand("%s/quantification.{sample}.log"% __QUANT_OUTDIR__, sample=dconfig["samples"].keys()),
    kohits       = expand("%s/ko_hits.{sample}.fastq"% __QUANT_OUTDIR__, sample=dconfig["samples"].keys())


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


rule fragLengthMeans:
  output:
    fragLengthMeans = "%s/fragLengthMean.tsv" % __QUANT_OUTDIR__
  params:
    samples = dconfig["samples"].keys(),
    k       = dconfig["k"]
  run:
    with open(outut.fragLengthMeans, "w") as ofd:
      for sample in samples:
        ofd.write("%s\t%d\n" % (sample, params.k))
      #efor
    #ewith

###############################################################################

include: "%s/deseq.Snakefile" % __PC_DIR__

rule deseq:
  input:
    tests = rules.deseqTests.output.tests if len(dconfig["tests"]) > 0 else [],
    norm  = rules.normalize.output.norm,
    map   = rules.deSeqMapTargetNames.output if dconfig["targetMap"] != "" else []

###############################################################################


