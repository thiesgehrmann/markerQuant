configfile: "config.json"

__INSTALL_DIR__ = "./"
__RUN_DIR__ = config["outdir"] + "run"


__MARKER_OUTDIR__ = "%s/markers" % __RUN_DIR__
rule generateMarkers:
  input:
    targets = config.targets,
    genomes = config.genomes
  output:
    markers = "%s/markers.fa" % __MARKER_OUTDIR__
    gaps    = "%s/gaps.fa"    % __MARKER_OUTDIR__
  params:
    install_dir = __INSTALL_DIR__
  shell: """
    outdir=`dirname {output.markers}`
    mkdir -p "$outdir"
    java -jar {params.install_dir}/markerQuant.jar markers -g {input.genomes} -t {input.targets} -o "$outdir"
  """

__QUANT_OUTDIR__ = "%s/quantification" % __RUN_DIR__
rule quantifyMarkers:
  input:
    markers = rules.generateMarkers.output.markers,
    gaps    = rules.generateMarkers.output.gaps,
    fastq   = lambda wildcards: config["fastq"][wildcards.sample]
  output:
    quant = "%s/sample.%s.tsv" % __QUANT_OUTDIR__
  shell: """
    outdir=`dirname {output.quant}`
    mkdir -p "$outdir"
    java -jar {params.install_dir}/markerQuant.jar quant -m {input.markers} -g {input.gaps} -f {input.fastq} -o {output.quant}
  """
