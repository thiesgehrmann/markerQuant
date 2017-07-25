
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
  params:
    deseq_adjust = "%s/deseq_adjust.R" % __PC_DIR__
  conda : "%s/pipeline_components/env.yaml"% __INSTALL_DIR__
  shell: """
    echo -e '"target"\t"baseMean"\t"log2FoldChange"\t"lfcSE"\t"stat"\t"pvalue"\t"padj"\t"condition"' > {output.tests}.combined
    cat {input.tests} \
     | grep -v '^"baseMean"' \
     >> {output.tests}.combined
    Rscript {params.deseq_adjust} {output.tests}.combined {output.tests}
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


rule deSeqMapTargetNames:
  input:
    targetMap = dconfig["targetMap"],
    quant = rules.deseqNorm.output.norm,
    tests = rules.deseqTests.output.tests
  output:
    quant = "%s/quantification.normalized.map.tsv" % __DIFF_OUTDIR__,
    tests = "%s/tests.map.tsv" % __DIFF_OUTDIR__,
  run:
    import utils
    utils.mapTargetNames(input.targetMap, input.quant, output.quant)
    utils.mapTargetNames(input.targetMap, input.tests, output.tests)

