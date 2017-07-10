#!/usr/bin/bash

genes="tdh_genes.fasta"
outdir="fastq"


mkdir -p $outdir/unstranded
mkdir -p $outdir/stranded
Rscript gen_expr_data.R $genes $outdir/unstranded "false" 0.005
Rscript gen_expr_data.R $genes $outdir/stranded "true" 0.005
Rscript gen_expr_data.R $genes $outdir/errorfree "false" 0

find $outdir \
 | grep -e '[.]fasta$' \
 | while read fn; do 
  newname=`echo $fn | sed -e 's/[.]fasta$/.fastq/'`
  cat $fn | ./fasta2fastq.sh > $newname
  rm $fn
done
