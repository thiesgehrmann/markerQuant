#!/usr/bin/bash

genes="tdh_genes.fasta"
outdir="simulated_data"

mkdir -p $outdir
Rscript gen_expr_data.R $genes $outdir

find $outdir \
 | grep -e '[.]fasta$' \
 | while read fn; do 
  newname=`echo $fn | sed -e 's/[.]fasta$/.fastq/'`
  cat $fn | ./fasta2fastq.sh > $newname
  rm $fn
done
