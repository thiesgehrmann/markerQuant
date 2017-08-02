#!/usr/bin/env bash

genes="tdh_genes.fasta"
outdir="fastq"

./generate_expr_data.sh $genes $outdir/unstranded "false" 0.005
./generate_expr_data.sh $genes $outdir/stranded "true" 0.005
./generate_expr_data.sh $genes $outdir/errorfree "false" 0


