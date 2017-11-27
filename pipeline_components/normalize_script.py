#!/usr/bin/env python

import utils
import gff3utils as gutils

import os
import sys
import numpy as np
import pandas as pd

###############################################################################

def usage():
  print("Normalize read counts to FPKM scores")
  print("")
  print(" fpkm_norm.py <gffFile> <quantFile> <sampleInfoFile> <fragLengthFile> <attrField> <featureType> <normType> <outFile>")
  print("")
  print("  gffFile: a GFF34 file")
  print("  quantFile: A tab separated file with quantification. Top row takes form: 'target\\tsample1\\tsample2...'")
  print("  sampleInfoFile: A file describing the relationships between samples. Each row has the following fields: 'conditionid\\tsampleid'.") 
  print("  attrField: Which field in the GFF3 file to aggregate lengths.")
  print("  featureType: Which features to aggregate in the GFF3 file.")
  print("  normType: Which normalization to use. (FPKM is only implemented type)")
  print("  outFile: The file to output to")
  print("")

###############################################################################

def fpkm(Q, conditions, samples, SI, sampleSizes, targetLengths, SFL):

  norms = []
  for i, row in Q.iterrows():
    target = row['target']
    fpkms = [ target ]
    for condition in conditions:
      conditionFPKMs = []
      for sample in SI[condition]:
        fpkm = (10**9) * float(row[sample])/( (targetLengths[target] - SFL[sample] + 1) * sampleSizes[sample])
        conditionFPKMs.append(fpkm)
      #efor
      fpkms.append('%f' % np.mean(conditionFPKMs))
    #efor
    norms.append(fpkms)
  #efor

  return norms
#edef

###############################################################################

def tpm(Q, conditions, samples, SI, sampleSizes, targetLengths, SFL):

  sampleFractions = { sample: [] for sample in samples  }
  for sample in samples:
    for i, (_, row) in enumerate(Q.iterrows()):
      sampleFractions[sample].append( row[sample] / (targetLengths[row['target']] - SFL[sample] - 1))
    #efor
  #efor

  sampleSum = { sample : sum(sampleFractions[sample]) for sample in samples }

  norm = []
  for i, (_, row) in enumerate(Q.iterrows()):
    tpms = [ row["target"] ]
    for condition in conditions:
      conditionTPMs = []
      for sample in SI[condition]:
        conditionTPMs.append( (10**6 * sampleFractions[sample][i]) / sampleSum[sample])
      #efor
      tpms.append( '%f' % np.mean(conditionTPMs))
    #efor
    norm.append(tpms)
  #efor

  return norm

#edef

normMethods = { "tmm" : fpkm,
                "fpkm" : fpkm,
                "tpm" : tpm }

###############################################################################

if __name__ == "__main__":

  if len(sys.argv) != 9:
    usage()
    sys.exit(1)
  #fi

  ###############################################################################

  gffFile     = sys.argv[1]
  quantFile   = sys.argv[2]
  sampleInfoFile  = sys.argv[3]
  fragLengthsFile = sys.argv[4]
  attrField   = sys.argv[5]
  featureType = sys.argv[6]
  normType    = sys.argv[7]
  outFile     = sys.argv[8]

  ###############################################################################

  # Prepare the data

  G = gutils.readGFF3File(gffFile)
  GE = G.indexByAttr(attrField)
  SI = { condition : [ x.sample for x in samples ] for (condition, samples) in utils.indexListBy(utils.readColumnFile(sampleInfoFile, "sample condition", skip=1), lambda x: x.condition).items() }  
  SFL = { x.sample : float(x.mfl) for x in utils.readColumnFile(fragLengthsFile, "sample mfl") } 
  
  Q = pd.read_csv(quantFile, "\t")
  Q = Q[np.logical_not(Q['target'].str.startswith('__'))]
  
  conditions = sorted(SI.keys())
  samples = [ s for sg in SI.values() for s in sg ]
  targets = Q['target']
  sampleSizes = { sample : Q[sample].sum() for sample in samples }
  
  targetLengths = { target : sum([ x.end - x.start for x in GE[target] if x.type == featureType]) for target in targets if target in GE  }


  ###############################################################################

  # Calculate the normalization, averaged per condition

  norms = normMethods[normType](Q, conditions, samples, SI, sampleSizes, targetLengths, SFL)

  ###############################################################################

  # Write output

  with open(outFile, "w") as ofd:
    ofd.write("target\t%s\n" % ('\t'.join(conditions)))
    for line in norms:
      ofd.write("%s\n" % '\t'.join(line))
    #efor
  #ewith

#fi
