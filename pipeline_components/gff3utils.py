import errno
import os
import csv
import gzip
from collections import namedtuple
import gzip

import utils

###############################################################################

gff3Entry = namedtuple("gff3Entry", "seqid, source, type, start, end, score, strand, phase, attr")

def parseGFF3entry(fields):
   (seqid, source, feature, start, end, score, strand, phase, attr) = fields
   attr = dict( tuple(x) if len(x) == 2 else (x[0], '='.join(x[1:])) for x in [ y.strip().split('=') for y in attr.split(";") ] )
   return gff3Entry(seqid, source, feature, int(start), int(end), score, strand, phase, attr)
#edef

###############################################################################

class GFF3(object):
  def __init__(self, entries):
    self.entries = entries
    self.seqids  = set([ e.seqid for e in self.entries])
  #edef

  def indexByAttr(self, attr="Parent"):
    return utils.indexListBy(self.entries, lambda x: x.attr.get(attr, ""))
  #edef

#eclass

###############################################################################

def readGFF3File(filename):
  G = []
  with (gzip.open(filename, "rt") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
    gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
    for row in gffReader:
      if len(row) != 9:
        continue
      #fi
      G.append(parseGFF3entry(row))
    #efor
  #ewith
  return GFF3(G)
#edef

###############################################################################
