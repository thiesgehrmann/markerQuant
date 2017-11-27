import csv
from collections import namedtuple

###############################################################################

def mapTargetNames(targetMapFile, inputFile, outputFile):
  import csv
  mapN = {}
  with open(targetMapFile, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t")
    for row in reader:
      if len(row) < 2:
        continue
      #fi
      mapN[row[0]] = row[1]
    #efor
  #ewith

  with open(inputFile, "r") as ifd:
    with open(outputFile, "w") as ofd:
      reader = csv.reader(ifd, delimiter="\t")
      for i, row in enumerate(reader):
        if i == 0:
          ofd.write("%s\tmapID\t%s\n" % (row[0], '\t'.join(row[1:])))
          continue
        #fi

        if row[0] in mapN:
          ofd.write("%s\t%s\t%s\n" % (row[0], mapN[row[0]], "\t".join(row[1:])))
        else:
          ofd.write("%s\t-\t%s\n" % (row[0], '\t'.join(row[1:])))
        #fi
      #efor
    #ewith
  #ewith
#edef

###############################################################################

def readColumnFile(filename, columnNames, delimiter='\t', types="", skip=0):
  import csv
  L = []
  typeFunctions = { "str" : lambda x: str(x),
                    "int" : lambda x: int(x),
                    "float" : lambda x: float(x) }

  if types != "":
    types = [ typeFunctions[c] for c in types.split(" ") ]
  #fi

  nColumns = len(columnNames.split(" "))

  lineType = namedtuple("lineType", columnNames)
  skipped = 0
  with open(filename, "r") as ifd:
    reader = csv.reader(ifd, delimiter=delimiter)
    for row in reader:
      if (row[0][0] == '#') or (skipped < skip):
        skipped += 1
        continue
      #fi
      if len(types) == len(row):
        row = [ tf(v) for (tf, v) in zip(types, row) ]
      #fi
      rowLen = len(row)
      if rowLen < nColumns:
        row = [ row[i] if i < rowLen else ""  for i in range(nColumns) ]
      #fi
      L.append(lineType(*row))
    #efor
  #ewith
  return L
#edef

###############################################################################

def readCSV(fileName, delimiter=',', quotechar='"', skiprows=0):
  D = []
  with open(fileName, 'r', newline='') as ifd:
    reader = csv.reader(ifd, delimiter=delimiter, quotechar=quotechar)
    for row in reader:
      D.append(row)
    #efor
  #ewith
  return D
#edef

###############################################################################

def indexListBy(L, key=lambda x: x[0]):
  G = {}
  for item in L:
    k = key(item)
    if k not in G:
      G[k] = []
    #fi
    G[k].append(item)
  #efor
  return G
#edef
