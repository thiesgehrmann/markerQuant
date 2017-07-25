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
      for row in reader:
        if row[0] in mapN:
          ofd.write("%s\t%s\t%s\n" % (row[0], mapN[row[0]], "\t".join(row[1:])))
        else:
          ofd.write("%s\t-\t%s\n" % (row[0], '\t'.join(row[1:])))
        #fi
      #efor
    #ewith
  #ewith
#edef
