#!/usr/bin/env python

import json
import sys
import inspect, os

import utils as utils

###############################################################################

__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

if len(sys.argv) < 3:
  print("Error, incomplete arguments to checkInput.py")
  sys.exit(1)
#fi

configFile = sys.argv[1];
action     = sys.argv[2];

config = {}
if not(os.path.isfile(sys.argv[1])):
  errors.append("ConfigFile '%s'doesn't exist!"% sys.argv[1])
else:
  config = json.load(open(sys.argv[1],"r"))
#fi
dconfig = json.load(open("%s/defaults.json" % __INSTALL_DIR__,"r"))

###############################################################################

errors = []
warnings = []

filegroups = [ ("genomes", "genome", True) ]

if action != "align":
  filegroups += [ ("targets", "target", True), ("transcriptome", "transcriptome", False) ]
#fi

for (key, title, required) in filegroups:
  if key not in config:
    if required:
      errors.append("No %s sequences defined." % title)
    #fi
  else:
    files = config[key].split(",")
    if len(files) == 0:
      errors.append("No %s sequences defined." % title)
    else:
      for fileName in files:
        if not(os.path.isfile(fileName)):
          errors.append("File (%s) '%s' does not exist." % (title, fileName))
        #fi
      #efor
    #fi
  #fi
#efor


if action == "align":
  if "genes" not in config:
    errors.append("When aligning, you must define a GFF file with genes under 'genes' in the JSON file.")
  else:
    if not(os.path.isfile(config["genomes"])):
      errors.append("Genes file '%s' doesn't exist." % config["genes"])
    #fi
  #fi
#fi

replicate_groups = {}
if ("samples" not in config) or (len(config["samples"]) == 0):
  errors.append("No samples defined")
else:
  for sample in config["samples"]:
    sampleData = config["samples"][sample]
    if "replicate_group" not in sampleData:
      errors.append("Undefined replicate group for sample '%s'." % sample)
    else:
      replicate_groups[sampleData["replicate_group"]] = replicate_groups.get(sampleData["replicate_group"], 0) + 1
    #fi

    if ("fastq" not in sampleData) or (len(sampleData["fastq"]) == 0):
      errors.append("No FASTQ files defined for sample '%s'." % sample)
    else:
      for fastq in sampleData["fastq"]:
        if not(os.path.isfile(fastq)):
          errors.append("FASTQ file '%s' for sample '%s' does not exist." % (fastq, sample))
        #fi
      #efor
    #fi
  #efor
#fi

if len(replicate_groups) < 2:
  warnings.append("We can't normalize properly with only one condition. We will calculate a mean of the FPKMs for each of the replicates.")
#fi 

if "tests" not in config:
  warnings.append("No diff. ex. tests defined.")
else:
  for test in config["tests"]:
    if len(test) != 2:
      errors.append("Test '%s' is invalid, they must always contain two conditions." % ','.join(test))
    elif test[0] == test[1]:
      errors.append("Test '%s' is invalid, the conditions must be different." % ','.join(test))
    #fi

    for condition in test:
      if condition not in replicate_groups:
        errors.append("Replicate group (condition) '%s' is unknown." % condition)
      elif replicate_groups[condition] < 2:
        warnings.append("There is only one replicate for replicate group '%s'. Diff. Ex. might not work properly." % condition)
      #fi
    #efor
  #efor
#fi
        

if "targetsMap" in config:
  if not(os.path.isfile(config["targetsMap"])):
    errors.append("TargetsMap file '%s' does not exist." % config["targetsMap"])
  #fi
  # CHECK FORMAT OF TARGETSMAP
#fi

if "knockouts" in config:
  if not(os.path.isfile(config["knockouts"])):
    errors.append("Knockouts file '%s' does not exist." % config["knockouts"])
  #fi
  # CHECK FORMAT OF TARGETSMAP
#fi

    
###############################################################################
    
for error in errors:
  print("ERROR: %s" % error)
#efor

for warning in warnings:
  print("WARNING: %s" % warning)
#efor

if len(errors) > 0:
  sys.exit(1)
#fi

sys.exit(0)
