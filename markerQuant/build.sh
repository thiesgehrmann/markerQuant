#!/bin/bash

scala_classpath="./:./lib:./lib/scala-parser-combinators_2.12-1.0.4.jar"

scalac="`which scalac` -classpath $scala_classpath -deprecation -feature"
scala="`which scala` -classpath $scala_classpath"

echo "Building markerQuant..."
mkdir -p panalysis
#$scalac -d ./ `grep 'package[ ]\+markerQuant' *.scala | cut -d: -f1`

echo "Constructing jar..."
cat > markerQuant.mf << EOF
Main-Class: markerQuant.Main
Class-Path: $( echo $scala_classpath | cut --complement -d: -f1 | tr ':' '\n' | sed -e 's/^.*$/ & /')
 lib/scala-library.jar
EOF

jar -cmf markerQuant.mf markerQuant.jar $(find markerQuant | grep -e '.*class$')
