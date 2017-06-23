#!/bin/bash

package="markerQuant"

scala_classpath="./:./lib:./lib/scala-parser-combinators_2.12-1.0.4.jar:./lib/atk-90.jar"

scalac="`which scalac` -classpath $scala_classpath -deprecation -feature"
scala="`which scala` -classpath $scala_classpath"

echo "Building $package..."
mkdir -p $package
$scalac -d ./ `grep -i "package[ ]\+$package" *.scala | cut -d: -f1`

echo "Constructing jar..."
cat > $package.mf << EOF
Main-Class: $package.Main
Class-Path: $( echo $scala_classpath | cut --complement -d: -f1 | tr ':' '\n' | sed -e 's/^.*$/ & /')
 lib/scala-library.jar
EOF

jar -cmf $package.mf $package.jar $(find $package | grep -e '.*class$')
