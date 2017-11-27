#!/usr/bin/env bash

sam="$1"
sampleSize="$2"
sample="$3"
out="$4"

mfl=`samtools view "$sam" \
      | head -n$sampleSize \
      |  awk '
        function abs(v) {{
          return v < 0 ? -v : v
        }}
        BEGIN{{
          OFS=FS;
          sum=0;
          count=0
        }}
        {{
          if($9 != ""){{
            sum+=abs($9);
            count++
          }}
        }}
        END{{
          print  sum/count
        }}'`
echo -e "$sample\t$mfl" > "$out" 
