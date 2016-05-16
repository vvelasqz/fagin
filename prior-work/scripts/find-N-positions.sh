#!/bin/bash
set -u

target="$1"
output=$(sed -r 's/.*\/([^.]+).*/\1.N-string.tab/' <<< "$target")
smof grep -qP '(N+)' --gff < $target |
    awk 'BEGIN{OFS="\t"} {print $1, $4, $5 - $4 + 1}' > $output
