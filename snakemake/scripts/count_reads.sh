#!/usr/bin/env bash

# Count reads in a sam file
samtools view -F 260 ${snakemake_input[0]} | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > ${snakemake_output[0]}
echo "Done counting reads"

