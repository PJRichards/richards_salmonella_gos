#!/bin/bash

for i in $(seq 756 1 871); do
    if [[ ! -f "resources/reads/SRR24547${i}_1.fastq.gz" || ! -f "resources/reads/SRR24547${i}_2.fastq.gz" ]]
    then
    fastq-dump --split-files SRR24547$i -O resources/reads
    gzip resources/reads/SRR24547${i}_1.fastq
    gzip resources/reads/SRR24547${i}_2.fastq
    fi
done

