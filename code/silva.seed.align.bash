#!/usr/bin/env bash

# Download the SILVA reference file (v.128). We will pull out the bacteria-specific sequences and
# clean up the directories to remove the extra files
wget https://mothur.org/w/images/a/a4/Silva.seed_v128.tgz
tar xvzf Silva.seed_v128.tgz silva.seed_v128.align silva.seed_v128.tax
code/mothur/mothur "#get.lineage(fasta=silva.seed_v128.align, taxonomy=silva.seed_v128.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v128.pick.align)"
mv silva.seed_v128.pick.align data/references/silva.seed.align
rm Silva.seed_v128.tgz silva.seed_v128.*
rm mothur.*.logfile

# Download the RDP taxonomy references (v14), put the necessary files in data/references, and
# clean up the directories to remove the extra files
wget -N https://mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz
tar xvzf Trainset16_022016.pds.tgz
mv trainset16_022016.pds/train* data/references/
rm -rf trainset16_022016.pds
rm Trainset16_022016.pds.tgz

# Generate a customized version of the SILVA reference database that targets the V4 region
code/mothur/mothur "#pcr.seqs(fasta=data/references/silva.seed.align, start=11894, end=25319, keepdots=F, processors=16)"
mv data/references/silva.seed.pcr.align data/references/silva.v4.align