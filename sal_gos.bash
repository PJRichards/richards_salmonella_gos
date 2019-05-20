#!/usr/bin/env bash

# download and unpack mother v1.41.1
bash code/get_mothur_linux_64_1_41_1.bash

# get references
# generate a customized version of the SILVA reference database that targets the V4 region
# using Silva.seed_v128 and Trainset16_022016
bash code/silva.seed.align.bash

# run mothur through quality control steps
code/mothur/mothur code/get_good_seqs.batch

# get zymo mock seqs
bash code/get_mock_community_data.bash