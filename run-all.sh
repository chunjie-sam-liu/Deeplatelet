#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-29 22:10:52
# @DESCRIPTION:

# Number of input parameters
root_dir=/home/liucj/github/TEP-prognosis
data_dir=/home/liucj/github/TEP-prognosis/data

cd ${root_dir}

rm -rf logs/0*
rm data
ln -s ${data_dir} data

Rscript ${root_dir}/analysis/01-counts2se.R 1> logs/01-counts2se.R.log 2> logs/01-counts2se.R.err