#!/bin/bash

wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 \
  https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
  https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2

tools=("samtools-1.19.2" "bcftools-1.19" "htslib-1.19.1")
for tool in "${tools[@]}"; do tar xjf "${tool}".tar.bz2 ; done
for tool in "${tools[@]}"; do cd "${tool}" && ./configure && make && make install && cd .. && rm -rf "${tool}" ; done