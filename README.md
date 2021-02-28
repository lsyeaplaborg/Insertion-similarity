# Insertion-similarity
***

## Introduction

### This pipeline can be used to calculate similarity score between insertions and adjacent genomic segments. Insertions were then classfied based on insertion length and similarity score.

> 1: >3bp, same ratio>=60%

> 2: 1bp, same ratio=100%

> 3: 3bp, same ratio>=60%

> 11: 2bp, same ratio=100%

> 7: >=3bp, same ratio<60%

> 12: 2bp, same ratio<100%

> 13: 1bp, same ratio=0

## Getting Started

### Specify the input, output, scripts and reference path in the `submit.count.sh` (mouse) and`submit.cate.count.cell_line.sh` (cell line) and `submit.cate.count.no_12d.sh` file

### Run pipeline with command line:

`bash submit.count.sh` for mouse data

`bash submit.cate.count.cell_line.sh` for cell line data

`bash submit.cate.count.no_12d.sh` for no_12d data
