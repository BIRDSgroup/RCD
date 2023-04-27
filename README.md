# Robust discovery of causal gene networks via measurement error estimation and correction.

This is the official repository of the Manuscript "Robust discovery of causal gene networks via measurement
error estimation and correction" by Rahul Biswas, Brintha V P, Amol Dumrewal and Manikandan Narayanan. 

The code was developed by Rahul Biswas and Manikandan Narayanan (in consultation with the other co-authors of the paper listed above), and these developers are jointly referred to as the "BIRDS Group, IIT Madras" in the preamble of all code files in this RCD project.  


## License preamble 

Copyright 2020 BIRDS Group, IIT Madras

Robust Causal Discovery (RCD) is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

RCD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with BioEmbedS.  If not, see <https://www.gnu.org/licenses/>.

## Installation

See ```requirements.txt``` file for the list of dependencies.

## Usage

### Step I. Measurement error Prediction

Determine the normalized count and sizefactors of the input count matrix (rows correspond to genes and columns to the samples)

```
mrna.count <- readRDS('data/yeast/mrna.count.rds')     #load the count data
mrna.count[, 1:ncol(mrna.count)] <- sapply(mrna.count[, 1:ncol(mrna.count)], as.integer)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = mrna.count, colData = as.data.frame(colnames(mrna.count)), design = ~ 1)
dds <- estimateSizeFactors(dds)
mrna.count.normalized <- counts(dds, normalized=TRUE)
  
saveRDS(sizeFactors(dds), file = 'data/yeast/mrna.count.DESeq.sizeFactors.rds')
saveRDS(mrna.count.normalized, file = "data/yeast/mrna.count.DESeqnormalized.rds")
```

Run the measurement_error_estimation.R script in scripts folder to estimate the variance. Please note that the file paths have to be set correctly.

### Step II. RCD causal tests

The R functions to get the p-values and direction is in the scripts/modelstage2/adj_cit1_4.R file.
```
source("scripts/modelstage2/adj_cit1_4.yeast.R")

pvaladj_AB <- get_adj_cit_pvals(L, Gp, Tp, v_eG, v_eT, bootstrap = 300 ,resampl = 50)
pvaladj_BA <- get_adj_cit_pvals(L, Gp, Tp, v_eG, v_eT, bootstrap = 300 ,resampl = 50)
  
adj_cit_res <- get_cit_direction(pvaladj_AB[1], pvaladj_BA[1], thresh = 0.05)
```

L = vector representing the instrumental variable(s).
&nbsp;

L1, L2 = one hot encoded vector corresponding to each element in L.
&nbsp;

Gp  = continuous vector representing the causal mediating variable G.
&nbsp;

Tp  = continuous vector representing the outcome variable T
&nbsp;

v_eG = variance of G determined in step 1
&nbsp;

v_eT = variance of G determined in step 2
