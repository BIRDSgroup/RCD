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

- The R functions to get the p-values and direction is in the scripts/modelstage2/adj_cit1_4.R file.
```
source("scripts/modelstage2/adj_cit1_4.yeast.R")

pvaladj_AB <- get_adj_cit_pvals(L, Gp, Tp, v_eG, v_eT, bootstrap = 300 ,resampl = 50)
pvaladj_BA <- get_adj_cit_pvals(L, Gp, Tp, v_eG, v_eT, bootstrap = 300 ,resampl = 50)
  
adj_cit_res <- get_cit_direction(pvaladj_AB[1], pvaladj_BA[1], thresh = 0.05)
```

L = vector representing the instrumental variable(s).

L1, L2 = one hot encoded vector corresponding each element in L.
&nbsp;

Gp  = continuous vector representing the causal mediating variables.
&nbsp;
Tp  = continuous vector representing the outcome variables

