# confidence-multiple-alternatives

This repository accompanies the following paper:

Li, H.H. and Ma, WJ (2020) Confidence reports in decision-making with multiple alternatives violate the Bayesian confidence hypothesis. Nature Communications.
DOI: 10.1038/s41467-020-15581-6

The repository contains data and MATLAB codes for fitting three main models to the data.

## Data
The data of each participant are saved in .mat files in the 'data' folder.

Data are saved in a structure variable named 'trl' that contains trial-by-trial information.

trl.config: the configuration (stimulus distribution)

trl.estC: confidence reports 

trl.estD : category decisions

trl.target_x: horizontal locations of the target dot

trl.target_y: vertical locations of the target dot

trl.theta: the color (angle in unit of radian on a color wheel) of the leftmost category.

## Model fitting
Execute wrapper_fitmodel.m to run the optimization. In this script, one can determine which model, experiment and subject to fit. 
wrapper_fitmodel.m calls fitmodel_on_individual.m to set up and start the optimization. 

In fitmodel_on_individual.m, codes for optimization using either BADS or patternsearch are provided.
nll_diff_is.m, nll_map_is.m and, nll_en_is.m compute the values of the objective functions for the Difference, Max and Entropy models respectively.

In these three functions, response distributions of category decisions and confidence reports are estimated by Monte-Carlo simulations, 
and negative log-likelihood of the model is then computed.

