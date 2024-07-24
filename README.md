# Approximating Unidentifiable Causal Effects
*Code repository for a Bachelor Thesis written by Lisa Schulze-Bergkamen*

This repository contains the following:
- Code provided by the authors of both papers in the respective folders
- Code written by me to run the experiments in Chapter 4 of my thesis

## Non-Continuous Paper (An Automated Approach to Causal Inference in Discrete Settings)
Two types of experiments were done for this paper:
- Analysing the length of the run time when the number of nodes in the graph increase. The code 
  for this can be found in `Non-Continuous Paper/Experiments/Runtime Nodes`
- Analysing the calculated causal effect for both identifiable and non-identifiable causal 
  effects when the sample size of the observational data increases. The code for this can be 
  found in `Non-Continuous Paper/Experiments/Bounds Samples`

The folder `Source Code` contains the code and demos provided by the authors. To run the program 
please refer to the `README.md` file within this folder.

## Continuous Paper (Bounding Causal Effects on Continuous Outcome)
Due to the limited amount of code the author provided, I was not able to run experiments 
regarding the calculations of the causal bounds. Due to this paper being older, the author 
could only recover the code for the kl-UCB algorithm, which is the algorithm in the file 
`Continuous Paper/Source Code/b_kl_ucb.m`. The code for this algorithm was written in `Matlab`.
