# Antonio Cozzolino Portfolio

This GitHub repository contains personal sample codes written in various programming languages, primarily including Python, Matlab, R, and Julia. The repository is still a work in progress, so some sections may not be complete yet.

Here is a brief description of the folders within the repository:

## 1. Bachelor Thesis

Brief translation of the work, WIP.

## 2. Master Thesis
**Title**: Skill Heterogeneity, Idiosyncratic Risk, and Labor Market Fluctuations. 

**Abstract**: This thesis proposes a theoretical framework to study the relationship between skill heterogeneity, income uncertainty, and non-pecuniary motivations influencing occupational choices. I develop a dynamic discrete choice model of occupational choice where workers differ by productivity, and are subject to income uncertainty and taste shocks. I use this model to study the behavior of individuals and run some numerical simulations to assess the effects of such heterogeneity on the economy. The main finding is that financially disadvantaged households prioritize jobs that maximize income and align with their comparative advantage. In contrast, more rich households are willing to sacrifice potential earnings for jobs that provide greater satisfaction. At an aggregate level, heterogeneity in preference for occupations leads to movements in effective labor supply that affect aggregate labor productivity.

The Matlab codes are the original one used for the thesis. Julia codes instead are a re-implementation that allows for Parallelization and Interpolation of the Value Function of the model. The Julia Code is a work in progress, I'm moving the code there since I'm planning to do an estimation exercise of the model. 

## 3. Robust Identification and Inference in Repeated Games

**Abstract**: We study identification in repeated games with imperfect public monitoring under no assumptions on equilibrium selection. We leverage the restrictions on equilibrium payoffs in Awaya and Krishna (2019) to construct a law of large numbers that implies bounds on payoff parameters, the discount factor, and monitoring quality. Our method is computationally tractable in generic games. We illustrate the identified sets of parameters in simulations, highlighting the trade-off between maintaining strong assumptions and the tightness of the identified set.

This work is a collaboration with Cristina Gualdani, Niccolò Lomys, and Lorenzo Magnolfi. It features a sample code extracted and simplified from the Python library I am currently developing for the estimation procedure. For illustrative purposes, I have chosen to showcase the estimation of a repeated Prisoner's Dilemma. The essential functions are included in the "RIIRG_lib.py" file, while the usage of the method is demonstrated in the "main.ipynb" notebook.
 
## 4. Data Analysis
This section contains some data tasks example from "predoc.org". It contains the exercise and the code. In the "Data Manipulation Medicare Advantage" exercise I show some code in R that involves a simple data cleaning and dataset creation.

The "Economic Analysis" exercise instead is more involved, and requires to use a subsample from "U.S. Current Population Survey" to produce an answer to the question of how hourly wages (“wage”) and labor force participation (“lfp”) evolved for skilled and unskilled workers since 1976. WIP

## 5. ANN and Function Approximation

This folder contains some example of code I produced about Neural Networks for Function Approximation. 

- The file "Ramsey_Model_ANN_Approx.ipynb" is a replication of the paper "Deep learning for solving dynamic economic models, Maliar et al. 2021" in a simpler setting. 

- The file "Static_Unconstrained_Optimization.ipynb" is instead an example of how one can solve for a Value function of an optimization problem by means of ANN. The advantage of doing this is that, instead of solving multiple optimiziation problem for each parameter, one can sample points from the parameter space and solve the entire optimization problem in that space.

## 6. Computational Economics

This folder contains some other example of code in Julia. There is the standard Cake Eating Problem solved via VFI. Instead the file "Brock_Mirman_Interp_Quadrat_Julia.ipynb" is a routine that solves the Brock-Mirman model via VFI with Cubic Spline Interpolation and Quadrature to approximate expectations.

## 7. Other Projects

This last folder contains other projects that I think may be interesting. It includes sample codes from my previous exams:

- Derivatives/Option Pricing: This folder contains a Matlab file created for the take-home exam, in which I had to perform numerical simulation and pricing of a 'down-and-out' put option. The folder includes the Matlab script and the exercise's text."

- Econometrics_Problem_Set: Is a collection of R script produced for various Problem Sets. WIP