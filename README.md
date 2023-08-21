This repository contains public code for the manuscript:

Zou, Heng-Xing, SJ Schreiber, and Volker HW Rudolf. "Stage-mediated priority effects and season lengths shape long-term competition dynamics" *bioRXiv* (2020). [https://doi.org/10.1101/2020.08.14.251926](https://doi.org/10.1101/2020.08.14.251926)

For any questions, please [contact Heng-Xing Zou](hengxingzou@rice.edu)

# Abstract

The relative arrival time of species can affect their interactions and thus determine which species persist in a community. Although this phenomenon, called priority effect, is widespread in natural communities, it is unclear how it depends on the length of growing season. Using a seasonal stage-structured model, we show that differences in stages of interacting species could generate priority effects by altering the strength of stabilizing and equalizing coexistence mechanisms, changing outcomes between exclusion, coexistence, and positive frequency dependence. However, these priority effects are strongest in systems with just one or a few generations per season and diminish in systems where many overlapping generations per season dilute the importance of stage-specific interactions. Our model reveals a novel link between the number of generations in a season and the consequences of priority effects, suggesting that consequences of phenological shifts driven by climate change should depend on specific life histories of organisms.

# Code Files

All simulations were run in R version 4.2.1.

- `Functions.R` contains all functions used for contructing the model and running simulations. 
- `Analysis.Rmd` contains all code for generating figures in the main text and the appendix. Running any portion of this file will automatically run `Functins_SpatialBH.R`. Note that running all simulations takes a long time; we therefore do not recommend knitting this file. Instead, we have provided `SavedEnviron.RData` that contains all generated data.

# Data

`SavedEnviron.RData` contains generated data from all simulations in the main text. This is a pre-packaged R dataset that includes all variables in the environment. To load the dataset, simply double click while `Functions.R` and `Analysis.Rmd` are opened in RStudio. Loading `SavedEnviron.RData` will load results from invasion analyses using the following model formulation:
- **Stage-mediated interspecific competition (3 data frames)**: `lowintra`, `midintra` and `highintra`. 
- **No stage-mediated interspecific competition (3 data frames)**: `xlowintra`, `xmidintra` and `xhighintra`.
- **Both stage-mediated intra- and interspecific competition (3 data frames)**: `vlowintra`, `vmidintra` and `vhighintra`.
- **Periodic arrival time with stage-mediated interspecific competition (1 data frame)**: `var2`.
- **Periodic arrival time without stage-mediated interspecific competition (1 data frame)**: `xvar2`.
- **Simulation with $\Delta s>5$, with stage-mediated interspecific competition (1 data frame) (1 data frame)**: `long_midintra`.
- **Simulation with $\Delta s>5$, without stage-mediated interspecific competition (1 data frame) (1 data frame)**: `xlong_midintra`.
- **Simulation with adult density dependence, with stage-mediated interspecific competition (1 data frame) (1 data frame)**: `adu_dd_midintra`.
- **Simulation with adult density dependence, without stage-mediated interspecific competition (1 data frame) (1 data frame)**: `xadu_dd_midintra`.

# Reproducing Figures

- Make sure all files are under the same directory.
- Open `Analysis.Rmd`, then load `SavedEnviron.RData`.
- Follow specific instructions in `Analysis.Rmd` to generate figures in the main text, or change model parameters to produce figures in the appendix.

Latest release: version Jan 2023; uploaded Aug 21 2023
