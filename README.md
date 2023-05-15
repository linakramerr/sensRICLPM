# sensRICLPM
package to conduct sensitivity analysis of RI-CLPM results to unmodeled measurement error.

To install package type: 

install.packages("devtools")

devtools::install_github("linakramerr/sensRICLPM")

# example:

The package includes the amotivation.Rdata data set. Type:

load(~/data/amotivation.RData) # to load data

example <- sensRICLPM(sensRICLPM::amotivation) # to run sensitvity analysis on data

example$p2 # to obtain sensitvity plots
example$outputtable # to obtain all sensitivity analysis results
