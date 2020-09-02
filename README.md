# bstats
Bootstrapped t-values for MHT robust p-values in Stata

Romano & Wolf R code for calculating adjusted p-values requires the input of bootstrapped t-values from the models being estimated. This code creates these files ready for upload to R to run the MHT correction code

This code is 99% based on the rwolf command by Damian Clarke (@damiancclarke) but stops short of actually calculating MHT adjusted p-values. This is because it is meant for use in non-standard examples, where the outcome variable is the same one but what changes in each hypotheses tested are the independent variables.
