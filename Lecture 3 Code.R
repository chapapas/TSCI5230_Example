###### Holly Chapapas - 8Sept2021 - TSCI5230 Lecture 3 Code #####

library('rio');
library('synthpop');

######################   import your synthetic data
orig0 <- import('C:/Users/chapa/Documents/UTHSCSA/Course Work/Spring 2021/TSCI5230_Example/over_expressed_genes2.csv')
orig0cb <- codebook.syn(orig0);
sim0 <- syn(orig0);
compare(sim0,orig0);
export(sim0$syn,"sim_data.csv");

inputdata <- c(dat0='over_expressed_genes_sim.csv');
if(file.exists('local.config.R')) source('local.config.R', local = TRUE, echo = FALSE);
dat0 <- import(inputdata['dat0'])