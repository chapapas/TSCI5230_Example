###### Holly Chapapas - 8Sept2021 - TSCI5230 Sim Data File Code #####

library('rio');
library('synthpop');

######################   import your original data
orig0 <- import('C:/Users/chapa/Documents/UTHSCSA/Course Work/Fall 2021/TSCI5230/over_expressed_genes.csv')
orig0cb <- codebook.syn(orig0);
sim0 <- syn(orig0);
compare(sim0,orig0);
# Export for later use
export(sim0$syn,"sim_data.csv")