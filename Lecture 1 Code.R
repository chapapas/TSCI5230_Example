###### Holly Chapapas - 25Aug2021 - TSCI5230 Lecture 1 Code #####

#Generate R report with Cntrl+Shift+K or click little notebook button

#+ init, echo=FLASE, message=FALSE, warning=FALSE, results='hide'
debug <- 0;
knitr::opts_chunk$set(echo=debug>0, warning=debug>0, message=debug>0);

#libraries R currently sees
#search() %>% pander()
#library(pander) will format the search command in the R report nicer for publication purposes
library(GGally)
library(rio)
library(dplyr)
library(synthpop)

#add data to environment
dat0 <- survival::veteran
export(dat0, 'example_file.xlsx')
dat0 <- import('example_file.xlsx')
# scatter plot matrix
ggpairs(dat0)
unique(dat0$status)
#assess with true false
dat0$status == 1
unique(dat0$prior)
#assess with true false
dat0$prior > 0
#replace status column to remove continuous variables
dat0$status <- dat0$status == 1
#look at the new results
ggpairs(dat0)
# set all the 2 value columns to be TRUE/FALSE
dat1 <- mutate((dat0), across(where( function(xx) length(unique(xx))==2), as.logical))
# scatterplot matrix again...
ggpairs(dat1)

