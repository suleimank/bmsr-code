#clean-up console
rm(list=ls())

#load packages
library(rstan)

#load model and utility libraries
source('demo_utils.R')
source("bmsr.R")

#set random seed
set.seed(101)

#run demo for bmsmtr
demo_bmsr(file = "bmsmtr.stan",dY = 3)
