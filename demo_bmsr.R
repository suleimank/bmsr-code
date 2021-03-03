#clean-up console
rm(list=ls())

#load packages
library(rstan)

#load model and utility libraries
source('demo_utils.R')
source("bmsr.R")

#set random seed
set.seed(101)

#run demo for bmsr
demo_bmsr(file = "bmsr.stan",dY = 1)
