# NULL Hypothesis 


setwd("~/Dropbox/AVONET")
require(tidyr)
require(stringr)
require(phytools)
require(dplyr)
require(caper)
require(mFD)
library(foreach)
library(doParallel)

#we read in  our tree and our species coordinates in morphospace, and our presence absence information
tree <- read.tree("Weeks_et_al_2022/singe_bird_phylo.tre")
load('sp_faxes_coord_by_bodysize.Rdata')

comm <- read.csv('presence_absence_Weeks_et_al.csv')
#lets randomly sample comm

head(comm)[1:10]
dim(comm)

#load("sp_faxes_coord_by_bodysize.Rdata")

get_MPFD <- function(as){
  MPFD <- mean(dist(sp_faxes_coord[which(rownames(sp_faxes_coord) %in% names(as)[(which(as == 1))]),], method = 'euclidean'))
}

#create a sample of or presence absence matrix
dist_weights <- as.matrix(comm[,-c(1,2,3)])
dist_weights <- dist_weights[,which(colnames(dist_weights) %in% rownames(sp_faxes_coord) == T)]
rownames(dist_weights) <- paste0('site', seq(1, nrow(dist_weights), 1))
MPFD <- apply(dist_weights, 1, get_MPFD)

#dist_weights <- dist_weights[sample(seq(1, nrow(dist_weights), 1), 1000),]
lat_long <- comm[as.numeric(str_remove(rownames(dist_weights), 'site')),1:2]
rm(comm)
head(dist_weights)[,1:5]

#function to randomise species composition in each assemblage 
make_weights <- function(site){
  n <- sum(dist_weights[site,])
  assemb <- sample(colnames(dist_weights), n)
  weight <- rep(0,ncol(dist_weights))
  names(weight) <- colnames(dist_weights)
  weight[which(names(weight) %in% assemb)] <- 1
  return(weight)
}

simulate_weights <- function(n){
  new_data <- lapply(rownames(dist_weights), FUN = make_weights)
  new_data <- do.call(rbind, new_data)
  saveRDS(new_data, file = paste0('NULL/newdata_run', n, '.rds'))
  rm(new_data)
}

#simulate_weights(1)

require(pbmcapply)
Sys.setenv('R_MAX_VSIZE'=32000000000)
new_data <- pbmclapply(X = seq(5,999,1), FUN = simulate_weights, mc.cores = detectCores()-3)


### Calculate FD and PD for Synthetic Assemblages ###
site_names <- rownames(dist_weights)
rm(dist_weights)
rm(distances)
rm(lat_long)

clade_matrix <- clade.matrix(tree)
get_pd <- function(as) {
  t <- keep.tip(phy = tree, tip = names(as)[(which(as == 1))])
  pd <- pd.calc(cm = clade_matrix, tip.subset = names(as)[(which(as == 1))])
  co <- mean(cophenetic(t))
  return(setNames(c(pd, co), nm = c("PD", "MPD")))
  #return(setNames(c(pd), nm = c("PD")))
}


get_stats <- function(run){
  new_data <- readRDS(paste0('NULL/newdata_run', run, '.rds'))
  PD <- apply(X = new_data, 1, get_pd)
  MPFD <- apply(X = new_data, 1, get_MPFD)
  stats <- tibble(t(PD)[,1], t(PD)[,2], MPFD, site = site_names)
  colnames(stats) <- c('PD', 'MPD', 'MPFD', 'site')
  #stats <- tibble(PD, MPFD, site = site_names)
  #colnames(stats) <- c('PD', 'MPFD', 'site')
  stats$sp_richn <- rowSums(new_data)
  stats$run <- run
  saveRDS(stats, file = paste0('NULL/summary_run', run, '.rds'))
  rm(new_data)
  return(stats)
}

system.time(get_stats(1))

require(pbmcapply)
pbmclapply(seq(201,999,1), get_stats, mc.cores = detectCores()-6)
stats <- do.call(rbind, stats)
stats

x <- readRDS('NULL/summary_run1.rds')
x
