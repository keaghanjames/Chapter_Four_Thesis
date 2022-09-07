# Script for undertaking principle coordinate analysis of assemblages and species and then calculating functional dispersion scores for each assemblage

### Weeks et al 2022 data

setwd("~/Dropbox/AVONET")
require(tidyr)
require(stringr)
require(phytools)
require(dplyr)
require(caper)
require(mFD)
library(foreach)
library(doParallel)

#just establishing some shortcuts for the file structure
TraitDataFolder <- "ELEData/TraitData/"
PhylogeneticDataFolder <- "ELEData/PhylogeneticData/"
codeFolder <- "ELECode/"
figFolder <- "figFolder/"
SpatialDataFolder <- 'ELEData/SpatialData/'

#read in copy of Jetzz et al. 2018 tree as used Weeks et al. 2020
tree <- read.tree("Weeks_et_al_2022/singe_bird_phylo.tre")
Jetz.PAM2 <- read.csv('Weeks_et_al_2022/Presence-Absence-Matrix_7-19-21.csv')
tibble(Jetz.PAM2)
Jetz.PAM2 <- Jetz.PAM2[which(Jetz.PAM2$Tip_Label %in% tree$tip.label),]

#read in Weeks et al. (2020) assemblage data
comm <- read.csv('presence_absence_Weeks_et_al.csv')



#load in AVONET trait dataset
df <- as_tibble(read.csv("ELEData/TraitData/AVONET1_BirdLife.csv"))

df$Species1 <- gsub(" ", "_", df$Species1)

#next we need to add the Jetz taxonomy to our data set
crossWalk<-read.csv(paste0(PhylogeneticDataFolder,'BirdLife-BirdTree crosswalk.csv')) # taxonomy crosswalk
crossWalk$Species3<-gsub(" ","_",crossWalk$Species3) # add underscore to species names in the crosswalk to match those in the phylogeny
crossWalk$Species1 <- gsub(" ", '_', crossWalk$Species1)
Species3 <- vector()

for(i in 1:length(df$Species1)){
  Species3 <- c(Species3,
                crossWalk$Species3[which(df$Species1[i] == crossWalk$Species1)[1]]     
  )
}

length(Species3) - length(unique(Species3))

df$Species3 <- Species3
#now we need to drop double ups
df <- df[-(which(duplicated(df$Species3) == T)),]
#now lets drop any species that are not in our presence absence data
df <- df[which(df$Species3 %in% colnames(comm) == T),]
df <- df[which(df$Species3 %in% tree$tip.label == T),]
dim(df)



#next we need to calculate functional dispersion for each of our assemblages
df_traits <-  as.data.frame(df[,c(11:15, 20:21)]) #we select only the traits used by Weeks et al. (2020) and body mass
rownames(df_traits) <- df$Species3
write.csv(df_traits, 'df_traits.csv')
df_traits <- read.csv('df_traits.csv', row.names = 1)
df_traits <- df_traits[,1:6]/df_traits[,7] #correct all traits by body mass
tr_cat <- tibble(trait_name = colnames(df_traits[1:6]), trait_type = 'Q')
tr_cat

#implement principle coordinate analysis for trait and assemblage data 
#first calculate functional distance matrix
dist <- funct.dist(sp_tr = df_traits, tr_cat = tr_cat, metric = 'euclidean', scale_euclid = 'scale_center')
#Sys.setenv('R_MAX_VSIZE'=32000000000) #may need to increase maximim vector size
#run principle coordinate analysis
system.time(quality <- quality.fspaces(sp_dist = dist, maxdim_pcoa = 6, fdendro = NULL))
save(quality, file = 'quality_fspaces_by_bodysize.Rdata')

q <- quality$quality_fspaces
format(q, scientific = F)

#we save our species principle coordinates because we will need them later
sp_faxes_coord <- quality$details_fspaces$sp_pc_coord
save(sp_faxes_coord, file = 'sp_faxes_coord_by_bodysize.Rdata')


#we need to convert our presence absence data into a matrix so it can be parsed as the species weight information 
#it is a bit fidly to get it in the right format
dist_weights <- as.matrix(comm[,-c(1,2,3)])
dist_weights <- dist_weights[,which(colnames(dist_weights) %in% rownames(sp_faxes_coord) == T)]
rownames(dist_weights) <- paste0('site', seq(1, nrow(dist_weights), 1))

#okay finally we can estimate the functional dispersion at each of the assemblages 
alpha.fd.indices <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')],
  asb_sp_w = dist_weights,
  ind_vect = c('fdis'), 
  scaling = T, 
  check_input = T, 
  details_returned = T
)
save(alpha.fd.indices, file = 'alpha.fd.indices_by_mass.Rdata')


# next we will calculate the phylogenetic diversity of each assemblage
comm <- read.csv('presence_absence_Weeks_et_al.csv')


require(caper)
cores=detectCores()

#function for getting PD and mean pairwise distance
get_pd <- function(as) {
  t <- keep.tip(phy = tree, tip = colnames(comm)[(which(comm[as,] == 1))])
  pd <- pd.calc(t)[[1]]
  co <- mean(cophenetic(t))
  return(setNames(c(pd, co), nm = c("PD", "MPD")))
}

assemb <- seq(1, nrow(dist_weights), 1)

require(pbmcapply)
pd <- pbmclapply(assemb, get_pd, mc.cores = detectCores()-1)
pd <- do.call(rbind, pd)


#we will add this to output of our functional dipsersion analysis
alpha.fd.indices$functional_diversity_indices$pd <- pd[,1]
alpha.fd.indices$functional_diversity_indices$mpd <- pd[,2]

alpha.fd.indices$functional_diversity_indices$pd <- pd
alpha.fd.indices$functional_diversity_indices$lat <- comm$lat[assemb]
alpha.fd.indices$functional_diversity_indices$long <- comm$long[assemb]
write.csv(alpha.fd.indices$functional_diversity_indices, 'functional_diversity_indicesFULL.csv', row.names = T, col.names = T)
