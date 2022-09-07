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

tree <- read.tree("Weeks_et_al_2022/singe_bird_phylo.tre")
#plot(tree, show.tip.label = F, type = 'fan')
Jetz.PAM2 <- read.csv('Weeks_et_al_2022/Presence-Absence-Matrix_7-19-21.csv')
tibble(Jetz.PAM2)

Jetz.PAM2 <- Jetz.PAM2[which(Jetz.PAM2$Tip_Label %in% tree$tip.label),]
# 
# #Convertthe format so that each row is a cell, and lat/long columns along with species presence/absence are included
# ids <- unique(Jetz.PAM2$WorldID)
# species <- unique(Jetz.PAM2$Tip_Label)
# id <- unique(Jetz.PAM2$WorldID)
# library(foreach)
# library(doParallel)
# 
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# data <- foreach(i = 1:length(species), .combine = cbind) %dopar% {
#   pres <- numeric(length=length(unique(Jetz.PAM2$WorldID)))
#   pres[ids %in% Jetz.PAM2$WorldID[Jetz.PAM2$Tip_Label==species[i]]]<-1
#   pres
# }
# stopCluster(cl)
# 
# 
# data <- cbind(ids, data)
# colnames(data) <- c("id", species)
# tibble(data)
# 
# 
# 
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# 
# latlong <- foreach(i = 1:nrow(data), .combine = rbind) %dopar% {
#  pres <- c(unique(Jetz.PAM2$y[Jetz.PAM2$WorldID==data[i,1]]),
#                   unique(Jetz.PAM2$x[Jetz.PAM2$WorldID==data[i,1]]) )
#  pres
# }
# 
# head(latlong)
# 
# data <- cbind(latlong, data)
# colnames(data)[1:2] <- c('lat', 'long')
# 
# 
# stopCluster(cl)
# 
# ## Trim the matrix to communities with 7 or more species, so that functional richness can be calculated (see Weeks et al. 2020 Methods), and remove any species that have no presence (i.e. they were only in the communties with fewer than 7 species)
# Jetz.PAM.trimmed<-data[which(rowSums(data[,4:ncol(data)])>6),]
# comm<-Jetz.PAM.trimmed[,which(colSums(Jetz.PAM.trimmed)>0)]
# 
#write.csv(comm, file = 'presence_absence_Weeks_et_al.csv', row.names = F, col.names = T)

comm <- read.csv('presence_absence_Weeks_et_al.csv')



#load in our dataset
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
#we may need to fix up the presence absence data

dim(df)



#next let us get our function dispersion

df_traits <-  as.data.frame(df[,c(11:15, 20:21)])
rownames(df_traits) <- df$Species3

head(df_traits)

write.csv(df_traits, 'df_traits.csv')
df_traits <- read.csv('df_traits.csv', row.names = 1)

head(df_traits)

df_traits <- df_traits[,1:6]/df_traits[,7]


# 
tr_cat <- tibble(trait_name = colnames(df_traits[1:6]), trait_type = 'Q')
tr_cat
dist <- funct.dist(sp_tr = df_traits, tr_cat = tr_cat, metric = 'euclidean', scale_euclid = 'scale_center')
rm(comm)
rm(df)
rm(df_traits)
rm(tree)
Sys.setenv('R_MAX_VSIZE'=32000000000)
system.time(quality <- quality.fspaces(sp_dist = dist, maxdim_pcoa = 6, fdendro = NULL))
save(quality, file = 'quality_fspaces_by_bodysize.Rdata')

#load('quality_fspaces.Rdata')

q <- quality$quality_fspaces
format(q, scientific = F)
#quality.fspaces.plot(quality, quality_metric = 'mad', fspaces_plot = rownames(quality$quality_fspaces))


sp_faxes_coord <- quality$details_fspaces$sp_pc_coord
save(sp_faxes_coord, file = 'sp_faxes_coord_by_bodysize.Rdata')


#here
df_traits <- read.csv('df_traits.csv', row.names = 1)
comm <- read.csv('presence_absence_Weeks_et_al.csv')
#load('sp_faxes_coord.Rdata')
#tree <- read.tree('Weeks_et_al_2022/singe_bird_phylo.tre')

#funct.space.plot(sp_faxes_coord = sp_faxes_coord)
traits.faxes.cor(sp_tr = df_traits, sp_faxes_coord = sp_faxes_coord)

dist_weights <- as.matrix(comm[,-c(1,2,3)])
head(dist_weights)

dist_weights <- dist_weights[,which(colnames(dist_weights) %in% rownames(sp_faxes_coord) == T)]

rownames(dist_weights) <- paste0('site', seq(1, nrow(dist_weights), 1))
head(dist_weights)[,1:10]

Sys.setenv('R_MAX_VSIZE'=32000000000)

#assemb <- sample(1:nrow(dist_weights), 1000)

alpha.fd.indices <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6')],
  asb_sp_w = dist_weights,
  ind_vect = c('fdis'), 
  scaling = T, 
  check_input = T, 
  details_returned = T
)

save(alpha.fd.indices, file = 'alpha.fd.indices_by_mass.Rdata')
load('alpha.fd.indices_backup.Rdata')

alpha.fd.indices$functional_diversity_indices
tree <- read.tree("Weeks_et_al_2022/singe_bird_phylo.tre")


#lets get pd

rm(alpha.fd.indices)

comm <- read.csv('presence_absence_Weeks_et_al.csv')


require(caper)
cores=detectCores()

get_pd <- function(as) {
  t <- keep.tip(phy = tree, tip = colnames(comm)[(which(comm[as,] == 1))])
  pd <- pd.calc(t)[[1]]
  co <- mean(cophenetic(t))
  return(setNames(c(pd, co), nm = c("PD", "MPD")))
}

assemb <- seq(1, nrow(dist_weights), 1)
length(assemb)

require(pbmcapply)
pd <- pbmclapply(assemb, get_pd, mc.cores = detectCores()-1)
#pd <- unlist(pd)
pd <- do.call(rbind, pd)
pd

#load('alpha.fd.indices_backup.Rdata')

alpha.fd.indices$functional_diversity_indices$pd <- pd[,1]
alpha.fd.indices$functional_diversity_indices$mpd <- pd[,2]

alpha.fd.indices$functional_diversity_indices$pd <- pd
alpha.fd.indices$functional_diversity_indices$lat <- comm$lat[assemb]
alpha.fd.indices$functional_diversity_indices$long <- comm$long[assemb]
write.csv(alpha.fd.indices$functional_diversity_indices, 'functional_diversity_indicesFULL.csv', row.names = T, col.names = T)
