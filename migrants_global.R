setwd('~/Dropbox/AVONET')

require(ggplot2)
require(dplyr)
comm <- read.csv('presence_absence_Weeks_et_al.csv')
head(comm)[,1:10]

migrants <- read.csv('Dufour_et_al_2019.csv')
head(migrants)


length(which(migrants$Sp.Scien.jetz %in% colnames(comm) == F))

mig_data <- tibble(Sp = colnames(comm)[-c(1:3)], migrant = NA)

mig_data$migrant <- migrants$strategy_3[match(mig_data$Sp, migrants$Sp.Scien.jetz)]

mig_data$migrant <- as.factor(mig_data$migrant)


# proportion_strict <- vector()
# proportion_migrant <- vector()
strict_migrants <- vector()
all_migrants <- vector()
richness <- vector()

df <- comm[4:ncol(comm)]

for(i in 1:nrow(df)){
  assem <- colnames(df)[which(df[i,] == 1)]
  strict <- length(which(mig_data$Sp %in% assem == T & mig_data$migrant == 'strict_mig'))
  all <-    length(which(mig_data$Sp %in% assem == T & mig_data$migrant == 'strict_mig' |
                             mig_data$Sp %in% assem == T & mig_data$migrant ==  'partial_mig'))
  strict_migrants <- c(strict_migrants, strict)
  all_migrants <- c(all_migrants, all)
  richness <- c(richness, length(which(mig_data$Sp %in% assem == T)))
  print(i)
}

migration_data <- tibble(site = comm$id, lat = comm$lat, long = comm$long,
                         richness = richness, strict_migrants = strict_migrants,
                         all_migrants = all_migrants)

migration_data$prop_migrants <- migration_data$all_migrants/migration_data$richness
migration_data$prop_strict_migrants <- migration_data$strict_migrants/migration_data$richness
migration_data$Region <- NA

migration_data$Region[which(migration_data$long  >= 65 & migration_data$long < 180)] <- 'Asia and Oceania'
migration_data$Region[which(migration_data$long  >= -30 & migration_data$long < 65)] <- 'Africa and Europe'
migration_data$Region[which(migration_data$long  >= -180 & migration_data$long < -30)] <- 'Americas'

migration_data$Region <- factor(migration_data$Region, ordered = T, levels = c('Americas', 'Africa and Europe', 'Asia and Oceania'))
colours <- c( "#FE5D26", "#340068", "#EBDDFF", 'darkgrey')

functional_indices <- read.csv('functional_diversity_indicesFULL.csv')
migration_data$FD <- functional_indices$fdis
migration_data$PD <- functional_indices$pd
prop_mig <- ggplot(data = migration_data, aes(x = lat, y = prop_migrants, col = prop_migrants)) +
        geom_point(alpha = .4) +
        geom_smooth(se = F, col = 'black', method = 'gam', aes(linetype = Region), name = 'Region', lwd = .7) +
        theme_classic() +
        scale_x_discrete(limits = seq(-40,80,20)) +
        scale_linetype_manual(values=c("solid", 'dashed', "dotdash"))+
        scale_color_gradient(high = '#86A8E7',  low = '#D16BA5', name = 'Proportion migratory species') +
        ylab('Proportion of migratory species') +
        xlab('Latitude')
prop_mig

prop_strict_mig <- ggplot(data = migration_data, aes(x = lat, y = prop_strict_migrants, col = abs(lat))) +
  geom_point(alpha = .4) +
  geom_smooth(se = F, col = colours[4], method = 'loess', aes(linetype = Region), name = 'Region', lwd = .7) +
  theme_classic() +
  scale_x_discrete(limits = seq(-40,80,20)) +
  scale_linetype_manual(values=c("solid", 'dashed', "dotdash"))+
  scale_color_gradient(high = '#86A8E7',  low = '#D16BA5', name = 'Absolute latitude') +
  ylab('Proportion of strictly migratory species') +
  xlab('Latitude') 

prop_strict_mig


world_map = map_data("world") %>% 
  filter(! long > 180)

library(ggthemes)
map_prop_migrants <- 
  ggplot() +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  theme_map() +
  geom_point(data = migration_data, aes(x = long, y = lat, col = prop_migrants), size = .6, alpha = 1) +
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  scale_color_gradient(high = '#86A8E7',  low = '#D16BA5', name = '') +
  geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none")     

map_prop_migrants

map_prop_strict_migrants <- 
  ggplot() +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  theme_map() +
  geom_point(data = migration_data, aes(x = long, y = lat, col = prop_strict_migrants), size = .6, alpha = 1) +
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  scale_color_gradient(low = '#86A8E7',  high = '#D16BA5', name = '') +
  geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "bottom")  +
  ggtitle('Proportion of strict migrants')     

map_prop_strict_migrants


require(egg)
grid_plot <- ggarrange(map_prop_migrants, prop_mig,
                        nrow = 2, widths = c(3,1.5))
grid_plot
ggsave(plot = grid_plot, filename = 'migrants_proportion.pdf', device = 'pdf', scale = .8)


load('sp_faxes_coord_by_bodysize.Rdata')
head(sp_faxes_coord)

sp_faxes_coord <- as_tibble(sp_faxes_coord, rownames = 'Species')
sp_faxes_coord

sp_faxes_coord$migrant <- mig_data$migrant[match(sp_faxes_coord$Species, mig_data$Sp)]
sp_faxes_coord$migrant_binary <- 'resident'
sp_faxes_coord$migrant_binary[which(sp_faxes_coord$migrant == 'partial_mig' | sp_faxes_coord$migrant ==  'strict_mig')] <- 'migratory'
sp_faxes_coord$migrant_binary <- as.factor(sp_faxes_coord$migrant_binary)

require(ggalt)
require(ggpubr)
pc1_2 <- ggplot(data = na.omit(sp_faxes_coord), aes(x = PC1, y = PC2, group = migrant_binary))+
  geom_encircle(aes(fill = migrant_binary), s_shape = 1, expand = 0,
               alpha = .1, color = "black", show.legend = F) +  
  #stat_chull(aes(fill = migrant_binary)) +
  geom_point(aes(col = migrant_binary, shape = migrant_binary))+
  scale_color_manual(values = alpha(c('black',  'grey'), 1), name = '', labels = c('Migratory', 'Resident')) +
  scale_fill_manual(values = alpha(c('black',  'lightgrey'), 1), name = '') +
  scale_shape_manual(values=c(3, 4), name = '', labels = c('Migratory', 'Resident'))+
  theme_classic() +
  theme(legend.position = 'bottom')
pc1_2

pc3_4 <- ggplot(data = na.omit(sp_faxes_coord), aes(x = PC3, y = PC4, group = migrant_binary))+
  geom_encircle(aes(fill = migrant_binary), s_shape = 1, expand = 0,
                alpha = .1, color = "black", show.legend = F) +  
  #stat_chull(aes(fill = migrant_binary)) +
  geom_point(aes(col = migrant_binary, shape = migrant_binary))+
  scale_color_manual(values = alpha(c('black',  'grey'), 1), name = '', labels = c('Migratory', 'Resident')) +
  scale_fill_manual(values = alpha(c('black',  'lightgrey'), 1), name = '') +
  scale_shape_manual(values=c(3, 4), name = '', labels = c('Migratory', 'Resident'))+
  theme_classic() +
  theme(legend.position = 'bottom')
pc3_4

pc5_6 <- ggplot(data = na.omit(sp_faxes_coord), aes(x = PC5, y = PC6, group = migrant_binary))+
  geom_encircle(aes(fill = migrant_binary), s_shape = 1, expand = 0,
                alpha = .1, color = "black", show.legend = F) +  
  #stat_chull(aes(fill = migrant_binary)) +
  geom_point(aes(col = migrant_binary, shape = migrant_binary))+
  scale_color_manual(values = alpha(c('black',  'grey'), 1), name = '', labels = c('Migratory', 'Resident')) +
  scale_fill_manual(values = alpha(c('black',  'lightgrey'), 1), name = '') +
  scale_shape_manual(values=c(3, 4), name = '', labels = c('Migratory', 'Resident'))+
  theme_classic() +
  theme(legend.position = 'bottom')
pc5_6


grid_plot2 <- ggarrange(pc1_2, pc3_4, pc5_6,
                       nrow = 1, widths = c(1,1,1))
grid_plot2
ggsave(plot = grid_plot2, filename = 'PCoA_migrantory.pdf', height = 5, width = 12, device = 'pdf', scale = .8)
