#script for calculating standerdised effect sizes of PD and FD from distance-weighted null and undertaking structual equation analysis
#the script also includes code for producing Figures 1-5 and Figures S1-S2.


setwd("~/Dropbox/AVONET")
require(tidyr)
require(stringr)
require(phytools)
require(dplyr)
require(caper)
require(mFD)
library(foreach)
library(doParallel)
require(fields)
require(picante)
require(metricTester)
require(ggpubr)
col <- c('#0A2342', '#0B7A75', '#FA8334', '#7B2D26', '#9593D9')
colours <- c('#340068', 'darkgrey', '#FE5D26')

#############################################
### Calculating Standerdised Effect Sizes ###
#############################################
setwd("D_NULL")
df_null <- readRDS('full_null.rds')
df_emp <- read.csv('~/Dropbox/AVONET/ses_empircal.csv')

#we need to get the mean of each indicator index for each site
df_null$site <- as.factor(df_null$site)
df_null$run <- as.factor(df_null$run)


df <- df_null %>%
  group_by(site) %>%
  summarise(PDn = mean(PD), PDsd = sd(PD),
            MPDn = mean(MPD), MPDsd = sd(MPD),
            MPFDn = mean(MPFD), MPFDsd = sd(MPFD),)
#we just need to do some reworking 
df$site_num <- as.character(df$site)
df$site_num <- as.numeric(str_remove(df$site_num, 'site'))
df <- df[order(df$site_num),]

#now we add the empirical observations
df$PDe <- df_emp$pd
df$MPDe <- df_emp$mpd
df$MPFDe <- df_emp$MPFD
df$lat <- df_emp$lat
df$abs_lat <- df_emp$abs_lat
df$long <- df_emp$long

#and now we can calculate the ses
df$PDses <- (df$PDe - df$PDn)/df$PDsd
df$MPDses <- (df$MPDe - df$MPDn)/df$MPDsd
df$MPFDses <- (df$MPFDe - df$MPFDn)/df$MPFDsd
df$sp_richn <- df_emp$sp_richn
write.csv(df, 'standerdised_effect_sizes.csv', row.names = F)

plot(df$PDses, df$MPFDses)

df$lat_cat <- cut(df$abs_lat, 4)
df$long_cats <- df_emp$long_cats
require(ggplot2)

#Plot for SES MPD and MPFD

p0 <- ggplot(df, aes(x = MPDses, y = MPFDses, col = abs_lat)) +
  geom_point() +
  geom_smooth(method = 'lm', color = '#78C0E0') +
  scale_color_gradient2(low = col[1], mid = col[2], high = col[3]) +
  #scale_color_gradient(low = '#FE5D26', high = '#340068') + 
  ylab('SES of Mean Pairwise Functional Distance') +
  xlab('SES of Mean Pairwise Phylogenetic Distance') +
  guides(col=guide_legend(title= 'Ab. Latitude')) +
  xlim(c(-11, 6.5)) +
  ylim(c(-5, 5)) +
  theme_classic() +
  theme(legend.position='bottom')

p0 + facet_wrap(~ lat_cat, nrow = 4)

#ggsave('sesMPD_sesMFPD.pdf', scale = 0.8)

#plot for SES PD and MPFD
p0.1 <- ggplot(df, aes(x = PDses, y = MPFDses, col = abs_lat)) +
  geom_point() +
  geom_smooth(method = 'lm', color = '#78C0E0') +
  scale_color_gradient2(low = col[1], mid = col[2], high = col[3]) +
  ylab('SES of Mean Pairwise Functional Distance') +
  xlab('SES of Phylogenetic Diversity') +
  guides(col=guide_legend(title= 'Ab. Latitude')) +
  xlim(c(-11, 6.5)) +
  ylim(c(-5, 5)) +
  theme_classic() +
  theme(legend.position='bottom')

p0.1 + facet_wrap(~ lat_cat, nrow = 4)



#now we are going to look at the strength of the correlation between our variables
#of interest but stratify over latitude and longitude

df_null$abs_lat <- rep(df$abs_lat, max(as.numeric(df_null$run)))

df_emp$lat_cat <- cut(df_emp$abs_lat, breaks = c(0, 20, 40, 60, 90))
df_null$lat_cat <- cut(df_null$abs_lat, breaks = c(0, 20, 40, 60, 90))
df$lat_cat <- cut(df$abs_lat, breaks = c(0,20,40,60,90))

df_null$long <- rep(df$long, max(as.numeric(df_null$run)))

df_emp$long_cats <- cut(df_emp$long, breaks = c(-180, -30, 65, 180))
df_null$long_cats <- cut(df_null$long, breaks = c(-180, -30, 65, 180))
levels(df_null$long_cats) <- c('Americas', 'Africa and Europe', 'Oceania and Asia')
levels(df_emp$long_cats) <- c('Americas', 'Africa and Europe', 'Oceania and Asia')

df_cor <- df_null %>%
  group_by(lat_cat, run) %>%
  summarise(PD_MPFD = cor(PD, MPFD), 
            MPD_MPFD = cor(MPD, MPFD),
            PD_SR = cor(PD, sp_richn),
            MPD_SR = cor(MPD, sp_richn), 
            MPFD_SR = cor(MPFD, sp_richn))
df_cor

require(ggplot2)
p1 <- ggplot(df_cor, aes(x = MPD_MPFD))+
  geom_histogram() + theme_classic()
p1 <- p1 +facet_wrap(~ lat_cat)
p1

emp_cor <- df_emp %>%
  group_by(lat_cat) %>%
  summarise(PD_MPFD = cor(MPFD, pd),
            MPD_MPFD = cor(mpd, MPFD))

p1 + geom_vline(
  data    = emp_cor,
  mapping = aes(xintercept = MPD_MPFD),
  linetype="dotted"
)

hist(df_cor$MPD_MPFD, xlim = c(0, 0.16))
abline(v = cor(df_emp$mpd, df_emp$MPFD))


df_cor_long <- df_null %>%
  group_by(lat_cat, run, long_cats) %>%
  summarise(PD_MPFD = cor(PD, MPFD),
            MPD_MPFD = cor(MPD, MPFD),
            PD_SR = cor(PD, sp_richn),
            MPD_SR = cor(MPD, sp_richn), 
            MPFD_SR = cor(MPFD, sp_richn))

colours <- c('#FE5D26', '#340068', '#161717')

#Correlation of PD and MPFD
require(ggridges)
p2 <- ggplot(df_cor_long, aes(x = PD_MPFD, y = long_cats, fill = long_cats))+
  #geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_density_ridges(aes(scale = 7), alpha = 0.7, color = 'white') + 
  scale_fill_manual(values = rev(colours)) +
  xlab(label = beta~of~PD~and~MPFD)+
  ylab('Absolute latitude') +
  theme_classic() +  
  xlim(c(-0.5, 1)) +
  theme(legend.position="bottom", axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.title =element_blank())
p2
labels <- setNames(c('0-20', '20-40', '40-60', '60-90'), nm = levels(df_cor_long$lat_cat))

p2 <- p2 +facet_wrap(.~ lat_cat, nrow = 4, strip.position = 'left', labeller = labeller(lat_cat = labels))
p2



emp_cor_long <- df_emp %>%
  group_by(lat_cat, long_cats) %>%
  summarise(PD_MPFD = cor(MPFD, pd),
            MPD_MPFD = cor(mpd, MPFD))


p2 <- p2 + geom_vline(
  data    = emp_cor_long,
  mapping = aes(xintercept = PD_MPFD, col = long_cats),
  linetype="dashed", alpha = 0.7
) + scale_color_manual(values = rev(colours))
p2


#correlation of MPD and MFPD
p2.1 <- ggplot(df_cor_long, aes(x = MPD_MPFD, y = long_cats, fill = long_cats))+
 # geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_density_ridges(aes(scale = 7), alpha = 0.7) + 
  scale_fill_manual(values = colours) +
  xlab(expression(paste(beta, ' MPD and MFPD')))+
  ylab(element_blank()) +
  theme_classic() +
  xlim(c(-1, 1)) +
  theme(legend.position="none", axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.title =element_blank())
p2.1

new_labels <- c('0-20', '20-40', '40-60', '60-90')
names(new_labels) <- levels(df_cor_long$lat_cat)
new_labels
p2.1 <- p2.1 +facet_wrap(~ lat_cat, nrow = 4, strip.position = 'right', labeller = labeller(lat_cat = new_labels) )

p2.1
p2.1 <- p2.1 + geom_vline(
  data    = emp_cor_long,
  mapping = aes(xintercept = MPD_MPFD, col = long_cats),
  linetype="dashed"
) + scale_color_manual(values = colours)
p2.1

dev.off()
p1.slopes <- ggplot(data = df_cor_long, aes(line = lat_cat))+
  geom_abline( data = df_cor_long,
  mapping = aes(intercept = 0, slope = df_cor_long$PD_MPFD, col = long_cats), alpha = 0.01) +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  scale_color_manual(values = colours) +
  theme_classic()

p1.slopes <- p1.slopes + facet_wrap(~lat_cat, labeller = labeller(lat_cat = new_labels) )
p1.slopes
p1.slopes + geom_abline(data = emp_cor_long,
                        mapping = aes(intercept = 0, slope = PD_MPFD, alpha = 1, col = long_cats)
)

lm <- 


require(egg)
col_plot <- ggarrange(p2, p2.1,
                       widths = c(2,2), nrow = 1, ncol = 2)
col_plot


world <- map_data("world")



ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region)
  ) +   coord_map("moll") +
  geom_point(data = df, aes(x = long, y = lat, col = PDses), size = 1) +
  scale_color_viridis()


library(tidyverse)
library(ggthemes)

world_map = map_data("world") %>% 
  filter(! long > 180)

colours <- c('#340068', rgb(235,221,255, maxColorValue = 255), '#FE5D26')
colour_breaks <- c(-8, 0, 8)


mapPDses <- 
  ggplot() +
              expand_limits(x = world_map$long, y = world_map$lat) +
              coord_map("moll") +
              theme_map() +
              geom_point(data = df, aes(x = long, y = lat, col = PDses), size = .6, alpha = 1) +
  scale_colour_gradientn(
    limits  = range(df$PDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$PDses)), 1),
    name = expression(PD[SES])
  )+
    geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
    geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none")     +
  ggtitle(bquote(PD[SES]))     
    #scale_color_viridis(alpha = 0.3, option = 'D', limits = c(-7, 7))
mapPDses
ggsave('map_PDses.pdf', device = 'pdf')

functional_indices2 <- read.csv('~/Dropbox/AVONET/functional_diversity_indicesFULL.csv')

mapPD <- 
  ggplot() +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  theme_map() +
  geom_point(data = functional_indices2, aes(x = long, y = lat, col = (pd)), size = .6, alpha = 1) +
  scale_colour_gradientn(
    limits  = range((functional_indices2$pd)),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range((functional_indices2$pd))), 1),
    name = "PD")+
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none")     +
  ggtitle('PD')     
mapPD

df$long_cats <- df_emp$long_cats

hist_PDses <- ggplot() +
  geom_histogram(data = df, aes(x = PDses, fill = long_cats), alpha = 1, col = 'black') +
  scale_fill_manual(values = col) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 'dashed')
hist_PDses
ggsave('hist_PDses.pdf', device = 'pdf')

cor_SR_PDses <- 
  ggplot() +
  geom_point(data = df, aes(y = PDses, x = sp_richn, col = PDses), size = .4) +
  #scale_color_gradient2(low = col[1], mid = col[2], high = col[3]) +
  scale_colour_gradientn(
    limits  = range(df$PDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$PDses)), 1),
    name = expression(PD[SES])
  )+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2)+
  xlab('Species richness') +
  ylab(element_blank())+
  ylim(c(-10,10))
cor_SR_PDses <- cor_SR_PDses + facet_wrap(~long_cats, nrow = 3)
cor_SR_PDses
ggsave('cor_SR_PDses.pdf', device = 'pdf')


mapMPFDses <- 
  ggplot() +
    expand_limits(x = world_map$long, y = world_map$lat) +
    coord_map("moll") +
    theme_map() +
    geom_point(data = df, aes(x = long, y = lat, col = MPFDses), size = .6, alpha = 1) +
  scale_colour_gradientn(
    limits  = range(df$MPFDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$MPFDses)), 1),
    name = expression(MPFD[SES])
    )+
    geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
    geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none") +
  ggtitle(bquote(MPFD[SES]))          
mapMPFDses
ggsave('map_MPFDses.pdf', device = 'pdf')



mapFD <- 
  ggplot() +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  theme_map() +
  geom_point(data = functional_indices2, aes(x = long, y = lat, col = (fdis)), size = .6, alpha = 1) +
  scale_colour_gradientn(
    limits  = range((functional_indices2$fdis)),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range((functional_indices2$fdis))), 1),
    name = "FD")+
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none") +
  ggtitle('FD')          
mapFD

hist_MPFDses <- ggplot() +
  geom_histogram(data = df, aes(x = MPFDses, fill = long_cats), alpha = 1, col = 'black') +
  scale_fill_manual(values = col) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 'dashed')
hist_MPFDses
ggsave('hist_MPFDses.pdf', device = 'pdf')


cor_SR_MPFDses <- 
  ggplot() +
  geom_point(data = df, aes(y = MPFDses, x = sp_richn, col = MPFDses), size = .4) +
  scale_colour_gradientn(
    limits  = range(df$MPFDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$MPFDses)), 1),
    name = expression(MPFD[SES])
    
  )+
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2)+
 # theme(legend.position = "none")+
  xlab('Species richness') +
  ylab(element_blank())+
  ylim(c(-10,10))

cor_SR_MPFDses <- cor_SR_MPFDses + facet_wrap(~long_cats, nrow = 3)
cor_SR_MPFDses
ggsave('cor_SR_MPFDses.pdf', device = 'pdf')



functional_indices2$lat_cats <- df$lat_cat
functional_indices2$long_cats <- df$long_cats

cor_SR_FD <- 
  ggplot() +
  geom_point(data = functional_indices2, aes(y = (fdis), x = sp_richn, col = (fdis)), size = .4) +
  scale_colour_gradientn(
    limits  = range((functional_indices2$fdis)),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range((functional_indices2$fdis))), 1),
    name = "FD")+
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2)+
  # theme(legend.position = "none")+
  xlab('Species richness') +
  ylab(element_blank())#+
  #ylim(c(-6,6))

cor_SR_FD <- cor_SR_FD + facet_wrap(~long_cats, nrow = 3)
cor_SR_FD


mapMPDses <-   
  ggplot() +
    expand_limits(x = world_map$long, y = world_map$lat) +
    coord_map("moll") +
    theme_map() +
    geom_point(data = df, aes(x = long, y = lat, col = MPDses), size = .6, alpha = 1) +
  scale_colour_gradientn(
    limits  = range(df$MPDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$MPDses)), 1),
  )+
    geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
    geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none") +
  ggtitle(bquote(MPD[SES]))
mapMPDses
ggsave('map_MPDses.pdf', device = 'pdf')
  
hist_MPDses <- ggplot() +
  geom_histogram(data = df, aes(x = MPDses, fill = long_cats), alpha = 1, col = 'black') +
  scale_fill_manual(values = col) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 'dashed')
hist_MPDses  
ggsave('hist_MPDses.pdf', device = 'pdf')


cor_SR_MPDses <- 
  ggplot() +
  geom_point(data = df, aes(y = MPDses, x = sp_richn, col = MPDses), size = .4) +
  scale_colour_gradientn(
    limits  = range(df$MPDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$MPDses)), 1),
  )+  theme_classic() + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2) +
#  theme(legend.position = "none") +
  xlab('Species richness') +
  ylab(element_blank()) +
  ylim(c(-10,10))
cor_SR_MPDses <- cor_SR_MPDses + facet_wrap(~long_cats, nrow = 3)
cor_SR_MPDses
ggsave('cor_SR_MPDses.pdf', device = 'pdf')

cor_SR_PD <- 
  ggplot() +
  geom_point(data = functional_indices2, aes(y = (pd), x = sp_richn, col = (pd)), size = .4) +
  scale_colour_gradientn(
    limits  = range((functional_indices2$pd)),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range((functional_indices2$pd))), 1),
    name = "PD")+
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2)+
  # theme(legend.position = "none")+
  xlab('Species richness') +
   ylab(element_blank())
  # ylim(c(-4,4))

cor_SR_PD <- cor_SR_PD + facet_wrap(~long_cats, nrow = 3)
cor_SR_PD

alt <- read.csv('~/Dropbox/AVONET/bioclim_data.csv')
df$altitude <- alt$altitude_mean

colour_breaks <- c(-50, 0, 5500)
p7 <- ggplot(data = df, aes(x = PDses, y = MPFDses, col = altitude)) +
  geom_point( alpha = 0.8, size = .3) +
  theme_classic() +
  geom_smooth(method = 'lm', col = 'darkgrey',  lwd = .5) +
  scale_colour_gradientn( colours = colours, name = 'Altitude'
  ) +
  geom_rug(alpha = .2) +
  xlab(bquote(PD[SES])) +
  ylab(bquote(MPFD[SES])) 
p7
p7 <- p7 + facet_wrap(~lat_cat, labeller = labeller(lat_cat = labels))
p7

alt_breaks <- c(-50, 500, 1000, 2000, 3000, 6000)
alt_labels <- c("-50-500masl", "500-1000masl", "1000-2000masl", "2000-3000masl", ">3000masl")
df$alt_cat <- cut(df$altitude, breaks = alt_breaks)
names(alt_labels) <- levels(df$alt_cat)

p8 <- ggplot(data = df[which(is.na(df$altitude) == F),], aes(x = PDses, y = MPFDses, col = abs_lat)) +
  geom_point( alpha = 0.8, size = .3) +
  theme_classic() +
  geom_smooth(method = 'lm', col = 'darkgrey',  lwd = .5) +
  scale_colour_gradientn( colours = colours, name = 'Latitude'
  ) +
  geom_rug(alpha = .2) +
  xlab(bquote(PD[SES])) +
  ylab(bquote(MPFD[SES])) 
p8
p9 <- p8 + facet_wrap(~long_cats + alt_cat, labeller = labeller(alt_cat = alt_labels), ncol = 5, drop = F)
p8 <- p8 + facet_wrap(~alt_cat, labeller = labeller(alt_cat = alt_labels), nrow = 1)
p8
ggsave(plot = p8, filename = 'SESscores_AltLat.pdf', device = 'pdf', height = 5, width = 7)
ggsave(plot = p9, filename = 'SESscores_AltLatLong.pdf', device = 'pdf', height = 5, width = 7)

functional_indices2 <- read.csv('~/Dropbox/AVONET/functional_diversity_indicesFULL.csv', row.names = 1)

transformations <- 
  list(mean_pd = mean(functional_indices2$pd),
       sd_pd = sd(functional_indices2$pd), 
       mean_fd = mean(functional_indices2$fdis),
       sd_fd = sd(functional_indices2$fdis))
transformations

trans_pd <- function(x){(x*transformations$sd_pd)+transformations$mean_pd}
trans_fd <- function(x){(x*transformations$sd_fd)+transformations$mean_fd}

p10 <- ggplot(data = df[which(is.na(df$altitude) == F),], aes(x = trans_pd(functional_indices2$pd[which(is.na(df$altitude) == F)]),
                                                              y = trans_fd(functional_indices2$fdis[which(is.na(df$altitude) == F)]), col = abs_lat)) +
  geom_point( alpha = 0.8, size = .3) +
  theme_classic() +
  geom_smooth(method = 'lm', col = 'darkgrey',  lwd = .5) +
  scale_colour_gradientn( colours = colours, name = 'Latitude'
  ) +
 # scale_fill_grey(start = 0.1, end = 0.9, name = 'Latitude') +
  geom_rug(alpha = .2) +
  xlab('PD') +
  ylab('FD') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1))
p10

p10 <- p10 + facet_wrap(~alt_cat, labeller = labeller(alt_cat = alt_labels), nrow = 1)
p10
ggsave(plot = p10, filename = 'FDPDscores_AltLat.pdf', device = 'pdf', height = 5, width = 7)

require(egg)
grid_plot0 <- ggarrange(p10, p8,
          nrow = 2, labels = c('a', 'b'))

grid_plot0
ggsave(plot = grid_plot0, filename = 'scatters_by_alt_lat.pdf', device = 'pdf')


grid_plot <- ggarrange(mapMPFDses, cor_SR_MPFDses,
                       mapPDses, cor_SR_PDses,
                       widths = c(3,1), nrow = 2, ncol = 2, labels = c('a', '', 'b', ''))
grid_plot

ggsave(plot = grid_plot, filename = 'maps_and_dists.pdf', device = 'pdf')


grid_plot1 <- ggarrange(mapFD, cor_SR_FD,
                       mapPD, cor_SR_PD,
                       widths = c(3,1), nrow = 2, ncol = 2, labels = c('a', '', 'b', ''))
grid_plot1
ggsave(plot = grid_plot1, filename = 'FDPDmaps_and_dists.pdf', device = 'pdf')


####################################
### Structual Equation Modelling ###
####################################
require(lavaan)

df_standardise <- df[c(13,15:18)]
df_standardise$sp_richn2 <- df_standardise$sp_richn^2

df_standardise <- standardize(df_standardise)

model1 <- 'MPFDses ~ sp_richn + PDses + abs_lat + abs_lat:PDses
PDses ~ sp_richn + abs_lat + sp_richn2
sp_richn ~ abs_lat'
fit1 = cfa(model1, data = df_standardise, estimator = 'MLM', se = 'robust')

summary(fit1, fit.measures = T, standardized = T, rsquare = T)
path1 <- semPaths(fit1, 'std', layout = 'tree3', residuals = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4)
m1 <- lm(MPFDses ~ sp_richn + PDses + abs_lat + abs_lat:PDses + sp_richn2, data = df_standardise)
summary(m1)

breaks <- df %>%
  group_by(lat_cat) %>%
  summarise(min = min(abs_lat), max = max(abs_lat), breaks = mean(c(max(abs_lat), min(abs_lat))))
breaks

df_emp$altitude <- df$altitude

p8 <- ggplot(data = df_emp, aes(y = fdis, x = pd, col = abs_lat)) +
  geom_point() +
 # scale_colour_manual(values = colours) +
  geom_smooth(method = 'lm') +
  theme_classic()
p8 + facet_wrap(~ lat_cat)


require(sesem)
require(fields)
require(jtools)
distance_matrix <- round(rdist.earth(x1 = as.matrix(df[,c(12,14)]), miles = F)[lower.tri(matrix(nrow = nrow(df), ncol = nrow(df)), diag=T)])
save(distance_matrix, file = 'distance_matrix.Rdata')
load('distance_matrix.Rdata')


system.time(bins<-make.bin(distance_matrix,type='n.bins',p.dist=10, n.bins = 3))
save(bins, file = 'bins.Rdata')
load('bins.Rdata')
binsize<-bins[1][[1]] #truelove lowland bin sizes
binname<-bins[2][[1]] #truelove lowland bin names
plotbin(distance_matrix,binsize)

colnames(df)

data <- df[,c(12,14,13,15,17,18,21)]
data
data$PDses_int_abs_lat <- data$PDses*data$abs_lat
data$sp_richn2 <- data$sp_richn^2
head(data)
data <- standardize(data, vars = colnames(data)[3:9])


system.time(covariances <- make.covar(data, distance_matrix, binsize, binname))
covariances
save(covariances, file = 'covariances_ses.Rdata')
load('covariances_ses.Rdata')
model1 <- 'MPFDses ~ sp_richn + PDses + abs_lat + PDses_int_abs_lat + altitude
PDses ~ sp_richn + abs_lat + sp_richn2 + altitude
sp_richn ~ abs_lat + altitude'

results <- runModels(model1, covariances)
plotmodelfit(results, rmsea_err=F)
modelsummary(results)
save(results, file = 'results.Rdata')
load('results.Rdata')
results

require(lavaan)
fit1 = cfa(model1, data = data, estimator = 'MLM', se = 'robust')
summary(fit1, standardized = T)
path1 <- semPaths(fit1, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))
fit1.Americas <- cfa(model1, data = subset(data, long >= -180 & long < -30), estimator = 'MLM', se = 'robust')
semPaths(fit1.Americas, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
         pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))
fit1.Africa_Europe <- cfa(model1, data = subset(data, long >= -30 & long < 65), estimator = 'MLM', se = 'robust')
fit1.Asia_Oceania <- cfa(model1, data = subset(data, long >= 65 & long < 180), estimator = 'MLM', se = 'robust')

require(tidyr)
coefficients <- as_tibble(bind_rows(standardizedSolution(fit1)[1:11,],
                                    standardizedSolution(fit1.Americas)[1:11,],
                                    standardizedSolution(fit1.Africa_Europe)[1:11,],
                                    standardizedSolution(fit1.Asia_Oceania)[1:11,])
)

coefficients$section <- as.factor(c(rep('Global', 11), rep('Americas', 11), rep('Europe and Africa', 11), rep('Asia and Oceania', 11)))
coefficients$op <- NULL
#coefficients$exo <- NULL
coefficients$lhs <- as.factor(coefficients$lhs)
coefficients$rhs <- as.factor(coefficients$rhs)
coefficients$rhs
coefficients

colours <- c( "#340068","#FE5D26", "#EBDDFF",  'darkgrey')
labels <- setNames(c(bquote(MPFD[SES]), bquote(PD[SES]), bquote(SR)), nm = levels(coefficients$lhs))
jitter <- position_jitter(width = 0.2, height = 0)

p1 <- ggplot(data = coefficients, aes(x = rhs, y = est.std, col = section, group = rhs)) +
  geom_pointrange(aes(ymin = ci.lower, ymax = ci.upper, xmin = rhs, xmax = rhs, col = section), position = jitter) +
   scale_color_manual(values = colours)+
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  grids(linetype = "dashed") +
  ylab('Coefficient')+
  xlab(element_blank()) +
  scale_x_discrete(labels = c('Lat', 'Alt', bquote(PD[SES]), bquote(PD[SES]:Lat), 'SR', bquote('SR'^2)))
p1
p1 <- p1 +facet_wrap(~lhs, nrow = 1, labeller = labeller(lhs = labels)) + theme(legend.position="bottom") +theme(legend.title = element_blank())
p1

ggsave('continental_coef_ses.pdf', device = 'pdf', width = 12, height = 6)


### Code to plot the slope of MPFDses and PDses against absolute latitude 

colours <- c( "#340068","#FE5D26", "#EBDDFF",  'darkgrey')

models <- list(fit1, fit1.Africa_Europe, fit1.Asia_Oceania, fit1.Americas)
model_names <- c('Global', 'Africa and Europe', 'Asia and Oceania', 'Americas')
DF <- matrix(nrow = 0, ncol = 5)
colnames(DF) <- c('b1', 'sd1', 'UB', 'LB', 'Region')
for(i in 1:length(models)) {
  parm<-as.data.frame(parameterEstimates(models[[i]]))
  parm[,1:3]
  
  COV<-vcov(fit1)
  COVd<-as.data.frame(COV)
  head(COVd)
  
  b1<-parm[2,4]#slope for PD
  b3<-parm[4,4]#slope for PD x abs_lat
  s11<-COV[2,2]#variance for PD
  s13<-COV[4,2]#covariance for PD x abs_lat
  s33<-COV[4,4]#variance for abs_lat
  
  sd1<-c(-2, -1.75, -1.5, -1.25, -1, -.75, -.5, -.25, 0,
         .25,.5,.75,1,1.25,1.5,1.75,2, 2.25, 2.5)
  
  sd2<-sd1*sd(data$MPFDses)
  se<-sqrt(s11+2*(sd2)*s13+(sd2)*(sd2)*s33)
  se
  b1<-b1+b3*(sd2)
  
  b1
  
  UB<-b1+se*qnorm(.975)
  LB<-b1+se*qnorm(.025)
  
  DF.temp<-data.frame(b1, sd1, UB, LB)
  DF.temp$Region <- model_names[i]
  
  DF <- rbind(DF, DF.temp)
}

DF$Region <- factor(DF$Region, ordered = T, levels = c('Americas', 'Asia and Oceania', 'Africa and Europe', 'Global'))

transformations <- 
  list(mean_lat = mean(df$abs_lat),
       sd_lat = sd(df$abs_lat)
  )
transformations
trans_lat <- function(x){(x*transformations$sd_lat)+transformations$mean_lat}

g1<-ggplot(DF)+
  geom_ribbon(aes(ymin=(LB), ymax=(UB), x=trans_lat(sd1), fill = Region), alpha=.7)+
  geom_line(aes(y=b1, x=trans_lat(sd1), col = Region))+
  geom_hline(aes(yintercept=0), lty='dashed')+
  xlab('Absolute latitude')+
  ylab(Standardised~slope~of~MPFD[SES]~-~PD[SES])+
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme_classic() +
  xlim(c(0,90)) +
  theme(legend.position = 'none')
g1

p1g1 <- ggarrange(p1, g1, nrow = 1, widths = c(3,1), labels = c('b', ''))
p1g1
ggsave(p1g1, filename = 'continental_coef_int_slope_ses.pdf', device = 'pdf', width = 12, height =4)

#######################################
### Southern v Northern Hempisphere ###
#######################################

fit1.North <- cfa(model1, data = subset(data, lat >= 0), estimator = 'MLM', se = 'robust')
summary(fit1.North, standardized = T)


fit1.South <- cfa(model1, data = subset(data, lat <= 0), estimator = 'MLM', se = 'robust')
summary(fit1.South, standardized = T)



coefficients <- as_tibble(bind_rows(standardizedSolution(fit1)[1:11,],
                                    standardizedSolution(fit1.North)[1:11,],
                                    standardizedSolution(fit1.South)[1:11,]
))



coefficients$section <- as.factor(c(rep('Global', 11), rep('North', 11), rep('South', 11)))
coefficients$op <- NULL
#coefficients$exo <- NULL
coefficients$lhs <- as.factor(coefficients$lhs)
coefficients$rhs <- as.factor(coefficients$rhs)
coefficients$rhs
tail(coefficients)


colours <- c("darkgrey", "#340068",  "#FE5D26")
pp2 <- ggplot(data = coefficients, aes(x = rhs, y = est.std, col = section, group = rhs)) +
  geom_pointrange(aes(ymin = ci.lower, ymax = ci.upper, xmin = rhs, xmax = rhs, col = section), position = jitter) +
  scale_color_manual(values = colours)+
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  grids(linetype = "dashed") +
  ylab('Coefficient')+
  xlab(element_blank()) +
  scale_x_discrete(labels = c('Lat', 'Alt', bquote(PD[SES]), bquote(PD[SES]:Lat), 'SR', bquote('SR'^2)))
pp2
pp2 <- pp2 +facet_wrap(~lhs, nrow = 1, labeller = labeller(lhs = labels)) + theme(legend.position="bottom") +theme(legend.title = element_blank())
pp2
ggsave('hemisphere_coef_ses.pdf', device = 'pdf', width = 12, height = 6)


models <- list(fit1, fit1.North, fit1.South)
model_names <- c('Global', 'North', 'South')
DF <- matrix(nrow = 0, ncol = 5)
colnames(DF) <- c('b1', 'sd1', 'UB', 'LB', 'Region')
for(i in 1:length(models)) {
  parm<-as.data.frame(parameterEstimates(models[[i]]))
  parm[,1:3]
  
  COV<-vcov(fit1)
  COVd<-as.data.frame(COV)
  head(COVd)
  
  b1<-parm[2,4]#slope for PD
  b3<-parm[4,4]#slope for PD x abs_lat
  s11<-COV[2,2]#variance for PD
  s13<-COV[4,2]#covariance for PD x abs_lat
  s33<-COV[4,4]#variance for abs_lat
  
  sd1<-c(-2, -1.75, -1.5, -1.25, -1, -.75, -.5, -.25, 0,
         .25,.5,.75,1,1.25,1.5,1.75,2, 2.25, 2.5)
  
  sd2<-sd1*sd(data$MPFDses)
  se<-sqrt(s11+2*(sd2)*s13+(sd2)*(sd2)*s33)
  se
  b1<-b1+b3*(sd2)
  
  b1
  
  UB<-b1+se*qnorm(.975)
  LB<-b1+se*qnorm(.025)
  
  DF.temp<-data.frame(b1, sd1, UB, LB)
  DF.temp$Region <- model_names[i]
  
  DF <- rbind(DF, DF.temp)
}

DF$Region <- factor(DF$Region, ordered = T, levels = c('Global', 'North', 'South'))

g2<-ggplot(DF)+
  geom_ribbon(aes(ymin=LB, ymax=UB, x=trans_lat(sd1), fill = Region), alpha=.7)+
  geom_line(aes(y=b1, x=trans_lat(sd1), col = Region))+
  geom_hline(aes(yintercept=0), lty='dashed')+
  xlab('Absolute latitude')+
  xlim(c(0, 90)) +
  ylab(Standardised~slope~of~MPFD[SES]~-~PD[SES])+
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme_classic() +
  theme(legend.position = 'none')
g2

pp2g2 <- ggarrange(pp2, g2, nrow = 1, widths = c(3,1), labels = c('b', ''))
pp2g2
ggsave(pp2g2, filename = 'hemisphere_coef_int_slope_ses.pdf', device = 'pdf', width = 12, height =4)


regions <- ggplot()+
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  theme_map() +
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  geom_vline(xintercept = c(-30, 65), lty = 2, lwd = .3) +
  geom_vline(xintercept = c(-180, 180), lwd = .1)
regions  

ggsave('regions.pdf', device = 'pdf', width = 8, height = 6)
