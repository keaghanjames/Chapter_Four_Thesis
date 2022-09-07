## NULL ANALYSIS

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
require(fields)
require(picante)
require(metricTester)
col <- c('#0A2342', '#0B7A75', '#FA8334', '#7B2D26', '#9593D9')
colours <- c('#340068', 'darkgrey', '#FE5D26')

setwd("NULL")

df_null <- readRDS('summary_run1.rds')

for(i in 2:999){
  df_null <- bind_rows(df_null, readRDS(paste0('summary_run', i, '.rds')))
  print(i)
}
tail(df_null)

saveRDS(df_null, 'full_null.rds')
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

ggsave('sesMPD_sesMFPD.pdf', scale = 0.8)

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

ggsave('sesPD_sesMFPD.pdf', scale = 0.8)


#now we are going to look at the strength of the correlation between our variables
#of interest but stratify over latitude and longitude

df_null$abs_lat <- rep(df$abs_lat, max(as.numeric(df_null$run)))

df_emp$lat_cat <- cut(df_emp$abs_lat, breaks = c(0, 20, 40, 60, 90))
df_null$lat_cat <- cut(df_null$abs_lat, breaks = c(0, 20, 40, 60, 90))

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
write.csv(df_cor, 'simulated_correlations.csv', row.names = F)  

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
# you should subset by absolute latitude

df_cor_long <- df_null %>%
  group_by(lat_cat, run, long_cats) %>%
  summarise(PD_MPFD = cor(PD, MPFD), 
            MPD_MPFD = cor(MPD, MPFD),
            PD_SR = cor(PD, sp_richn),
            MPD_SR = cor(MPD, sp_richn), 
            MPFD_SR = cor(MPFD, sp_richn))
df_cor_long


#Correlation of PD and MPFD
require(ggridges)
p2 <- ggplot(df_cor_long, aes(x = PD_MPFD, y = long_cats, fill = long_cats))+
  #geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_density_ridges(aes(scale = 7), alpha = 0.7) + 
  scale_fill_manual(values = colours) +
  xlab(expression(paste(beta, ' PD and MFPD')))+
  ylab(element_blank()) +
  theme_classic() +  
  xlim(c(-1, 1)) +
  theme(legend.position="bottom", axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.title =element_blank())
p2

p2 <- p2 +facet_wrap(~ lat_cat, nrow = 4, strip.position = 'right') +theme(
  strip.background = element_blank(),
  strip.text.y = element_blank()
)
p2



emp_cor_long <- df_emp %>%
  group_by(lat_cat, long_cats) %>%
  summarise(PD_MPFD = cor(MPFD, pd),
            MPD_MPFD = cor(mpd, MPFD))


p2 <- p2 + geom_vline(
  data    = emp_cor_long,
  mapping = aes(xintercept = PD_MPFD, col = long_cats),
  linetype="dashed"
) + scale_color_manual(values = colours)

ggsave('corPD_MFPD.pdf', scale = 0.8)


#cor MPD and MFPD
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
ggsave('corMPD_MFPD.pdf', scale = 0.8)

require(egg)
col_plot <- ggarrange(p2, p2.1,
                      widths = c(2,2), nrow = 1, ncol = 2)
col_plot

ggsave('simulated_correlations.pdf', plot = col_plot, device = 'pdf')

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

#colours <- c(rgb(74,150,236, maxColorValue = 255), rgb(231,239,242, maxColorValue = 255), rgb(241,166,127, maxColorValue = 255))
colours <- c('#340068', rgb(235,221,255, maxColorValue = 255), '#FE5D26')

colour_breaks <- c(-8, 0, 8)
#scale_color_gradient(low = '#FE5D26', high = '#340068') + 
#colours[c(1,3)] <- c('#340068', '#FE5D26')

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
  )+
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none")     +
  ggtitle(bquote(PD[SES]))     
#scale_color_viridis(alpha = 0.3, option = 'D', limits = c(-7, 7))
mapPDses
ggsave('map_PDses.pdf', device = 'pdf')

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
  )+
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2)+
  theme(legend.position = "none")+
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
    limits  = range(df$PDses),
    colours = colours[c(1, seq_along(colours[1:3]), length(colours[1:3]))],
    values  = c(0, scales::rescale(colour_breaks, from = range(df$PDses)), 1),
  )+
  geom_map(map = world_map, data = world_map, aes(long, lat, map_id = region), fill = rgb(240,240,240, 1,maxColorValue = 255), col = 'black', lwd = .2) +
  geom_hline(yintercept = 0, linetype="dashed", size = .2)+
  theme(legend.position = "none") +
  ggtitle(bquote(MPFD[SES]))          
mapMPFDses
ggsave('map_MPFDses.pdf', device = 'pdf')


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
  )+
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2)+
  theme(legend.position = "none")+
  xlab('Species richness') +
  ylab(element_blank())+
  ylim(c(-10,10))

cor_SR_MPFDses <- cor_SR_MPFDses + facet_wrap(~long_cats, nrow = 3)
cor_SR_MPFDses
ggsave('cor_SR_MPFDses.pdf', device = 'pdf')

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
  theme(legend.position = "none") +
  xlab('Species richness') +
  ylab(element_blank()) +
  ylim(c(-10,10))
cor_SR_MPDses <- cor_SR_MPDses + facet_wrap(~long_cats, nrow = 3)
cor_SR_MPDses
ggsave('cor_SR_MPDses.pdf', device = 'pdf')

require(egg)
grid_plot <- ggarrange(mapMPDses, cor_SR_MPDses,
                       mapPDses, cor_SR_PDses,
                       mapMPFDses, cor_SR_MPFDses,
                       widths = c(2,1), nrow = 3, ncol = 2)
grid_plot

ggsave(plot = grid_plot, filename = 'maps_and_dists.pdf', device = 'pdf')

