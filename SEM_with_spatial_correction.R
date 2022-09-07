#this code implements the structual equation models but includes versions which correct for spatial autocorrelation


library(lavaan)
library(semPlot)
library(OpenMx)
library(GGally)
library(corrplot)



setwd("~/Dropbox/AVONET")
require(tidyr)
require(stringr)
require(phytools)
require(dplyr)
require(caper)
require(mFD)
library(foreach)
library(doParallel)
require(interactions)
require(jtools)

require(ggdag)
require(ggplot2)
theme_set(theme_dag())
colours <- c( "#FE5D26", "#340068", "#EBDDFF", 'darkgrey')

#lets quickly make a dag to establish what our hypothesized causal chain is 
dagify(Sp.R ~ Lat,
       PD ~ Sp.R,
       FDis ~ PD,
       FDis ~ Sp.R,
       FDis ~ Lat,
       PD ~ Lat) %>% 
  ggdag()

#so we need to find a way to isolate the effects of our three predictors, even though they all impact each other

#lets investigate whether PD is capable of predicting FD in birds
require(nlme)
require(visreg)
require(dplyr)
require(car)
require(performance)
require(GGally)

functional_indices <- read.csv('functional_diversity_indicesFULL.csv', row.names = 1)
head(functional_indices)
functional_indices <- functional_indices[,-c(3:8)]
functional_indices$abs_lat <- abs(functional_indices$lat)
functional_indices$lat_cat <- cut(functional_indices$abs_lat, breaks = c(0, 20, 40, 60, 90))
functional_indices$long_cats <- cut(functional_indices$long, breaks = c(-180, -30, 65, 180))
levels(functional_indices$long_cats) <- c('Americas', 'Africa and Europe', 'Oceania and Asia')

bioclim <- read.csv('bioclim_data.csv')
functional_indices$altitude <- bioclim$altitude_mean

head(functional_indices)


p1 <- ggplot(data = functional_indices, aes(x = sp_richn, y = pd, col = pd)) +
  geom_point() +
  geom_smooth() +
  theme_classic()
p1

ggsave('pd_rich.pdf')

p2 <- ggplot(data = functional_indices, aes(x = pd, y = fdis, col = abs_lat)) +
  geom_point() +
  geom_smooth(method  = 'lm') +
  theme_classic()
p2

ggsave('fdis_pd_lat.pdf')

p3 <- ggplot(data = functional_indices, aes(x = sp_richn, y = fdis, col = abs_lat)) +
  geom_point() +
  geom_smooth(method  = 'lm') +
  theme_classic()
p3
ggsave('fdis_rich_lat.pdf')



p4 <- ggplot(data = functional_indices, aes(x = sp_richn, y = fdis,alpha = 0.2)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic()
p4
ggsave('pd_rich_split.pdf')


p5 <- ggplot(data = functional_indices, aes(x = altitude, y = fdis, col = abs_lat)) +
  geom_point() +
  geom_smooth(method  = 'lm') +
  theme_classic()
p5

#lets drop the pincriple components

functional_indices$sp_richn2 <- functional_indices$sp_richn^2
functional_indices$abs_lat2 <- functional_indices$abs_lat^2

functional_indices <- standardize(functional_indices, vars = c("fdis","sp_richn", "pd", "mpd", "abs_lat", 
                                                               "sp_richn2", "abs_lat2", "altitude"))


model0 <- 'fdis ~ 1
pd ~ 1
sp_richn ~ 1'
fit0 = cfa(model0, data = functional_indices, estimator = 'MLM', se = 'robust')
summary(fit0, fit.measures = T, standardized = T, rsquare = T)
path0 <- semPaths(fit0, 'std', layout = 'tree3', residuals = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4)

model1 <- 'fdis ~ sp_richn + pd + abs_lat + abs_lat:pd
pd ~ sp_richn + abs_lat + sp_richn2 + altitude
sp_richn ~ abs_lat + altitude'
fit1 = cfa(model1, data = functional_indices, estimator = 'MLM', se = 'robust')
summary(fit1, fit.measures = T, standardized = T, rsquare = T)
path1 <- semPaths(fit1, 'std', layout = 'tree3', residuals = F, curve = 2, curvature = 1, nCharNodes = 0,
         pastel = F, style = 'lisrel', rotation = 4, edge.color = rev(c('#FE5D26', '#340068')))


model2 <- 'fdis ~ sp_richn + pd + abs_lat + abs_lat:pd + altitude
pd ~ sp_richn + abs_lat + sp_richn2 + altitude
sp_richn ~ abs_lat + altitude'


fit2 = cfa(model2, data = functional_indices, estimator = 'MLM', se = 'robust')
summary(fit2, standardized = T)
path2 <- semPaths(fit2, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))



model3 <- 'fdis ~ sp_richn + pd + abs_lat + abs_lat:pd + altitude + altitude:pd
pd ~ sp_richn + abs_lat + sp_richn2 + altitude
sp_richn ~ abs_lat + altitude'


fit3 = cfa(model3, data = functional_indices, estimator = 'MLM', se = 'robust')
summary(fit3, standardized = T)
path3 <- semPaths(fit3, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))


model4 <- 'fdis ~ sp_richn + pd + abs_lat + abs_lat:pd + altitude + altitude:pd + altitude:sp_richn
pd ~ sp_richn + abs_lat + sp_richn2 + altitude + altitude:sp_richn
sp_richn ~ abs_lat + altitude'


fit4 = cfa(model4, data = functional_indices, estimator = 'MLM', se = 'robust')
summary(fit4, standardized = T)
path4 <- semPaths(fit4, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))


model5 <- 'fdis ~ sp_richn + pd + abs_lat + abs_lat:pd + altitude:pd + altitude:sp_richn
pd ~ sp_richn + abs_lat + sp_richn2 + altitude
sp_richn ~ abs_lat + altitude'


fit5 = cfa(model5, data = functional_indices, estimator = 'MLM', se = 'robust')
summary(fit5, standardized = T)
path5 <- semPaths(fit5, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))



AIC(fit1, fit2, fit3, fit4, fit5)

m1 <- lm(fdis ~ sp_richn + sp_richn2 + pd + abs_lat + abs_lat:pd + altitude, data = functional_indices)
summary(m1)

#m1.spatial <- gls(fdis ~ sp_richn + sp_richn2 + pd + abs_lat + abs_lat:pd, data = functional_indices, 
#                  correlation = corExp(form=~lat+long))



# dist <- as.matrix(dist(cbind(functional_indices$lat, functional_indices$long)))
# dist.inv <- 1/dist
# diag(dist.inv) <- 0
# Moran.I(functional_indices$fdis, dist.inv)

fit2.Americas <- cfa(model2, data = subset(functional_indices, long >= -180 & long < -30), estimator = 'MLM', se = 'robust')
summary(fit2.Americas, standardized = T)
semPaths(fit2.Americas, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
         pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F, edge.color = rev(c('#FE5D26', '#340068')))
fit2.Africa_Europe <- cfa(model2, data = subset(functional_indices, long >= -30 & long < 65), estimator = 'MLM', se = 'robust')
summary(fit2.Africa_Europe, standardized = T)

fit2.Asia_Oceania <- cfa(model2, data = subset(functional_indices, long >= 65 & long < 180), estimator = 'MLM', se = 'robust')
summary(fit2.Asia_Oceania, standardized = T)

require(tidyr)
coefficients <- as_tibble(bind_rows(standardizedSolution(fit2)[1:11,],
                                    standardizedSolution(fit2.Americas)[1:11,],
                                    standardizedSolution(fit2.Africa_Europe)[1:11,],
                                    standardizedSolution(fit2.Asia_Oceania)[1:11,])
)

coefficients$section <- as.factor(c(rep('Global', 11), rep('Americas', 11), rep('Europe and Africa', 11), rep('Asia and Oceania', 11)))
coefficients$op <- NULL
#coefficients$exo <- NULL
coefficients$lhs <- as.factor(coefficients$lhs)
coefficients$rhs <- as.factor(coefficients$rhs)
coefficients$rhs
coefficients
m1.Americas <- update(m1, .~., data = subset(functional_indices, long >= -180 & long < -30))
m1.Europe_Africa <- update(m1, .~., data = subset(functional_indices, long >= -30 & long < 65))
m1.Asia_Oceania <- update(m1, .~., data = subset(functional_indices, long >= 65 & long < 180))
require(broom)
m1.1 <- rbind(tidy(m1, conf.int = T),
              tidy(m1.Americas, conf.int = T),
              tidy(m1.Europe_Africa, conf.int = T),
              tidy(m1.Asia_Oceania, conf.int = T)
              
              
              )
colnames(m1.1) <- c('rhs', 'est.std', 'se', 'z', 'pvalue', 'ci.lower', 'ci.upper')
m1.1$lhs <- 'MR(fdis)'
m1.1 <- m1.1[which(m1.1$rhs != '(Intercept)'),]
m1.1$rhs[which(m1.1$rhs == 'pd:abs_lat')] <- 'abs_lat:pd'
m1.1$section <- as.factor(c(rep('Global', 6), rep('Americas', 6), rep('Europe and Africa', 6), rep('Asia and Oceania', 6)))
m1.1
# coefficients <- bind_rows(coefficients, m1.1)
# coefficients
# coefficients$lhs <- factor(coefficients$lhs, ordered = T, levels = c('fdis', 'pd', 'sp_richn'))
# coefficients$rhs <- factor(coefficients$rhs, ordered = T, levels = c('abs_lat', 'altitude', 'pd', 'abs_lat:pd', 'sp_richn', 'sp_richn2'))
col <- c('#0A2342', '#0B7A75', '#FA8334', '#7B2D26', '#9593D9')

labels <- setNames(c('FD', 'PD', 'SR'), nm = levels(coefficients$lhs))
require(ggpubr)
require(egg)
jitter <- position_jitter(width = 0.2, height = 0)

#labels <- setNames(c('Lat', 'Lat:PD', 'Alt', 'PD', 'SR', bquote('SR'^2)), nm = levels(coefficients$rhs))
#labels
pp1 <- ggplot(data = coefficients, aes(x = rhs, y = est.std, col = section, group = rhs)) +
  geom_pointrange(aes(ymin = ci.lower, ymax = ci.upper, xmin = rhs, xmax = rhs, col = section), position = jitter) +
  scale_color_manual(values = colours)+
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  grids(linetype = "dashed") +
  ylab('Coefficient')+
  xlab(element_blank()) +
  scale_x_discrete(labels = c('Lat', 'Lat:PD','Alt', 'PD', 'SR', bquote('SR'^2)))
pp1
pp1 <- pp1 +facet_wrap(~lhs, nrow = 1, labeller = labeller(lhs = labels)) + theme(legend.position="bottom") +theme(legend.title = element_blank())
pp1
ggsave('continental_coef.pdf', device = 'pdf', width = 12, height = 6)


semPaths(fit2.Americas, 'std', layout = 'tree3', residuals = T, intercepts = F, curve = 2, curvature = 1, nCharNodes = 0,
         pastel = F, style = 'lisrel', rotation = 4, fade = F, exoCov = F)

head(functional_indices)

breaks <- functional_indices %>%
  group_by(lat_cat) %>%
  summarise(min = min(abs_lat), max = max(abs_lat), breaks = mean(c(max(abs_lat), min(abs_lat))))
breaks


#we want to untransform our data 
functional_indices2 <- read.csv('functional_diversity_indicesFULL.csv', row.names = 1)

transformations <- 
  list(mean_pd = mean(functional_indices2$pd),
            sd_pd = sd(functional_indices2$pd), 
            mean_fd = mean(functional_indices2$fdis),
            sd_fd = sd(functional_indices2$fdis),
            mean_lat = mean(df$abs_lat),
            sd_lat = sd(df$abs_lat)
)
transformations
trans_lat <- function(x){(x*transformations$sd_lat)+transformations$mean_lat}

trans_pd <- function(x){(x*transformations$sd_pd)+transformations$mean_pd}
trans_fd <- function(x){(x*transformations$sd_fd)+transformations$mean_fd}

colours <- scales::seq_gradient_pal("#FE5D26", "lightblue", "Lab")(seq(0,1,length.out=4))
p1 = visreg(m1.Americas, 'pd', 'abs_lat', data = subset(functional_indices, long >= -180 & long < -30), 
            breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
            line = list(col = alpha(colours, 0.7), lty = 1:4), 
            fill = list(col = alpha(colours, 0.4)),
            points = list(col = alpha(colours, 0.4)),
            xlab = 'PD', ylab = 'FD',
            strip.names=c('0-20', '20-40', '40-60', '60-90'), 
            xlim = c(0, 14200),
            ylim = c(-0.05, .2)
)

dev.off()



p.AM <- ggplot(p1$fit, aes(pd, visregFit, linetype=factor(abs_lat), fill=factor(abs_lat))) +
  
      geom_point(data = p1$res, aes(pd, visregRes, col = factor(abs_lat)), alpha = 0.5, size = .5) +
      geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, 
              linetype=1, size=0.2, outline.type = 'both') +
      scale_fill_manual(values = col)+
      scale_color_manual(values = col)+
      geom_line() +
      ylab("FD") +
      xlab("PD") +
      ylim(c(-10, 10))+
  xlim(c(-2, 4)) +
  theme_classic() +theme(legend.title = element_blank())+
      ggtitle('Americas')
p.AM

p2 <- visreg(m1.Europe_Africa, 'pd', 'abs_lat', data = subset(functional_indices, long >= -30 & long < 65), 
             breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
             line = list(col = alpha(colours, 0.4)), 
             fill = list(col = alpha(colours, 0.4)),
             points = list(col = alpha(colours, 0.4)),
             xlab = 'PD', ylab = 'FD',
             strip.names=c('0-20', '20-40', '40-60', '60-90'), 
             xlim = c(0, 14200),
             ylim = c(-0.05, .2)
)



p.EA <- ggplot(p2$fit, aes(pd, visregFit, linetype=factor(abs_lat), fill=factor(abs_lat))) +
  
  geom_point(data = p2$res, aes(pd, visregRes, col = factor(abs_lat)), alpha = 0.5, size = .5) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, 
              linetype=1, size=0.2) +
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)+
  geom_line() +
  ylab("FD") +
  xlab("PD") +
  ylim(c(-10, 10))+
  xlim(c(-2, 4)) +
  theme_classic() +theme(legend.title = element_blank())+
  ggtitle('Africa and Europe')
p.EA


p3 <-  visreg(m1.Asia_Oceania, 'pd', 'abs_lat', data = subset(functional_indices, long >= 65 & long < 180), 
              breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
              line = list(col = alpha(colours, 0.4)), 
              fill = list(col = alpha(colours, 0.4)),
              points = list(col = alpha(colours, 0.4)),
              xlab = 'PD', ylab = 'FD',
              strip.names=c('0-20', '20-40', '40-60', '60-90'), 
              xlim = c(0, 14200),
              ylim = c(-0.05, .2),
              main = 'Oceania and Asia', 
              legend = F
)



p.OA <- ggplot(p3$fit, aes(pd, visregFit, linetype=factor(abs_lat), fill=factor(abs_lat))) +
  
  geom_point(data = p3$res, aes(pd, visregRes, col = factor(abs_lat)), alpha = 0.5, size = .5) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, 
              linetype=1, size=0.2) +
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)+
  geom_line() +
  ylab("FD") +
  xlab("PD") +
  ylim(c(-10, 10))+
  xlim(c(-2, 4)) +
  theme_classic() +theme(legend.title = element_blank())+
  ggtitle('Oceania and Asia')
p.OA


colours <- c( "#FE5D26", "#340068", "#EBDDFF", 'darkgrey')

pdf(file = "partial_resid.pdf",   
    width = 12, # The width of the plot in inches
    height = 4, 
    pointsize = 12) # The height of the plot in inches
par(mfrow = c(1,4))
p0 <- visreg(m1, 'pd', 'abs_lat', data = functional_indices, 
             breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
             line = list(col = alpha(colours, 0.7), lwd = 1), 
             fill = list(col = alpha(colours, 0.4)),
             points = list(col = alpha(colours, 0.4), cex = 0.2),
             xlab = 'PD', ylab = 'FD',
             strip.names=c('0-20', '20-40', '40-60', '60-90'), 
             xlim = c(0, 14200),
             ylim = c(-0.05, .2),
             main = 'Global', 
             legend = T
)
p1 = visreg(m1.Americas, 'pd', 'abs_lat', data = subset(functional_indices, long >= -180 & long < -30), 
            breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
            line = list(col = alpha(colours, 0.4), lwd = 1), 
            fill = list(col = alpha(colours, 0.4)),
            points = list(col = alpha(colours, 0.4), cex = 0.2),
            xlab = 'PD', ylab = 'FD',
            strip.names=c('0-20', '20-40', '40-60', '60-90'), 
            xlim = c(0, 14200),
            ylim = c(-0.05, .2),
            main = 'Americas', 
            legend = F
)
p2 <- visreg(m1.Europe_Africa, 'pd', 'abs_lat', data = subset(functional_indices, long >= -30 & long < 65), 
             breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
             line = list(col = alpha(colours, 0.4), lwd = 1), 
             fill = list(col = alpha(colours, 0.4)),
             points = list(col = alpha(colours, 0.4), cex = 0.2),
             xlab = 'PD', ylab = 'FD',
             strip.names=c('0-20', '20-40', '40-60', '60-90'), 
             xlim = c(0, 14200),
             ylim = c(-0.05, .2),
             main = 'Africa and Europe', 
             legend = F
)
p3 <-  visreg(m1.Asia_Oceania, 'pd', 'abs_lat', data = subset(functional_indices, long >= 65 & long < 180), 
              breaks = breaks$breaks, overlay = T, xtrans = trans_pd, trans = trans_fd, gg = F, partial = T,
              line = list(col = alpha(colours, 0.4), lwd = 1), 
              fill = list(col = alpha(colours, 0.4)),
              points = list(col = alpha(colours, 0.4), cex = 0.2),
              xlab = 'PD', ylab = 'FD',
              strip.names=c('0-20', '20-40', '40-60', '60-90'), 
              xlim = c(0, 14200),
              ylim = c(-0.05, .2),
              main = 'Oceania and Asia', 
              legend = F
)




### Code to plot the slope of MPFDses and PDses against absolute latitude 

colours <- c( "#340068","#FE5D26", "#EBDDFF",  'darkgrey')

models <- list(fit2, fit2.Africa_Europe, fit2.Asia_Oceania, fit2.Americas)
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
  
  sd2<-sd1*sd(data$fdis)
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

g1<-ggplot(DF)+
  geom_ribbon(aes(ymin=LB, ymax=UB, x=trans_lat(sd1), fill = Region), alpha=.7)+
  geom_line(aes(y=b1, x=trans_lat(sd1), col = Region))+
  geom_hline(aes(yintercept=0), lty='dashed')+
  xlab('Absolute latitude')+
  xlim(c(0, 90)) +
  ylab(Standardised~slope~of~FD~-~PD)+
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme_classic() +
  theme(legend.position = 'none')
g1


pp1g1 <- ggarrange(pp1, g1, nrow = 1, widths = c(3,1), labels = c('a', ''))

ggsave(pp1g1, filename = 'continental_coef_int_slope.pdf', device = 'pdf', width = 12, height =4)


dev.off()

require(egg)
grid_plot <- ggarrange(p0, p1,
                       p2, p3,
                       widths = c(1,1,1,1), nrow = 1, ncol = 4)
grid_plot


require(gridExtra)
require(grid)
grid.arrange(p.AM, p.EA, p.OA, nrow = 1)
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }


gp <- grid_arrange_shared_legend(p1, p2, p3, nrow = 1)
gp
ggsave(plot = gp, 'FD_PD_LAT_by_LONG.pdf', width = 12, height = 6)


options(digits = 3)
ols_table <- tidy(m1, conf.int = T)


write.csv(ols_table, file = "olstab.csv",  quote = FALSE, row.names = F)

ols_table <- all_models %>%
  select(-statistic, -p.value) %>%
  mutate_each(funs(round(., 2)), -term) %>% 
  gather(key, value, estimate:std.error) %>%
  spread(model, value) 


parTable(fit2)

fitmeasures(fit2)
inspect(fit2, what = 'cor.all')
lavCor(fit2)
resid(fit2, "cor")


modificationindices(fit2, minimum.value = 20)

functional_indices$abs_lat2 <- functional_indices$abs_lat^2
model3 <- 'fdis ~ sp_richn + pd + abs_lat + abs_lat:pd
pd ~ sp_richn + sp_richn2 + abs_lat
sp_richn ~ abs_lat + abs_lat2'
fit3 = cfa(model3, data = functional_indices)
summary(fit3, standardized = T)
path3 <- semPaths(fit3, 'std', layout = 'tree3', residuals = T, curve = 2, curvature = 1, nCharNodes = 0,
                  pastel = F, style = 'ram', rotation = 4)

AIC(fit1, fit2, fit3)


plot(path1)
summary(fit1)
plot(fit1@Model)

AIC(fit1)


m1 <- lm(fdis ~ sp_richn + pd + abs_lat + pd:abs_lat + sp_richn:abs_lat, dat = functional_indices)
m2 <- lm(fdis ~ sp_richn + pd + abs_lat + pd:abs_lat, dat = functional_indices)
m2 <- lm(fdis ~ sp_richn + pd + abs_lat + sp_richn:abs_lat, dat = functional_indices)

AIC(m1, m2, m3)
summary(m1)
semPaths(m1)

visreg(m1., 'sp_richn', 'abs_lat', breaks = 5)

v <- visreg(m1, 'sp_richn', 'abs_lat', breaks = c(-2,-1,0,1,2), xlab = 'PD', ylab = 'FD', overlay = T)
v <- visreg(m1, 'pd', 'abs_lat', breaks = c(-2,-1,0,1,2), xlab = 'PD', ylab = 'FD', overlay = T)


gm1 <- lm(fdis ~ sp_richn + pd + abs_lat + pd:abs_lat + sp_richn2 + long_cats, data = functional_indices)
summary(gm1)

ggplot(v$fit, aes(pd, visregFit, linetype=factor(abs_lat), fill=factor(abs_lat))) +
  
  geom_point(data = v$res, aes(pd, visregRes, col = factor(abs_lat)), alpha = 0.5, size = .5) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, 
              linetype=1, size=0.2) +
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)+
  geom_line() +
  ylab("FD") +
  xlab("PD") +
  ylim(c(-10, 10))+
  #xlim(c(-1.5, 4)) +
  theme_classic() +theme(legend.title = element_blank()) + theme(legend.position = 'bottom')
ggsave('FD_PD_Lat.pdf', height = 5, width = 6)


require(sesem)
require(fields)
distance_matrix <- round(rdist.earth(x1 = functional_indices[,5:6], miles = F)[lower.tri(matrix(nrow = nrow(functional_indices), ncol = nrow(functional_indices)), diag=T)])
system.time(bins<-make.bin(distance_matrix,type='n.bins',p.dist=10, n.bins = 3))
binsize<-bins[1][[1]] #truelove lowland bin sizes
binname<-bins[2][[1]] #truelove lowland bin names
plotbin(distance_matrix,binsize)

colnames(functional_indices)

data <- functional_indices[, c(5,6,1,2,3,7,10,11)]

head(functional_indices)
head(data)
data$pd_int_abs_lat <- data$pd*data$abs_lat
#data$abs_lat2 <- NULL
#data$sp_richn2 <- data$sp_richn^2
head(data)
#data <- standardize(data, vars = colnames(data))

system.time(covariances <- make.covar(data, distance_matrix, binsize, binname))
covariances
save(covariances, file = 'covariances.Rdata')
load('~/Dropbox/AVONET/covariances.Rdata')

model1 <- 'fdis ~ sp_richn + pd + abs_lat + pd_int_abs_lat + altitude
pd ~ sp_richn + abs_lat + sp_richn2 + altitude
sp_richn ~ abs_lat + altitude'

results <- runModels(model1, covariances)
plotmodelfit(results, rmsea_err=F)
modelsummary(results)
save(results, file = 'results.Rdata')
results


new_data <- functional_indices
new_data$sp_richn <- median(new_data$sp_richn)
new_data$sp_richn2 <- new_data$sp_richn^2
new_data$abs_lat <- runif(nrow(functional_indices), range(new_data$abs_lat[which(new_data$lat_cat == '(0,20]')]))

head(functional_indices)

pred <- lavPredict(fit1, newdata = new_data)


### Southern v Northern

fit2.North <- cfa(model2, data = subset(functional_indices, lat >= 0), estimator = 'MLM', se = 'robust')
summary(fit2.North, standardized = T)


fit2.South <- cfa(model2, data = subset(functional_indices, lat < 0), estimator = 'MLM', se = 'robust')
summary(fit2.South, standardized = T)



coefficients <- as_tibble(bind_rows(standardizedSolution(fit2)[1:11,],
                                    standardizedSolution(fit2.North)[1:11,],
                                    standardizedSolution(fit2.South)[1:11,]
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
  scale_x_discrete(labels = c('Lat', 'Lat:PD','Alt', 'PD', 'SR', bquote('SR'^2)))
pp2
pp2 <- pp2 +facet_wrap(~lhs, nrow = 1, labeller = labeller(lhs = labels)) + theme(legend.position="bottom") +theme(legend.title = element_blank())
pp2
ggsave('hemisphere_coef.pdf', device = 'pdf', width = 12, height = 6)


models <- list(fit2, fit2.North, fit2.South)
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
  
  sd2<-sd1*sd(data$fdis)
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
  ylab(Standardised~slope~of~FD~-~PD)+
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme_classic() +
  theme(legend.position = 'none')
g2

pp2g2 <- ggarrange(pp2, g2, nrow = 1, widths = c(3,1), labels = c('a', ''))
pp2g2
ggsave(pp2g2, filename = 'hemisphere_coef_int_slope.pdf', device = 'pdf', width = 12, height =4)

