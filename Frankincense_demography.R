
### ************************************** ###
### ****** IPM Boswellia papyrifera ****** ###

#### Load packages and main data ----
#
#install.packages("lme4")

library(reshape2)
library(plyr)
library(MASS)
library(stats)
library(mgcv)
library(Matrix)
library(IPMpack)
library(fields)
library(ggplot2)
library(ggpmisc)
library(popbio)
library(gtable)
library(extrafont)
library(devtools)
require(lme4)
library(tidyr)
library(dplyr)
library(broom)

rm(list=ls(all=TRUE))

font_import(pattern="[A/a]rial")
sessionInfo()

wd_txt <- "" # add working directory folder

setwd(wd_txt)
dff_all_ini      <- read.table("data_in.txt",header=TRUE,sep="\t")
overview_data    <- read.table("overview_data_in.txt",header=TRUE,sep="\t")
yield            <- read.table("yield_data.txt", header = T, sep = "\t")
diam_data        <- read.table("diam_data.txt",header=TRUE,sep="\t")
#

#### Load functions: future distributions + multiplot ####

futureDist <- function (startingSizes, IPM, n.time.steps, startingEnv = 1) 
{
  breakpoints <- c(IPM@meshpoints - (IPM@meshpoints[2] - IPM@meshpoints[1]), 
                   IPM@meshpoints[length(IPM@meshpoints)] + (IPM@meshpoints[2] - IPM@meshpoints[1]))
  if (IPM@nEnvClass > 1) {
    if (length(startingEnv) == 1) 
      startingEnv <- rep(startingEnv, length(startingSizes))
    compound <- TRUE
    env.index <- IPM@env.index
    n.new.dist <- rep(0, length(IPM[1, ]))
    for (ev in 1:IPM@nEnvClass) {
      index.new.dist <- findInterval(startingSizes[startingEnv == ev], breakpoints, all.inside = TRUE)
      loc.sizes <- table(index.new.dist)
      n.new.dist[ev == IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
    }
    n.new.dist0 <- n.new.dist
  }
  else {
    compound <- FALSE
    index.new.dist <- findInterval(startingSizes, breakpoints, all.inside = TRUE)
    loc.sizes <- table(index.new.dist)
    env.index <- rep(1, length(IPM@meshpoints))
    n.new.dist <- rep(0, length(IPM@meshpoints))
    n.new.dist[as.numeric(names(loc.sizes))] <- loc.sizes
    n.new.dist0 <- n.new.dist
  }
  for (t in 1:n.time.steps) n.new.dist <- IPM@.Data %*% n.new.dist
  #plot(IPM@meshpoints, n.new.dist0[env.index == 1], type = "l", xlab = "size", ylab = "n in each size class", ylim = range(c(n.new.dist0, n.new.dist)))
  #points(IPM@meshpoints, n.new.dist[env.index == 1], type = "l", col = 2)
  if (compound) {
    for (j in 1:max(env.index)) {
      #points(IPM@meshpoints, n.new.dist0[env.index == j], type = "l", col = 1, lty = j)
      #points(IPM@meshpoints, n.new.dist[env.index == j], type = "l", col = 2, lty = j)
    }
  }
  #legend("topright", legend = c("current", "future"), col = 1:2, lty = 1, bty = "n")
  return(list(n.new.dist0 = n.new.dist0, n.new.dist = n.new.dist))
} 

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) 
{
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#



#### Boostrap to calculate confidence intervals ####

runs <- 500                # runs to bootstrap
yrs  <- 50                 # years to project transient dynamics
mesh <- 500                # size of matrix mesh
mnsz <- 1                  # minimum size (1 cm dbh)
mxsz <- 1.05*(52.5)        # maximum size (5% above maxium observed size)
bopt <- mnsz+c(0:mesh)*(mxsz-mnsz)/mesh     # boundary points
mept <- 0.5*(bopt[1:mesh]+bopt[2:(mesh+1)]) # mesh points
stsz <- mept[2]-mept[1]                     # step size

dff_all_ini      <- dff_all_ini[c(1:17,24)]
dff_all_ini$incr <- (dff_all_ini$sizeNext-dff_all_ini$size)
diam_data        <- diam_data[c(2,3,6)]

## IPM model and choices 
# choose hilghland or lowland data, regeneration yes or no and whether to use ring data or the combination of plot and ring data 

lo_hi    <- "highland"   # "lowland" or "highland"
reg_no   <- "no_reg"       # "reg" or "no_reg"
rings_no <- "rings"     # "rings" or "rings_plots"

dff_all  <- droplevels(subset(dff_all_ini, dff_all_ini$low_high == lo_hi))

lm_yield <- lm(yield$yield[yield$spots>6 & yield$low_high == lo_hi] ~ yield$dbh[yield$spots>6 & yield$low_high == lo_hi])
overview_lo_hi <- droplevels(subset(overview_data, overview_data$low_high == lo_hi& overview_data$reg == reg_no))

pred_pop_temp           <- data.frame(matrix(NA, ncol = 8, nrow = 0))
colnames(pred_pop_temp) <- c("site", "site_nr", "rel_pop", "rel_yield", "run", "year", "abs_pop", "abs_yield")
pred_pop                <- pred_pop_temp

for (ran in 1:runs) { # bootstrap of the IMP model
  tryCatch({
    dff      <- dff_all
    rand_ind <- c(sample(dff$ID, floor(max(dff$ID)*0.01))) # select 1% of individuals
    dff      <- subset(dff, !(dff$ID %in% rand_ind))       # remove 1% of ind for bootstrap
    dff_gro_rng      <- droplevels(subset(dff, dff$data_type == "rings"))
    dff_gro_plo      <- droplevels(subset(dff, dff$data_type == "plots" & dff$stage > 4 & dff$incr < 4 & dff$incr > - 1.5))
    
    if (rings_no == "rings") {
      dff_gronly <- dff_gro_rng
    } else { 
      dff_gronly <- merge(dff_gro_rng, dff_gro_plo, all = T)}
    
    dff_gronly$size2 <- dff_gronly$size^2
    dff_gronly$size3 <- dff_gronly$size^3
    
    grh <- makeGrowthObj(dataf = dff_gronly, Formula = incr ~ size + size2, regType = "changingVar")
    sv  <- makeSurvObj(dff, Formula=surv~size)
    dff_suronly  <- subset(dff, dff$stage >4) # Adapt survival curve to yearly value
    xx           <- seq(min(dff_suronly$dbh1,na.rm=T),ceiling(max(dff_suronly$dbh3,na.rm=T)),by=.1)
    
    if (sum(dff_suronly$surv5,na.rm = T) > 0) {
      surv.reg1 <- glm(surv1~dbh1,data=dff_suronly,family=binomial()) # calculate curve year 1
      surv.reg2 <- glm(surv3~dbh3,data=dff_suronly,family=binomial()) # calculate curve year 2
      surv.reg3 <- glm(surv5~dbh5,data=dff_suronly,family=binomial()) # calculate curve year 3
      pred.s    <- data.frame(size = xx, reg1 = predict(surv.reg1, data.frame(dbh1=xx),type="response"), reg2 = predict(surv.reg2, data.frame(dbh3=xx),type="response"), reg3 = predict(surv.reg3, data.frame(dbh5=xx),type="response"))
      
    } else {
      surv.reg1 <- glm(surv1~dbh1,data=dff_suronly,family=binomial()) # calculate curve year 1
      surv.reg2 <- glm(surv3~dbh3,data=dff_suronly,family=binomial()) # calculate curve year 2
      pred.s    <- data.frame(size = xx, reg1 = predict(surv.reg1, data.frame(dbh1=xx),type="response"), reg2 = predict(surv.reg2, data.frame(dbh3=xx),type="response"))}
    
    pred.s$mean <- rowMeans(pred.s[,2:3]) # mean of predicted yearly curves
    surv.reg    <- glm(mean~size,data=pred.s,family= quasibinomial(logit)) # calculate avg curve
    sv@fit$coefficients <- coef(surv.reg) ## change coefficients survival object
    dff_feconly <- droplevels(subset(dff, dff$data_type == "fecundity")[,c("size","fec", "sizeNext")])
    fec_recr    <- ((9/42)+(3/43))/2     # average values for year 1-2 and 2-3
    mnoffspring <- 1.6126  # offspr size = based on DBH-ht relationship (else on DBH-RCD = 3.7258)
    stoffspring <- 0.7486  # offspr stdv = based on DBH-ht relationship (else on DBH-RCD = 0.3011)
    
    fv  <- makeFecObj(dff_feconly, Formula=list(fec~size), Family ="binomial", meanOffspringSize = mnoffspring, sdOffspringSize = stoffspring, fecConstants = data.frame(fec1Recruits = fec_recr))
    Pmatrix <- makeIPMPmatrix(growObj = grh, survObj = sv, nBigMatrix = mesh, minSize = mnsz, maxSize = mxsz, nEnvClass = 1, integrateType = "midpoint", correction="constant")
    Fmatrix <- makeIPMFmatrix(nEnvClass = 1, fecObj = fv, nBigMatrix = mesh, minSize = mnsz, maxSize = mxsz, integrateType="midpoint", correction="constant", preCensus = F)
    
    if (reg_no == "reg") {
      sc_non_rep <- length(Fmatrix@meshpoints[Fmatrix@meshpoints<10]) # meshpoint @10cm DBH
      Fmatrix@.Data[, 1:sc_non_rep] <- 0 # set classes <10cm DBH as non-reproductive
    } else {
      Fmatrix@.Data <- Fmatrix@.Data*0}
    
    ### --- ### --- ### --- ### --- ### --- ### 
    IPM <- makeIPMmatrix(Pmatrix,Fmatrix)
    ### --- ### --- ### --- ### --- ### --- ###
    
    for (si_te in 1:nrow(overview_lo_hi)){ # calculate size distributions for the next 50 years
    
    site_name         <- "Asgede-Tsimbla" #paste(overview_lo_hi$site[si_te])
    site_nr           <- 17 #paste(overview_lo_hi$site_nr[si_te])
    print(paste(site_name, "- round", ran))
    start_size        <- droplevels(subset(diam_data,diam_data$site == site_name))
    
    for (dis in 1:yrs+1){
      dist_yr           <- futureDist(startingSizes = start_size$dbh, IPM = IPM, n.time.steps = dis-1)
      current_pop       <- data.frame(size = mept, freq = dist_yr$n.new.dist0)
      pop_yr            <- data.frame(size = mept, freq = dist_yr$n.new.dist)
      
      current_pop$yield <- ifelse(current_pop$size >= 10, (current_pop$size*summary(lm_yield)$coefficients[2] + summary(lm_yield)$coefficients[1])*current_pop$freq, 0) # calculate yield from current population distribution
      pop_yr$yield      <- ifelse(pop_yr$size >= 10,  (pop_yr$size*summary(lm_yield)$coefficients[2] + summary(lm_yield)$coefficients[1])*pop_yr$freq, 0)             # calculate yield from future population distribution
      
      pred_pop_temp[dis,1]   <- site_name
      pred_pop_temp[dis,2]   <- site_nr
      pred_pop_temp[dis,3]   <- sum(pop_yr$freq)/sum(current_pop$freq)
      pred_pop_temp[dis,4]   <- sum(pop_yr$yield)/sum(current_pop$yield)
      pred_pop_temp[dis,5]   <- ran
      pred_pop_temp[dis,6]   <- dis-1
      pred_pop_temp[dis,7]   <- sum(pop_yr$freq)
      pred_pop_temp[dis,8]   <- sum(pop_yr$yield)
    }
    pred_pop_temp[1,] <- c(site_name, site_nr, 1, 1, ran, 0, sum(current_pop$freq), sum(current_pop$yield))
    pred_pop <- rbind(pred_pop,pred_pop_temp)
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

pred_pop$ring_plot <- paste(rings_no)
write.table(pred_pop, paste("pred_pop_", lo_hi, "_", reg_no, "_", rings_no, ".txt", sep = ""), sep = '\t', col.names = T, row.names = F)

#### Figures main text ----

my_theme <- theme(plot.background  = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), panel.spacing = unit(0.75, "lines"), text = element_text(family= "sans", size= 8), axis.text = element_text(family = "sans", size = 8))

age_diam  <- read.table("age_diam.txt", sep = "\t", header = T)
age_diam$Site_name <- ifelse(age_diam$area_nr < 10, paste(0,age_diam$area_nr, ".", age_diam$area, sep=""), paste(age_diam$area_nr, ".", age_diam$area, sep="")) 
my.formula <- y ~ x

age_hist           <- read.table("age_hist.txt",header=TRUE,sep="\t")
age_hist$Site_name <- ifelse(age_hist$area_nr < 10, paste(0,age_hist$area_nr, ".", age_hist$site, sep=""), paste(age_hist$area_nr, ".", age_hist$site, sep="")) 
age_hist$hi_lo     <- ifelse(age_hist$site == "Abergelle", "Highland", "Lowland")
age_hist           <- droplevels(subset(age_hist, age_hist$site != "Kuara"))

## Figure 3 - Age distribution histograms
age_hist_new  <- read.table("age_hist_new.txt",header=TRUE,sep="\t")

hist(age_hist_new$rec_yr_lo[age_hist_new$Site_name == "09.Lemlem Terara"]) #lower
hist(age_hist_new$rec_yr_hi[age_hist_new$Site_name == "09.Lemlem Terara"], breaks = 59)

age_span              <- ceiling((max(age_hist_new$rec_yr) - min(age_hist_new$rec_yr))/30)
age_hist_new$age_aggr <- floor(age_hist_new$ rec_yr/age_span)*age_span
pivot_recr            <- ddply(age_hist_new, c("site", "area_nr", "Site_name", "age_aggr"), summarise, n_trees = length(rec_yr), total_area = mean(total_area), n_tree_ha = (length(rec_yr)/mean(total_area)))
pivot_recr$hi_lo      <- ifelse(pivot_recr$site == "Abergelle", "Highland", "Lowland")

summ_hist  <- ddply(age_hist_new, c("Site_name"), summarise, ln = length(dbh), rec_yr = 1860, y = 0.05, hi_lo = "Highland")
summ_pivot <- ddply(pivot_recr, c("Site_name"), summarise, ln = length(area_nr), sum = sum(n_trees), age_aggr = 1860, y = 1.1*max(pivot_recr$n_tree_ha), hi_lo = "Highland")

pdf("Figure 3 - histogram age distributions density plot.pdf", h=4.6, w=4.6) 
ggplot(age_hist_new, aes(x=rec_yr, color = hi_lo, fill = hi_lo)) +
  geom_density(aes(x=rec_yr_hi), color = "grey85", fill = "#00000022")+
  geom_density(aes(x=rec_yr_lo), color = "grey85", fill = "#00000022")+
  geom_histogram(aes(y=..density..), bins = 30, position="identity", colour = "transparent")+
  geom_density(fill = "transparent")+
  #scale_y_continuous(limits = c(0,0.06))+ 
  scale_x_continuous(breaks = c(1850,1900,1950,2000), expand=c(0,5))+
  facet_wrap(~Site_name) +
  labs(x = "Recruitment year", y= "Age structure (Kernel density estimation)") +
  theme(axis.title = element_text(size = 9), 
        axis.text =  element_text(size = 7), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing.x=unit(0.3, "lines"), 
        panel.spacing.y=unit(0.4,"lines"),        
        legend.position = c(.91,.91),
        legend.background = element_blank(),
        legend.text  = element_text(size = 6),
        legend.key = element_blank(),
        legend.direction = "vertical")+
  geom_text(data = summ_hist, aes(label = paste(Site_name, "  (n =", ln, ")", sep=""), y=0.058), colour = "black", hjust=0, vjust=0, size = 2.4)+
  scale_colour_manual(name = "", values = c("dodgerblue3", "brown4"), labels = c('Highland','Lowland'))+
  scale_fill_manual(name = "", values = c("#1C86EE60", "#CD333360"), labels = c('Highland','Lowland'))
dev.off()

## Figure 4 - predicted populations
pred_pop_all <- read.table(paste(wd_txt, "Predicted pop", "pred_pop_all.txt", sep= "\\"),header=TRUE,sep="\t")
pivot_pred   <- ddply(pred_pop_all, c("ring_plot", "site", "site_nr", "year", "hi_lo", "header"), summarise, mean_pop = mean(rel_pop), mean_yield = mean(rel_yield))

subs_pop50 <- subset(pivot_pred, pivot_pred$mean_pop <0.5)
min_yr_pop50 <- aggregate(year ~ site, subs_pop50, function(x) min(x))
mean_pop50 = mean(min_yr_pop50$year); sd_pop50 = sd(min_yr_pop50$year)

subs_yld50 <- subset(pivot_pred, pivot_pred$mean_yield <0.5)
min_yr_yld50 <- aggregate(year ~ site, subs_yld50, function(x) min(x))
mean_yld50 = mean(min_yr_yld50$year); sd_yld50 = sd(min_yr_yld50$year)

avg_pred               <- droplevels(subset(pivot_pred, pivot_pred$ring_plot == "rings"))
avg_pred_yld           <- avg_pred[,-7]
avg_pred_yld$pop_yld   <- "Frankicense yield"
avg_pred_pop           <- avg_pred[,-8]
avg_pred_pop$pop_yld   <- "Population size"
colnames(avg_pred_yld) <- c("ring_plot", "site", "site_nr", "year", "hi_lo", "header", "value", "pop_yld")
colnames(avg_pred_pop) <- colnames(avg_pred_yld)

avg_pop_yld         <- rbind(avg_pred_yld,avg_pred_pop)
avg_pop_yld$reg     <- ifelse(avg_pop_yld$site_nr %in% c(3,6,7,19,21), "Regeneration", "No Regeneration")
avg_pop_yld$reg     <- factor(avg_pop_yld$reg, levels = c("No Regeneration","Regeneration"))
avg_pop_yld$pop_yld <- factor(avg_pop_yld$pop_yld, levels = c("Population size", "Frankicense yield"))

dummy <- data.frame(pop_yld = c("Population size", "Frankicense yield"), X = c(c(mean_pop50,mean_pop50+sd_pop50,mean_pop50-sd_pop50), c(mean_yld50,mean_yld50+sd_yld50,mean_yld50-sd_yld50)))

fig3 <- 
  ggplot(avg_pop_yld, aes(x=year, y=value, color = hi_lo)) +  
  geom_hline(yintercept = 1, colour = "grey70",linetype = 2) + 
  geom_vline(data = dummy, aes(xintercept = X), colour = "grey50", linetype = 1) + 
  geom_line(aes(group = site))+ 
  facet_wrap(~pop_yld, as.table = T, scales = "free_y", ncol = 1 )+ 
  scale_y_continuous(breaks = c(0,.25,.5,.75,1,1.25,1.5),labels = scales::percent)+ 
  scale_x_continuous(limits=c(0,52),expand=c(0,0))+ 
  labs(x = "Time (years)", y= "Projected size (relative to t=0)")+
  scale_colour_manual(name = "", values = c("#1874CD65", "#8B232365"), labels = c('Highland','Lowland'))+
  theme(axis.title = element_text(size = 10), 
        axis.text =  element_text(size = 9), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        strip.text.x = element_blank(),
        #strip.text.y = element_text(size = 11, face = "bold"),
        strip.background = element_blank(), 
        panel.spacing.x=unit(0.3, "lines"), 
        panel.spacing.y=unit(0.4,"lines"),
        legend.background = element_blank(),
        legend.position = c(.7, .965),
        legend.text  = element_text(size = 9),
        legend.key = element_blank(),
        legend.direction = "horizontal")

pdf("Figure 4 - Predicted populations and yields NEW.pdf", paper = "a4", h=120/25.4, w= 93/25.4)
fig3
dev.off()

#### Figures Supporting information ----

my_theme <- theme(plot.background  = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), panel.spacing = unit(0.75, "lines"), text = element_text(family= "sans", size= 7), axis.text = element_text(family = "sans", size = 7))

## Figure S1a - Age-diameter lines (waaier)
dff_rings <- droplevels(subset(dff_all_ini, dff_all_ini$data_type == "rings"))
dff_rings <- merge(dff_rings, (droplevels(subset(overview_data, overview_data$site %in% unique(dff_rings$site)))[1:2]), all.x = T)
dff_rings$Site_name <- ifelse(dff_rings$site_nr < 10, paste(0,dff_rings$site_nr, ".", dff_rings$site, sep=""), paste(dff_rings$site_nr, ".", dff_rings$site, sep=""))
age_dia_rings <- ddply(dff_rings, c("site", "ID", "low_high"), summarise, max_age = max(age, na.rm = T), max_size = max(size, na.rm = T))
age_dia_rings <- merge(age_dia_rings, (droplevels(subset(overview_data, overview_data$site %in% unique(age_dia_rings$site)))[1:2]), all.x = T)
age_dia_rings$Site_name <- ifelse(age_dia_rings$site_nr < 10, paste(0,age_dia_rings$site_nr, ".", age_dia_rings$site, sep=""), paste(age_dia_rings$site_nr, ".", age_dia_rings$site, sep=""))

## Figure S1 - Diameter-age dots
fig_s1a <- 
  ggplot(dff_rings, aes(x=age , y=size, color = low_high)) + 
  geom_line(size = 0.75, aes(group = ID), alpha = 0.15) + 
  ggtitle("a") +
  geom_point(data= age_dia_rings, aes(x= max_age, y= max_size, color= low_high), 
             shape= 20, size = 1.5, alpha=0.85, stroke= 0) + 
  labs(x = "Age (years)", y= "Diameter (cm)")+
  facet_wrap(~Site_name) +
  scale_fill_discrete(name = "Category", guide="none")+
  scale_color_manual(values=c("dodgerblue3","brown4")) +  
  scale_x_continuous(breaks=seq(0, 150, 25)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10.6, b = 0, l = 0)), 
        plot.title = element_text(family= "sans", face="bold", size = 8, hjust = -0.08), 
        strip.text.x = element_text(family= "sans", size= 7)) +
  my_theme

fig_s1b <- 
  ggplot(age_diam, aes(x=dbh, y=age, color = Site_name)) + 
  geom_point(alpha = 0.5)  + 
  ggtitle("b")+ 
  geom_smooth(method = "lm", fill="grey80") +
  facet_wrap(~Site_name) +
  labs(x = "Diameter (cm)", y= "Age (years)")+
  scale_fill_discrete(name = "Category", guide="none") + 
  coord_cartesian(xlim = c(0,50), ylim = c(0,175))+
  scale_color_manual(values=c("brown4", "brown4", "brown4","dodgerblue3"))+
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), size = rel(2.5), colour = "black", parse = TRUE, label.x = 17, label.y = 163)+ 
  theme(plot.title = element_text(family= "sans", face="bold", size = 8, hjust = -0.08),
        strip.text.x = element_text(family= "sans", size= 7)) +
  my_theme

pdf("Sup Mat Fig S1 - Waaier + size age.pdf", paper = "a4", h= 9, w= 169/25.41, useDingbats = F)
multiplot(fig_s1a,fig_s1b)
dev.off()

## Figure S2 - dbh-growth rate
dff_all_ini    <- read.table("data_in.txt",header=TRUE,sep="\t")
dff_all_ini$type_lo_hi <- paste(dff_all_ini$data_type,dff_all_ini$low_high)
dff_plots      <- droplevels(subset(dff_all_ini, dff_all_ini$data_type == "plots"))
dff_rings      <- droplevels(subset(dff_all_ini, dff_all_ini$data_type == "rings"))
dff_fecun      <- droplevels(subset(dff_all_ini, dff_all_ini$data_type == "fecundity"))
dff_rings$site <- as.character(dff_rings$site)
dff_ploadu     <- droplevels(subset(dff_all_ini, dff_all_ini$data_type == "plots" & dff_all_ini$adult == 1))
dff_riplo      <- rbind(dff_rings, dff_ploadu)

dff_rings$pred_line <- ifelse(dff_rings$low_high == "highland", (0.2338853882-0.0072039955*(dff_rings$size)+0.0003157039*(dff_rings$size)^2), (0.3101773932-0.0069247572*(dff_rings$size)+0.0003202773*(dff_rings$size)^2))
dff_riplo$pred_line <- ifelse(dff_riplo$low_high == "highland", 
                              (0.2659381038-0.0143587704*(dff_riplo$size)+0.0002701484*(dff_riplo$size)^2), 
                              (0.2983835177-0.0055675924*(dff_riplo$size)+0.0003963726*(dff_riplo$size)^2))

ann_text <- data.frame(size = c(8,8,8,8), incr = c(1,1,1.5,1.5), lab=c("form_hi_ri","form_hi_rp","form_lo_ri","form_lo_rp"), type_lo_hi = factor(c("rings highland", "plots highland", "plots lowland","rings lowland"),levels = c("rings highland", "plots highland", "plots lowland","rings lowland")))

## Figure S2a - Yields
yield  <- read.table("yield_data.txt", header = T, sep = "\t")
yield$low_high  <- ifelse(yield$low_high == "highland", "highland", "lowland")
yield6 <- yield[yield$spots >= 6,]
my.formula <- y ~ x

fig_s2a <- 
  ggplot(dff_ploadu, aes(x=size , y=incr, color = type_lo_hi)) + 
  ggtitle("a")+
  geom_hline(yintercept = 0, lty = 1, colour = "grey") +
  geom_point(shape=20, size = 0.75, alpha = 0.075, na.rm = T, stroke = 0) + 
  coord_cartesian(xlim = c(0,50)) + 
  geom_point(data= dff_rings, aes(x=size , y=incr, color= type_lo_hi), size = 0.5, alpha=0.04, na.rm = T) + 
  geom_line(data= dff_rings, aes(x=size , y=pred_line, color= type_lo_hi), na.rm = T) +
  geom_line(data= dff_riplo, aes(x=size , y=pred_line, color= low_high), na.rm = T) +
  labs(x = "Diameter (cm)", y= bquote('Growth rate ('*'cm.' ~year^-1*')'))+
  facet_wrap(~low_high) +
  scale_fill_discrete(name = "Category", guide="none") + 
  scale_color_manual(values=c("dodgerblue3","brown4","dodgerblue3","brown4","cyan4","orange"), name = "")+  
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15.7, b = 0, l = 0)), 
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        plot.title = element_text(family= "sans", face="bold", size = 8, hjust = -0.08)) +
  my_theme

## Figure S2b - Fecundity
fig_s2b <- 
  ggplot(dff_fecun, aes(x=size , y=fec, color = low_high)) + 
  ggtitle("b")+ 
  geom_point(shape=20, alpha=0.2, position=position_jitter(height=0.01), stroke = 0) + 
  geom_smooth(method="glm", method.args = list(family="binomial"), size = 0.5) +
  facet_wrap(~low_high) +
  coord_cartesian(xlim = c(0,50)) + 
  labs(x = "Diameter (cm)", y= "Fecundity probability")+
  scale_fill_discrete(name = "Category", guide="none") + 
   scale_color_manual(values=c("dodgerblue3","brown4")) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 11, b = 0, l = 0)), 
        strip.text.x = element_blank(), 
        plot.title = element_text(family= "sans", face="bold", size = 8, hjust = -0.08)) +
  my_theme

## Figure S2c - Yield
fig_s2c <- 
  ggplot(yield6, aes(x=dbh, y=yield, color = low_high)) + 
  ggtitle("c")+ 
  geom_point(shape = 20, alpha = 0.5, stroke = 0) + 
  geom_smooth(method = "lm", fill="grey80", size = 0.5) +
  facet_wrap(~low_high) +
  coord_cartesian(xlim = c(0,50), ylim = c(0,2860))+  
  labs(x = "Diameter (cm)", y= bquote('Frankincense yield ('*'g.'~year^-1*')'))+
  scale_color_manual(values=c("dodgerblue3", "brown4"))+
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), size = rel(2.2), colour = "black", parse = TRUE, label.x = 25, label.y = 2760)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 8.7, b = 0, l = 0)),
        strip.text.x = element_blank(), 
        plot.title = element_text(family= "sans", face="bold", size = 8, hjust = -0.08)) +
  my_theme

pdf("Sup Mat Fig S2 - Growth Fecundity Yield.pdf", paper = "a4", h= 180/25.41, w= 169/25.41, useDingbats = F)
multiplot(fig_s2a, fig_s2b, fig_s2c)
dev.off()

## Figure S3 - Survival
surv_plot1 <- 
  ggplot(dff_ploadu, aes(x=dbh1 , y=surv1, color = low_high)) + 
  geom_point(shape=20, alpha=0.1, position=position_jitter(height=0.01), stroke = 0) + 
  geom_smooth(method="glm", method.args = list(family = "binomial")) +
  facet_wrap(~low_high) +
  scale_fill_discrete(name = "Category", guide="none")+
  labs(x = "", y= "Survival probability (year 1)") +
  theme(axis.text.x=element_blank()) +
  coord_cartesian(xlim = c(0,55), ylim = c(-0.1,1.1))+ scale_x_continuous(limits = c(0,55), breaks=seq(0, 50, 10)) +
  scale_y_continuous(limits = c(-0.1,1.1), breaks=seq(0, 1, 0.25)) + 
  scale_color_manual(values=c("dodgerblue3","brown4")) + 
  theme(axis.text.x= element_blank(), strip.text.x = element_text(family= "sans", size= 8)) +
  my_theme

surv_plot2 <- 
  ggplot(dff_ploadu, aes(x=dbh3 , y=surv3, color = low_high)) + 
  geom_point(shape=20, alpha=0.1, position=position_jitter(height=0.01), stroke = 0) + 
  geom_smooth(method="glm", method.args = list(family = "binomial")) +
  facet_wrap(~low_high) +
  scale_fill_discrete(name = "Category", guide="none") +
  scale_x_continuous(limits = c(0,55), breaks=seq(0, 50, 10)) +
  scale_y_continuous(limits = c(-0.1,1.1), breaks=seq(0, 1, 0.25)) + 
  scale_color_manual(values=c("dodgerblue3","brown4")) + 
  labs(x = "", y= "Survival probability (year 2)") + 
  theme(axis.text.x= element_blank(), strip.text = element_blank()) +
  my_theme

surv_plot3 <- ggplot(dff_ploadu, aes(x=dbh5 , y=surv5, color = low_high)) + 
  geom_point(shape=20, alpha=0.1, position=position_jitter(height=0.01), stroke = 0) + 
  geom_smooth(method="glm", method.args = list(family = "binomial")) +
  facet_wrap(~low_high) +
  scale_fill_discrete(name = "Category", guide="none")+
  scale_x_continuous(limits = c(0,55), breaks=seq(0, 50, 10)) +
  scale_y_continuous(limits = c(-0.1,1.1), breaks=seq(0, 1, 0.25)) + 
  scale_color_manual(values=c("dodgerblue3","brown4")) + 
  labs(x = "Diameter (cm)", y= "Survival probability (year 3)") + 
  theme(strip.text = element_blank()) +
  my_theme

pdf("Sup Mat Fig S3 - Survival all sites_years.pdf", paper = "a4", h= 180/25.41, w= 169/25.41, useDingbats = F)
multiplot(surv_plot1, surv_plot2, surv_plot3)
dev.off()


#### Additional figures ####
## Figure Predicted populations
pred_pop_all <- read.table(paste(wd_txt, "Predicted pop", "pred_pop_all.txt", sep= "\\"),header=TRUE,sep="\t")
pivot_pred   <- ddply(pred_pop_all, c("ring_plot", "site", "site_nr", "year", "hi_lo", "header"), summarise, mean_pop = mean(rel_pop), mean_yield = mean(rel_yield))

y_lim_pop <- round(max(pivot_pred$mean_pop), digits=1)
y_lim_yld <- round(max(pivot_pred$mean_yield), digits=1)

df_out <- data.frame(matrix(NA, ncol = ncol(pred_pop_all), nrow = 0))
colnames(df_out) <- colnames(pred_pop_all)

for (i in 1:length(unique(pred_pop_all$site))) {
  df_temp_r   <- droplevels(subset(pred_pop_all, pred_pop_all$site_nr == i & pred_pop_all$ring_plot == "rings"))
  list_runs_r <- sample(unique(df_temp_r$run), 250)
  df_temp_r2  <- droplevels(subset(df_temp_r, df_temp_r$run %in% list_runs_r))
  
  df_temp_p   <- droplevels(subset(pred_pop_all, pred_pop_all$site_nr == i & pred_pop_all$ring_plot == "rings_plots"))
  list_runs_p <- sample(unique(df_temp_p$run), 250)
  df_temp_p2  <- droplevels(subset(df_temp_p, df_temp_p$run %in% list_runs_p))
  
  df_out    <- rbind(df_out, df_temp_r2, df_temp_p2)
}

summ_pred_pop <- ddply(df_out, c("header"), summarise, ln = length(run), year = 0, mean_pop = 1.15*(y_lim_pop), mean_yield = 1.15*(y_lim_pop), ring_plot = "rings")
summ_pred_yld <- ddply(df_out, c("header"), summarise, ln = length(run), year = 0, mean_pop = 1.15*(y_lim_yld), mean_yield = 1.15*(y_lim_yld), ring_plot = "rings")


pdf("Sup Mat Fig S7 - pred_pop_all ring and plot data.pdf", h=130/25.4, w= 169/25.4)
ggplot(data = pivot_pred, aes(x=year, y=mean_pop, colour = ring_plot))+ 
  geom_hline(yintercept = 1, colour = "grey70",linetype = 2)+ 
  facet_wrap(~header)+
  geom_line()+
  scale_y_continuous(limits = c(0,y_lim_pop+.3), breaks = c(0,.5,1,1.5),labels = scales::percent)+ 
  labs(x = "Time (years)", y= "Projected population size (relative to t=0)")+
  geom_text(data = summ_pred_pop, aes(label = header), colour = "black", hjust=0, vjust=0, size = 2.7)+
    scale_colour_manual(name = "Source of growth data:", values = c("dodgerblue3", "brown4"), labels = c("Tree rings", "Tree rings and plots"))+
theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 9), 
        axis.text =  element_text(size = 8), 
        panel.background =  element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        strip.text.x = element_blank(),
        strip.background = element_blank(), 
        panel.spacing.x=unit(0.1, "lines"), 
        panel.spacing.y=unit(0.2,"lines"),
        legend.position = c(.77, .06),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background= element_rect(fill=NA, colour = "grey", size =.5, linetype= "solid"))
dev.off()


## Figure Predicted yields
pdf("Sup Mat Fig S8 - pred_yield_all ring and plot data.pdf", h=130/25.4, w= 169/25.4)
ggplot(data = pivot_pred, aes(x=year, y=mean_yield, colour = ring_plot))+ 
  geom_hline(yintercept = 1, colour = "grey70",linetype = 2)+ 
  facet_wrap(~header)+
  geom_line()+
  scale_y_continuous(limits = c(0,y_lim_yld+.4), breaks = c(0,.5,1,1.5),labels = scales::percent)+ 
  labs(x = "Time (years)", y= "Projected frankincense yield (relative to t=0)")+
  geom_text(data = summ_pred_yld, aes(label = header), colour = "black", hjust=0, vjust=0, size = 2.7)+
    scale_colour_manual(name = "Source of growth data:", values = c("dodgerblue3", "brown4"), labels = c("Tree rings", "Tree rings and plots")) +
theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 9), 
        axis.text =  element_text(size = 8), 
        panel.background =  element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        strip.text.x = element_blank(),
        strip.background = element_blank(), 
        panel.spacing.x=unit(0.1, "lines"), 
        panel.spacing.y=unit(0.2,"lines"),
        legend.position = c(.77, .06),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background= element_rect(fill=NA, colour = "grey", size =.5, linetype= "solid"))
  dev.off()
  
  ## Figure new Overview sites
  
  overview_data    <- read.table("overview_sites.txt",header=TRUE,sep="\t")
  head(overview_data)
  
  pdf("Sup Mat Fig new S8 - elevation temp.pdf", paper = "a4", h=75/25.4, w= 100/25.4)
  ggplot(data = overview_data, aes(x=elevation, y=temp, colour = low_high))+ 
    geom_point()+ 
    geom_smooth(method = "lm", se = F)+
    labs(x = "Elevation (m a.s.l.)", y= "Temperature (ºC)")+
    scale_colour_manual(name = "", values = c("dodgerblue3", "brown4"), labels = c("Highland", "Lowland")) +
    theme(plot.title = element_text(size = 10), 
          axis.title = element_text(size = 9), 
          axis.text =  element_text(size = 8), 
          panel.background =  element_blank(), 
          panel.border = element_rect(colour = "grey", fill = NA), 
          strip.text.x = element_blank(),
          strip.background = element_blank(), 
          panel.spacing.x=unit(0.1, "lines"), 
          panel.spacing.y=unit(0.2,"lines"),
          legend.position = c(.65, .9),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.background= element_rect(fill=NA, colour = NA))
  dev.off()

#### Data summaries and statistics ----

dff_suronly  <- subset(dff_all_ini, dff_all_ini$stage >4 & dff_all_ini$data_type == "plots")

surv_summ <- ddply(dff_suronly, c("site", "low_high"), summarise, 
                   surv1 = sum(surv1, na.rm = T)/length(dbh1[!is.na(dbh1)]), 
                   surv2 = sum(surv3, na.rm = T)/length(dbh3[!is.na(dbh3)]),
                   surv3 = sum(surv5, na.rm = T)/length(dbh5[!is.na(dbh5)]))

surv_summ2 <- ddply(dff_suronly, c("site", "low_high"), summarise, 
                    surv1 = 1 - (length(dbh3[!is.na(dbh3)])/length(dbh1[!is.na(dbh1)])), 
                    surv2 = 1 - (length(dbh5[!is.na(dbh5)])/length(dbh3[!is.na(dbh3)])),
                    surv3 = 1 - (length(dbh7[!is.na(dbh7)])/length(dbh5[!is.na(dbh5)])))

dff_feconly <- dff_all_ini %>% 
  filter(data_type == "fecundity") %>% 
  select(site, site_nr, ID, size, fec, low_high, site_old)

fec_summ    <- ddply(dff_feconly, c("low_high"), summarise, fec = sum(fec, na.rm = T)/length(fec[!is.na(fec)]))


### Statistics highland vs. lowland ###
## Regeneration index
head(overview_data)
reg_test4 <- t.test(perc.4cmDBH ~ low_high, data = overview_data)
reg_test7 <- t.test(perc.7cmDBH ~ low_high, data = overview_data)

## Survival
dff_suronly  <- subset(dff_all_ini, dff_all_ini$stage >4 & dff_all_ini$data_type == "plots")

dbh_long <- dff_suronly %>% 
  gather(key="dbh_year", value= "dbh", 11:13) %>% 
  select(site, site_nr, ID, dbh_year, dbh,low_high, site_old)
surv_long <- dff_suronly %>% 
  gather(key="surv_year", value= "survival", 15:17) %>% 
  select(site, site_nr, ID,surv_year,survival,low_high,site_old)

surv_long <- left_join(dbh_long, surv_long)
surv_long$dummy_sel <- paste(surv_long$dbh_year,surv_long$surv_year, sep = "_")
surv_long <- filter(surv_long, dummy_sel %in% c("dbh1_surv1","dbh3_surv3","dbh5_surv5"))

surv_long$year <- ifelse(surv_long$surv_year=="surv1",1,ifelse(surv_long$surv_year=="surv2",2,3))

## Mixed models survival
m_full <- glmer(survival ~ dbh * low_high + (1 | site_old) + (1 | year), data = surv_long, family = binomial)
m_1    <- glmer(survival ~ dbh            + (1 | site_old) + (1 | year), data = surv_long, family = binomial)

AIC(m_full) < AIC(m_1)-2
anova(m_full, m_1)
summary(m_full)
summ_surv <- tidy(m_full)

surv.reg <- glm(survival ~ dbh * low_high, data = surv_long, family=binomial())
summary(surv.reg)

## Fecundity
dff_feconly <- dff_all_ini %>% 
  filter(data_type == "fecundity") %>% 
  select(site, site_nr, ID, size, fec, low_high, site_old)

head(dff_feconly)

fecu.reg <- glm(fec ~ size * low_high, data = dff_feconly, family=binomial())
summary(fecu.reg)
summ_fecu <- tidy(fecu.reg)

## Growth
dff_gro_rng      <- dff_all_ini %>% 
  filter(data_type == "rings", incr) %>% 
  select(site, site_nr, ID, year, age, size, incr, low_high, site_old)%>% 
  droplevels() 

dff_gro_plo      <- dff_all_ini %>% 
  filter(data_type == "plots" & stage > 4 & incr < 4 & incr > - 1.5) %>% 
  select(site, site_nr, ID, year, age, size, incr, low_high, site_old) %>% 
  droplevels()

gro_mod <- lm(incr ~ size * low_high, data = dff_gro_plo) # test for differences in growth using plot data
summary(gro_mod)
plot(resid(gro_mod))
qqnorm(resid(gro_mod))
qqline(resid(gro_mod))

gro_mod2 <- lm(log(incr) ~ size * low_high, data = dff_gro_rng) # test for differences in growth using ring data
summary(gro_mod2)
plot(resid(gro_mod2))
qqnorm(resid(gro_mod2))
qqline(resid(gro_mod2))

summ_lm2 <- tidy(gro_mod2)

dff_gro_rng$log_incr <- log(dff_gro_rng$incr) #log transform data to build plot

ggplot(dff_gro_rng, aes(x=size , y=log_incr, color = low_high)) + 
  geom_point(shape=20, size = 1.5, alpha = 0.075, na.rm = T, stroke = 0) + 
  coord_cartesian(xlim = c(0,50)) + 
  geom_smooth(method = "lm", se = F)+
  labs(x = "Diameter (cm)", y= bquote('Growth rate ('*'cm.' ~year^-1*')'))+
  scale_color_manual(values=c("#1874CD", "#8B2323"), labels = c('Highland','Lowland'), name = "")+  
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15.7, b = 0, l = 0)), 
        plot.background  = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        legend.position = c(0.9,0.9), 
        legend.background = element_blank(), 
        legend.key = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        panel.spacing = unit(0.75, "lines"), 
        text = element_text(family= "sans", size= 8), 
        axis.text = element_text(family = "sans", size = 8))

#### Population SWAP  ----

pred_pop_hi_reg        <- read.table("pred_pop_data_type_highland.txt",header=TRUE,sep="\t")
pred_pop_hi_reg$reg_no <- "reg" 
pred_pop_hi_no_reg     <- read.table("pred_pop_data_type_highland_no_reg.txt",header=TRUE,sep="\t")
pred_pop_hi_no_reg$reg_no <- "no_reg" 
pred_pop_lo_reg        <- read.table("pred_pop_data_type_lowland.txt", header=TRUE,sep="\t")
pred_pop_lo_reg$reg_no <- "reg" 
pred_pop_lo_no_reg     <- read.table("pred_pop_data_type_lowland_no_reg.txt", header=TRUE,sep="\t")
pred_pop_lo_no_reg$reg_no <- "no_reg" 

pred_pop_hi <- rbind(subset(pred_pop_hi_reg,     pred_pop_hi_reg$site    %in% c("Adi Arkay", "Kisha")), 
                     subset(pred_pop_hi_no_reg, !pred_pop_hi_no_reg$site %in% c("Adi Arkay", "Kisha")))

pred_pop_lo <- rbind(subset(pred_pop_lo_reg,     pred_pop_lo_reg$site    %in% c("Kurmuk", "Mengue", "Guba")), 
                     subset(pred_pop_lo_no_reg, !pred_pop_lo_no_reg$site %in% c("Kurmuk", "Mengue", "Guba")))

unique(pred_pop_lo$reg_no)
pred_pop_swap <- rbind(pred_pop_hi,pred_pop_lo)

pivot_pred_hi   <- ddply(pred_pop_hi, c("site", "site_nr", "year", "data_source"), summarise, 
                         mean_pop = mean(rel_pop), 
                         mean_yield = mean(rel_yield))

subs_pop50_hi   <- subset(pivot_pred_hi, pivot_pred_hi$mean_pop <0.5)
min_yr_pop50_hi <- aggregate(year ~ site, subs_pop50_hi, function(x) min(x))
mean_pop50_hi   <- mean(min_yr_pop50_hi$year)
sd_pop50_hi     <- sd(min_yr_pop50_hi$year)

subs_yld50_hi   <- subset(pivot_pred_hi, pivot_pred_hi$mean_yield <0.5)
min_yr_yld50_hi <- aggregate(year ~ site, subs_yld50_hi, function(x) min(x))
mean_yld50_hi   <- mean(min_yr_yld50_hi$year)
sd_yld50_hi     <- sd(min_yr_yld50_hi$year)

pivot_pred_lo   <- ddply(pred_pop_lo, c("site", "site_nr", "year", "data_source"), summarise, mean_pop = mean(rel_pop), mean_yield = mean(rel_yield))

subs_pop50_lo   <- subset(pivot_pred_lo, pivot_pred_lo$mean_pop <0.5)
min_yr_pop50_lo <- aggregate(year ~ site, subs_pop50_lo, function(x) min(x))
mean_pop50_lo   <- mean(min_yr_pop50_lo$year)
sd_pop50_lo     <- sd(min_yr_pop50_lo$year)

subs_yld50_lo   <- subset(pivot_pred_lo, pivot_pred_lo$mean_yield <0.5)
min_yr_yld50_lo <- aggregate(year ~ site, subs_yld50_lo, function(x) min(x))
mean_yld50_lo   <- mean(min_yr_yld50_lo$year)
sd_yld50_lo     <- sd(min_yr_yld50_lo$year)

dummy <- data.frame(pop_yld = c("Population size", "Frankicense yield"), 
                    data_source = rep(c("highland", "lowland"), each = 6), 
                    X = c(c(mean_pop50_hi, mean_pop50_hi + sd_pop50_hi, mean_pop50_hi - sd_pop50_hi), 
                          c(mean_yld50_hi, mean_yld50_hi + sd_yld50_hi, mean_yld50_hi - sd_yld50_hi),
                          c(mean_pop50_lo, mean_pop50_lo + sd_pop50_lo, mean_pop50_lo - sd_pop50_lo), 
                          c(mean_yld50_lo, mean_yld50_lo + sd_yld50_lo, mean_yld50_lo - sd_yld50_lo)))

dummy_pop <- data.frame(data_source = rep(c("highland", "lowland"), each = 3), 
                        X = c(c(mean_pop50_hi,mean_pop50_hi+sd_pop50_hi,mean_pop50_hi-sd_pop50_hi), 
                              c(mean_pop50_lo,mean_pop50_lo+sd_pop50_lo,mean_pop50_lo-sd_pop50_lo)))

dummy_yld <- data.frame(data_source = rep(c("highland", "lowland"), each = 3), 
                        X = c(c(mean_yld50_hi,mean_yld50_hi+sd_yld50_hi,mean_yld50_hi-sd_yld50_hi), 
                              c(mean_yld50_lo,mean_yld50_lo+sd_yld50_lo,mean_yld50_lo-sd_yld50_lo)))

mean_pop_yield <- data.frame(data_source = c("highland", "lowland"), 
                             mean_pop = c(mean_pop50_hi,mean_pop50_lo), 
                             sd_pop = c(sd_pop50_hi, sd_pop50_lo),
                             mean_yld = c(mean_yld50_hi,mean_yld50_lo), 
                             sd_yld = c(sd_yld50_hi, sd_yld50_lo))

swap_pop <- 
  ggplot(pred_pop_swap, aes(x=year, y=rel_pop, color = data_source)) +  
  geom_hline(yintercept = 1, colour = "grey70",linetype = 2) + 
  geom_vline(data = dummy_pop, aes(xintercept = X), colour = "grey50", linetype = 1) + 
  geom_line(aes(group = site))+ 
  facet_wrap(~data_source, as.table = F, scales = "free_y", ncol = 2 )+ 
  scale_y_continuous(limits = c(0,1.55), breaks = c(0,.25,.5,.75,1,1.25,1.5),labels = scales::percent)+ 
  scale_x_continuous(limits=c(0,52),expand=c(0,0))+ 
  labs(x = "Time (years)", y= "Projected size (relative to t=0)")+
  scale_colour_manual(name = "", values = c("#1874CD65", "#8B232365"), labels = c('Highland data','Lowland data'))+
  theme(axis.title = element_text(size = 10), 
        axis.text =  element_text(size = 9), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        strip.text.x = element_blank(),
        #strip.text.y = element_text(size = 11, face = "bold"),
        strip.background = element_blank(), 
        panel.spacing.x=unit(0.3, "lines"), 
        panel.spacing.y=unit(0.4,"lines"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(.05, .965),
        legend.text  = element_text(size = 9),
        legend.key = element_blank(),
        legend.direction = "horizontal")


swap_yld <- ggplot(pred_pop_swap, aes(x=year, y=rel_yield, color = data_source)) +  
  geom_hline(yintercept = 1, colour = "grey70",linetype = 2) + 
  geom_vline(data = dummy_yld, aes(xintercept = X), colour = "grey50", linetype = 1) + 
  geom_line(aes(group = site))+ 
  facet_wrap(~data_source, as.table = F, scales = "free_y", ncol = 2 )+ 
  scale_y_continuous(limits = c(0,1.55), breaks = c(0,.25,.5,.75,1,1.25,1.5),labels = scales::percent)+ 
  scale_x_continuous(limits=c(0,52),expand=c(0,0))+ 
  labs(x = "Time (years)", y= "Projected yield (relative to t=0)")+
  scale_colour_manual(name = "", values = c("#1874CD65", "#8B232365"), labels = c('Highland data','Lowland data'))+
  theme(axis.title = element_text(size = 10), 
        axis.text =  element_text(size = 9), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey", fill = NA), 
        strip.text.x = element_blank(),
        #strip.text.y = element_text(size = 11, face = "bold"),
        strip.background = element_blank(), 
        panel.spacing.x=unit(0.3, "lines"), 
        panel.spacing.y=unit(0.4,"lines"),
        legend.background = element_blank(),
        legend.position = "none")

pdf("Sup Mat Fig - Predicted populations SWAP.pdf", paper = "a4", h= 6, w= 169/25.41)
plot_grid(swap_pop,swap_yld, nrow = 2)
dev.off()