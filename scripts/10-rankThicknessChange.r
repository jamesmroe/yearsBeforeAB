#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Compute the rank order of thickness changes across all converter MRIs and correlate with the rank order of regional Aβ deposition. 
#          Uses linear time-to-Aβ trajectories or SILA input, and computes time-to-Aβ thickness trajectories in regions of high-v-low Aβ.
#          Script requires individual-level data as input and is not executable
#========================================================================================#


#---load packages
loadPackages = function() {
  packages = c("here", "tidyverse","magrittr","gamm4","itsadug","numDeriv","gratia","mgcv","viridis","wesanderson","asbio","broom","cowplot","data.table","stringi","tictoc","lmerTest","effects","ggpubr")
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  print(sapply(packages, require, character.only = T))
  print(sapply(packages, function(p) as.character(packageVersion(p))))
}
loadPackages()
here()



#---set dir
b = "/cluster/projects/p274/projects/p040-ad_change/Berkeley"
# b = here()
setwd(b)
#---set dir


#---make dirstruct
plotdir = "plots"; if (! dir.exists(plotdir)) { dir.create(plotdir)}
resdir = "results"; if (! dir.exists(resdir)) { dir.create(resdir)}


#---load data
savefigs=F
nTime=2
agecut=30
saveres=T


load(file.path(b, "reproduce/data/DF_LONG_4570.Rda"))
dim(DF); length(unique(DF$subject_id))



# load converter/nonconverter data
load(file.path(b, "reproduce/data/converters_all_negfirst_ADNINC_UPDATE.Rda"))
load(file.path(b, "reproduce/data/converters_all_negfirst_BACS_UPDATE.Rda"))
load(file.path(b, "reproduce/data/converters_all_negfirst_LCBC_UPDATE_REPRO_corthresh.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_BACS_UPDATE_REPRO.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_ADNINC_UPDATE_REPRO.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_LCBC_UPDATE_REPRO_corthresh.Rda"))


converters_allMRI_negfirst_LCBC = left_join(converters_allMRI_negfirst_LCBC,
                                            converters_allMRI_negfirst_LCBC_plotdat %>% select(imageLink, slope, intercept, age_at_threshold, contains("CL_at_thresh")) %>% 
                                              rename(slope_centiloid = slope,
                                                     intercept_centiloid = intercept)
)

converters_allMRI_negfirst_ADNINC = left_join(converters_allMRI_negfirst_ADNINC,
                                              converters_allMRI_negfirst_ADNINC_plotdat %>% select(imageLink, slope, intercept, age_at_threshold, contains("CL_at_thresh")) %>% 
                                                rename(slope_centiloid = slope,
                                                       intercept_centiloid = intercept)
)

converters_allMRI_negfirst_BACS = left_join(converters_allMRI_negfirst_BACS,
                                            converters_allMRI_negfirst_BACS_plotdat %>% select(imageLink, slope, intercept, age_at_threshold, contains("CL_at_thresh")) %>% 
                                              rename(slope_centiloid = slope,
                                                     intercept_centiloid = intercept)
)



allFeat = readLines(file.path(b, "reproduce/data/allFeatures364.txt"))
adnioutlier = "029_S_0845"
rois=allFeat[grepl("thickness", allFeat)]
nrois=length(rois)
subset.size=nrois; jj = 1
N = ceiling(nrois/subset.size)
print(N)
start = (jj*subset.size)-subset.size+1

if (jj == N) {
  end = nrois
  loopend = length(rois[start:end])
} else {
  end = jj*subset.size
  loopend = subset.size
}
print(paste("subsetting cols", start, "-", end))
rois = rois[start:end]
print(rois)
ROIs = rois


pb = txtProgressBar(min=2, max=end, style=3)
Usubs = length(unique(DF$subject_id))


# load sila outputs
osila_bacs = fread("/cluster/projects/p274/projects/p040-ad_change/Berkeley/scripts/SILA-AD-Biomarker/demo/output/testBACS.csv")
osila_adni = fread("/cluster/projects/p274/projects/p040-ad_change/Berkeley/scripts/SILA-AD-Biomarker/demo/output/testADNINC.csv")
osila_lcbc = fread("/cluster/projects/p274/projects/p040-ad_change/Berkeley/scripts/SILA-AD-Biomarker/demo/output2/testLCBC.csv")


# load sila inputs (ids get changed in sila modelling)
isila_all = fread("/cluster/projects/p274/projects/p040-ad_change/Berkeley/scripts/SILA-AD-Biomarker/demo/df_silo_amyloidtimeCorrect.csv")
isila_adni = isila_all[isila_all$cohort == "ADNINC",] %>% rename(subject_id = subid) %>% rename(subid = subjid)
isila_bacs = isila_all[isila_all$cohort == "BACS",] %>% rename(subject_id = subid) %>% rename(subid = subjid)
isila_lcbc = isila_all[isila_all$cohort == "LCBC",] %>% rename(subject_id = subid) %>% rename(subid = subjid)


range(osila_adni$subid)
range(isila_adni$subid)
range(osila_bacs$subid)
range(isila_bacs$subid)
range(isila_lcbc$subid)
range(osila_lcbc$subid)



osila_adni = left_join(osila_adni, isila_adni)
osila_bacs = left_join(osila_bacs, isila_bacs)
osila_lcbc = left_join(osila_lcbc, isila_lcbc)



mytheme = theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  title = element_text(size=17),
  text = element_text(color = "black", size = 18, family="Nimbus Sans Narrow"),
  plot.title = element_text(hjust = 0.5),
  # axis.ticks = element_blank(),
  axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
  axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,20,0)),
  axis.text = element_text(color = "black", size = 18),
  legend.key.size = unit(1,"cm"))


pal = wes_palette("FantasticFox1", n = 5)
plotSila = function(dat, cohort) {
  
  # dat = osila_adni
  (p_sila1 = 
     dat %>% 
     ggplot(.) +
     geom_line(data=dat,aes(x=age,val,group=subid, col = factor(conv)),alpha=0.6, size=0.5) +
     geom_point(data=dat,aes(x=age,val,group=subid, col = factor(conv)),stat="identity",alpha=1, size=0.5) +
     scale_color_manual(values = c(pal[2], pal[5])) +
     # geom_smooth(method = "gam", col = "black", se = F) +
     ggtitle(cohort) +
     labs(x = "Age") +
     theme_classic() + mytheme)
  
  (p_sila2 = 
      dat %>% 
      ggplot(.) +
      geom_line(data=dat,aes(x=estdtt0,val,group=subid, col = factor(conv)),alpha=0.6, size=0.5) +
      geom_point(data=dat,aes(x=estdtt0,val,group=subid, col = factor(conv)),stat="identity",alpha=1, size=0.5) +
      scale_color_manual(values = c(pal[2], pal[5])) +
      geom_hline(yintercept = dat$valt0, linetype = 2, col = "black") +
      ggtitle(cohort) +
      labs(x = "Years to Aβ+ (SILA)") +
      theme_classic() + mytheme)
  
  return(list(p_sila1 = p_sila1, p_sila2 = p_sila2, threshold = dat$valt0[1]))
  
}


p_sila_adni = plotSila(osila_adni, "ADNI")
p_sila_bacs = plotSila(osila_bacs, "BACS")
p_sila_lcbc = plotSila(osila_lcbc, "LCBC")
p_sila_adni$p_sila2
p_sila_bacs$p_sila2
p_sila_lcbc$p_sila2


if (savefigs) {
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_sila_adni.pdf",
         plot = p_sila_adni$p_sila2 + theme(legend.position = "none"),
         width = 13,
         height = 13,
         dpi = 600,
         units = "cm",
         device = cairo_pdf
  )
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_sila_bacs.pdf",
         plot = p_sila_bacs$p_sila2 + theme(legend.position = "none"),
         width = 13,
         height = 13,
         dpi = 600,
         units = "cm",
         device = cairo_pdf
  )
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_sila_lcbc_threshcorrect.pdf",
         plot = p_sila_lcbc$p_sila2 + theme(legend.position = "none"),
         width = 13,
         height = 13,
         dpi = 600,
         units = "cm",
         device = cairo_pdf
  )
}


DF_SILA = rbind(osila_adni,
                osila_bacs,
                osila_lcbc)
DF_SILA$visit_age = DF_SILA$age
DF_SILA$subject_id %in% DF$subject_id
DF_SILA$subject_id[!DF_SILA$subject_id %in% DF$subject_id]
DF_SILA$age_at_sila_threshold = DF_SILA$age - DF_SILA$estdtt0


DF_SILA %>%
  filter(subject_id == "002_S_4213") %>%
  pull(age_at_sila_threshold) %>%
  dput()

# fix four subjects that have very slightly different age at threshold across their observations
(checkSubs = DF_SILA %>%
  group_by(subject_id) %>%
  summarise(n_unique = n_distinct(round(age_at_sila_threshold, 4))) %>%
  filter(n_unique != 1)
)
DF_SILA[DF_SILA$subject_id == "021_S_4276",]
DF_SILA[DF_SILA$subject_id == "031_S_4021",]
DF_SILA[DF_SILA$subject_id == "036_S_4491",]
DF_SILA[DF_SILA$subject_id == "1100591",]

# by taking their last estimate
DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "021_S_4276"] = DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "021_S_4276"][4]
DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "031_S_4021"] = DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "031_S_4021"][3]
DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "036_S_4491"] = DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "036_S_4491"][4]
DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "1100591"] = DF_SILA$age_at_sila_threshold[DF_SILA$subject_id == "1100591"][2]

DF_SILA$age_at_sila_threshold = round(DF_SILA$age_at_sila_threshold, 3)
DF_SILA %>%
  group_by(subject_id) %>%
  summarise(n_unique = n_distinct(age_at_sila_threshold)) %>%
  filter(n_unique != 1)


DF = left_join(DF, 
          DF_SILA %>% select(subject_id, age_at_sila_threshold, conv) %>% distinct()
          )

DF$visit_age
DF$time_from_sila_threshold = DF$age_at_sila_threshold - DF$visit_age
DF$time_from_sila_threshold_flip = DF$time_from_sila_threshold*-1
head(DF %>% select(matches("age|time", ignore.case = TRUE)))



# NB! this is correctly 76 (not 77) due to the negative centiloid slope BACS person -------
# and this person had 10 MRI scans
# hence the difference between 477 scans (in this analysis) and max 487 MRI scans in converters (in paper) is correct
load("/cluster/projects/p274/projects/p040-ad_change/Berkeley/reproduce/data/DF.convallMRI.rda")
length(unique(DF.convallMRI$subject_id)); dim(DF.convallMRI)
DF.convallMRI$diff_mriAge_predABpos


# choose to estimate via original method (linear estimates)
# or SILA
estType = "SILA"
convonly = 1
estType = "ORIG"

if (estType != "SILA") {
  
  # if not SILA analysis (original with converters only)
  DF = DF.convallMRI
  dim(DF)
  length(unique(DF$subject_id))
  converters_allMRI_negfirst_BACS[converters_allMRI_negfirst_BACS$SID == "B16-220",] %>% dim()
  convonly = 0
  
} else if (estType == "SILA") {
  
  # if SILA analysis (review)
  DF = DF %>% filter(!is.na(time_from_sila_threshold))
  dim(DF)
  length(unique(DF$subject_id))
  length(unique(DF$subject_id[DF$conv == 1]))
  
  # if testing SILA across only converter group
  if (convonly) {
    dim(DF %>% filter(conv == 1))
    DF %<>% filter(conv == 1)
  }
}




# high ab region map
frontal_regions_lh <- c(
  "lh_superiorfrontal_thickness.aparcnative71",
  "lh_rostralmiddlefrontal_thickness.aparcnative71",
  "lh_caudalmiddlefrontal_thickness.aparcnative71",
  "lh_parsopercularis_thickness.aparcnative71",
  "lh_parstriangularis_thickness.aparcnative71",
  "lh_parsorbitalis_thickness.aparcnative71",
  "lh_lateralorbitofrontal_thickness.aparcnative71",
  "lh_medialorbitofrontal_thickness.aparcnative71",
  "lh_frontalpole_thickness.aparcnative71"
)
frontal_regions_rh <- c(
  "rh_superiorfrontal_thickness.aparcnative71",
  "rh_rostralmiddlefrontal_thickness.aparcnative71",
  "rh_caudalmiddlefrontal_thickness.aparcnative71",
  "rh_parsopercularis_thickness.aparcnative71",
  "rh_parstriangularis_thickness.aparcnative71",
  "rh_parsorbitalis_thickness.aparcnative71",
  "rh_lateralorbitofrontal_thickness.aparcnative71",
  "rh_medialorbitofrontal_thickness.aparcnative71",
  "rh_frontalpole_thickness.aparcnative71"
)
parietalregions = c(
  "lh_precuneus_thickness.aparcnative71",
  "rh_precuneus_thickness.aparcnative71",
  "lh_inferiorparietal_thickness.aparcnative71",
  "rh_inferiorparietal_thickness.aparcnative71",
  "lh_supramarginal_thickness.aparcnative71",
  "rh_supramarginal_thickness.aparcnative71"
)
cingulate_regions = c(
  "lh_posteriorcingulate_thickness.aparcnative71",
  "rh_posteriorcingulate_thickness.aparcnative71",
  "lh_isthmuscingulate_thickness.aparcnative71",
  "rh_isthmuscingulate_thickness.aparcnative71",
  "lh_rostralanteriorcingulate_thickness.aparcnative71",
  "rh_rostralanteriorcingulate_thickness.aparcnative71",
  "lh_caudalanteriorcingulate_thickness.aparcnative71",
  "rh_caudalanteriorcingulate_thickness.aparcnative71"
)
temporal_regions = c(
  "lh_middletemporal_thickness.aparcnative71",
  "rh_middletemporal_thickness.aparcnative71"
)


ab_regions_high = c(frontal_regions_lh, frontal_regions_rh, parietalregions, temporal_regions, cingulate_regions)


# low ab regions - everything else
ROIs[!ROIs %in% ab_regions_high]
ab_regions_low = ROIs[!ROIs %in% ab_regions_high]
ab_regions_low = ab_regions_low[!grepl("Mean", ab_regions_low)]


# FDR regions in thickness analysis + homologues
FDR_mask = c(frontal_regions_lh, frontal_regions_rh,
             "lh_precentral_thickness.aparcnative71",
             "rh_precentral_thickness.aparcnative71",
             "lh_paracentral_thickness.aparcnative71",
             "rh_paracentral_thickness.aparcnative71",
             "lh_insula_thickness.aparcnative71",
             "rh_insula_thickness.aparcnative71",
             "lh_supramarginal_thickness.aparcnative71",
             "rh_supramarginal_thickness.aparcnative71"
)


# minus frontal poles which were not FDR sig
FDR_mask = FDR_mask[!grepl("pole", FDR_mask)]
ROIs = FDR_mask


# reverse years to AB to be correct (-years to AB)
if (estType != "SILA") {
  DF$diff_mriAge_predABpos_flip = DF$diff_mriAge_predABpos*-1
} else if (estType == "SILA") {
  DF$diff_mriAge_predABpos_flip = DF$time_from_sila_threshold_flip
}
DF$brainvarLow = rowMeans(DF[,ab_regions_low])
DF$brainvarHigh = rowMeans(DF[,ab_regions_high])
DF$subject_id = as.factor(DF$subject_id)
fullDF = DF


if (estType != "SILA") {
  facSmooth = T
} else {
  facSmooth = F
}


if (facSmooth) {

  # ordered factor approach to get test statistics for GAMM interaction
  dat = DF
  dat$brainvarfac = dat$brainvarHigh
  stackdat = rbind(
    dat %>% mutate(highlow = "high"),
    dat %>% mutate(highlow = "low")
  )
  stackdat$brainvarfac[stackdat$highlow == "low"] = dat$brainvarLow
  stackdat = mutate(stackdat,
                    ohighlow = factor(highlow, levels = c("low","high"),ordered = T))


  gamm.trajectories = gamm4(brainvarfac ~ s(diff_mriAge_predABpos_flip, by = as.factor(highlow)) + as.factor(highlow) + visit_age + subject_sex + cohort + scanStrength + ICV,
                            data = stackdat, random = ~ (1 |subject_id))
  gamm.sum = summary(gamm.trajectories$gam)
  g = gamm.trajectories$gam
  plot.gam(g)
   
  
  # estimate smooth for set reflevel and a smoothed difference between ref and other levels
  ogamm.trajectories = gamm4(brainvarfac ~ as.factor(highlow) + s(diff_mriAge_predABpos_flip) + s(diff_mriAge_predABpos_flip, by = ohighlow) + visit_age + subject_sex + cohort + scanStrength + ICV,
                             data = stackdat, random = ~ (1 |subject_id))
  ogamm.sum = summary(ogamm.trajectories$gam)
  plot.gam(ogamm.trajectories$gam)
 
  
  # thickness trajectory in high Aβ regions was significantly different than low Aβ regions
  ogamm.sum
}



# three analyses to run through - select which here
regionTest = "high"
regionTest = "low"
regionTest = "cortex"



if (regionTest == "high") {
  loopend = 1
} else if (regionTest == "low") {
  loopend = 1
} else if (regionTest == "cortex") {
  loopend = length(ROIs)
}



for (i in 1:loopend) {
  # tic()
  print(paste(i,"/",length(ROIs)))
  if (i == 1) {
    derivMat = gratiaderivMat = gratiaderivCIlwrMat = gratiaderivCIuprMat = matrix(NA, nrow = 1000, ncol = length(ROIs))
    fitMat = matrix(NA, nrow = 100, ncol = length(ROIs))
    RR = list()
    p_derivs = p_derivs_se = list()
    yearsBeforepredAB_on_ci_exclusion = c()
    yearsBeforepredAB_on_se_exclusion = c()
    maxaccels = c()
  }
  DF = fullDF
  setTxtProgressBar(pb,i)
  set.seed(123)
  
  
  # ab regions low / high
  if (regionTest == "high") {
    ROI = "high AB composite"
    DF$brainvar = rowMeans(DF[,ab_regions_high])
    palcol = pal[5]
    anaTitle = "Aβ high"
  }
  
  if (regionTest == "low") {
    ROI = "low AB composite"
    DF$brainvar = rowMeans(DF[,ab_regions_low])
    palcol = pal[2]
    anaTitle = "Aβ low"
  }
  
  if (regionTest == "cortex") {
    ROI = ROIs[i]
    DF$brainvar = DF[[ROI]]
    palcol = "darkgrey"
    anaTitle = ROI
  }
  
  
  g_convall = gamm4(brainvar ~ s(diff_mriAge_predABpos_flip) + visit_age + subject_sex + cohort + scanStrength + ICV, data = DF, random = ~ (1 | subject_id))
  g_convall_mgcv = gam(
    brainvar ~ 
      s(diff_mriAge_predABpos_flip) +
      s(subject_id, bs = "re") +
      visit_age + subject_sex + cohort + scanStrength + ICV,
    data = DF,
    method = "REML"
  )
  summary(g_convall$gam)
  g_sum = summary(g_convall_mgcv)
  RR[[i]] = g_sum
  
  
  lmm_convall = lmer(
    brainvar ~ diff_mriAge_predABpos_flip + visit_age + subject_sex + cohort + scanStrength + ICV + (1 | subject_id),
    data = DF
  )
  # plot.gam(g_convall$gam, residuals = T)
  summary(lmm_convall)
  
  
  # clearly the covariates are important
  ggplot(DF, aes(y=brainvar, x = diff_mriAge_predABpos_flip)) + geom_point(aes(group = subject_id)) + geom_smooth(method = "gam")
  plot(effect("diff_mriAge_predABpos_flip", lmm_convall, residuals=TRUE)) #warning is due to scaling of Y and X
  
  
  predictions <- DF %>% 
    mutate(visit_age = mean(DF$visit_age), subject_sex = "Female",  cohort = "ADNINC", scanStrength = "1-5T", ICV=0) %>% 
    select(diff_mriAge_predABpos_flip, visit_age,subject_sex,ICV,mri_info_site_name,scanStrength, cohort) %>% 
    predict(g_convall$gam, newdata = ., se.fit = T)
  
  
  residualsg <- residuals(g_convall$mer)
  DF$partial_residuals = predictions$fit + residualsg
  DF$fit = predictions$fit
  DF$sefit = predictions$se.fit
  DF$cifit = predictions$se.fit*1.96
  
  
  
  #colour palette ---
  pal = wesanderson::wes_palettes$FantasticFox1
  pointcol = "#6faca8"
  #colour palette ---
  
  if (estType != "SILA") {
    tmpDF = DF %>% filter(diff_mriAge_predABpos >= 1)
  } else {
    tmpDF = DF %>% filter(time_from_sila_threshold >= 1)
  }
  (fig1 = 
      DF %>% #filter(subject_id != adnioutlier) %>%
      ggplot(.) +
      geom_line(data=tmpDF,aes(x=diff_mriAge_predABpos_flip,partial_residuals,group=subject_id),color=pointcol,alpha=0.6, size=0.5) +
      geom_point(data=tmpDF,aes(x=diff_mriAge_predABpos_flip,partial_residuals,group=subject_id),color=pointcol,stat="identity",alpha=1, size=0.5) +
      # geom_ribbon(data=dug,aes(x=diff_mriAge_predABpos_flip,ymin=fit-CI,ymax=fit+CI),alpha=.7,show.legend=F,fill="dark grey") +
      geom_line(data=tmpDF,aes(x=diff_mriAge_predABpos_flip,y=fit),col="black") +
      ggtitle(ROI) +
      labs(x = "Age") +
      theme_classic() + mytheme)
  
  
  
  # Extract plotting data without rendering the plot
  plot_data = plot.gam(g_convall$gam, select = 1, se = TRUE, rug = FALSE, shade = T, pages = 0)
  
  
  # get gam intercept
  g_sum = summary(g_convall$gam)
  g_sum$p.coeff
  
  
  # The output is a list, extract x, fit, and se
  df = data.frame(
    x = plot_data[[1]]$x,
    fit = plot_data[[1]]$fit + g_sum$p.coeff[1],
    se = plot_data[[1]]$se
  )
  
  
  # Calculate upper and lower confidence intervals
  df$upper = df$fit + 1 * df$se #NB! help(plot.gam) shows it is already using CI
  df$lower = df$fit - 1 * df$se
  
  
  ggplot(df, aes(x = x, y = fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(
    ) +
    theme_minimal()
  
  
  # Approximate first derivative
  df$deriv = c(NA, diff(df$fit) / diff(df$x))
  
  # Find index of minimum (i.e., where uptick starts)
  min_idx = which.min(df$fit)
  
  # Optionally, find first point where derivative turns positive
  first_uptick_idx = which(df$deriv > 0 & seq_along(df$deriv) > min_idx)[1]
  
  # Get corresponding x-value
  uptick_point = df$x[first_uptick_idx]
  uptick_pointy = df$fit[first_uptick_idx]
  
  
  fitMat[,i] = df$fit
  
  
  (pderiv0 = ggplot(df, aes(x = x, y = fit)) +
    geom_line() +
    labs(
      x = "Years to Aβ+",
      y = "Thickness"
    ) +
    # geom_point(aes(x = uptick_point, y = uptick_pointy), col = "black", size = 5) +
    theme_minimal())
  
  
  (pderiv0_ci = ggplot(df, aes(x = x, y = fit)) +
    geom_line(col = palcol) +
    ggtitle(ROI) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, col = palcol, fill = palcol) +
    labs(
      x = "Years to Aβ+",
      y = "Thickness"
    ) +
    # geom_point(aes(x = uptick_point, y = uptick_pointy), col = "black", size = 5) +
    theme_minimal())
  
  
  # first derivative
  df$deriv = c(NA, diff(df$fit) / diff(df$x))
  df$deriv_lower = c(NA, diff(df$lower) / diff(df$x))
  df$deriv_upper = c(NA, diff(df$upper) / diff(df$x))
  
  
  (pderiv1 = ggplot(df, aes(x = x, y = deriv)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Years to Aβ+", y = "Rate of change") +
    theme_minimal())
  
  
  (pderiv1_se = ggplot(df, aes(x = x, y = deriv)) +
    geom_line() +
    geom_ribbon(aes(ymin = deriv_lower, ymax = deriv_upper), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Years to Aβ+", y = "Rate of change") +
    theme_minimal())
  
  
  derivMat[,i] = df$deriv
  
  
  # Compute second derivative (numerical)
  second_derivative = diff(diff(df$fit)) / diff(df$x[-1])
  
  # Align x-axis (midpoints between x values)
  second_x = df$x[-c(1, 2)] + diff(df$x)[-1] / 2
  
  #point of maximum accerelated change
  maxaccel_y = second_derivative[which(second_derivative == max(second_derivative))]
  maxaccel_x = second_x[which(second_derivative == max(second_derivative))]
  
  if (length(maxaccel_x) > 1) {
    maxaccel_y = NA
    maxaccel_x = NA
  }
  if (is.na(maxaccel_x)) {
    maxaccels[i] = NA
  } else {
    maxaccels[i] = maxaccel_x
    
    pderiv2 = ggplot(data.frame(x = second_x, second_derivative = second_derivative), aes(x = x, y = second_derivative)) +
        geom_line() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(y = "Acceleration", x = "Years to Aβ+") +
        geom_point(aes(x = maxaccel_x[i], y = maxaccel_y[i]), col = "black", size = 2) +
        theme_minimal()
    
    # check all derivatives
    (pcheckderiv = cowplot::plot_grid(pderiv0, pderiv1, pderiv2, nrow = 3))
    
  }
  
  
  
  # NB! se of calculated derivative not ideal - using gratia instead
  (psumse1 = cowplot::plot_grid(pderiv0_ci, pderiv1_se, nrow = 3))
  
  
  # Compute first derivative of the smooth term
  d1 = gratia::derivatives(g_convall_mgcv, term = "s(diff_mriAge_predABpos_flip)", interval = "confidence", n = 1000)
  d1 = as.data.frame(d1)
  d1$.lower_ci[which(d1$.lower_ci > 0)]
  
  
  # point at which the derivative CI excluded zero used to estimate the rank order of thickness changes
  cross_point = which(d1$.lower_ci > 0)[1]
  
  
  data.frame(d1$diff_mriAge_predABpos_flip, d1$.lower_ci, logical = d1$.lower_ci > 0)
  
  
  if (is.na(cross_point)) {
    print("CIs do not exclude zero")
    yearsBeforepredAB_on_ci_exclusion[i] = "never"
    exclude0 = 0
  } else if (cross_point == 1) {
    print("CIs always exclude zero")
    yearsBeforepredAB_on_ci_exclusion[i] = "always"
    cross_y = 0
    cross_x = d1$diff_mriAge_predABpos_flip[cross_point]
    exclude0 = 1
  } else {
    print("CIs exclude zero")
    cross_y = d1$.lower_ci[cross_point]
    cross_x = d1$diff_mriAge_predABpos_flip[cross_point]
    yearsBeforepredAB_on_ci_exclusion[i] = cross_x
    exclude0 = 1
  }
  
  
  #repeat for SE (since not all cross)
  d1$upper_se = d1$.derivative + d1$.se
  d1$lower_se = d1$.derivative - d1$.se
  
  
  
  # point at which the derivative SE excluded zero
  cross_point_se = which(d1$lower_se > 0)[1]
  
  # two instances where the deriv excluded zero on the downward trajectory - fixed to first point on upward trajectory
  data.frame(d1$diff_mriAge_predABpos_flip, d1$lower_se, logical = d1$lower_se > 0)
  # if (estType == "SILA" & convonly == 0) { if (i == 12) { cross_point_se = 420 } }  # fix the one across full group
  # if (estType == "SILA" & convonly == 0) { if (i == 2) { cross_point_se = 320 } }  # fix the one across full group
  # if (estType == "SILA" & convonly == 0) { if (i == 9) { cross_point_se = 537 } }  # fix the one across full group
  # if (estType == "SILA" & convonly == 0) { if (i == 10) { cross_point_se = 322 } }  # fix the one across full group
  # if (estType == "SILA" & convonly == 1) { if (i == 4) { cross_point_se = 262 } }  # fix if only converters
  # if (estType == "SILA" & convonly == 1) { if (i == 2) { cross_point_se = 281 } }  # fix the one across full group
  # if (estType == "SILA" & convonly == 1) { if (i == 15) { cross_point_se = 296 } }  # fix the one across full group
  # if (estType == "SILA" & convonly == 1) { if (i == 21) { cross_point_se = NA } }  # fix the one across full group
  
  if (is.na(cross_point_se)) {
    print("SEs do not exclude zero")
    yearsBeforepredAB_on_se_exclusion[i] = "never"
  } else if (cross_point_se == 1) {
    print("SEs always exclude zero")
    cross_yse = 0
    cross_xse = d1$diff_mriAge_predABpos_flip[cross_point_se]
    yearsBeforepredAB_on_se_exclusion[i] = "always"
  } else {
    print("SEs exclude zero")
    cross_yse = d1$lower_se[cross_point_se]
    cross_xse = d1$diff_mriAge_predABpos_flip[cross_point_se]
    yearsBeforepredAB_on_se_exclusion[i] = cross_xse
  }
  
  
  # check crossing points
  if (!is.na(cross_point)) {
    ggplot(d1, aes(x = diff_mriAge_predABpos_flip, y = .derivative)) +
      geom_line() +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(x = cross_x, y = cross_y), col = "black", size = 5)
  }
  if (!is.na(cross_point_se) & cross_point_se != 1) {
    ggplot(d1, aes(x = diff_mriAge_predABpos_flip, y = .derivative)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower_se, ymax = upper_se), alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(x = cross_xse, y = cross_yse), col = "black", size = 5)
  }
  
  
  
  (pderiv1_gratia_ci = ggplot(d1, aes(x = diff_mriAge_predABpos_flip, y = .derivative)) +
      geom_line(col = palcol) +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.1, colour = palcol, fill = palcol) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        x = "Years to Aβ+",
        y = "First derivative"
      ) +
      # geom_point(aes(x = cross_x, y = cross_y), col = "black", size = 1) +
      theme_classic() + mytheme
  )
  
  
  (pderiv1_gratia_se = ggplot(d1, aes(x = diff_mriAge_predABpos_flip, y = .derivative)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower_se, ymax = upper_se), alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        x = "Years to Aβ+",
        y = "First derivative"
      ) +
      # geom_point(aes(x = cross_xse, y = cross_yse), col = "black", size = 1) +
      theme_classic() + mytheme
  )
  
  
  
  if (regionTest == "cortex") {
    if (!is.na(cross_point)) {
      # add crossing point to plot
      pderiv1_gratia_ci = pderiv1_gratia_ci + geom_point(aes(x = cross_x, y = cross_y), size = 6, fill = "#cbac09", col = "#cbac09")
    }
    if (!is.na(cross_point_se)) {
      pderiv1_gratia_se = pderiv1_gratia_se + geom_point(aes(x = cross_xse, y = cross_yse), size = 6, fill = "#cbac09", col = "#cbac09")
    }
  }
  
  
  # save gratia outputs
  gratiaderivMat[,i] = d1$.derivative
  gratiaderivCIlwrMat[,i] = d1$.lower_ci
  gratiaderivCIuprMat[,i] = d1$.upper_ci
  
  
  # main plot
  (
    p_combine_fitderiv_ci = ggpubr::ggarrange(
      pderiv0_ci + theme_classic() + mytheme + ggtitle(anaTitle) +
        pderiv1_gratia_ci,
      nrow = 1,
      align = "hv"
    )
  )
  
  (
    p_combine_fitderiv_se = ggpubr::ggarrange(
      pderiv0_ci + theme_classic() + mytheme + ggtitle(anaTitle) +
        pderiv1_gratia_se,
      nrow = 1,
      align = "hv"
    )
  )  
  
  
  if (regionTest == "cortex") {
    (
      p_combine_fitderiv_ci = ggpubr::ggarrange(
        pderiv0_ci + theme_classic() + mytheme + ggtitle(anaTitle) +  theme(axis.ticks = element_blank(), axis.line = element_blank()),
        pderiv1_gratia_ci + theme(axis.ticks = element_blank(), axis.line = element_blank()),
        nrow = 1,
        align = "hv"
      )
    )

  }
  p_derivs[[i]] = p_combine_fitderiv_ci
  p_derivs_se[[i]] = p_combine_fitderiv_se
  
  
  if (savefigs == 1) {
    if (regionTest == "high") {
      print("saving high plot")
      # ggsave(plot = p_combine_fitderiv_ci,
      #               filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_thickTraj_ABhigh.pdf"),
      #               width=20, height=9, units="cm", dpi=600, device = cairo_pdf
      #        )
    } else if (regionTest == "low") {
      print("saving low plot")
      # ggsave(plot = p_combine_fitderiv_ci,
      #        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_thickTraj_ABlow.pdf"),
      #        width=20, height=9, units="cm", dpi=600, device = cairo_pdf
      # )
    } else {
      # ggsave(plot = p_combine_fitderiv_ci,
      #        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_thickTraj_", ROI, ".pdf"),
      #        width=20, height=9, units="cm", dpi=600, device = cairo_pdf
      # )
    }
  }
  
}

# check cross point on all 24 roi derivs
p_derivs[[1]]
p_derivs[[2]]
p_derivs[[3]]
p_derivs[[4]]
p_derivs[[5]]
p_derivs[[6]]
p_derivs[[7]]
p_derivs[[8]]
p_derivs[[9]]
p_derivs[[10]]
p_derivs[[11]]
p_derivs[[12]]
p_derivs[[13]]
p_derivs[[14]]
p_derivs[[15]]
p_derivs[[16]]
p_derivs[[17]]
p_derivs[[18]]
p_derivs[[19]]
p_derivs[[20]]
p_derivs[[21]]
p_derivs[[22]]
p_derivs[[23]]
p_derivs[[24]]


p_derivs_se[[1]]
p_derivs_se[[2]]
p_derivs_se[[3]]
p_derivs_se[[4]]
p_derivs_se[[5]]
p_derivs_se[[6]]
p_derivs_se[[7]]
p_derivs_se[[8]]
p_derivs_se[[9]]
p_derivs_se[[10]]
p_derivs_se[[11]]
p_derivs_se[[12]]
p_derivs_se[[13]]
p_derivs_se[[14]]
p_derivs_se[[15]]
p_derivs_se[[16]]
p_derivs_se[[17]]
p_derivs_se[[18]]
p_derivs_se[[19]]
p_derivs_se[[20]]
p_derivs_se[[21]]
p_derivs_se[[22]]
p_derivs_se[[23]]
p_derivs_se[[24]]
# all confirmed correct


# rank order on thickness change
data.frame(yearsBeforepredAB_on_ci_exclusion,
           yearsBeforepredAB_on_se_exclusion
)

rankThickTraj = data.frame(rois = ROIs, 
                           rankedThickTraj = yearsBeforepredAB_on_ci_exclusion,
                           rankedThickTrajSE = yearsBeforepredAB_on_se_exclusion,
                           rankedThickTrajMaxAccel = maxaccels)

# set as 0 if CI always excludes 0
rankThickTraj$rankedThickTraj[rankThickTraj$rankedThickTraj == "always"] = 0
rankThickTraj$rankedThickTraj[rankThickTraj$rankedThickTraj == "never"] = NA
rankThickTraj$rankedThickTrajSE[rankThickTraj$rankedThickTrajSE == "always"] = 0
rankThickTraj$rankedThickTrajSE[rankThickTraj$rankedThickTrajSE == "never"] = NA
rankThickTraj$rankedThickTraj = as.numeric(rankThickTraj$rankedThickTraj)
rankThickTraj$rankedThickTrajSE = as.numeric(rankThickTraj$rankedThickTrajSE)



# rank and reorder ROIs based on the point at which the CI / SE of the derivative crosses 0
# NB! no need to reverse as yearsBeforepredAB was flipped in model (diff_mriAge_predABpos_flip)
# CI
rankThickTraj$rankedThickTrajRev <- ifelse(
  is.na(rankThickTraj$rankedThickTraj),
  NA,
  ifelse(rankThickTraj$rankedThickTraj == 0,
         0,
         rank(rankThickTraj$rankedThickTraj, ties.method = "first"))
)
rankThickTraj %<>% arrange(rankedThickTrajRev)


# SE
rankThickTraj$rankedThickTrajSERev <- ifelse(
  is.na(rankThickTraj$rankedThickTrajSE),
  NA,
  ifelse(rankThickTraj$rankedThickTrajSE == 0,
         0,
         rank(rankThickTraj$rankedThickTrajSE, ties.method = "first"))
)
rankThickTraj %<>% arrange(rankedThickTrajSERev)



# load in amyloid order
rankAB = fread(file.path(b, "reproduce/data/rankABpredTraj.csv"))
rankAB %<>% rename(rankedABTraj = rankedTraj,
                   rankedABTrajRev = revRankedTraj)
rankThickTraj$rois = paste0("CTX_", toupper(gsub("_thickness.aparcnative71", "", rankThickTraj$rois)), "_SUVR")
rankThickTraj$rois %in% rankAB$rois

rankedBoth = merge(
  rankThickTraj, rankAB
)

# filter data where there is no rank order for thickness (i.e. ROI derivative did not exclude 0 and thus could not be ranked)
rankedBoth_cut = rankedBoth %>% filter(!is.na(rankedThickTrajRev))


# make ranking plots
rankedBoth$roi = gsub("CTX_", "", rankedBoth$rois)
rankedBoth$roi = gsub("_SUVR", "", rankedBoth$roi)
rankedBoth$roi = tolower(rankedBoth$roi)
orderplot = arrange(rankedBoth, rankedABTrajRev) %>% select(roi)
rankedBoth$roi = factor(rankedBoth$roi, levels = rev(orderplot$roi))
rankRange = range(rankedBoth$rankedABTrajRev, na.rm= T)

table(rankedBoth$rankedABTrajRev)
rankedBoth$rankedABTrajRev = rankedBoth$rankedABTrajRev-1 # make 0 indexed (so colours match up with thickness rank plot)
(p_rankAB_FDRregions = ggplot(rankedBoth %>%
                                arrange(rankedABTrajRev),
                              aes(x = rankedABTrajRev, y = factor(roi, levels = rev(unique(roi))))) +
    geom_bar(aes(x=rankedABTrajRev, y =roi), stat = "identity", colour = "grey", alpha = 0.5, width = 0.01, size=0.5) +
    geom_point(aes(colour = rankedABTrajRev), fill = "white", shape = 21, stroke = 2, size = 5) +
    scale_colour_viridis(option = "E", direction = -1, name = "Rank", limits = c(0, rankRange[2]), oob = scales::squish) +
    theme_classic() +
    mytheme +
    labs(x = "Rank (Aβ)",y = NULL) + theme(legend.position = "none")
)


# ggsave(plot = p_rankAB_FDRregions,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankAB_FDRregions.png"),
#        width=7, height=26, units="cm", dpi=600)


# ensure rank plot for thickness is ordered the same
(p_rankThick_frontal = ggplot(rankedBoth,
                              aes(x = rankedThickTrajRev, y = factor(roi, levels = rev(unique(roi))))) +
    geom_bar(aes(x=rankedThickTrajRev, y =roi), stat = "identity", colour = "grey", alpha = 0.5, width = 0.01, size=0.5) +
    geom_point(aes(colour = rankedThickTrajRev), fill = "white", shape = 21, stroke = 2, size = 5) +
    scale_colour_viridis(option = "E", direction = -1, name = "Rank", limits = c(0, rankRange[2]), oob = scales::squish) +
    theme_classic() +
    mytheme +
    labs(x = "Rank (Thickness)",y = NULL) + theme(legend.position = "none")
)

table(rankedBoth$rankedThickTrajSERev)
(p_rankThick_frontal_bySE = ggplot(rankedBoth,
                                   aes(x = rankedThickTrajSERev, y = factor(roi, levels = rev(unique(roi))))) +
    geom_bar(aes(x=rankedThickTrajSERev, y =roi), stat = "identity", colour = "grey", alpha = 0.5, width = 0.01, size=0.5) +
    geom_point(aes(colour = rankedThickTrajSERev), fill = "white", shape = 21, stroke = 2, size = 5) +
    scale_colour_viridis(option = "E", direction = -1, name = "Rank", limits = c(0, rankRange[2]), oob = scales::squish) +
    theme_classic() +
    mytheme +
    labs(x = "Rank (Thickness)",y = NULL) + theme(legend.position = "none")
)

cp1 = cowplot::plot_grid(p_rankAB_FDRregions, p_rankThick_frontal + theme(axis.text.y = element_text(color = "transparent"), axis.line.y = element_blank(), axis.ticks.y = element_blank()))
cp2 = cowplot::plot_grid(p_rankAB_FDRregions, p_rankThick_frontal_bySE + theme(axis.text.y = element_text(color = "transparent"), axis.line.y = element_blank(), axis.ticks.y = element_blank()))

# ggsave(plot = p_rankAB_FDRregions,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankAB_FDRregions.pdf"),
#        width=12, height=26, units="cm", dpi=600, device = cairo_pdf)

# ggsave(plot = p_rankThick_frontal,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankRThick_FDRregions.pdf"),
#        width=12, height=26, units="cm", dpi=600, device = cairo_pdf)



table(rankedBoth_cut$rankedThickTrajRev)
pear1 = cor.test(rankedBoth_cut$rankedThickTrajRev, rankedBoth_cut$rankedABTrajRev)
pear2 = cor.test(rankedBoth_cut$rankedThickTrajSERev, rankedBoth_cut$rankedABTrajRev)

pval1 = format(
  tidy(pear1)$p.value[1],
  scientific = TRUE, digits = 2)
bval1 = round(
  tidy(pear1)$estimate[1],
  digits = 2)
pval2 = format(
  tidy(pear2)$p.value[1],
  scientific = TRUE, digits = 2)
bval2 = round(
  tidy(pear2)$estimate[1],
  digits = 2)

table(rankedBoth_cut$rankedThickTrajRev)


# correlation plot based on CI (fig 5g)
(p_rankCI = ggplot(rankedBoth_cut, aes(x = rankedThickTrajRev, y = rankedABTrajRev)) +
    geom_point(aes(col = rankedThickTrajRev), size = 3) +
    scale_colour_viridis(option = "E", direction = -1, name = "Rank", limits = c(0, rankRange[2]), oob = scales::squish) +
    theme_classic() +
    mytheme +
    # coord_fixed() +
    geom_smooth(method = "lm", se = T, col = "black", alpha = .1) +
    labs(x = "Rank (Thickness)",y = "Rank (Aβ)") + theme(legend.position = "none") +
    ggtitle("order of CIs excluding 0") +
    annotate("text", x = 1, y = max(rankedBoth_cut$rankedABTrajRev)+2, label = paste("p =", pval1), color = "black", hjust = 0) +
    annotate("text", x = 1, y = max(rankedBoth_cut$rankedABTrajRev)+3, label = paste("B =", bval1), color = "black", hjust = 0)
)  

# ggsave(plot = p_rankCI,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankCI.pdf"),
       # width=9, height=12, units="cm", dpi=600, device = cairo_pdf)

# ggsave(plot = p_rankCI,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankCI_", estType, "convonly", convonly, ".pdf"),
# width=9, height=12, units="cm", dpi=600, device = cairo_pdf)


table(rankedBoth_cut$rankedThickTrajSERev)


# correlation plot based on SE (fig 5h)
(p_rankSE = ggplot(rankedBoth_cut, aes(x = rankedThickTrajSERev, y = rankedABTrajRev)) +
    geom_point(aes(col = rankedThickTrajSERev), size = 3) +
    scale_colour_viridis(option = "E", direction = -1, name = "Rank", limits = c(0, rankRange[2]), oob = scales::squish) +
    theme_classic() +
    mytheme +
    # coord_fixed() +
    geom_smooth(method = "lm", se = T, col = "black", alpha = .2) +
    labs(x = "Rank (Thickness)",y = "Rank (Aβ)") + theme(legend.position = "none") +
    ggtitle("order of SEs excluding 0") +
    annotate("text", x = 1, y = max(rankedBoth_cut$rankedABTrajRev)+2, label = paste("p =", pval2), color = "black", hjust = 0) +
    annotate("text", x = 1, y = max(rankedBoth_cut$rankedABTrajRev)+3, label = paste("B =", bval2), color = "black", hjust = 0)
)

# ggsave(plot = p_rankSE,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankSE.pdf"),
#        width=9, height=12, units="cm", dpi=600, device = cairo_pdf)

# ggsave(plot = p_rankSE,
#        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_rankSE_", estType, "convonly", convonly, ".pdf"),
#        width=9, height=12, units="cm", dpi=600, device = cairo_pdf)

