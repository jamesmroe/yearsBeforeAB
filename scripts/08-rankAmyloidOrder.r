#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Calculate AB deposition maps and rank order of regional deposition
#========================================================================================#


#---load packages
loadPackages = function() {
  packages = c("here", "data.table","tidyverse","magrittr","viridis","ggpubr","tidyr","broom","ggseg","broom.mixed","lme4","lmerTest")
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  sapply(packages, require, character.only = T)
}
loadPackages()


rm(list=ls())
b = "/cluster/projects/p274/projects/p040-ad_change/Berkeley"
load(file.path(b, "reproduce/data/df_amyloidOrder.Rda"))

df_AB = rbind(
  df_lcbc_convab,
  df_lcbc_nonconvab,
  df_adni_convab,
  df_adni_nonconvab,
  df_bacs_convab,
  df_bacs_nonconvab
)
length(unique(df_lcbc_convab$PTID))
length(unique(df_adni_convab$PTID))
length(unique(df_bacs_convab$PTID)) # NB! one missing BACS converter is negative centiloid slope

df_AB$TRACER <- factor(df_AB$TRACER)
df_AB$TRACER <- relevel(df_AB$TRACER, ref = "PiB")

tmp_CONV = df_AB %>% filter(convgroup == 1)
tmp_NONCONV = df_AB %>% filter(convgroup == 0)

lhrois = names(df_AB)[grepl(glob2rx("CTX*_LH_*SUVR"), names(df_AB))]
rhrois = names(df_AB)[grepl(glob2rx("CTX*_RH_*SUVR"), names(df_AB))]



# calculate average AB SUVR map -----
calcAvgMap = function(df) {
  
  df_submeans = df %>%
    select(newlink, PTID, all_of(c(lhrois, rhrois))) %>%
    group_by(PTID) %>%
    summarise(across(all_of(c(lhrois, rhrois)), ~mean(.x, na.rm = TRUE)))
  
  length(unique(df_submeans$PTID)) == nrow(df_submeans)
  average_map = colMeans(df_submeans[2:ncol(df_submeans)], na.rm = T)
  average_map = data.frame(names(average_map), unlist(unname(average_map)))
  return(average_map)
}


getDim = function(df) {
  print(paste(nrow(df), "obs"))
  print(paste(length(unique(df$PTID)), "subs"))
}


# wrangle for ggseg -----
ggWrangle = function(imap) {
  
  omap = imap
  names(omap) = c("ROI","SUVR")
  omap$plotrois = tolower(omap$ROI)
  omap$plotrois = gsub("ctx_", "", omap$plotrois)
  omap$plotrois = gsub("_suvr", "", omap$plotrois)
  omap$hemi = substr(omap$plotrois, 1,2)
  
  omap = omap[omap$hemi == "lh" | omap$hemi == "rh",]
  omap$hemi = ifelse(omap$hemi == "lh", "left",
                     ifelse(omap$hemi == "rh", "right", NA))
  omap$plotrois[!omap$plotrois %in% dk$data$label]
  omap$plotrois[omap$plotrois %in% dk$data$label]
  omap$label = omap$plotrois
  
  return(omap)
}


# average persons AB across ABpos timepoints
# mean across 1315 AB positive obs from 734 independent subs
mapABpos_noNC = calcAvgMap(AB_posnoNC)
mapABpos_noNC = ggWrangle(mapABpos_noNC)
getDim(AB_posnoNC)


# pre-correct for tracer
nrois = c(lhrois, rhrois)
for (ii in 1:length(nrois)) {
  if (ii == 1) {
    modlist = list(); suvr_corrMat = matrix(ncol = length(nrois), nrow = nrow(df_AB))
  }
  roi = c(lhrois, rhrois)[ii]
  print(paste(roi, ii))
  tmp = df_AB
  tmp$suvr = tmp[[roi]]
  lmer_model <- lmer(
    suvr ~ TRACER + (1 | PTID),
    data = tmp
  )
  # ggplot(tmp, aes(TRACER, suvr, col = factor(convgroup))) + geom_jitter() + geom_smooth(method = "lm")
  modlist[[ii]] = summary(lmer_model)$coefficients
  tmp$suvr_cor = tmp$suvr
  tmp %<>%
    mutate(
      suvr_cor = case_when(
        TRACER == "F18" ~ suvr_cor - modlist[[1]]["TRACERF18", 1],
        TRACER == "FBB" ~ suvr_cor - modlist[[1]]["TRACERFBB", 1],
        TRACER == "FBP" ~ suvr_cor - modlist[[1]]["TRACERFBP", 1],
        TRUE ~ suvr_cor
      )
    )
  # ggplot(tmp, aes(TRACER, suvr)) + geom_point() + geom_smooth(method = "lm") +
  #   geom_point(aes(TRACER, suvr_cor), col = "red")
  suvr_corrMat[,ii] = tmp$suvr_cor
}
df_suvr_corrMat = data.frame(suvr_corrMat)
names(df_suvr_corrMat) = c(lhrois, rhrois)


df_suvr_corr = data.frame(
  df_AB %>% select(-one_of(c(lhrois, rhrois))), 
  df_suvr_corrMat)


abmap_convCOR = ggWrangle(
  calcAvgMap(df_suvr_corr[df_suvr_corr$convgroup == 1,])
)
abmap_nonconvCOR = ggWrangle(
  calcAvgMap(df_suvr_corr[df_suvr_corr$convgroup == 0,])
)
abmap_conv = ggWrangle(
  calcAvgMap(df_AB[df_AB$convgroup == 1,])
)
abmap_nonconv = ggWrangle(
  calcAvgMap(df_AB[df_AB$convgroup == 0,])
)


mytheme = theme(plot.background = element_rect(fill = "white"),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                title = element_text(size=12),
                text = element_text(color = "black", size = 12, family="Nimbus Sans Narrow"),
                plot.title = element_text(hjust = 0.5),
                # axis.ticks = element_blank(),
                axis.title.y = element_text(color = "black", size = 12, vjust =-1),
                axis.title.x = element_text(color = "black", size = 12, vjust = -2),
                axis.text = element_text(color = "black", size = 12),
                legend.key.size = unit(1,"cm"))


common_text = theme(panel.background = element_rect(fill = "white"),
                    plot.background = element_rect(fill = "white"),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    # plot.subtitle = element_text(hjust = 0.5),
                    plot.subtitle = element_blank(),
                    plot.title = element_text(size = 10),
                    legend.title = element_text(size = 10),
                    legend.text = element_text(size = 10),
                    axis.line = element_blank(),
                    axis.ticks = element_blank())


# plot with ggseg
ggSegIt = function(imap, map_text, lwr_thresh, upr_thresh) {
  
  newAtlas = dk %>% 
    as_tibble() %>% 
    left_join(as.data.frame(imap))
  
  newBrainAtlas = newAtlas %>%
    as_brain_atlas()
  
  newbrain = ggseg(.data = newAtlas, mapping = aes(fill = SUVR), col = "dark grey", size= 0.5, atlas = "dk")
  
  # SUVR cortical map -----
  (ggSuvr = newbrain + theme_classic() + mytheme +
     scale_fill_viridis(option = "E", na.value = "white", limits = c(lwr_thresh, upr_thresh), oob = scales::squish) +
     common_text +
     ggtitle(label = map_text)
  )
  
  return(list(ggSuvr = ggSuvr))
}


maplist_abconv = ggSegIt(imap = abmap_conv,
                         map_text = "all TP average of AB conv obs (sample)",
                         lwr_thresh = 1, upr_thresh = 1.3)
maplist_abconvCOR = ggSegIt(imap = abmap_convCOR,
                            map_text = "all TP average of AB conv obs (sample, tracer-corrected)",
                            lwr_thresh = 1, upr_thresh = 1.3)
cowplot::plot_grid(maplist_abconv$ggSuvr, 
                   maplist_abconvCOR$ggSuvr,
                   nrow=2)


getCohens = function(imod) {
  iT = tidy(imod)[c("statistic")]
  iDF = df.residual(imod)
  iEffect = data.frame(tidy(imod)[c("term")], effectsize::t_to_d(iT$statistic, df = iDF))
  return(iEffect)
}


# cortex-wide tests ----
testCortex =function(df) {
  for (ii in 1:length(nrois)) {
    if (ii == 1) {
      modlist = effectlist = combinelist = list()
    }
    roi = c(lhrois, rhrois)[ii]
    print(paste(roi, ii))
    tmp = df
    tmp$suvr = tmp[[roi]]
    lmer_model <- lmer(
      suvr ~ convgroup + TRACER + (1 | PTID),
      data = tmp
    )
    
    # ggplot(tmp, aes(TRACER, suvr, col = factor(convgroup))) + geom_jitter() + geom_smooth(method = "lm")
    modlist[[ii]] = tidy(lmer_model)
    
    # get cohens d
    effectlist[[ii]] = getCohens(lmer_model)
    
    # combine stats
    combinelist[[ii]] = left_join(modlist[[ii]], effectlist[[ii]])
    
  }
  return(list(modlist = modlist, effectlist = effectlist, combinelist = combinelist))
}


# add FDR and significance indicators
ggStatWrangle = function(imap, statistical_test = T) {
  
  omap = imap
  # names(omap) = c("ROI","SUVR")
  omap$plotrois = tolower(omap$ROI)
  omap$plotrois = gsub("ctx_", "", omap$plotrois)
  omap$plotrois = gsub("_suvr", "", omap$plotrois)
  omap$hemi = substr(omap$plotrois, 1,2)
  
  omap = omap[omap$hemi == "lh" | omap$hemi == "rh",]
  omap$hemi = ifelse(omap$hemi == "lh", "left",
                     ifelse(omap$hemi == "rh", "right", NA))
  omap$plotrois[!omap$plotrois %in% dk$data$label]
  omap$plotrois[omap$plotrois %in% dk$data$label]
  omap$label = omap$plotrois
  
  
  if (statistical_test) {
    tmpbh = sgof::BH(c(omap$p.value
    ))
    fdthresh = max(tmpbh$data[tmpbh$Adjusted.pvalues < .05])
    
    omap$pCol = "white"
    omap$pInd = 0
    omap$pIndFDR = 0
    omap$FDRind = 0
    
    omap$FDRind = ifelse(omap$p.value <= fdthresh, 1, 0)
    omap$pInd[omap$p.value<.05] = 1
    omap$pCol[omap$p.value<.05] = "#aad6e7"
    
    omap$pIndFDR[omap$pInd == 1 & omap$FDRind != 1] = 1
    omap$pIndFDR[omap$pInd == 1 & omap$FDRind == 1] = 2
  }
  
  return(omap)
}



ggStatIt = function(imap, map_text, statistical_test = T, lwr_thr, upr_thr, lwr_d = 0, upr_d = 0.4) {
  metric="thickness"
  
  newAtlas = dk %>% 
    as_tibble() %>% 
    left_join(as.data.frame(imap))
  
  
  if (statistical_test) {
    
    newbrain = ggseg(.data = newAtlas, mapping = aes(fill = factor(pIndFDR)), col = "dark grey", size= 0.5, atlas = "dk")
    newbrainT = ggseg(.data = newAtlas, mapping = aes(fill = statistic), col = "dark grey", size= 0.5, atlas = "dk")
    newbrainP = ggseg(.data = newAtlas, mapping = aes(fill = -log10(p.value)), col = "dark grey", size= 0.5, atlas = "dk")
    newbrainD = ggseg(.data = newAtlas, mapping = aes(fill = d), col = "dark grey", size= 0.5, atlas = "dk")
    
    sighits = sum(imap$pIndFDR == 1)
    fdrhits = sum(imap$pIndFDR == 2)
    
    # cortical map pIndFDR -----
    (ggFacmap = newbrain + scale_fill_manual(values = c(
      "0" = "white",
      "1" = "#83b8c9",
      "2" = "black"
    ),
    na.value = "white"
    ) + theme_classic() + mytheme +
      common_text +
      ggtitle(label = paste0(map_text, " ", metric, " ",
                             sighits, "sig, ", fdrhits, "rejections")))
    
    (ggPmap = newbrainP + theme_classic() + mytheme +
        common_text +
        ggtitle(label = paste0(map_text, " ", metric, " ",
                               sighits, "sig, ", fdrhits, "rejections")))
    
    # cortical map T -----
    (ggStatmap = newbrainT + theme_classic() + mytheme +
       scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0, limits = c(lwr_thr, upr_thr), oob = scales::squish, na.value = "white") +
       common_text +
       ggtitle(label = paste0(map_text, " ", metric, " ",
                              sighits, "sig, ", fdrhits, "rejections")) +
       theme(legend.position = "top")
    )
    
    # cortical map Cohens D -----
    (ggDmap = newbrainD + theme_classic() + mytheme +
       scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0, limits = c(lwr_d, upr_d), oob = scales::squish, na.value = "white") +
       common_text +
       ggtitle(label = paste0(map_text, " ", metric, " ",
                              sighits, "sig, ", fdrhits, "rejections")) +
       theme(legend.position = "top")
    )
    
    return(list(ggFacmap = ggFacmap, ggStatmap = ggStatmap, ggPmap = ggPmap, ggDmap = ggDmap))
    
  } else {
    
    # if only correlation
    newbrain = ggseg(.data = newAtlas, mapping = aes(fill = r), col = "dark grey", size= 0.5, atlas = "dk")
    
    (ggStatmap = newbrain + theme_classic() + mytheme +
        scale_fill_viridis(option = "D", na.value = "white", limits = c(-1,1)) +
        common_text +
        ggtitle(label = paste0(map_text)) +
        theme(legend.position = "top")
    )
    
    return(list(ggStatmap = ggStatmap))
    
  }
  
}


# load in years2ab data -----
load(file.path(b, "reproduce/data/converters_data_for_plot_BACS_UPDATE_REPRO.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_ADNINC_UPDATE_REPRO.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_LCBC_UPDATE_REPRO.Rda"))


DF_CONV_COMBINE = rbind(DF_CONV_adni,
                        DF_CONV_bacs,
                        DF_CONV_lcbc)

DF_CONV_COMBINE$yearsBeforePredABpos = DF_CONV_COMBINE$age_at_threshold - DF_CONV_COMBINE$agePet


dim(DF_CONV_COMBINE) #260 regional
length(unique(DF_CONV_COMBINE$PTID)) #76 converters


# pred data
newdata = data.frame(
  yearsBeforePredABpos = seq(0, 10, by = 0.01),
  TRACER = "FBP",
  PTID = NA  # Dummy subject to avoid random intercepts
)


# linear mixed models of time-to-predicted Aβ positivity
# 68 regional Aβ SUVRs adjusted for tracer
rois = c(lhrois, rhrois)
for (ii in 1:length(rois)) {
  if (ii == 1) {
    modlist = list()
    predtraj = matrix(data = NA, nrow = nrow(newdata), ncol = length(rois))
  }
  print(ii)
  roi = rois[ii]
  lmer_model <- lmer(
    DF_CONV_COMBINE[[roi]] ~ yearsBeforePredABpos + TRACER + (1 | PTID),
    data = DF_CONV_COMBINE
  )
  # ggplot(ABlong_noNC, aes(yearsBeforePredABpos, ABlong_noNC[[rois[2]]])) + geom_point() + geom_smooth(method = "lm")
  modlist[[ii]] = summary(lmer_model)$coefficients
  predtraj[,ii] = predict(lmer_model, newdata = newdata, re.form = NA)
}
predtraj = as.data.frame(predtraj)
names(predtraj) = rois
predtraj = data.frame(yearsBeforePredABpos = newdata$yearsBeforePredABpos, 
                      predtraj)


# get predicted trajectory at given yearsBeforeAB (predicted)
predWrangle = function(years) {
  predtraj_out = predtraj[predtraj$yearsBeforePredABpos == years,]
  predtraj_out = data.frame(names(predtraj_out), unlist(unname(predtraj_out)))
  predtraj_out = predtraj_out[-1,] 
  row.names(predtraj_out) = NULL
  return(predtraj_out)
}
predtraj1 = ggWrangle(predWrangle(1))
predtraj2 = ggWrangle(predWrangle(2))
predtraj3 = ggWrangle(predWrangle(3))
predtraj4 = ggWrangle(predWrangle(4))
predtraj5 = ggWrangle(predWrangle(5))
predtraj6 = ggWrangle(predWrangle(6))
predtraj7 = ggWrangle(predWrangle(7))
predtraj8 = ggWrangle(predWrangle(8))
predtraj9 = ggWrangle(predWrangle(9))
predtraj10 = ggWrangle(predWrangle(10))


# check worked
roi_num = 1
ggplot(DF_CONV_COMBINE, aes(yearsBeforePredABpos, DF_CONV_COMBINE[[rois[roi_num]]])) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  geom_point(aes(x=1, y=predtraj1[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=2, y=predtraj2[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=3, y=predtraj3[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=4, y=predtraj4[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=5, y=predtraj5[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=6, y=predtraj6[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=7, y=predtraj7[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=8, y=predtraj8[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=9, y=predtraj9[roi_num,2]), size = 5, col = "red") +
  geom_point(aes(x=10, y=predtraj10[roi_num,2]), size = 5, col = "red") +
  theme_classic() + mytheme +
  ylab("SUVR") +
  ggtitle(rois[roi_num])



# now make the maps
maplistpred1 = ggSegIt(imap = predtraj1 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 1 yearsBeforePredAB+",
                       1, 1.2)
maplistpred2 = ggSegIt(imap = predtraj2 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 2 yearsBeforePredAB+",
                       1, 1.2)
maplistpred3 = ggSegIt(imap = predtraj3 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 3 yearsBeforePredAB+",
                       1, 1.2)
maplistpred4 = ggSegIt(imap = predtraj4 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 4 yearsBeforePredAB+",
                       1, 1.2)
maplistpred5 = ggSegIt(imap = predtraj5 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 5 yearsBeforePredAB+",
                       1, 1.2)
maplistpred6 = ggSegIt(imap = predtraj6 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 6 yearsBeforePredAB+",
                       1, 1.2)
maplistpred7 = ggSegIt(imap = predtraj7 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 7 yearsBeforePredAB+",
                       1, 1.2)
maplistpred8 = ggSegIt(imap = predtraj8 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 8 yearsBeforePredAB+",
                       1, 1.2)
maplistpred9 = ggSegIt(imap = predtraj9 %>% 
                         mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                ),
                       map_text = "Predicted AB at 9 yearsBeforePredAB+",
                       1, 1.2)
maplistpred10 = ggSegIt(imap = predtraj10 %>% 
                          mutate(SUVR = ifelse(SUVR < 1, NA, SUVR)
                                 ),
                        map_text = "Predicted AB at 10 yearsBeforePredAB+",
                        1, 1.2)


predtrajOut = rbind(
  predtraj1 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 1),
  predtraj2 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 2),
  predtraj3 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 3),
  predtraj4 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 4),
  predtraj5 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 5),
  predtraj6 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 6),
  predtraj7 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 7),
  predtraj8 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 8),
  predtraj9 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 9),
  predtraj10 %>%
    mutate(SUVR_above_reference = ifelse(SUVR < 1, NA, SUVR),
           predictedYearsBeforepredAB = 10)
)
# write.table(predtrajOut, file = file.path(b, "reproduce/data/cortmapPredTraj.csv"),
#             quote = F, row.names = F, col.names = T)
# write.table(predtraj, file = file.path(b, "reproduce/data/predTraj.csv"),
#             quote = F, row.names = F, col.names = T)


# Fig 5c combine maps predicted -----------
(ggPredTraj = ggarrange(
  maplistpred1$ggSuvr + ggtitle(NULL),
  maplistpred2$ggSuvr + ggtitle(NULL),
  maplistpred3$ggSuvr + ggtitle(NULL),
  maplistpred4$ggSuvr + ggtitle(NULL),
  maplistpred5$ggSuvr + ggtitle(NULL),
  maplistpred6$ggSuvr + ggtitle(NULL),
  maplistpred7$ggSuvr + ggtitle(NULL),
  maplistpred8$ggSuvr + ggtitle(NULL),
  maplistpred9$ggSuvr + ggtitle(NULL),
  maplistpred10$ggSuvr + ggtitle(NULL),
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  align = "hv"
))

# ggsave(plot = ggPredTraj,
#        filename = paste0(file.path(b, "paper2/figs_yearsBeforeAB/p_Seg_predictedABtrajConverters.png")),
#        width=20, height=30, units="cm", dpi=600,
#        background = "white")



# Rank order of regions crossing the 1.0 threshold ----
# (uptake greater than the reference region)
trajectory_cols <- setdiff(names(predtraj), "yearsBeforePredABpos")
for (jjj in 2:length(predtraj)) {
  
  if (jjj == 2) {
    rankCols = c()
    age_crossing = c()
  }
  
  tmp_col = predtraj[,jjj]
  idx_crossing_1 = max(which(tmp_col > 1))
  age_crossing_1 = predtraj$yearsBeforePredABpos[idx_crossing_1]

  if (idx_crossing_1 == length(tmp_col)) {
    rankCols[jjj] = "always"
    age_crossing[jjj] = "always"
  } else if (is.infinite(idx_crossing_1)) {
    rankCols[jjj] = "never"
    age_crossing[jjj] = "never"
  } else {
    print(tmp_col[(idx_crossing_1-1):(idx_crossing_1+1)])
    rankCols[jjj] = as.character(idx_crossing_1)
    age_crossing[jjj] = age_crossing_1
  }
}
names(predtraj)[which(rankCols == "always")]


# NB! if there are competing positions you need to expand the prediction grid
table(rankCols)


# you have an amyloid ordering
rankTraj = data.frame(rois = trajectory_cols,
                      rankedTraj = rankCols[2:length(rankCols)],
                      age_crossing = age_crossing[2:length(rankCols)])

main_ROIS = c(
  rankTraj$rois[grepl("FRON", rankTraj$rois)],
  rankTraj$rois[grepl("PARS", rankTraj$rois)],
  rankTraj$rois[grepl("PRECENTRAL", rankTraj$rois)],
  rankTraj$rois[grepl("PARACENTRAL", rankTraj$rois)],
  rankTraj$rois[grepl("INSULA", rankTraj$rois)],
  rankTraj$rois[grepl("SUPRAMA", rankTraj$rois)],
  rankTraj$rois[grepl("BANKSS", rankTraj$rois)]
)
rankedROIs = rankTraj %>% filter(rois %in% main_ROIS)
rankedROIs$rankedTraj[rankedROIs$rankedTraj == "never"] = NA
rankedROIs$rankedTraj = as.numeric(rankedROIs$rankedTraj)


# reverse the rank ordering of indices (because yearsBeforeAB was positive in model)
rankedROIs$revRankedTraj <- ifelse(
  is.na(rankedROIs$rankedTraj),
  NA,
  rank(-rankedROIs$rankedTraj, ties.method = "first")
)
rankedROIs %<>% arrange(revRankedTraj)
head(rankedROIs)
# write ranked amyloid list
# write.table(rankedROIs, file.path(b, "reproduce/data/rankABpredTraj.csv"), col.names = T, row.names = F, quote = F)



# group difference in SUVR between groups in your sample ----
modlist = testCortex(df_AB)[[3]]
cortmap_suvrdiff = bind_rows(modlist, .id = "roi") %>% filter(term == "convgroup")
cortmap_suvrdiff$ROI = c(lhrois, rhrois)
cortmap_suvrdiff$group = NA
cortmap_suvrdiff = ggStatWrangle(cortmap_suvrdiff)
names(cortmap_suvrdiff)[1] = "model" #NB! ggseg not working because "roi" in df


maplist_suvrdiff = ggStatIt(cortmap_suvrdiff, "test", lwr_thr = min(cortmap_suvrdiff$statistic), upr_thr = max(cortmap_suvrdiff$statistic),
                            lwr_d = min(cortmap_suvrdiff$d), upr_d = max(cortmap_suvrdiff$d))

maplist_suvrdiff$ggDmap
maplist_suvrdiff$ggFacmap
maplist_suvrdiff$ggPmap
maplist_suvrdiff$ggStatmap



# group difference in SUVR between ab- (your sample) and ab+ (diagnosed adni sample) ----
tmp_ABdiagnosed = AB_posnoNC %>% mutate(convgroup = "diagnosed", cohort = "ADNIdiagnosed") %>% select(all_of(names(tmp_NONCONV)))
tmp_CONV
tmp_NONCONV
tmp_ABdiagnosedVnonconv = rbind(tmp_ABdiagnosed %>% mutate(convgroup = 1),
                                tmp_NONCONV %>% mutate(convgroup = 0))

modlist = testCortex(tmp_ABdiagnosedVnonconv)[[3]]
cortmap_suvrdiff_ABdiagnosed = bind_rows(modlist, .id = "roi") %>% filter(term == "convgroup")
cortmap_suvrdiff_ABdiagnosed$ROI = c(lhrois, rhrois)
cortmap_suvrdiff_ABdiagnosed$group = NA
cortmap_suvrdiff_ABdiagnosed = ggStatWrangle(cortmap_suvrdiff_ABdiagnosed)
names(cortmap_suvrdiff_ABdiagnosed)[1] = "model" #NB! ggseg not working because "roi" in df


maplist_suvrdiff_ABdiagnosed = ggStatIt(cortmap_suvrdiff_ABdiagnosed, "test", lwr_thr = min(cortmap_suvrdiff_ABdiagnosed$statistic), upr_thr = max(cortmap_suvrdiff_ABdiagnosed$statistic),
                                       lwr_d = min(cortmap_suvrdiff_ABdiagnosed$d), upr_d = 1.5)
maplist_suvrdiff_ABdiagnosed$ggDmap
maplist_suvrdiff_ABdiagnosed$ggFacmap
maplist_suvrdiff_ABdiagnosed$ggPmap
maplist_suvrdiff_ABdiagnosed$ggStatmap


plot(abmap_conv$SUVR, mapABpos_noNC$SUVR)
plot(abmap_nonconv$SUVR, mapABpos_noNC$SUVR)
cor(abmap_nonconv$SUVR, mapABpos_noNC$SUVR)
cor.test(abmap_nonconv$SUVR, mapABpos_noNC$SUVR, method = "spearman")


# writeMaps = F
# if (writeMaps) {
  # path_normmaps = file.path(b, "results/normMaps")
  # write.table(cortmap_suvrdiff, file.path(path_normmaps, "amyloid_map_convVabminus_suvrdiff_d.csv"), row.names = F, col.names = T, quote = F)
  # write.table(cortmap_suvrdiff_adnifucked, file.path(path_normmaps, "amyloid_map_ABdiagnosedVnonconv_suvrdiff_d.csv"), row.names = F, col.names = T, quote = F)
  # write.table(mapABpos_noNC, file.path(path_normmaps, "amyloid_map_ABdiagnosed_suvrmean.csv"), row.names = F, col.names = T, quote = F)
  # write.table(abmap_conv, file.path(path_normmaps, "amyloid_map_conv_suvrMean.csv"), row.names = F, col.names = T, quote = F)
  # write.table(abmap_nonconv, file.path(path_normmaps, "amyloid_map_nonconv_suvrMean.csv"), row.names = F, col.names = T, quote = F)
  
  
  
  # ggsave(plot = maplist_suvrdiff_ABdiagnosed$ggStatmap,
  #        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_Seg_Tsuvrdiff_ABpos_noNC-v-nonconv.png"),
  #        width=20, height=7, units="cm", dpi=600)
  # 
  # ggsave(plot = maplist_suvrdiff$ggStatmap,
  #        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_Seg_Tsuvrdiff_insample_convnonconv.png"),
  #        width=20, height=7, units="cm", dpi=600)
  # 
  # ggsave(plot = maplist_suvrdiff$ggD,
  #        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_Seg_Dsuvrdiff_insample_convnonconv.png"),
  #        width=20, height=7, units="cm", dpi=600)
  # 
  # 
  # ggsave(plot = maplist_suvrdiff_ABdiagnosed$ggDmap,
  #        filename = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_Seg_Dsuvrdiff_ABpos_noNC-v-nonconv.png"),
  #        width=20, height=7, units="cm", dpi=600)
# }

