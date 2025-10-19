#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose:
#========================================================================================#


#================ INPUTS =================#
parse_args <- function(args) {
  if (length(args) != 8) {
    stop("Wrong number of arguments provided. Expected 8, received ", length(args))
  }
  
  list(
    analysnum = as.character(args[1]),
    analysname = as.character(args[2]),
    cohortSelect = as.character(args[3]),
    yearsBeforeAB = as.integer(args[4]),
    matchFollow = as.logical(as.integer(args[5])),
    procType = as.character(args[6]),
    predAB = as.logical(as.integer(args[7])),
    singleT = as.logical(as.integer(args[8]))
  )
}


args <- commandArgs(TRUE)
print(args)
input_params <- parse_args(args)
list2env(input_params, envir = .GlobalEnv)
cat(sprintf("analysnum = %s\n", analysnum))
cat(sprintf("analysname = %s\n", analysname))
cat(sprintf("cohortSelect = %s\n", cohortSelect))
cat(sprintf("yearsBeforeAB = %d\n", yearsBeforeAB))
cat(sprintf("matchFollow = %s\n", matchFollow))
cat(sprintf("procType = %s\n", procType))
cat(sprintf("predAB = %s\n", predAB))
cat(sprintf("singleT = %s\n", singleT))
#=========================================#


# # testing ---------
# rm(list=ls())
# print("---------_")
# print(cohortSelect)
# print("---------_")

# # input tests
analysnum=1
cohortSelect = "adnibacslcbc" #"adnibacs" "bacslcbc" "adnilcbc"
yearsBeforeAB = 7
matchFollow = 0
procType = "LONG"
predAB = 0
singleT = 0
analysname = paste0(
  "analysnum",
  analysnum,
  "matchF",
  matchFollow,
  "proc",
  procType,
  "predAB",
  predAB,
  "singleT",
  singleT,
  "reproduce1"
)


analysisString = "tesla"
savefigs=T
nTime=2
agecut=30
saveres=T
scaleRes = F
writeMaps = 1


#---load packages
loadPackages = function() {
  packages = c("here", "tidyverse","magrittr","gamm4","itsadug","numDeriv","gratia","mgcv","viridis","wesanderson","asbio","broom","cowplot","data.table","stringi","ggseg")
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  print(sapply(packages, require, character.only = T))
  print(sapply(packages, function(p) as.character(packageVersion(p))))
}
loadPackages()
# here()



#---set dir
b = "/cluster/projects/p274/projects/p040-ad_change/Berkeley"
#---set dir


#---make dirstruct
setwd(b)
resdir = file.path(b,
                   paste0("results/prepSlopes/", analysname))
if (! dir.exists(resdir)) {
  stop(paste("ERROR:", resdir, "does not exist. Quitting"))
}

ifile =
  file.path(
    resdir,
    paste0(
      "results.prepSlopes",
      cohortSelect,
      agecut,
      "_",
      analysisString,
      "yearsBeforeAB",
      yearsBeforeAB,
      "_",
      analysname,
      ".Rda"
    )
  )



#load and combine data
load(file.path(b,"reproduce/data/DF_BACSUPDATE_PREPPED_AB.Rda"))
load(file.path(b,"DF_ADNIUPDATE_PREPPED_AB.Rda"))
ADNI = ADNIEXPORT
BACS = BACSEXPORT
allFeat = readLines(file.path(b, "reproduce/data/allFeatures364.txt"))
adnioutlier = "029_S_0845"


if (procType == "LONG") {
  load(file.path(b, "reproduce/data/DF_LONG_4570.Rda"))
} else if (procType == "CROSS") {
  load(file.path(b, "reproduce/data/DF_CROSS_4570.Rda"))
}
dim(DF); length(unique(DF$subject_id))



# load converter/nonconverter data
load(file.path(b, "reproduce/data/converters_all_negfirst_ADNINC_UPDATE.Rda"))
load(file.path(b, "reproduce/data/converters_all_negfirst_BACS_UPDATE.Rda"))
load(file.path(b, "reproduce/data/converters_all_negfirst_LCBC_UPDATE.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_BACS_UPDATE_REPRO.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_ADNINC_UPDATE_REPRO.Rda"))
load(file.path(b, "reproduce/data/converters_data_for_plot_LCBC_UPDATE_REPRO.Rda"))
converters_both_negfirst_ADNINC_plotdat$SID_scan = converters_both_negfirst_ADNINC_plotdat$PTID_scan
nonconverters_both_LCBC_plotdat$SID_scan = nonconverters_both_LCBC_plotdat$subject_id_scan
AB_nonconv_LCBC_plotdat$SID_scan = AB_nonconv_LCBC_plotdat$subject_id_scan
identical(converters_allAB_negfirst_ADNINC$newlink, converters_allAB_negfirst_ADNINC_plotdat$newlink)
identical(converters_allAB_negfirst_LCBC$imageLink, converters_allAB_negfirst_LCBC_plotdat$Folder)
identical(converters_allAB_negfirst_BACS$imageLink, converters_allAB_negfirst_BACS_plotdat$imageLink)



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


# ------ load prepSlopes ------
print("----------------------")
print(paste("loading", ifile))
load(ifile)
U.all$subject_id[U.all$cohort=="ADNINC"] = gsub("-", "_", U.all$subject_id[U.all$cohort=="ADNINC"])
ls(pattern = "\\.all")
print("----------------------")




# load LCBC amyloid data ------
lcbcABnew = 1
if (lcbcABnew == 1) {

  df_AB_LCBC = read.csv(file.path(b, "amyloid_to_p274/AB_classification.tsv"), sep = "\t")
  df_AB_LCBC$subject_id = substr(df_AB_LCBC$raw_name, 1, 7)
  df_AB_LCBC %<>% group_by(subject_id) %>% mutate(nTimesAB = n(), nTimeAB = row_number(), converter = ifelse(length(unique(AB_class))>1, "converter", "nonconverter") )
  length(unique(df_AB_LCBC$subject_id))
  df_AB_LCBC %<>% group_by(subject_id) %>% mutate(ABlong = ifelse(any(AB_class == 1), 1, 0))
  names(df_AB_LCBC)[1] = "Folder"

  # pet timestamps
  pet_times = fread(file.path(b, "pet_timestamps.txt")) %>% rename(Folder = SUBJECTID)
  pet_times$subject_id = as.character(pet_times$subject_id)
  df_AB_LCBC = left_join(df_AB_LCBC, pet_times)
  
  # add centiloids
  df_AB_LCBC$eq1 = (df_AB_LCBC$mean_uptake_values - 0.006) / 0.969
  df_AB_LCBC$centiloids = 100*(df_AB_LCBC$eq1 - 1.176)/(2.439-1.176)
  plot(x = df_AB_LCBC$mean_uptake_values, y = df_AB_LCBC$centiloids)
  
  # mean centiloids
  df_AB_LCBC %<>% group_by(subject_id) %>% mutate(meanCentiloids = mean(centiloids)) %>% ungroup()
  U_AB_LCBC = df_AB_LCBC %>% filter(nTimeAB == 1)
  U_AB_LCBC %<>% mutate(convertgroup =
                         case_when(
                           converter == "nonconverter" & ABlong == 1 ~ "nonconverter_pos",
                           converter == "nonconverter" & ABlong == 0 ~ "nonconverter_neg",
                           converter == "converter" & ABlong == 1 ~ "converter",
                           converter == "converter" & ABlong == 0 ~ "check",
                         ))
}


# prep brainslopes by converter group -------
converterlist = c(converters_allAB_negfirst_ADNINC$PTID, converters_allAB_negfirst_BACS$SID, converters_allAB_negfirst_LCBC$subject_id) %>% unique()
converterlist %in% U.all$subject_id %>% sum() #- length(converterlist) #lost due to long constraint
groupconv = U.all[U.all$subject_id %in% converterlist,]

neverABlist = c(
  gsub("-","_",U_AB_ADNI$PTID[U_AB_ADNI$longABGroup == "nonconverter_neg"]), 
  U_AB_BACS$SID[U_AB_BACS$longABGroup == "nonconverter_neg"],
  U_AB_LCBC$subject_id[U_AB_LCBC$convertgroup == "nonconverter_neg"]
)

alwaysABlist = c(
  gsub("-","_",U_AB_ADNI$PTID[U_AB_ADNI$longABGroup == "nonconverter_pos"]), 
  U_AB_BACS$SID[U_AB_BACS$longABGroup == "nonconverter_pos"],
  U_AB_LCBC$subject_id[U_AB_LCBC$convertgroup == "nonconverter_pos"]
)

groupnever = U.all[U.all$subject_id %in% neverABlist,]
groupconv = U.all[U.all$subject_id %in% groupconv$subject_id,]
groupalways = U.all[U.all$subject_id %in% alwaysABlist,]


# sanity check
# group differences in hippocampal vol 
t.test(groupnever$vol.Left.Hippocampus, groupconv$vol.Left.Hippocampus)
t.test(groupnever$vol.Left.Hippocampus, groupalways$vol.Left.Hippocampus)
t.test(groupconv$vol.Left.Hippocampus, groupalways$vol.Left.Hippocampus)
t.test(groupnever$vol.Right.Hippocampus, groupconv$vol.Right.Hippocampus)
t.test(groupnever$vol.Right.Hippocampus, groupalways$vol.Right.Hippocampus)
t.test(groupconv$vol.Right.Hippocampus, groupalways$vol.Right.Hippocampus)


# ensure adnioutlier is removed for all tests
neverABlist = neverABlist[neverABlist != adnioutlier]
groupconv = U.all[U.all$subject_id %in% converterlist,]
groupnever = U.all[U.all$subject_id %in% neverABlist,]


# brain slopes by group
nrow(U.all) == nrow(rSlopeMat.all) #TRUE
groupconv_slopes = rSlopeMat.all[U.all$subject_id %in% converterlist,]
groupnever_slopes = rSlopeMat.all[U.all$subject_id %in% neverABlist,]
groupalways_slopes = rSlopeMat.all[U.all$subject_id %in% alwaysABlist,]


# recalculate meanCentiloids across years before AB+ -------
converters_allAB_negfirst_ADNINC_beforeAB <- converters_allAB_negfirst_ADNINC %>%
  group_by(PTID) %>%
  filter(AgeAtABPET < first(ageAtFirstABpos)) %>%
  mutate(meanCentiloidsBeforeAB = mean(CENTILOIDS)) %>% 
  ungroup()


converters_allAB_negfirst_BACS_beforeAB <- converters_allAB_negfirst_BACS %>%
  group_by(SID) %>%
  filter(agePet < first(ageAtFirstABpos)) %>%
  mutate(meanCentiloidsBeforeAB = mean(CENTILOID)) %>% 
  ungroup()


converters_allAB_negfirst_LCBC_beforeAB <- converters_allAB_negfirst_LCBC %>%
  group_by(subject_id) %>%
  filter(ageAtPet < first(ageAtFirstABpos)) %>%
  mutate(meanCentiloidsBeforeAB = first(centiloids)) %>%
  ungroup()



link_meancent_beforeAB = rbind(
  converters_allAB_negfirst_ADNINC_beforeAB %>% dplyr::select(PTID, "meanCentiloidsBeforeAB", "ageAtFirstABpos") %>% distinct() %>% rename(subject_id = PTID),
  converters_allAB_negfirst_BACS_beforeAB %>% dplyr::select(SID, "meanCentiloidsBeforeAB", "ageAtFirstABpos") %>% distinct() %>% rename(subject_id = SID),
  converters_allAB_negfirst_LCBC_beforeAB %>% dplyr::select(subject_id, "meanCentiloidsBeforeAB", "ageAtFirstABpos") %>% distinct()
)


link_meancent_all = rbind(
  U_AB_ADNI %>% dplyr::select(PTID, contains("meanCentiloids")) %>% distinct() %>% rename(subject_id = PTID) %>% mutate(subject_id = gsub("-", "_", subject_id)),
  U_AB_BACS %>% dplyr::select(SID, contains("meanCentiloids")) %>% distinct() %>% rename(subject_id = SID),
  U_AB_LCBC %>% dplyr::select(subject_id, contains("meanCentiloids")) %>% distinct()
)


link_meancent = left_join(link_meancent_all, link_meancent_beforeAB)


# ROIs.all is order of vars in rSlopeMat.all
corthickrois = ROIs.all[grepl("thickness", ROIs.all, ignore.case = T)]
rois = corthickrois
  

# add in PC1 (removing mean)
tmprois = rois[!grepl("mean", rois, ignore.case = T)]
tmproisL = tmprois[grep("lh_", tmprois)]
tmproisR = tmprois[grep("rh_", tmprois)]


# whole brain
pca1 = prcomp(rSlopeMat.all[,which(ROIs.all %in% tmprois)],center=T,scale.=T)
pcaL = prcomp(rSlopeMat.all[,which(ROIs.all %in% tmproisL)],center=T,scale.=T)
pcaR = prcomp(rSlopeMat.all[,which(ROIs.all %in% tmproisR)],center=T,scale.=T)
(spca1 = summary(pca1))
(spcaL = summary(pcaL))
(spcaR = summary(pcaR))
pca_Full = pca1$x[,1]
pca_Left = pcaL$x[,1]
pca_Right = pcaR$x[,1]


# LOOP ROIS ------
for (ii in 1:length(rois)) {
  
  if (ii == 1) {
    
    # Dynamically create variables
    loopVar = function(varstring, num, type = "l", pos = "suffix") {
      template = switch(type,
                        l = list(),
                        v = c(),
                        stop("Invalid type: must be 'l (for list)' or 'v (for vector)'")
      )
      for (n in seq_len(num)) {
        # Try to insert the number between base and suffix if "_" is present
        if (grepl("_", varstring)) {
          parts = strsplit(varstring, "_", fixed = TRUE)[[1]]
          var_name = paste0(parts[1], n, "_", parts[2])
        } else {
          var_name = paste0(varstring, n)  # fallback to suffix
        }
        assign(var_name, template, envir = .GlobalEnv)
        message(sprintf("Created variable '%s' as %s", var_name, type))
      }
    }
    loopVar("modlist", 7, 'l')
    loopVar("modlist_cent", 7, 'l')
    loopVar("cor", 3, 'v')
    tmpp = tmpp_cent = c()
    p1 = pAll1 = list()
    convplotmodYintonly = convplotmodYmanchange = convplotmodmain = convplotmodYmanchangenoout = convplotslopeoys = convplotslointooys = list()
  }
  
  
  print(ii)
  roi = rois[ii]
  
  
  # overwrite with full data each time
  U.all.tmp = U.all %>% select(1:379)
  U.all.tmp = left_join(U.all.tmp, link_meancent)
  
  
  # add in brainslopes
  U.all.tmp$brainvarSlope = rSlopeMat.all[,ROIs.all == roi]
  U.all.tmp$brainvarAbs = absChangeMat.all[,ROIs.all == roi]
  
  
  tmp_conv = U.all.tmp %>% filter(subject_id %in% DF.conv.all$subject_id) %>% mutate(convgroup = 1)
  tmp_conv$meanCentiloids
  tmp_conv$meanCentiloidsBeforeAB
  
  
  tmp_never = U.all.tmp %>% filter(subject_id %in% neverABlist) %>% mutate(convgroup = 0)
  tmp_never$meanCentiloids
  tmp_always = U.all.tmp %>% filter(subject_id %in% alwaysABlist) %>% mutate(convgroup = 1)
  tmp_always$meanCentiloids
  
  
  # overwrite meanCentiloids with those calculated before AB+ (conv group)
  # setup group contrast
  U.all.tmp = rbind(
    tmp_conv %>% mutate(meanCentiloids = meanCentiloidsBeforeAB), 
    tmp_never)
  table(U.all.tmp$convgroup)
  
  
  # MAIN LINEAR MODELS ------
  mod = lm((brainvarSlope) ~ (convgroup) + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = U.all.tmp)
  mod_cent = lm((brainvarSlope) ~ (convgroup) + meanCentiloids + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = U.all.tmp)
  
  
  # track order of subs in model
  subs_order = U.all.tmp$subject_id
  modlist1[[ii]] = tidy(mod)
  modlist1_cent[[ii]] = tidy(mod_cent)
  tmpp[ii] = modlist1[[ii]]$p.value[2]
  tmpp_cent[ii] = modlist1_cent[[ii]]$p.value[2]
  
  
  # model predictions for plotting -------
  # raw
  predlm1 = predict(
    mod,
    newdata = data.frame(
      convgroup = U.all.tmp$convgroup,
      meanAge = mean(U.all.tmp$meanAge),
      cohort = "ADNINC",
      # cohort = "LCBC",
      subject_sex = "Female",
      nTimepoints = median(U.all.tmp$nTimepoints),
      intervalFirstlast = mean(U.all.tmp$intervalFirstlast)
    ),
    se.fit = T
  )
  
  # centiloid-corrected
  predlm1_cent = predict(
    mod_cent,
    newdata = data.frame(
      convgroup = U.all.tmp$convgroup,
      meanAge = mean(U.all.tmp$meanAge),
      meanCentiloids = mean(U.all.tmp$meanCentiloids),
      cohort = "ADNINC",
      # cohort = "LCBC",
      subject_sex = "Female",
      nTimepoints = median(U.all.tmp$nTimepoints),
      intervalFirstlast = mean(U.all.tmp$intervalFirstlast)
    ),
    se.fit = T
  )
  
  # residualized effects for plotting
  U.all.tmp %<>% mutate(predY1 = predlm1$fit,
                        predSE1 = predlm1$se.fit,
                        predCI1 = predlm1$se.fit*1.96,
                        resid1 = residuals(mod),
                        presid1 = predY1 + resid1,
                   
                   predY1_cent = predlm1_cent$fit,
                   predSE1_cent = predlm1_cent$se.fit,
                   predCI1_cent = predlm1_cent$se.fit*1.96,
                   resid1_cent = residuals(mod_cent),
                   presid1_cent = predY1_cent + resid1_cent,
                   X1 = U.all.tmp$convgroup)
  
  
  if (exists("means1")) {
    rm(means1); rm(means2)
  }
  means0 = mean(U.all.tmp$presid1[U.all.tmp$convgroup == 0])
  means1 = mean(U.all.tmp$presid1[U.all.tmp$convgroup == 1])
  
  
  pal = wesanderson::wes_palettes$FantasticFox1
  mytheme = theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    title = element_text(size=17),
    text = element_text(color = "black", size = 18, family="Nimbus Sans Narrow"),
    plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(),
    axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
    axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,20,0)),
    axis.text = element_text(color = "black", size = 18),
    legend.key.size = unit(1,"cm"))
  
  
  p1[[ii]] = ggplot(U.all.tmp,
                    aes(x=(factor(convgroup)), y=(presid1), col=factor(convgroup))) +
    geom_hline(yintercept = 0, linetype=2, col="light grey") +
    geom_jitter(alpha=0.7, position = position_jitter(seed = 123,width = 0.1)) +
    geom_segment(aes(x=0.75, y=means0, xend=1.25, yend=means0), size=.75, color="black") +
    geom_segment(aes(x=1.75, y=means1, xend=2.25, yend=means1), size=.75, color="black") +
    theme_classic() + mytheme +
    ylab(expression("Additional Change (mm/year)")) +
    labs(x = NULL) +
    scale_color_manual(values = c(pal[2], pal[5])) + #fantastic fox
    theme(legend.position = "top",
          axis.ticks= element_line())
  
  
  res0 = U.all.tmp[U.all.tmp$convgroup == 0,]
  res1 = U.all.tmp[U.all.tmp$convgroup == 1,]
  res0$zero=0
  res1$zero=0
  
  ggdens <- ggplot() +
      PupillometryR::geom_flat_violin(data=res0,aes(x=zero,y=presid1),alpha = 0.3, col = pal[2], fill = pal[2]) +
      PupillometryR::geom_flat_violin(data=res1,aes(x=zero,y=presid1),alpha = 0.3, col = pal[5], fill = pal[5]) +
      theme_void()
  
  
  (pAll1[[ii]] = p1[[ii]] + ggdens +
      patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
  )
  
  
  plotdir = file.path(b, paste0("plots/convplots/", analysname, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime))
  if (savefigs) {
    
    if (!dir.exists(plotdir)) { 
      system(paste("mkdir -p", plotdir))
    }
    
    sigflag = ifelse(tmpp[ii] < .05, 1, 0)
    filename = paste0(plotdir, "/pconv_cortrois-",ii,"-sig", sigflag,"-",rois[ii],"_scale", scaleRes, ".png")
    print(paste("saving", filename))
    print("-----------------------")
    ggsave(plot = pAll1[[ii]],
           filename = filename,
           width=9, height=9, units="cm", dpi=600)
  }
  
  
  
  # load in roi-specific manual slopes (not written to matrix)
  resdir = file.path(b, paste0("results/prepMods/", analysname, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime))
  Uconv = fread(list.files(resdir, pattern=paste0("^UconvNout.*", roi), full.names = T))
  Uconv = Uconv %>% select(subject_id, contains(c("meanbrain", "manual", "rInt", "rSlope")))
  Uconv %<>% filter(subject_id %in% U.all.tmp$subject_id)
  Uconv = Uconv[match(U.all.tmp$subject_id, Uconv$subject_id), ]
  identical(U.all.tmp$subject_id, Uconv$subject_id)
  
  
  # confirms it is the same and the values are roi-specific
  # note that U.all.tmp manual slope values will be the last ROI in loop
  cor(U.all.tmp$brainvarSlope, Uconv$rSlope)
  
  
  # var correlations across conv data
  cor1[ii] = cor(Uconv$rSlope, Uconv$rInt)
  cor2[ii] = cor(Uconv$rSlope, Uconv$rIntonly)
  cor3[ii] = cor(Uconv$rIntonly, Uconv$rSlopeOys)
  
  
  # remove predictions from data
  U.all.tmp = merge(U.all.tmp %>% select(-contains("meanbrain"),
                                         -starts_with("pred"),
                                         -starts_with("presid")), Uconv %>% select(-rSlope))
  
  
  # reorder data again to be same as in first model
  U.all.tmp <- U.all.tmp[match(subs_order, U.all.tmp$subject_id), ]
  U.all.tmp$subject_id == subs_order
  
  
  # outliers on manual_slope
  U.all.tmp$Q = as.numeric(scale(U.all.tmp$manual_slope))
  outlierThresh = 5
  range(U.all.tmp$Q)
  U.all.tmp.noout = U.all.tmp %>% filter( abs(Q) < outlierThresh)
  range(U.all.tmp.noout$Q)

    
  # check main model
  modmaincheck = lm((brainvarSlope) ~ (convgroup) + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = U.all.tmp)
  
  
  # SUBSEQUENT LINEAR MODELS ------
  
  # model 2 Yslopenoage (no additional age correction -confirmed results same) ---------
  modYslopenoage = lm((brainvarSlope) ~ (convgroup) + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  modYslopenoage_cent = lm((brainvarSlope) ~ (convgroup) + meanCentiloids + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  
  
  # model 3 Yintonly (intercept-only models) ---------
  modYintonly = lm((rIntonly) ~ (convgroup) + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = U.all.tmp)
  modYintonly_cent = lm((rIntonly) ~ (convgroup) + meanCentiloids + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = U.all.tmp)
  
  
  # model 4 Yslopeoys (T2 gamm slopes) ---------
  modYslopeoys = lm((rSlopeOys) ~ (convgroup) + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  modYslopeoys_cent = lm((rSlopeOys) ~ (convgroup) + meanCentiloids + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
   
  
  # model 5 Yslointooys (T2 gamm slopes intercept-corrected) ---------
  modYslointooys = lm((rSlopeOys) ~ (convgroup) + rIntonly + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  modYslointooys_cent = lm((rSlopeOys) ~ (convgroup) + rIntonly + meanCentiloids + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  
  
  # model 6 Ymanchange (manual linear change across GAMM residuals) ---------
  modYmanchange = lm((manual_slope) ~ (convgroup) + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  modYmanchange_cent = lm((manual_slope) ~ (convgroup) + meanCentiloids + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp)
  
  
  # model 7 Ymanchange (manual linear change across GAMM residuals - outliers removed)  ---------
  modYmanchangenoout = lm((manual_slope) ~ (convgroup) + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp.noout)
  modYmanchangenoout_cent = lm((manual_slope) ~ (convgroup) + meanCentiloids + subject_sex + cohort + nTimepoints + intervalFirstlast, data = U.all.tmp.noout)
  
  
  modlist2[[ii]] = tidy(modYslopenoage)
  modlist2_cent[[ii]] = tidy(modYslopenoage_cent)
  
  modlist3[[ii]] = tidy(modYintonly)
  modlist3_cent[[ii]] = tidy(modYintonly_cent)
  
  modlist4[[ii]] = tidy(modYslopeoys)
  modlist4_cent[[ii]] = tidy(modYslopeoys_cent)
  
  modlist5[[ii]] = tidy(modYslointooys)
  modlist5_cent[[ii]] = tidy(modYslointooys_cent)
  
  modlist6[[ii]] = tidy(modYmanchange)
  modlist6_cent[[ii]] = tidy(modYmanchange_cent)
  
  modlist7[[ii]] = tidy(modYmanchangenoout)
  modlist7_cent[[ii]] = tidy(modYmanchangenoout_cent)
  
  
  
  
  
  # calculate model plot for variable model inputs
  predictPlot = function(df, model, model_cent, plotstring, intcorrect = F) {
    
    tmp = df
    
    if (intcorrect == F) {
      print("not additionally correcting data for intercept")
      predlm = predict(
        model,
        newdata = data.frame(
          convgroup = tmp$convgroup,
          meanAge = mean(tmp$meanAge),
          cohort = "ADNINC",
          # cohort = "LCBC",
          subject_sex = "Female",
          nTimepoints = median(tmp$nTimepoints),
          intervalFirstlast = mean(tmp$intervalFirstlast)
        ),
        se.fit = T
      )
      
      predlm_cent = predict(
        model_cent,
        newdata = data.frame(
          convgroup = tmp$convgroup,
          meanAge = mean(tmp$meanAge),
          meanCentiloids = mean(tmp$meanCentiloids),
          # cohort = "LCBC",
          cohort = "ADNINC",
          subject_sex = "Female",
          nTimepoints = median(tmp$nTimepoints),
          intervalFirstlast = mean(tmp$intervalFirstlast)
        ),
        se.fit = T
      )
    } else {
      
      print("correcting additionally for intercept in oystein's T2 model (intercept-only model)")
      predlm = predict(
        model,
        newdata = data.frame(
          convgroup = tmp$convgroup,
          rIntonly = mean(tmp$rIntonly),
          meanAge = mean(tmp$meanAge),
          cohort = "ADNINC",
          # cohort = "LCBC",
          subject_sex = "Female",
          nTimepoints = median(tmp$nTimepoints),
          intervalFirstlast = mean(tmp$intervalFirstlast)
        ),
        se.fit = T
      )
      
      predlm_cent = predict(
        model_cent,
        newdata = data.frame(
          convgroup = tmp$convgroup,
          rIntonly = mean(tmp$rIntonly),
          meanAge = mean(tmp$meanAge),
          meanCentiloids = mean(tmp$meanCentiloids),
          # cohort = "LCBC",
          cohort = "ADNINC",
          subject_sex = "Female",
          nTimepoints = median(tmp$nTimepoints),
          intervalFirstlast = mean(tmp$intervalFirstlast)
        ),
        se.fit = T
      )
    }
    
    
    #overwrite data with residualized effects for plotting
    tmp %<>% mutate(predY = predlm$fit,
                          predSE = predlm$se.fit,
                          predCI = predlm$se.fit*1.96,
                          resid = residuals(model),
                          presid = predY + resid,
                          predY_cent = predlm_cent$fit,
                          predSE_cent = predlm_cent$se.fit,
                          predCI_cent = predlm_cent$se.fit*1.96,
                          resid_cent = residuals(model_cent),
                          presid_cent = predY_cent + resid_cent,
                          X1 = tmp$convgroup)
    
    
    #note that the data is in a different order from when first model was ran
    # tidy(mod); tidy(modmaincheck)
    # plot(residuals(model), residuals(mod))
    # plot(residuals(model), residuals(modmaincheck))
    
    
    means_nonconv = mean(tmp$presid[tmp$convgroup == 0])
    means_conv = mean(tmp$presid[tmp$convgroup == 1])
    
    meanscent_nonconv = mean(tmp$presid_cent[tmp$convgroup == 0])
    meanscent_conv = mean(tmp$presid_cent[tmp$convgroup == 1])
    
    (pconv = ggplot(tmp,
           aes(x=(factor(convgroup)), y=(presid), col=factor(convgroup))) +
      geom_hline(yintercept = 0, linetype=2, col="light grey") +
      geom_jitter(alpha=0.7, position = position_jitter(seed = 123,width = 0.1)) +
      geom_segment(aes(x=0.75, y=means_nonconv, xend=1.25, yend=means_nonconv), size=.75, color="black") +
      geom_segment(aes(x=1.75, y=means_conv, xend=2.25, yend=means_conv), size=.75, color="black") +
      theme_classic() + mytheme +
      
      ylab(plotstring) +
      labs(x = NULL) +
      scale_color_manual(values = c(pal[2], pal[5])) +
      theme(legend.position = "top",
            axis.ticks= element_line()) #+
      # ggtitle("raw"))
    )
    
    (pconv_cent = ggplot(tmp,
                   aes(x=(factor(convgroup)), y=(presid_cent), col=factor(convgroup))) +
      geom_hline(yintercept = 0, linetype=2, col="light grey") +
      geom_jitter(alpha=0.7, position = position_jitter(seed = 123,width = 0.1)) +
      geom_segment(aes(x=0.75, y=meanscent_nonconv, xend=1.25, yend=meanscent_nonconv), size=.75, color="black") +
      geom_segment(aes(x=1.75, y=meanscent_conv, xend=2.25, yend=meanscent_conv), size=.75, color="black") +
      theme_classic() + mytheme +
      
      ylab(plotstring) +
      labs(x = NULL) +
      scale_color_manual(values = c(pal[2], pal[5])) +
      theme(legend.position = "top",
            axis.ticks= element_line()) #+
      # ggtitle("centiloid-corrected"))
    )
    
    res0 = tmp[tmp$convgroup == 0,]
    res1 = tmp[tmp$convgroup == 1,]
    res0$zero=0
    res1$zero=0
    
    
    dens <- ggplot() +
        PupillometryR::geom_flat_violin(data=res0,aes(x=zero,y=presid),alpha = 0.3, col = pal[2], fill = pal[2]) +
        PupillometryR::geom_flat_violin(data=res1,aes(x=zero,y=presid),alpha = 0.3, col = pal[5], fill = pal[5]) +
        theme_void()
    
    dens_cent <- ggplot() +
        PupillometryR::geom_flat_violin(data=res0,aes(x=zero,y=presid_cent),alpha = 0.3, col = pal[2], fill = pal[2]) +
        PupillometryR::geom_flat_violin(data=res1,aes(x=zero,y=presid_cent),alpha = 0.3, col = pal[5], fill = pal[5]) +
        theme_void()

    
    (pConcat = pconv + dens +
      patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
    )
    (pConcat_cent = pconv_cent + dens_cent +
        patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
    )
    
    return(list(
      pConcat = pConcat, 
      pConcat_cent = pConcat_cent,
      
      smodel = tidy(model),
      smodel_cent = tidy(model_cent),
      
      means_nonconv = means_nonconv,
      means_conv = means_conv,
      
      meanscent_nonconv = meanscent_nonconv,
      meanscent_conv = meanscent_conv
    ))
  }
  
  # intercept only model plot
  convplotmodYintonly[[ii]] = predictPlot(U.all.tmp, modYintonly, modYintonly_cent, 
              plotstring = "Additional Thickness (mm)")
  
  convplotslopeoys[[ii]] = predictPlot(U.all.tmp, modYslopeoys, modYslopeoys_cent, 
                                        plotstring = "Additional Change (mm/year)")
  
  convplotslointooys[[ii]] = predictPlot(U.all.tmp, modYslointooys, modYslointooys_cent, 
                                          plotstring = "Additional Change (mm/year)", intcorrect = T)
  
  # manual change model plot
  convplotmodYmanchange[[ii]] = predictPlot(U.all.tmp, modYmanchange, modYmanchange_cent, 
              plotstring = "Additional Change (mm/year, Manual change)")
  
  # manual change model plot (outliers removed)
  convplotmodYmanchangenoout[[ii]] = predictPlot(U.all.tmp.noout, modYmanchangenoout, modYmanchangenoout_cent, 
                                            plotstring = "Additional Change (mm/year, Manual change)")
  
  # main model plot check
  convplotmodmain[[ii]] = predictPlot(U.all.tmp, mod, mod_cent, 
                                            plotstring = "Additional Change (mm/year)")
  
  
  
  # check function worked
  convplotmodmain[[ii]]$pConcat
  pAll1[[1]]
  convplotmodYintonly[[ii]]$pConcat
  convplotmodYmanchange[[ii]]$pConcat
  convplotmodYmanchangenoout[[ii]]$pConcat
  convplotslopeoys[[ii]]$pConcat
  convplotslointooys[[ii]]$pConcat
  
  
  # sigflags (for saving figs only)
  sigflag_Yintonly = ifelse(convplotmodYintonly[[ii]]$smodel$p.value[2] < .05, 1, 0)
  sigflag_Yslopoys = ifelse(convplotslopeoys[[ii]]$smodel$p.value[2] < .05, 1, 0)
  sigflag_Yslointooys = ifelse(convplotslointooys[[ii]]$smodel$p.value[2] < .05, 1, 0)
  sigflag_Ymanchange = ifelse(convplotmodYmanchange[[ii]]$smodel$p.value[2] < .05, 1, 0)
  sigflag_Ymanchangenoout = ifelse(convplotmodYmanchangenoout[[ii]]$smodel$p.value[2] < .05, 1, 0)
  
  
  sigflag_Yintonlycent = ifelse(convplotmodYintonly[[ii]]$smodel_cent$p.value[2] < .05, 1, 0)
  sigflag_Yslopoyscent = ifelse(convplotslopeoys[[ii]]$smodel_cent$p.value[2] < .05, 1, 0)
  sigflag_Yslointooyscent = ifelse(convplotslointooys[[ii]]$smodel_cent$p.value[2] < .05, 1, 0)
  sigflag_Ymanchangecent = ifelse(convplotmodYmanchange[[ii]]$smodel_cent$p.value[2] < .05, 1, 0)
  sigflag_Ymanchangenooutcent = ifelse(convplotmodYmanchangenoout[[ii]]$smodel_cent$p.value[2] < .05, 1, 0)
  
  
  convplotmodmain[[ii]]$smodel_cent; tidy(mod_cent)
  sigflag_maincent = ifelse(convplotmodmain[[ii]]$smodel_cent$p.value[2] < .05, 1, 0)
  
  
  if (savefigs) {
    
    plotdir_Yintonly = file.path(plotdir, "Yintonly")
    plotdir_Ymanchange = file.path(plotdir, "Ymanchange")
    plotdir_Ymanchangenoout = file.path(plotdir, paste0("Ymanchangenoout", outlierThresh))
    
    plotdir_maincent = file.path(plotdir, "centcor")
    plotdir_Yslopoys = file.path(plotdir, "Yslopoys")
    plotdir_Yslointooys = file.path(plotdir, "Yslointooys")
    
    mkdirFunc = function(pdir) {
      if (!dir.exists(pdir)) { 
        print(paste("making", pdir))
        system(paste("mkdir -p", pdir))
      }
    }
    mkdirFunc(plotdir_Yintonly)
    mkdirFunc(plotdir_Ymanchange)
    mkdirFunc(plotdir_Ymanchangenoout)
    mkdirFunc(plotdir_maincent)
    mkdirFunc(plotdir_Yslopoys)
    mkdirFunc(plotdir_Yslointooys)
    
    filename_Yintonly = paste0(plotdir_Yintonly, "/pconvYintonly_cortrois-",ii,"-sig", sigflag_Yintonly,"-",roi,"_scale", scaleRes, ".png")
    filename_Yintonlycent = paste0(plotdir_Yintonly, "/pconvYintonlycent_cortrois-",ii,"-sig", sigflag_Yintonlycent,"-",roi,"_scale", scaleRes, ".png")
    
    filename_Ymanchange = paste0(plotdir_Ymanchange, "/pconvYmanchange_cortrois-",ii,"-sig", sigflag_Ymanchange,"-",roi,"_scale", scaleRes, ".png")
    filename_Ymanchangecent = paste0(plotdir_Ymanchange, "/pconvYmanchangecent_cortrois-",ii,"-sig", sigflag_Ymanchangecent,"-",roi,"_scale", scaleRes, ".png")
    
    filename_Ymanchangenoout = paste0(plotdir_Ymanchangenoout, "/pconvYmanchangenoout", outlierThresh, "_cortrois-",ii,"-sig", sigflag_Ymanchangenoout,"-",roi,"_scale", scaleRes, ".png")
    filename_Ymanchangenooutcent = paste0(plotdir_Ymanchangenoout, "/pconvYmanchangenoout", outlierThresh, "cent_cortrois-",ii,"-sig", sigflag_Ymanchangenooutcent,"-",roi,"_scale", scaleRes, ".png")
    
    filename_Yslopoys = paste0(plotdir_Yslopoys, "/pconvYslopoys", outlierThresh, "_cortrois-",ii,"-sig", sigflag_Yslopoys,"-",roi,"_scale", scaleRes, ".png")
    filename_Yslopoyscent = paste0(plotdir_Yslopoys, "/pconvYslopoys", outlierThresh, "cent_cortrois-",ii,"-sig", sigflag_Yslopoyscent,"-",roi,"_scale", scaleRes, ".png")
    
    filename_Yslointooys = paste0(plotdir_Yslopoys, "/pconvYslointooys", outlierThresh, "_cortrois-",ii,"-sig", sigflag_Yslointooys,"-",roi,"_scale", scaleRes, ".png")
    filename_Yslointooyscent = paste0(plotdir_Yslopoys, "/pconvYslointooys", outlierThresh, "cent_cortrois-",ii,"-sig", sigflag_Yslointooyscent,"-",roi,"_scale", scaleRes, ".png")
    
    
    filename_maincent = paste0(plotdir_maincent, "/pconvcent_cortrois-",ii,"-sig", sigflag_maincent,"-",roi,"_scale", scaleRes, ".png")
    print(paste("saving Yintonly and manchange plots"))
    print("-----------------------")
    ggsave(plot = convplotmodYintonly[[ii]]$pConcat,
           filename = filename_Yintonly,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotmodYintonly[[ii]]$pConcat_cent,
           filename = filename_Yintonlycent,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotmodYmanchange[[ii]]$pConcat,
           filename = filename_Ymanchange,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotmodYmanchange[[ii]]$pConcat_cent,
           filename = filename_Ymanchangecent,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotmodYmanchangenoout[[ii]]$pConcat,
           filename = filename_Ymanchangenoout,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotmodYmanchangenoout[[ii]]$pConcat_cent,
           filename = filename_Ymanchangenooutcent,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotslopeoys[[ii]]$pConcat,
           filename = filename_Yslopoys,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotslopeoys[[ii]]$pConcat_cent,
           filename = filename_Yslopoyscent,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotslointooys[[ii]]$pConcat,
           filename = filename_Yslointooys,
           width=9, height=9, units="cm", dpi=600)
    
    ggsave(plot = convplotslointooys[[ii]]$pConcat_cent,
           filename = filename_Yslointooyscent,
           width=9, height=9, units="cm", dpi=600)
  }
  
}


modlist_to_df = function(modlist, Z = F) {
  
  if (Z) {
    df_modlist = data.frame(rois, bind_rows(modlist) %>% filter(term == "scale(convgroup)"))
  } else {
    df_modlist = data.frame(rois, bind_rows(modlist) %>% filter(term == "convgroup"))
  }
  return(df_modlist)
}
  
Z = scaleRes
cortmap = modlist_to_df(modlist1, Z)
cortmap_cent = modlist_to_df(modlist1_cent, Z)

cortmapnoage = modlist_to_df(modlist2, Z)
cortmapnoage_cent = modlist_to_df(modlist2_cent, Z)

cortmapYintonly = modlist_to_df(modlist3, Z)
cortmapYintonly_cent = modlist_to_df(modlist3_cent, Z)

cortmapYslopeoys = modlist_to_df(modlist4, Z)
cortmapYslopeoys_cent = modlist_to_df(modlist4_cent, Z)

cortmapYslointooys = modlist_to_df(modlist5, Z)
cortmapYslointooys_cent = modlist_to_df(modlist5_cent, Z)

cortmapYmanchange = modlist_to_df(modlist6, Z)
cortmapYmanchange_cent = modlist_to_df(modlist6_cent, Z)

cortmapYmanchangenoout = modlist_to_df(modlist7, Z)
cortmapYmanchangenoout_cent = modlist_to_df(modlist7_cent, Z)


df_cor1 = data.frame(rois, r = cor1) #slopeInt
df_cor2 = data.frame(rois, r = cor2) #rSlope rIntonly
df_cor3 = data.frame(rois, r = cor3) #slopeoys rIntonly


name_map = function(map_str, Z = F) {
  if (Z) {
    map_name = paste0(cortmapdir, "/", metric, "Z_cortmapAB", yearsBeforeAB, "_nTime", nTime, "_", analysname, map_str, ".csv")
  } else {
    map_name = paste0(cortmapdir, "/", metric, "_cortmapAB", yearsBeforeAB, "_nTime", nTime, "_", analysname, map_str, ".csv")
  }
  return(map_name)
}

metric = "thickness"
cortmapdir = file.path(b, paste0("results/cortmaps/", analysname, "/", metric))

map1 = name_map("", Z)
map1_cent = name_map("_centcor", Z)

map2 = name_map("_noage", Z)
map2_cent = name_map("_noagecentcor", Z)

map3 = name_map("_Yintonly", Z)
map3_cent = name_map("_Yintonlycentcor", Z)

map4 = name_map("_Yslopoys", Z) #NB! cannot name another map with Yslope* otherwise later script will not work
map4_cent = name_map("_Yslopoyscentcor", Z)

map5 = name_map("_Yslointooys", Z)
map5_cent = name_map("_Yslointooyscentcor", Z)

map6 = name_map("_Ymanchange", Z)
map6_cent = name_map("_Ymanchangecentcor", Z)

map7 = name_map("_Ymanchangenoout", Z)
map7_cent = name_map("_Ymanchangenooutcentcor", Z)


# dict of map names
map_names = c(
  map1      = "cortmap",
  map1_cent = "cortmap_centcor",
  
  map2      = "cortmapnoage",
  map2_cent = "cortmapnoage_centcor",
  
  map3     = "cortmapYintonly",
  map3_cent= "cortmapYintonlycentcor",
  
  map4     = "cortmapYslopeoys",
  map4_cent= "cortmapYslopeoyscentcor",
  
  map5     = "cortmapYslointooys",
  map5_cent= "cortmapYslopintooyscentcor",
  
  map6      = "cortmapYmanchange",
  map6_cent = "cortmapYmanchangecentcor",
  
  map7      = "cortmapYmanchangenoout",
  map7_cent = "cortmapYmanchangenooutcentcor"
)


# dict of variable names
map_vnames = c(
  vmap1      = "cortmap",
  vmap1_cent = "cortmap_cent",
  
  vmap2      = "cortmapnoage",
  vmap2_cent = "cortmapnoage_cent",
  
  vmap3     = "cortmapYintonly",
  vmap3_cent= "cortmapYintonly_cent",
  
  vmap4     = "cortmapYslopeoys",
  vmap4_cent= "cortmapYslopeoys_cent",
  
  vmap5     = "cortmapYslointooys",
  vmap5_cent= "cortmapYslointooys_cent",
  
  vmap6      = "cortmapYmanchange",
  vmap6_cent = "cortmapYmanchange_cent",
  
  vmap7      = "cortmapYmanchangenoout",
  vmap7_cent = "cortmapYmanchangenoout_cent"
  
)


if (writeMaps) {

  if (!dir.exists(cortmapdir)) { 
    system(paste("mkdir -p", cortmapdir))
  }
  print("writing cortical maps")
  
  saveMap = function(imap, filename) {
    write.table(imap, filename, row.names = F, col.names = T, quote = F)
  }
  
  saveMap(get(map_names['map1']), map1)
  saveMap(get(map_vnames['vmap1_cent']), map1_cent)
  
  saveMap(get(map_names['map2']), map2)
  saveMap(get(map_vnames['vmap2_cent']), map2_cent)
  
  saveMap(get(map_names['map3']), map3)
  saveMap(get(map_vnames['vmap3_cent']), map3_cent)
  
  saveMap(get(map_names['map4']), map4)
  saveMap(get(map_vnames['vmap4_cent']), map4_cent)
  
  saveMap(get(map_names['map5']), map5)
  saveMap(get(map_vnames['vmap5_cent']), map5_cent)

  saveMap(get(map_names['map6']), map6)
  saveMap(get(map_vnames['vmap6_cent']), map6_cent)
  
  saveMap(get(map_names['map7']), map7)
  saveMap(get(map_vnames['vmap7_cent']), map7_cent)

}


# prep data for ggseg
ggWrangle = function(imap, statistical_test = T) {
  
  omap = imap
  omap$plotrois = tolower(omap$rois)
  omap$plotrois = gsub(".aparcnative71", "", omap$plotrois)
  omap$plotrois = gsub(paste0("_", metric), "", omap$plotrois)
  omap$hemi = substr(omap$plotrois, 1,2)
  
  omap = omap[omap$hemi == "lh" | omap$hemi == "rh",]
  omap$hemi = ifelse(omap$hemi == "lh", "left",
                     ifelse(omap$hemi == "rh", "right", NA))
  omap$plotrois[!omap$plotrois %in% dk$data$label]
  omap$plotrois[omap$plotrois %in% dk$data$label]
  omap$label = omap$plotrois
  
  # add fdr correction indicators
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


common_text = theme(panel.background = element_rect(fill = "white"),
                    plot.background = element_rect(fill = "white"),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    plot.subtitle = element_blank(),
                    plot.title = element_text(size = 10),
                    legend.title = element_text(size = 10),
                    legend.text = element_text(size = 10),
                    axis.line = element_blank(),
                    axis.ticks = element_blank())



# make ggseg figs
ggSegIt = function(imap, map_text, statistical_test = T) {
  
  newAtlas = dk %>% 
    as_tibble() %>% 
    left_join(as.data.frame(imap))
  
  
  if (statistical_test) {
    
    newbrain = ggseg(.data = newAtlas, mapping = aes(fill = factor(pIndFDR)), col = "dark grey", size= 0.5, atlas = "dk")
    newbrainT = ggseg(.data = newAtlas, mapping = aes(fill = statistic), col = "dark grey", size= 0.5, atlas = "dk")
    
    sighits = sum(imap$pIndFDR == 1)
    fdrhits = sum(imap$pIndFDR == 2)
    
    # cortical map pIndFDR -----
    (ggFacmap = newbrain + scale_fill_manual(values = c("white", "#83b8c9", "black")) + theme_classic() + mytheme +
      common_text +
      ggtitle(label = paste0(map_text, " ", metric, " ",
                             sighits, "sig, ", fdrhits, "rejections")))
  
    
    # cortical map T -----
    (ggStatmap = newbrainT + theme_classic() + mytheme +
       scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0, limits = c(-2, 2), oob = scales::squish, na.value = "white") +
       common_text +
       ggtitle(label = paste0(map_text, " ", metric, " ",
                              sighits, "sig, ", fdrhits, "rejections")) +
       theme(legend.position = "top")
    )
  
    return(list(ggFacmap = ggFacmap, ggStatmap = ggStatmap))
    
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


# map1
cortmap = ggWrangle(cortmap)
cortmap_cent = ggWrangle(cortmap_cent)

# map2
cortmapnoage = ggWrangle(cortmapnoage)
cortmapnoage_cent = ggWrangle(cortmapnoage_cent)

# map3
cortmapYintonly = ggWrangle(cortmapYintonly)
cortmapYintonly_cent = ggWrangle(cortmapYintonly_cent)

# map4
cortmapYslopeoys = ggWrangle(cortmapYslopeoys)
cortmapYslopeoys_cent = ggWrangle(cortmapYslopeoys_cent)

# map5
cortmapYslointooys = ggWrangle(cortmapYslointooys)
cortmapYslointooys_cent = ggWrangle(cortmapYslointooys_cent)

# map6
cortmapYmanchange = ggWrangle(cortmapYmanchange)
cortmapYmanchange_cent = ggWrangle(cortmapYmanchange_cent)

# map7
cortmapYmanchangenoout = ggWrangle(cortmapYmanchangenoout)
cortmapYmanchangenoout_cent = ggWrangle(cortmapYmanchangenoout_cent)


# random effect correlations
df_cor1 = ggWrangle(df_cor1, F)
df_cor2 = ggWrangle(df_cor2, F)
df_cor3 = ggWrangle(df_cor3, F)



plotMaps = 0
if (plotMaps) {
  
  maplist1 = ggSegIt(cortmap, "main")
  maplist1_cent = ggSegIt(cortmap_cent, "main centcor")
  maplist1$ggFacmap; maplist1$ggStatmap
  maplist1_cent$ggFacmap; maplist1_cent$ggStatmap
  
  maplist2 = ggSegIt(cortmapnoage, "noage")
  maplist2_cent = ggSegIt(cortmapnoage_cent, "noage centcor")
  
  maplist3 = ggSegIt(cortmapYintonly, "Yintonly")
  maplist3_cent = ggSegIt(cortmapYintonly_cent, "Yintonly centcor")
  maplist3$ggFacmap; maplist3$ggStatmap
  
  maplist4 = ggSegIt(cortmapYslopeoys, "Yslopeoys")
  maplist4_cent = ggSegIt(cortmapYslopeoys_cent, "Yslopeoys centcor")
  
  maplist5 = ggSegIt(cortmapYslointooys, "Yslointooys")
  maplist5_cent = ggSegIt(cortmapYslointooys_cent, "Yslointooys centcor")
  
  maplist6 = ggSegIt(cortmapYmanchange, "Ymanchange")
  maplist6_cent = ggSegIt(cortmapYmanchange_cent, "Ymanchange centcor")
  
  maplist7 = ggSegIt(cortmapYmanchangenoout, "Ymanchangenoout")
  maplist7_cent = ggSegIt(cortmapYmanchangenoout_cent, "Ymanchangenoout centcor")
  
  
  ggarrange(
    maplist1$ggFacmap,
    maplist1$ggStatmap,
    maplist3$ggFacmap,
    maplist3$ggStatmap,
    ncol = 1,
    common.legend = TRUE,
    legend = "right",
    align = "hv"
  )
  ggarrange(
    maplist1_cent$ggFacmap,
    maplist1_cent$ggStatmap,
    maplist3_cent$ggFacmap,
    maplist3_cent$ggStatmap,
    ncol = 1,
    common.legend = TRUE,
    legend = "right",
    align = "hv"
  )
  
  
  #Yslope
  maplist1$ggFacmap
  maplist1$ggStatmap
  
  #Yslopenoage
  maplist2$ggFacmap
  maplist2$ggStatmap
  
  #Yintonly
  maplist3$ggFacmap
  maplist3$ggStatmap
  
  #Yslopoys (T2 gamm slopes)
  maplist4$ggFacmap
  maplist4$ggStatmap
  
  #Yslointoys (T2 gamm slopes intercept-corrected)
  maplist5$ggFacmap
  maplist5$ggStatmap

  #Ymanchange
  maplist6$ggFacmap
  maplist6$ggStatmap
  
  #Ymanchange (no outliers)
  maplist7$ggFacmap
  maplist7$ggStatmap
  
  
  maplistcor1 = ggSegIt(df_cor1, "rSlope rInt", F)
  maplistcor5 = ggSegIt(df_cor2, "rSlope rIntonly", F)
  maplistcor8 = ggSegIt(df_cor3, "rSlopoys rIntonly", F)
  
  
  plotDots = function(imap, x, y, fill_by = y) {
    
    if (x == "along_x") {
      x = "1:nrow(imap)"
    }
    
    (ggcor = ggplot(imap, aes_string(x = x, y = y)) +
       geom_point(aes_string(fill = fill_by, col = fill_by)) +
       scale_y_continuous(limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         labels = c("-1", "-0.5", "0", "0.5", "1")) +
      scale_colour_viridis(option = "D", na.value = "white", limits = c(-1,1)) +
      geom_hline(yintercept = 0, linetype = 2) +
      theme_classic() + mytheme +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(r = 10)),
            axis.title.y.left = element_text(margin = margin(r = 10))
      )
    )
    
    return(ggcor)
  }
  

  # correlations between random effect estimates
  cowplot::plot_grid(maplistcor1$ggStatmap + theme(legend.position = "right"), plotDots(df_cor1, "along_x", "r"))
  cowplot::plot_grid(maplistcor5$ggStatmap + theme(legend.position = "right"), plotDots(df_cor2, "along_x", "r"))
  cowplot::plot_grid(maplistcor8$ggStatmap + theme(legend.position = "right"), plotDots(df_cor3, "along_x", "r"))
  
  
  # SI Fig 12
  # ggsave(plot = maplistcor5$ggStatmap,
  #        filename = paste0(b, "/paper2/figs_yearsBeforeAB/p_Seg5_", analysname, ".png"),
  #        width=20, height=7, units="cm", dpi=600)
  #   
  # ptmp = plotDots(df_cor5, "along_x", "r")
  #   
  # ggsave(plot = ptmp,
  #        filename = paste0(b, "/paper2/figs_yearsBeforeAB/p_SegDots5_", analysname, ".png"),
  #        width=7, height=3.5, units="cm", dpi=600)
  # 
  # ggsave(plot = maplistcor8$ggStatmap,
  #        filename = paste0(b, "/paper2/figs_yearsBeforeAB/p_Seg8_", analysname, ".png"),
  #        width=20, height=7, units="cm", dpi=600)
  #   
  # ptmp = plotDots(df_cor8, "along_x", "r")
  #   
  # ggsave(plot = ptmp,
  #        filename = paste0(b, "/paper2/figs_yearsBeforeAB/p_SegDots8_", analysname, ".png"),
  #        width=7, height=3.5, units="cm", dpi=600)
    
    
}


# quick check of signifcant differences
sigidx = which(cortmap$p.value<.05)
fdrcor = sgof::BH(cortmap$p.value)
fdrthresh =  max(fdrcor$data[fdrcor$Adjusted.pvalues<.05])
fdridx = which(cortmap$p.value<=fdrthresh)
rois[sigidx]
rois[fdridx]


# resample groups -----
fullconv = U.all.tmp
table(fullconv$convgroup)


writeResamples = function(rep_num = 1000, rep_proportion = .90) {
  
  # rep_num = 1000; rep_proportion = .90
  rep_sample_size_group1 = floor(nrow(fullconv[fullconv$convgroup == 1,])*rep_proportion)
  rep_sample_size_group0 = floor(nrow(fullconv[fullconv$convgroup == 0,])*rep_proportion)
  
  table(fullconv$convgroup)
  
  set.seed(123)
  print("sampling datasets convgroup")
  sampled_datasets_group1 <- replicate(rep_num, {
    fullconv %>%
      filter(convgroup == 1) %>% 
      sample_n(rep_sample_size_group1)
  }, simplify = FALSE)
  
  print("sampling datasets AB- group")
  sampled_datasets_group0 <- replicate(rep_num, {
    fullconv %>%
      filter(convgroup == 0) %>% 
      sample_n(rep_sample_size_group0)
  }, simplify = FALSE)
  
  print("concatenating resampled datasets")
  sampled_datasets = vector("list", rep_num)
  for (i in 1:rep_num) {
    sampled_datasets[[i]] = rbind(sampled_datasets_group1[[i]],
                                  sampled_datasets_group0[[i]])
  }
  
  
  brainslopes = data.frame(rSlopeMat.all)
  names(brainslopes) = ROIs.all
  brainslopes$subject_id = U.all$subject_id
  brainslopes$yearsBeforeAB = yearsBeforeAB
  
  
  brainsintsonly = data.frame(rIntonlyMat.all)
  names(brainsintsonly) = ROIs.all
  brainsintsonly$subject_id = U.all$subject_id
  brainsintsonly$yearsBeforeAB = yearsBeforeAB
  
  
  odir = file.path(b, paste0("results/cortmaps_resample/", analysname, "/", metric, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime))
  if (!dir.exists(odir)) { 
    print(paste("writing resampled data to", odir))
    system(paste("mkdir -p", odir))
  }
  
  
  ofile1 = file.path(odir, paste0("U_brainslopes_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, ".csv"))
  ofile2 = file.path(odir, paste0("U_brainintsonly_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, ".csv"))
  if (! file.exists(ofile1)) { fwrite(brainslopes, file = ofile1) }
  if (! file.exists(ofile2)) { fwrite(brainsintsonly, file = ofile2) }
  
  
  for (ii in 1:length(sampled_datasets)) {
    print(paste("writing", ii))
    fwrite(sampled_datasets[[ii]], file = file.path(odir, paste0("resample_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-",sprintf("%04d", ii), ".csv")))
  }
}


if (writeMaps) {
  writeResamples()
} else if (writeMaps == 0) {
  quit()
}



# setup null model (permutation-based analyses - group differences in variance) ----------
set.seed(123)
samplerois = rois
samplerois = samplerois[1:70]
sample = fullconv


odir = file.path(b, paste0("results/cortmaps_resample/", analysname, "/", metric, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime))
ofile1 = file.path(odir, paste0("U_brainslopes_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, ".csv"))
slopes.all.tmp = fread(ofile1)


for (jj in 1:length(samplerois)) {

  print(jj)
  sampleroi = samplerois[jj]
  sample$brainvarSlope = NULL
  
  if (sampleroi != "pcaL" & sampleroi != "pcaR") {
    print("attaching ROI as brainvarSlope")
    sample = left_join(sample,
                       slopes.all.tmp %>% dplyr::select(subject_id,
                                                   all_of(sampleroi)) %>% rename(brainvarSlope = sampleroi)
    )
  }
  if (sampleroi == "pcaL") {
    sample$brainvarSlope = sample$pca_Left
  }
  if (sampleroi == "pcaR") {
    sample$brainvarSlope = sample$pca_Right
  }
  
  
  # null models (no group variable)
  nullmod = lm(scale(brainvarSlope) ~ subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = sample)
  nullmod_centcor = lm(scale(brainvarSlope) ~ meanCentiloids + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = sample)
  
  
  # predict
  nullpred = predict(nullmod, data = sample)
  nullresid = residuals(nullmod, data = fullconv)
  nullpred_centcor = predict(nullmod_centcor, data = sample)
  nullresid_centcor = residuals(nullmod_centcor, data = fullconv)
  
  
  # add to data
  sample$nullpred = nullpred
  sample$nullresid = nullresid
  sample$nullpred_centcor = nullpred_centcor
  sample$nullresid_centcor = nullresid_centcor
  
  
  # save in roi-specific dirs
  odir = paste0(b, "/results/cortmaps_resample/", analysname, "/", metric, "_nullmod/yearsBeforeAB", yearsBeforeAB, "nTime", nTime,
                "/", sampleroi)
  if (!dir.exists(odir)) { 
    system(paste("mkdir -p", odir))
  }
  
  # subset cols
  sample_dat = sample %>% select(-contains(c("w-g", "area", "volume", "vol.", "int.")))
  
  
  writeRegionsSample = 1
  if (writeMaps == 1) {
    if (writeRegionsSample) {
      fwrite(sample_dat, file = file.path(odir, paste0("region_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-",sampleroi, ".csv")))
    }
  }
}

quit()
