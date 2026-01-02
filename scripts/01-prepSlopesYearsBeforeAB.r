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


#### testing ------------
##input tests
# rm(list=ls())
#
# options(bitmapType='cairo')


# analysnum=1
# cohortSelect = "adnibacslcbc" #"adnibacs" "bacslcbc" "adnilcbc"
# yearsBeforeAB = 1
# matchFollow = 0
# procType = "LONG"
# predAB = 0
# singleT = 0
# analysname = paste0(
#   "analysnum",
#   analysnum,
#   "matchF",
#   matchFollow,
#   "proc",
#   procType,
#   "predAB",
#   predAB,
#   "singleT",
#   singleT,
#   "reproduce1"
# )


#---load packages
loadPackages = function() {
  packages = c("here", "tidyverse","magrittr","gamm4","itsadug","numDeriv","gratia","mgcv","viridis","wesanderson","asbio","broom","cowplot","data.table","stringi","tictoc","MatchIt","ggpubr")
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
savefigs=T
nTime=2
agecut=30
saveres=T



#load and combine data
load(file.path(b,"reproduce/data/DF_BACSUPDATE_PREPPED_AB.Rda"))
load(file.path(b,"DF_ADNIUPDATE_PREPPED_AB.Rda"))
ADNI = ADNIEXPORT
BACS = BACSEXPORT
allFeat = readLines(file.path(b, "reproduce/data/allFeatures364.txt"))
adnioutlier = "029_S_0845"



# SWITCH procType (cross / long) -----
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


# recalc difference between age at mri and first AB+
converters_allMRI_negfirst_ADNINC_plotdat$diff_mriAge_ABpos = converters_allMRI_negfirst_ADNINC_plotdat$ageAtFirstABpos - converters_allMRI_negfirst_ADNINC_plotdat$Age
converters_allMRI_negfirst_LCBC_plotdat$diff_mriAge_ABpos = converters_allMRI_negfirst_LCBC_plotdat$ageAtFirstABpos - converters_allMRI_negfirst_LCBC_plotdat$visit_age


# confirm same as loaded data
cor(converters_allMRI_negfirst_ADNINC$diff_mriAge_ABpos, converters_allMRI_negfirst_ADNINC_plotdat$diff_mriAge_ABpos)
cor(converters_allMRI_negfirst_BACS$diff_mriAge_ABpos, converters_allMRI_negfirst_BACS_plotdat$diff_mriAge_ABpos)
cor(converters_allMRI_negfirst_LCBC$diff_mriAge_ABpos, converters_allMRI_negfirst_LCBC_plotdat$diff_mriAge_ABpos)



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



changeImageLink = function(dat) {
  dat$imageLink = lapply(dat$imageLink, function(x) strsplit(x, "\\.")[[1]][1]) %>% unlist()
  return(dat$imageLink)
}

if (procType == "CROSS") {
  # change imagelinks
  converters_allMRI_negfirst_ADNINC$imageLink = changeImageLink(converters_allMRI_negfirst_ADNINC)
  converters_allMRI_negfirst_BACS$imageLink = changeImageLink(converters_allMRI_negfirst_BACS)
  converters_allMRI_negfirst_LCBC$imageLink = changeImageLink(converters_allMRI_negfirst_LCBC)
}
converters_allMRI_negfirst_ADNINC$ageAtFirstABpos
converters_allMRI_negfirst_BACS$ageAtFirstABpos
converters_allMRI_negfirst_LCBC$ageAtFirstABpos



# Get usable images from MRI dataset
getUsable = function(dat, yearsBeforeAB, predAB) {
  if (!"diff_mriAge_ABpos" %in% names(dat)) stop("Column diff_mriAge_ABpos is missing")
  if (!"imageLink" %in% names(dat)) stop("Column imageLink is missing")
  if (!"diff_mriAge_predABpos" %in% names(dat)) stop("Column diff_mriAge_predABpos is missing")
  
  if (predAB == F) {
    usable = dat %>% 
      filter(diff_mriAge_ABpos >= yearsBeforeAB) %>% 
      select(imageLink)
  } else if (predAB == T) {
    usable = dat %>% 
      filter(diff_mriAge_predABpos >= yearsBeforeAB) %>% 
      select(imageLink)
  }
  return(usable$imageLink)
}


# Get unusable image links from MRI dataset
getUnusable = function(dat, yearsBeforeAB, predAB) {
  if (!"diff_mriAge_ABpos" %in% names(dat)) stop("Column diff_mriAge_ABpos is missing")
  if (!"imageLink" %in% names(dat)) stop("Column imageLink is missing")
  if (!"diff_mriAge_predABpos" %in% names(dat)) stop("Column diff_mriAge_predABpos is missing")
  
  if (predAB == F) {
    unusable = dat %>% 
      filter(diff_mriAge_ABpos < yearsBeforeAB) %>% 
      select(imageLink)
  } else if (predAB == T) {
    unusable = dat %>% 
      filter(diff_mriAge_predABpos < yearsBeforeAB) %>% 
      select(imageLink)
  }
  return(unusable$imageLink)
}


#calculate difference from predicted AB pos
converters_allMRI_negfirst_ADNINC$diff_mriAge_predABpos = converters_allMRI_negfirst_ADNINC$age_at_threshold - converters_allMRI_negfirst_ADNINC$Age
converters_allMRI_negfirst_BACS$diff_mriAge_predABpos = converters_allMRI_negfirst_BACS$age_at_threshold - converters_allMRI_negfirst_BACS$ageMRI
converters_allMRI_negfirst_LCBC$diff_mriAge_predABpos = converters_allMRI_negfirst_LCBC$age_at_threshold - converters_allMRI_negfirst_LCBC$visit_age



# SWITCH get usable/unusable scans --------------------------------------------------------------------
# if predAB == 1 will return list of usable/unusable MRIs based on yearsBefore predABpos
# if predAB == 0 will return list of usable/unusable MRIs based on yearsBefore observedABpos
unusable =
  c(getUnusable(converters_allMRI_negfirst_ADNINC, yearsBeforeAB, predAB),
    getUnusable(converters_allMRI_negfirst_BACS, yearsBeforeAB, predAB),
    getUnusable(converters_allMRI_negfirst_LCBC, yearsBeforeAB, predAB))

usable =
  c(getUsable(converters_allMRI_negfirst_ADNINC, yearsBeforeAB, predAB),
    getUsable(converters_allMRI_negfirst_BACS, yearsBeforeAB, predAB),
    getUsable(converters_allMRI_negfirst_LCBC, yearsBeforeAB, predAB))


# check usable/unusable
checkUsable = function(dat) {
  checkDat = dat %>% filter(imageLink %in% usable) %>% select(imageLink, contains("age"), nTimepoints)
}
checkUnusable = function(dat) {
  checkDat = dat %>% filter(imageLink %in% unusable) %>% select(imageLink, contains("age"), nTimepoints)
}

check1 = checkUsable(converters_allMRI_negfirst_ADNINC)
check2 = checkUsable(converters_allMRI_negfirst_BACS)
check3 = checkUsable(converters_allMRI_negfirst_LCBC) %>% filter(imageLink %in% usable) %>% select(imageLink, contains("age"), nTimepoints)
ucheck1 = checkUnusable(converters_allMRI_negfirst_ADNINC)
ucheck2 = checkUnusable(converters_allMRI_negfirst_BACS)
ucheck3 = checkUnusable(converters_allMRI_negfirst_LCBC)

plotCheck = function(checkDat, ageVar, titleString) {
  ggplot(checkDat, aes(x = {{ ageVar }}, y = {{ ageVar }})) + 
    geom_point(col ="red") +
    geom_point(aes(y = {{ ageVar }}, x = ageAtFirstABpos), col = "blue") +
    labs(y = "Age", x = "ageAtFirstABpos", title = titleString) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1)
}
plotCheck(check1, Age, "usable ADNI")
plotCheck(check2, ageMRI, "usable BACS")
plotCheck(check3, visit_age, "usable LCBC")
plotCheck(ucheck1, Age, "unusable ADNI")
plotCheck(ucheck2, ageMRI, "unusable BACS")
plotCheck(ucheck3, visit_age, "unusable LCBC")


range(
  c(
    converters_allMRI_negfirst_BACS$diff_mriAge_predABpos,
    converters_allMRI_negfirst_LCBC$diff_mriAge_predABpos,
    converters_allMRI_negfirst_ADNINC$diff_mriAge_predABpos
  ),
  na.rm = T
)


# combine MRI data for plotting
# converters
convallMRI = rbind(
  converters_allMRI_negfirst_ADNINC_plotdat %>% mutate(imageLinkAge = paste0(imageLink, "-", round(Age, 2))) %>% dplyr::select(imageLink, imageLinkAge, PTID, Age, "ageAtFirstABpos",  "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(scanType = "MRI") %>% mutate(cohort = "ADNINC") %>% rename(SID = PTID),
  converters_allMRI_negfirst_BACS_plotdat %>% mutate(imageLinkAge = paste0(imageLink, "-", round(ageMRI, 2))) %>% dplyr::select(imageLink, imageLinkAge, SID, ageMRI, "ageAtFirstABpos",  "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(scanType = "MRI") %>% mutate(cohort = "BACS") %>% rename(Age = ageMRI),
  converters_allMRI_negfirst_LCBC_plotdat %>% mutate(imageLinkAge = paste0(imageLink, "-", round(visit_age, 2))) %>% dplyr::select(imageLink, imageLinkAge, subject_id, visit_age, "ageAtFirstABpos",  "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(scanType = "MRI") %>% rename(SID = subject_id) %>% mutate(cohort = "LCBC") %>% rename(Age = visit_age)
)
convallMRI %<>% mutate(usable = ifelse(imageLink %in% unusable, 0, 1))


# AB- group
nonconvallMRI = rbind(
  nonconverters_allMRI_ADNINC_plotdat %>% mutate(imageLinkAge = paste0(imageLink, "-", round(Age, 2))) %>% dplyr::select(imageLink, imageLinkAge, PTID, Age,  "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(scanType = "MRI") %>% mutate(cohort = "ADNINC") %>% rename(SID = PTID),
  nonconverters_allMRI_BACS_plotdat %>% mutate(imageLinkAge = paste0(imageLink, "-", round(ageMRI, 2))) %>% dplyr::select(imageLink, imageLinkAge, SID, ageMRI, "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(scanType = "MRI") %>% mutate(cohort = "BACS") %>% rename(Age = ageMRI),
  nonconverters_allMRI_LCBC_plotdat %>% mutate(imageLinkAge = paste0(imageLink, "-", round(visit_age, 2))) %>% dplyr::select(imageLink, imageLinkAge, subject_id, visit_age, "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(scanType = "MRI") %>% rename(SID = subject_id) %>% mutate(cohort = "LCBC") %>% rename(Age = visit_age)
)
nonconvallMRI %<>% filter(imageLink %in% DF$imageLink)


# combine AB data for plotting
# converters
AB_all_conv = rbind(
  AB_conv_ADNINC_plotdat %>% mutate(imageLinkAge = paste0(PTID, "-", round(Age, 2))) %>% dplyr::select(PTID, imageLinkAge, Age, "ageAtFirstABpos",  "ageAtFirstABscan", "ageAtLastABscan", "scanType") %>% rename(SID = PTID) %>% mutate(cohort = "ADNINC"),
  AB_conv_BACS_plotdat %>% mutate(imageLinkAge = paste0(SID, "-", round(Age, 2))) %>% dplyr::select(SID, imageLinkAge, Age, "ageAtFirstABpos",  "ageAtFirstABscan", "ageAtLastABscan", "scanType") %>% mutate(cohort = "BACS"),
  AB_conv_LCBC_plotdat %>% mutate(imageLinkAge = paste0(subject_id, "-", round(Age, 2))) %>% dplyr::select(subject_id, imageLinkAge, Age, "ageAtFirstABpos",  "ageAtFirstABscan", "ageAtLastABscan", "scanType") %>% mutate(cohort = "LCBC") %>% rename(SID = subject_id)
) %>% mutate(usable = 1)
setdiff(names(convallMRI), names(AB_all_conv))
convallMRI %<>% select(-imageLink)
convall = rbind(convallMRI, AB_all_conv)


# AB- group
AB_all_nonconv = rbind(
  AB_nonconv_ADNINC_plotdat %>% mutate(imageLinkAge = paste0(PTID, "-", round(Age, 2))) %>% dplyr::select(PTID, imageLinkAge, Age, "scanType", "ageAtFirstABscan", "ageAtLastABscan") %>% rename(SID = PTID) %>% mutate(cohort = "ADNINC"),
  AB_nonconv_BACS_plotdat %>% mutate(imageLinkAge = paste0(SID, "-", round(Age, 2))) %>% dplyr::select(SID, imageLinkAge, Age, "scanType", "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(cohort = "BACS"),
  AB_nonconv_LCBC_plotdat %>% mutate(imageLinkAge = paste0(subject_id, "-", round(Age, 2))) %>% dplyr::select(subject_id, imageLinkAge, Age, "scanType", "ageAtFirstABscan", "ageAtLastABscan") %>% mutate(cohort = "LCBC") %>% rename(SID = subject_id)
)

# ensure they are in the DF
AB_all_nonconv %<>% filter(SID %in% DF$subject_id) 
# table(AB_all_conv$SID %in% DF$subject_id)
# table(AB_all_nonconv$SID %in% DF$subject_id)


setdiff(names(nonconvallMRI), names(AB_all_nonconv))
nonconvallMRI %<>% select(-imageLink)
nonconvall = rbind(nonconvallMRI, AB_all_nonconv)


# order data for plotting
convall$cohort = factor(convall$cohort, levels = (c("LCBC", "BACS", "ADNINC")))
convall$scanInd = ifelse(convall$scanType == "AB", "AB", "MRI")
convall %<>% arrange(by = cohort, ageAtFirstABpos, SID, scanInd)
convall$SID_scan = paste0(convall$scanInd,"_",convall$SID)
convall$SID_scan = factor(convall$SID_scan, levels = rev(unique(convall$SID_scan)))


# converter group initially consisted of 77 individuals with 283 PET scans and a maximum of 487 MRI scans 
dim(convall)
convall %>% filter(scanType == "AB") %>% dim()
convall %>% filter(scanType == "MRI") %>% dim()
length(unique(convall$SID))

# Aβ− group initially consisted of 412 individuals with 954 Aβ PET scans and in total 1850 MRI scans
dim(nonconvall)
nonconvall %>% filter(scanType == "AB") %>% dim()
nonconvall %>% filter(scanType == "MRI") %>% dim()
length(unique(nonconvall$SID))



pal = wesanderson::wes_palettes$FantasticFox1
mytheme = theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  title = element_text(size=17),
  text = element_text(color = "black", size = 18, family="Nimbus Sans Narrow"),
  plot.title = element_text(hjust = 0.5),
  axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
  axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,20,0)),
  axis.text = element_text(color = "black", size = 18),
  legend.key.size = unit(1,"cm"))



# note that below needs to be commented out or script fails on cluster
# (pcovnall = convall %>% 
#     # filter(usable == 1) %>%
#     ggplot(.) +
#       geom_line(aes(x=Age,y=SID_scan,group=SID_scan), col="grey",alpha=0.7, size=1) +
#       geom_point(aes(x=Age,y=SID_scan,group=SID_scan, col=scanType),alpha=1, size=1) +
#       scale_color_manual(values = c(pal[2], pal[3])) +
#       geom_point(data = convall %>% filter(scanInd == "AB"), aes(x=ageAtFirstABpos,y=SID_scan,group=SID_scan),col=pal[5], alpha=1, shape=18, size=2) +
#       theme_classic() + mytheme + theme(
#         axis.text.y = element_text(size = 4.5),
#         axis.line.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.x = element_line(),
#         panel.grid.major.x = element_line(color = "grey85", linewidth = 0.25, linetype = 2)
#       ))
# pcovnall = pcovnall + theme(axis.text.y = element_blank(),
#                  axis.line.y = element_line())


# ggsave(filename = file.path(b, "paper2/figs_yearsBeforeAB/pabtrajALL_converters_middlesize_new770_AB1.pdf"),
#        plot = pcovnall + theme(axis.line.y = element_line(),
#                                        axis.text.y = element_blank(),
#                                        axis.ticks.y = element_blank()) + labs(x = "Age"),
#        width = 15,
#        height = 34,
#        dpi = 600,
#        units = "cm",
#        device = cairo_pdf
# )


nonconvall$cohort = factor(nonconvall$cohort, levels = (c("LCBC", "BACS", "ADNINC")))
nonconvall$scanInd = ifelse(nonconvall$scanType == "AB", "AB", "MRI")
nonconvall %<>% arrange(by = cohort, ageAtFirstABscan, SID, scanInd)
nonconvall$SID_scan = paste0(nonconvall$scanInd,"_",nonconvall$SID)
nonconvall$SID_scan = factor(nonconvall$SID_scan, levels = rev(unique(nonconvall$SID_scan)))


boundary_indices_nonconvall <- nonconvall %>%
  group_by(cohort) %>%
  summarize(boundary = max(as.numeric(SID_scan))) %>%
  pull(boundary)


# Calculate midpoint positions for geom_hline
hline_positions_nonconvall <- boundary_indices_nonconvall[-length(boundary_indices_nonconvall)] + 0.5
print(hline_positions_nonconvall)
hline_positions_nonconvall[1] = 426.5

# (pnonconvall = nonconvall %>%
#   ggplot(.) +
#   geom_line(aes(x=Age,y=SID_scan,group=SID_scan), col="grey",alpha=0.7, size=0.5) +
#   geom_point(aes(x=Age,y=SID_scan,group=SID_scan, col=scanType),alpha=1, size=0.5) +
#   geom_hline(yintercept = hline_positions_nonconvall, linetype = 1, color = "#E5E5E5") +
#   scale_color_manual(values = c(pal[2], pal[3])) +
#   theme_classic() + mytheme + theme(
#     axis.text.y = element_text(size = 4.5),
#     axis.line.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.x = element_line(),
#     panel.grid.major.x = element_line(color = "grey85", linewidth = 0.25, linetype = 2)
#   )
# )

# ggsave(filename = file.path("paper2/figs_yearsBeforeAB/pabtrajALL_nonconverters_middlesize_new2804_smallpoint.pdf"),
#        plot = pnonconvall + theme(axis.line.y = element_line(),
#                                axis.text.y = element_blank(),
#                                axis.ticks.y = element_blank()) + labs(x = "Age"),
#        width = 15,
#        height = 34,
#        dpi = 600,
#        units = "cm",
#        device = cairo_pdf
# )


dim(convall)
dim(convallMRI)
dim(AB_all_conv)

dim(nonconvall)
dim(nonconvallMRI)
dim(AB_all_nonconv)



# SWITCH singleT --------------
if (singleT) {
  print("-------------using only 1.5T scans------------------------")
  DF %<>% filter(scanStrength == "1-5T")
  dim(DF)
  # quit()
}


# overview of used scans in converter group -----
DF.conv = DF %>% filter(imageLink %in% usable) %>%
  # recalculate longitudinal
  group_by(subject_id) %>% mutate(TimeBsl = visit_age - min(visit_age),
                                  meanAge = mean(visit_age),
                                  nTimepoints= length(unique(visit_age)),
                                  intervalFirstlast = max(visit_age) - min(visit_age)) %>% filter(nTimepoints>=nTime)

dim(DF); length(unique(DF$subject_id))
dim(DF.conv); length(unique(DF.conv$subject_id))



# # load updated LCBC amyloid data
lcbcABnew=1
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
  U_AB_LCBC = df_AB_LCBC %>% filter(nTimeAB == 1)

  U_AB_LCBC %<>% mutate(convertgroup = 
                         case_when(
                           converter == "nonconverter" & ABlong == 1 ~ "nonconverter_pos",
                           converter == "nonconverter" & ABlong == 0 ~ "nonconverter_neg",
                           converter == "converter" & ABlong == 1 ~ "converter",
                           converter == "converter" & ABlong == 0 ~ "check",
                         ))
}


# total number of AB scans on MRI sample
all_ab_subs = c(
  U_AB_ADNI$PTID,
  U_AB_BACS$SID,
  U_AB_LCBC$subject_id) 

# number of AB scans total on DF subs
tmp1 = df_AB_LCBC %>% filter(subject_id %in% DF$subject_id)
dim(tmp1); length(unique(tmp1$subject_id))

tmp2 = df_ABbrain_BACS %>% filter(SID %in% DF$subject_id)
dim(tmp2); length(unique(tmp2$SID))

tmp3 = df_ABintersectNC_ADNI %>% filter(PTID %in% DF$subject_id)
dim(tmp3); length(unique(tmp3$PTID))

total_AB = nrow(tmp1) + nrow(tmp2) + nrow(tmp3)
total_AB_N = length(unique(tmp1$subject_id)) + length(unique(tmp2$SID)) + length(unique(tmp3$PTID))
total_AB; total_AB_N



range(DF$intervalFirstlast)


# CRITICAL remove unusable then redefine longitudinal -------------
DF %<>% filter(!imageLink %in% unusable) %>%
  
  # recalculate longitudinal
  group_by(subject_id) %>% mutate(TimeBsl = visit_age - min(visit_age),
                                  meanAge = mean(visit_age),
                                  nTimepoints= length(unique(visit_age)),
                                  intervalFirstlast = max(visit_age) - min(visit_age)) %>% filter(nTimepoints>=nTime)
    
  
dim(DF); length(unique(DF$subject_id))



converterlist = unique(c(converters_allAB_negfirst_ADNINC$PTID, converters_allAB_negfirst_BACS$SID, converters_allAB_negfirst_LCBC$subject_id))
nonconverterlist = unique(c(nonconverters_allAB_ADNINC$PTID, nonconverters_allAB_BACS$SID, nonconverters_allAB_LCBC$subject_id))
length(unique(nonconverterlist))
length(unique(converterlist))
sum(converterlist %in% DF$subject_id)
sum(nonconverterlist %in% DF$subject_id)



# SWITCH matchFollow --------------
# match on follow up data ----
if (matchFollow == 1) {
  
  print("matching converter / nonconverter groups on follow up time")  
  
  
  tmpU = DF %>% filter(TimeBsl == 0)
  U_convnonconv = rbind(
    tmpU %>% filter(subject_id %in% converterlist) %>% mutate(converter = 1),
    tmpU %>% filter(subject_id %in% nonconverterlist) %>% mutate(converter = 0)
    )
  
  
  # ensure adnioutlier is not drawn from the matched sample (found on first analysis, removed for all subsequent)
  U_convnonconv = U_convnonconv %>% filter(subject_id != adnioutlier)
  
  
  matchSample = function(dat) {

    # Create a matching object - NB! MASS package load masks dplyr select
    m.out <- matchit(as.factor(converter) ~ meanAge + nTimepoints + intervalFirstlast, data = dat,
                       method = "nearest"
      )
      
    matched_samples <- match.data(m.out)
    return(matched_samples)
  }
  n_unmatched = table(U_convnonconv$converter)
  
  
  if (yearsBeforeAB <= 7) {
    upper_lim = 0.9
  } else {
    upper_lim = 0.9
  }
  
  
  # plot unmatched variable densities
  p_unmatch1 = ggplot(U_convnonconv) + geom_density(aes(x=nTimepoints, fill=factor(converter)), alpha=0.2) + theme_classic() + mytheme +
    scale_fill_manual(values = rev(c(pal[5], pal[2]))) +
    scale_color_manual(values = rev(c(pal[5], pal[2]))) +
    scale_x_continuous(limits = c(2, 12), breaks = seq(2, 12, 2)) +
    scale_y_continuous(limits = c(0, upper_lim))
  
  
  p_unmatch2 = ggplot(U_convnonconv) + geom_density(aes(x=intervalFirstlast, fill=factor(converter)), alpha=0.2) + theme_classic() + mytheme +
    scale_fill_manual(values = rev(c(pal[5], pal[2]))) +
    scale_color_manual(values = rev(c(pal[5], pal[2]))) +
    labs(x = "Follow-up interval") +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, upper_lim))
  
  
  p_unmatch3 = ggplot(U_convnonconv) + geom_density(aes(x=meanAge, fill=factor(converter)), alpha=0.2) + ggtitle("Unmatched") + theme_classic() + mytheme +
    scale_fill_manual(values = rev(c(pal[5], pal[2]))) +
    scale_color_manual(values = rev(c(pal[5], pal[2]))) +
    labs(x = "Age (mean)") +
    scale_x_continuous(limits = c(40, 90)) +
    scale_y_continuous(limits = c(0, 0.075)) +
    labs(subtitle = paste0("conv = ", unname(n_unmatched[2]), "; Aβ- = ", unname(n_unmatched[1]))) +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 14, colour = "grey20"))
  
  
  
  # draw first matched sample
  matched_samples1 = matchSample(U_convnonconv)
  table(matched_samples1$converter)
  
  
  # non converter group is much larger
  table(U_convnonconv$converter)[1] / table(U_convnonconv$converter)[2]
  
  
  # draw second matched sample ensuring new draw from AB- group
  matched_samples2 = matchSample(
    U_convnonconv %>% filter(!subject_id %in% matched_samples1$subject_id[matched_samples1$converter == 0])
    )
  table(matched_samples2$converter)
  matched_samples2$subject_id[matched_samples2$converter == 0] %in% matched_samples1$subject_id
  
  
  tmpmatched = rbind(matched_samples1[matched_samples1$converter == 0,],
                     matched_samples2[matched_samples2$converter == 0,],
                     matched_samples1[matched_samples1$converter == 1,])
  
  n_matched = table(tmpmatched$converter)
  p_match1 = ggplot(tmpmatched) + geom_density(aes(x=nTimepoints, fill=factor(converter)), alpha=0.2)+ theme_classic() + mytheme +
    scale_fill_manual(values = rev(c(pal[5], pal[2]))) +
    scale_color_manual(values = rev(c(pal[5], pal[2]))) +
    scale_x_continuous(limits = c(2, 12), breaks = seq(2, 12, 2)) +
    scale_y_continuous(limits = c(0, upper_lim))
  
  p_match2 = ggplot(tmpmatched) + geom_density(aes(x=intervalFirstlast, fill=factor(converter)), alpha=0.2) + theme_classic() + mytheme +
    scale_fill_manual(values = rev(c(pal[5], pal[2]))) +
    scale_color_manual(values = rev(c(pal[5], pal[2]))) +
    labs(x = "Follow-up interval") +
    scale_x_continuous(limits = c(0, 16)) +
    scale_y_continuous(limits = c(0, upper_lim))
  
  p_match3 = ggplot(tmpmatched) + geom_density(aes(x=meanAge, fill=factor(converter)), alpha=0.2) + ggtitle("Matched") + theme_classic() + mytheme +
    scale_fill_manual(values = rev(c(pal[5], pal[2]))) +
    scale_color_manual(values = rev(c(pal[5], pal[2]))) +
    labs(x = "Age (mean)") +
    scale_x_continuous(limits = c(40, 90)) +
    scale_y_continuous(limits = c(0, 0.075)) +
    labs(subtitle = paste0("conv = ", unname(n_matched[2]), "; Aβ- = ", unname(n_matched[1]))) +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 14, colour = "grey20"))
  
  
  (p_matchoverview = ggarrange(
    ggarrange(p_unmatch3 + theme(legend.position = "none"), p_match3 + theme(legend.position = "none", axis.title.y = element_text(color ="white"))),
    ggarrange(p_unmatch2 + theme(legend.position = "none"), p_match2 + theme(legend.position = "none", axis.title.y = element_text(color ="white"))),
    ggarrange(p_unmatch1 + theme(legend.position = "none"), p_match1 + theme(legend.position = "none", axis.title.y = element_text(color ="white"))),
    nrow = 3)
  )
  
  
  # ggsave(filename = paste0(b, "/paper2/figs_yearsBeforeAB/pMatch_", analysname, "_AB", yearsBeforeAB,  ".png"),
  #        plot = p_matchoverview,
  #        width = 30,
  #        height = 17,
  #        dpi = 600,
  #        units = "cm"
  #        # device = cairo_pdf
  # )
  
  
  # make list of converter / nonconverter subjects included in the matching
  includelist = 
    c(
      matched_samples1$subject_id[matched_samples1$converter == 1], #converters (always same)
      matched_samples1$subject_id[matched_samples1$converter == 0],
      matched_samples2$subject_id[matched_samples2$converter == 0]
  )
  
  # make list of converter / nonconverter subjects NOT included in the matching
  notincludelist = U_convnonconv$subject_id[!U_convnonconv$subject_id %in% includelist]
  
  
  print("filtering out nonmatched subs")  
  DF %<>% filter(!subject_id %in% notincludelist)
  
}
dim(DF); length(unique(DF$subject_id))



rois=allFeat[grepl("thickness", allFeat)]
nrois=length(rois)
subset.size=nrois; jj = 1
N = ceiling(nrois/subset.size)
print(N)
start = (jj*subset.size)-subset.size+1
end = nrois
loopend = end
print(paste("subsetting cols", start, "-", end))
rois = rois[start:end]
print(rois)
ROIs = rois


pb = txtProgressBar(min=2, max=end, style=3)
Usubs = length(unique(DF$subject_id))


# log MRIs from AB- group
DF.nonconv = DF %>% filter(subject_id %in% nonconverterlist)
DF %>% filter(subject_id %in% nonconverterlist) %>% dim()


# log MRIs from converter group
# is already in DF.conv
DF %>% filter(subject_id %in% converterlist) %>% dim()
dim(DF.conv); length(unique(DF.conv$subject_id))


DF.both = rbind(DF.conv %>% mutate(convgroup = 1),
      DF.nonconv %>% mutate(convgroup = 0))


fullDF = DF


# LOOP ROIS  ------------------------------------------------
for (i in 1:length(ROIs)){
  
  print(paste(i,"/",subset.size))
  if (i == 1) {
    derivMat = absChangeMat = outlierMat1 = outlierMat2 = rIntonlyMat = rIntMat = rSlopeMat = matrix(NA, nrow = Usubs, ncol = loopend)
    predMat = seMat = matrix(NA, nrow = nrow(DF), ncol = loopend)
    RR = list()
  }
  DF = fullDF
  setTxtProgressBar(pb,i)
  set.seed(123)
  ROI = ROIs[i]
  # ROI = "lh_superiorfrontal_thickness.aparcnative71"
  # ROI = "lh_rostralmiddlefrontal_thickness.aparcnative71"
  DF$brainvar = DF[[ROI]]

  
  ggplot(DF, aes(y=brainvar, x = visit_age, col = scanStrength)) + geom_point(aes(group = subject_id)) + geom_smooth(method = "gam")
  
  
  # calc time_diff from mean age
  DF %<>% mutate(time_diff = visit_age - meanAge)
  
  
  analysisString = "tesla"
  print("running GAMMS")
  tic()
  
  
  # GAMM models ---------
  if (!singleT) {
    
    # ------ main GAMM model ------
    g = gamm4(brainvar ~ s(visit_age, k = 8) + subject_sex + cohort + scanStrength + ICV, data = DF, random = ~ (1+visit_age | subject_id))
    
    # ------ model with no slopes (intercept-only) ------
    g_noslope = gamm4(brainvar ~ s(visit_age, k = 8) + subject_sex + cohort + scanStrength + ICV, data = DF, random = ~ (1 | subject_id))
  
    # ------ T2 interaction models -----
    g_oys = gamm4(brainvar ~ t2(meanAge, time_diff, k = c(8, 3, 3)) + subject_sex + cohort + scanStrength + ICV, data = DF, random = ~ (1+time_diff | subject_id))
    g_oysnoslope = gamm4(brainvar ~ t2(meanAge, time_diff, k = c(8, 3, 3)) + subject_sex + cohort + scanStrength + ICV, data = DF, random = ~ (1 | subject_id))
    
  } else if (singleT) {
    
    # ------ main GAMM model ------
    g = gamm4(brainvar ~ s(visit_age, k = 8) + subject_sex + cohort + ICV, data = DF, random = ~ (1+visit_age | subject_id))
    
    # ------ model with no slopes (intercept-only) ------
    g_noslope = gamm4(brainvar ~ s(visit_age, k = 8) + subject_sex + cohort + ICV, data = DF, random = ~ (1 | subject_id))
    
    # ------ T2 interaction models ------
    g_oys = gamm4(brainvar ~ t2(meanAge, time_diff, k = c(8, 3, 3)) + subject_sex + cohort + ICV, data = DF, random = ~ (1+time_diff | subject_id))
    g_oysnoslope = gamm4(brainvar ~ t2(meanAge, time_diff, k = c(8, 3, 3)) + subject_sex + cohort + ICV, data = DF, random = ~ (1 | subject_id))
  }
  print("finished")
  toc()
  
  g.sum = summary(g$gam)
  dug = itsadug::get_predictions(g$gam,
                                 cond = list(visit_age = seq(
                                   min(DF$visit_age, na.rm = T),
                                   max(DF$visit_age, na.rm =
                                         T),
                                   length.out = nrow(DF)
                                 ), subject_sex = "Female",
                                 cohort = "ADNINC",
                                 scanStrength = "1-5T",
                                 ICV = 0,
                                 se = T))
  

  predictions <- DF %>% 
    mutate(subject_sex = "Female",  cohort = "ADNINC", scanStrength = "1-5T", ICV=0) %>% 
    select(visit_age,subject_sex,ICV,mri_info_site_name,scanStrength, cohort) %>% 
    predict(g$gam, newdata = ., se.fit = T)
 
  
  residualsg <- residuals(g$mer)
  DF$partial_residuals = predictions$fit + residualsg
  DF$fit = predictions$fit
  DF$sefit = predictions$se.fit
  DF$cifit = predictions$se.fit*1.96
  # plot(DF$visit_age, DF$partial_residuals)
  # points(DF$visit_age,DF$fit,col="blue")
  # plot(DF$visit_age,DF$brainvar)
  
   
  # correct for fixed effects of sex, scanner, and field strength (if not singleT)
  DF$brainvarcor = DF$brainvar
  DF$brainvarcor = ifelse(DF$subject_sex == "Male", DF$brainvarcor - g.sum$p.coeff["subject_sexMale"], DF$brainvarcor)
  DF$brainvarcor = ifelse(DF$cohort == "LCBC", DF$brainvarcor - g.sum$p.coeff["cohortLCBC"], DF$brainvarcor)
  DF$brainvarcor = ifelse(DF$cohort == "BACS", DF$brainvarcor - g.sum$p.coeff["cohortBACS"], DF$brainvarcor)
  if (!singleT) {
    DF$brainvarcor = ifelse(DF$scanStrength == "3T", DF$brainvarcor + g.sum$p.coeff["scanStrength3T"], DF$brainvarcor)
  }
  

  #colour palette ---
  pal = wesanderson::wes_palettes$FantasticFox1
  pointcol = "#6faca8"
  #colour palette ---
  
  
  #compute linear slopes across residuals
  DF$residuals_noslope = residuals(g_noslope$gam)
  manual_slopes <- DF %>%
    group_by(subject_id) %>%
    summarise(
      manual_slope = coef(lm(residuals_noslope ~ visit_age))[2], # main manual linear slope
      manual_lmslopetime = coef(lm(brainvarcor ~ TimeBsl))[2], # check linear slope across corrected data
      manual_lmslopeage = coef(lm(brainvarcor ~ visit_age))[2]
    )
  cor(manual_slopes[,2:4])
  
  
  # name and merge random effects -----
  rr = ranef(g$mer)$subject_id
  rrintonly = ranef(g_noslope$mer)$subject_id
  rroys = ranef(g_oys$mer)$subject_id
  rroysintonly = ranef(g_oysnoslope$mer)$subject_id
  
  length(unique(DF$subject_id))
  U = DF %>% group_by(subject_id) %>% mutate(meanbrainvarcor = mean(brainvarcor),
                                             meanbrainvar = mean(brainvar)) %>% ungroup() %>% filter(TimeBsl == 0)
  
  names(rr) = c("rInt", "rSlope")
  names(rrintonly) = c("rIntonly")
  names(rroys) = c("rIntOys", "rSlopeOys")
  names(rroysintonly) = c("rIntonlyOys")
  
  
  rr$subject_id = row.names(rr)
  rrintonly$subject_id = row.names(rrintonly)
  rroys$subject_id = row.names(rroys)
  rroysintonly$subject_id = row.names(rroysintonly)
  
  
  rr = merge(manual_slopes, rr, by = "subject_id")
  rr = merge(rr, rrintonly, by = "subject_id")
  rr = merge(rr, rroys, by = "subject_id")
  rr = merge(rr, rroysintonly, by = "subject_id")
  U = merge(U, rr, by = "subject_id")
  
  
  # SI trajectory figs
  (fig1 = 
      DF %>% filter(subject_id != adnioutlier) %>% 
      ggplot(.) +
      geom_line(data=DF,aes(x=visit_age,brainvarcor,group=subject_id),color=pointcol,alpha=0.6, size=0.5) +
      geom_point(data=DF,aes(x=visit_age,brainvarcor,group=subject_id),color=pointcol,stat="identity",alpha=1, size=0.5) +
      geom_ribbon(data=dug,aes(x=visit_age,ymin=fit-CI,ymax=fit+CI),alpha=.7,show.legend=F,fill="dark grey") +
      geom_line(data=DF,aes(x=visit_age,y=fit),col="black") +
      ggtitle(ROI) +
      labs(x = "Age") +
      theme_classic() + mytheme)
  
  
  
  # add in converter/nonconverter colours
  DF$convertgroup = NA
  DF$convertgroup[DF$subject_id %in% converterlist] = 1
  DF$convertgroup[DF$subject_id %in% nonconverterlist] = 0
  DF$convertgroup = as.factor(DF$convertgroup)
  DF %>% filter(convertgroup == 1) %>% dim()
  DF %>% filter(convertgroup == 0) %>% dim()
  
  
  #NB! this includes NA
  unique(DF$subject_id[DF$convertgroup == 1])
  DF$subject_id[DF$convertgroup == 1 & !is.na(DF$convertgroup)] %in% unique(DF.conv$subject_id)
  length(unique(DF$subject_id[DF$convertgroup == 1]))
  length(unique(DF.conv$subject_id))
  
  
  tmp = DF %>% filter(!is.na(convertgroup))
  tmp1 = tmp %>% filter(convertgroup == 1)
  tmp0 = tmp %>% filter(convertgroup == 0)
  
  
  tmp = DF
  tmp$newfac = as.character(tmp$convertgroup)
  tmp$newfac[is.na(tmp$newfac)] = "noData"
  tmpNo = tmp[tmp$newfac == "noData",]
  
  
  (p_allgroup = 
      tmp %>% filter(subject_id != adnioutlier) %>% 
    ggplot(., aes(x = visit_age, y = brainvarcor)) +
    geom_line(aes(col = newfac, group = subject_id), alpha = 0.6, size = 0.4) +
    scale_color_manual(values = rev(c("light grey", pal[5], pal[2]))) +
    geom_line(data = tmpNo, aes(group = subject_id), alpha = 0.3, size = 1, col = "light grey") +
    geom_line(data = tmp0, aes(group = subject_id), alpha = 0.3, size = 1, col = pal[2]) +
    geom_point(data = tmp0, aes(group = subject_id), alpha =1, size = 0.5, col = pal[2]) +
    geom_line(data = tmp1, aes(group = subject_id), alpha = 0.3, size = 1, col = pal[5]) +
    geom_point(data = tmp1, aes(group = subject_id), alpha = 1, size = 0.5, col = pal[5]) +
    geom_line(data=DF, aes(x=visit_age,y=fit),col="black") +
    theme_classic() + mytheme +
    ggtitle(ROI) +
    labs(x = "Age"))
    
  
  
  #check outliers
  outlierThresh = function(U, outlierthresh, ranef) {
    
    if (!singleT) {
      g1 = gam(get(ranef) ~ s(meanAge, k = 8) + subject_sex + cohort + scanStrength, data = U)
    } else if (singleT) {
      g1 = gam(get(ranef) ~ s(meanAge, k = 8) + subject_sex + cohort, data = U)
    }
    Opt.age = U$meanAge
    Q = scale(g1$residuals)
    db.out = data.frame(U$subject_id,Q) #residuals scaled Z
    db.out$QQ = NA
    db.out$QQ [ (abs(Q)>outlierthresh)] = Q[(abs(Q)>outlierthresh)] #places outliers in new col
    outliers = db.out$U.subject_id[which(!is.na(db.out$QQ))] #indexes
    cat("\noutlier thresh =", outlierthresh, "SD from GAMM model: highlighting", length(outliers), "cases for potential removal\n\n")
    N.outliers = length(outliers)
    
    if (ranef == "rSlope") {
      U$Qcol1="cyan2"
      U$Qcol1[(U$subject_id %in% outliers)]="black"
    } else {
      U$Qcol2="cyan2"
      U$Qcol2[(U$subject_id %in% outliers)]="black"
    }
    return(list(outliers = as.numeric(as.character(outliers)),
                U))
  }
  ranef="rSlope"
  thr = 6
  check = outlierThresh(U, thr, "rSlope")[[2]]
  # outlierMat1[,i] = outlierThresh(U, thr, "rSlope")[[2]] #rSlope
  
  
  
  (fig2 = ggplot() +
      U %>% filter(subject_id != adnioutlier) %>%
      geom_point(data=.,aes(x=meanAge,rSlope),color=pointcol,stat="identity",alpha=0.7) +
      geom_hline(yintercept = 0,linetype=2,size=1) +
      labs(y="rSlope",
           x="Age (mean)") +
      ggtitle("Age-relative change") + 
      ylab(expression("Additional Change (mm"^3*"/year)")) +
      theme_classic() + mytheme +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
      ) )
  
  
  predDat = U %>% 
    mutate(subject_sex = "Female",  cohort = "ADNINC", scanStrength = "1-5T", ICV=0) %>% 
    select(visit_age,subject_sex,ICV,mri_info_site_name,scanStrength, cohort)
  
  
  # calculate derivatives
  dd = derivatives(object=g$gam,
                 data=predDat,
                 eps = 1e-07, #default
                 order = 1L,
                 n = nrow(U),
                 interval = c("confidence")
  )
  
  
  #add random effect to curve derivative to get vals of slope
  plot(U$meanAge, dd$.derivative)
  U$absChange=U$rSlope+dd$.derivative
  plot(U$meanAge,U$absChange,col="blue", main = ROI)
  
  
  (fig3 = ggplot() +
      geom_point(data=U %>% filter(subject_id != adnioutlier),aes(x=meanAge,absChange),color=pointcol,stat="identity",alpha=0.7) +
      geom_smooth(data=U %>% filter(subject_id != adnioutlier),aes(x=meanAge,absChange),method="loess",col="black", size=0.5) +
      geom_hline(yintercept = 0,linetype=2,size=1) +
      labs(
        x="Age (mean)") +
      ggtitle("Absolute change") +
      ylab(expression("Change (mm"^3*"/year)")) +
      theme_classic() + mytheme +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
      ) )
  
  
  (cp1 = cowplot::plot_grid(fig1,fig3,fig2,ncol=3))
  
  
  derivMat[,i] = dd$.derivative
  absChangeMat[,i] = U$absChange
  rSlopeMat[,i] = U$rSlope
  rIntMat[,i] = U$rInt
  rIntonlyMat[,i] = U$rIntonly
  predMat[,i] = DF$fit
  seMat[,i] = DF$sefit
  RR$sTab[[i]] = g.sum$s.table
  RR$pTab[[i]] = g.sum$p.table
  
  
  # run/save quick models for checking against analysis script
  resdir = file.path(b, paste0("results/prepMods/", analysname, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime))
  
  
  Uconv = rbind(
    U %>% filter(subject_id %in% nonconverterlist) %>% mutate(convgroup = 0),
    U %>% filter(subject_id %in% converterlist) %>% mutate(convgroup = 1)
  )
  
  Uconvnoout = Uconv[Uconv$subject_id != adnioutlier,]
  tmpmod = tidy(summary(lm(scale(rSlope) ~ (convgroup) + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = Uconvnoout)))
  

  if (!dir.exists(resdir)) { 
    system(paste("mkdir -p", resdir))
  }
  
  
  sigflag = ifelse(tmpmod$p.value[2] < .05, 1, 0)
  filename1 = paste0(resdir, "/prepModNout-yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-sig", sigflag,"-",ROI,".csv")
  filename2 = paste0(resdir, "/UconvNout-yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-",ROI,".csv")
  print(paste("saving", filename1))
  print(paste("saving full data without adnioutlier", filename2))
  
  
  write.table(
    tmpmod, filename1
    )
  write.table(
    Uconvnoout, filename2
  )
  
  
  if (savefigs) {
    plotdir = file.path(b, paste0("plots/cowplots/", analysname, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "/", ROI, "-sig", sigflag))
    
    if (!dir.exists(plotdir)) { 
      system(paste("mkdir -p", plotdir))
    }
    
    
    filename3 = file.path(plotdir,(paste0("pAll-yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-", ROI,"-sig", sigflag, ".pdf")))
    print(paste("saving fig", filename3))
      
    ggsave(filename = filename3,
           plot = cp1,
           width = 30,
           height = 13,
           dpi = 600,
           units = "cm",
           device = cairo_pdf
    )
      
  }
  
}


# SAVERES ------
resdir = file.path(b,
                 paste0("results/prepSlopes/", analysname))
ofile =
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
if (! dir.exists(resdir)) {
  print(paste("saving results to", ofile))
  system(paste("mkdir -p", resdir))
}


ROIs.all=ROIs
U.all=U
Uconv.all = Uconv
Uconvnoout.all = Uconvnoout
derivMat.all=derivMat
absChangeMat.all=absChangeMat
rIntMat.all=rIntMat
rIntonlyMat.all=rIntonlyMat
rSlopeMat.all=rSlopeMat
predMat.all=predMat
seMat.all=seMat
RR.all=RR
# outlierMat1.all=outlierMat1
DF.all=DF
DF.conv.all = DF.conv
save(
  "DF.all",
  "ROIs.all",
  'U.all',
  "derivMat.all",
  "absChangeMat.all",
  "rIntMat.all",
  "rIntonlyMat.all",
  "rSlopeMat.all",
  "predMat.all",
  "seMat.all",
  "RR.all",
  # "outlierMat1.all",
  "yearsBeforeAB",
  "usable",
  "unusable",
  "DF.conv.all",
  file = ofile
)

quit()
# END ------
