#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Summarize results from resampling-based robustness check
#========================================================================================#


# rm(list=ls())
library("readr"); library("dplyr"); library("data.table"); library("magrittr")


#================ INPUTS =================#
args = commandArgs(TRUE)
analysname = as.character(args[1])
yearsBeforeAB = as.character(args[2])
metric = as.character(args[3])
nTime = as.character(args[4])
#=========================================#

# for (yearsBeforeAB in c(1:10)) {

# testing
# rm(list=ls())
# yearsBeforeAB=2
# metric = "thickness"
# analysname = "analysnum1matchF0procLONGpredAB0singleT0reproduce1"
# nTime = 2


print(analysname)
print(yearsBeforeAB)
print(metric)
print(nTime)


b = "/cluster/projects/p274/projects/p040-ad_change/Berkeley"
ufile = paste0("U_brainslopes_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, ".csv")

resdir = paste0(b, "/results/cortmaps_resample/", analysname, "/", metric, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime)
cortmapdir = paste0(b, "/results/cortmaps/", analysname, "/", metric)
print(paste("checking", cortmapdir))
files <- list.files(path = resdir, pattern = "results.*cortmap.*", full.names = TRUE)

if (length(files) < 2000) { print("the number of resampled results is < 2000. Check this."); quit() }
modelfiles1 = files[!grepl("centcor", files)]
modelfiles2 = files[grepl("centcor", files)]


# if (!exists("df_models")) {
  print("loading all model outputs")
  df_models = read_csv(modelfiles1, id = "src") |>
    mutate(src = basename(src))
# }
# if (!exists("df_modelscentcor")) {
  print("loading all model outputs centcor")
  df_modelscentcor = read_csv(modelfiles2, id = "src") |>
    mutate(src = basename(src))
# }


U.all.tmp = fread(file.path(resdir, ufile))
ROIs.all = names(U.all.tmp)[1:364]


pcaIn = 0
if (pcaIn == 1) {
  #sort out pca directions first
  signmeanL = sign(df_models$estimate[df_models$rois == "lh_MeanThickness_thickness.aparcnative71"])[1] 
  signpcaL = sign(df_models$estimate[df_models$rois == "pca_Left"])[1] 
  signmeanR = sign(df_models$estimate[df_models$rois == "rh_MeanThickness_thickness.aparcnative71"])[1] 
  signpcaR = sign(df_models$estimate[df_models$rois == "pca_Right"])[1] 
  if (signpcaL != signmeanL) {
    print("inversing pca direction Left to be same as mean")
    df_models$estimate[df_models$rois == "pca_Left"] = df_models$estimate[df_models$rois == "pca_Left"]*-1
  }
  if (signpcaL != signmeanL) {
    print("inversing pca direction Right to be same as mean")
    df_models$estimate[df_models$rois == "pca_Right"] = df_models$estimate[df_models$rois == "pca_Right"]*-1
  }
}

df_models %<>% mutate(sigdirection = case_when(
  sig == 1 & sign(estimate) == 1 ~ "positive",
  sig == 1 & sign(estimate) == -1 ~ "negative",
  TRUE ~ NA_character_
))
table(df_models$sigdirection)

df_modelscentcor %<>% mutate(sigdirection = case_when(
  sig == 1 & sign(estimate) == 1 ~ "positive",
  sig == 1 & sign(estimate) == -1 ~ "negative",
  TRUE ~ NA_character_
))
table(df_modelscentcor$sigdirection)


df_models %<>% group_by(rois) %>% 
  mutate(nSig = sum(sig)) %>% 
  mutate(nSigCheck = sum(p.value < .05)) %>% 
  mutate(
    nPositive = sum(sigdirection == "positive", na.rm = TRUE),
    nNegative = sum(sigdirection == "negative", na.rm = TRUE)
  ) %>% 
  ungroup()


df_modelscentcor %<>% group_by(rois) %>% 
  mutate(nSig = sum(sig)) %>% 
  mutate(nSigCheck = sum(p.value < .05)) %>% 
  mutate(
    nPositive = sum(sigdirection == "positive", na.rm = TRUE),
    nNegative = sum(sigdirection == "negative", na.rm = TRUE)
  ) %>% 
  ungroup()


df_models_summarised = df_models %>% filter(i == 1)
df_modelscentcor_summarised = df_modelscentcor %>% filter(i == 1)


#write summary map
map1 = paste0("cortmap_", metric, "_resamplesummary_AB", yearsBeforeAB, "nTime", nTime, "_", analysname, ".csv")
map2 = paste0("cortmap_", metric, "_resamplesummary_AB", yearsBeforeAB, "nTime", nTime, "_", analysname, "_centcor.csv")
print(paste("writing summary map:", map1))
print(paste("writing summary map:", map2))
fwrite(df_models_summarised,
       file = file.path(cortmapdir, map1)
)
fwrite(df_modelscentcor_summarised,
       file = file.path(cortmapdir, map2)
)
plot( 1:nrow(df_models_summarised), ((df_models_summarised$nSig)/1000)*100)
plot( 1:nrow(df_models_summarised), ((df_modelscentcor_summarised$nSig)/1000)*100)
                                     
# }
