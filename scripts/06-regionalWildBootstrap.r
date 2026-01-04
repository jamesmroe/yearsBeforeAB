#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: run regional wild bootstrap resampling (guard against group differences in variance)
#========================================================================================#


#================ INPUTS =================#
args = commandArgs(TRUE)
roi = as.character(args[1])
yearsBeforeAB = as.numeric(args[2])
analysname = as.character(args[3])
metric = as.character(args[4])
nTime = as.numeric(args[5])
#=========================================#


#testing
# rm(list=ls())
# i=1
# roi="lh_bankssts_thickness.aparcnative71"
# j = sprintf("%04d", i)
# yearsBeforeAB=1
# metric = "thickness"
# analysname = "analysnum1matchF0procLONGpredAB0singleT0reproduce1"
# nTime = 2


print("-----------------")
print("loading packages")
print("-----------------")
library("dplyr"); library("data.table"); library("broom")


idir = paste0(
  paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/results/cortmaps_resample/", analysname, "/", metric, "_nullmod/yearsBeforeAB", yearsBeforeAB, "nTime", nTime))

ifile = paste0(idir, "/", roi, "/region_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, "-", roi, ".csv")


print("-----------------")
print(paste("loading", ifile))
print("-----------------")
sample = fread(ifile)


iter = 1000
for (ii in 1:iter) {
  if (ii == 1) { 
    nullmodlist = nullmodlist_centcor = list() 
  }
  print(ii)
  
  wild <- sample(c(1, -1), size = nrow(sample), replace = TRUE)
  sample$nullSlope = sample$nullpred + (sample$nullresid*wild)
  sample$nullSlope_centcor = sample$nullpred_centcor + (sample$nullresid_centcor*wild)
  
  nullmodlist[[ii]] = tidy(lm(scale(nullSlope) ~ scale(convgroup) + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = sample))
  nullmodlist_centcor[[ii]] = tidy(lm(scale(nullSlope_centcor) ~ scale(convgroup) + meanCentiloids + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = sample))
  
}


# write maps
nulldistregion = data.frame(roi, bind_rows(nullmodlist) %>% filter(term == "scale(convgroup)")) %>% mutate(
  sig = ifelse(p.value <.05, 1, 0),
  N_conv = unname(table(sample$convgroup)[2]),
  N_never = unname(table(sample$convgroup)[1]),
  yearsBeforeAB = yearsBeforeAB,
  ifile = ifile
)
nulldistregion_centcor = data.frame(roi, bind_rows(nullmodlist_centcor) %>% filter(term == "scale(convgroup)")) %>% mutate(
  sig = ifelse(p.value <.05, 1, 0),
  N_conv = unname(table(sample$convgroup)[2]),
  N_never = unname(table(sample$convgroup)[1]),
  yearsBeforeAB = yearsBeforeAB,
  ifile = ifile
)

ofile1 = paste0(idir, "/", roi, "/nullresult_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, "-", roi, ".csv")
ofile2 = paste0(idir, "/", roi, "/nullresultcentcor_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, "-", roi, ".csv")

print("-----------------")
print(paste("saving", ofile1))
print("-----------------")
fwrite(nulldistregion,
       ofile1
)
print("-----------------")
print(paste("saving", ofile2))
print("-----------------")
fwrite(nulldistregion_centcor, 
       ofile2
)
quit()
