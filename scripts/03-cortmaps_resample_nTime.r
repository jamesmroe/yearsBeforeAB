#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Run resampling-based robustness analysis
#          Script requires individual-level data as input and is not executable
#========================================================================================#


#================ INPUTS =================#
args = commandArgs(TRUE)
start = as.numeric(args[1])
end = as.numeric(args[2])
yearsBeforeAB = as.numeric(args[3])
analysname = as.character(args[4])
metric = as.character(args[5])
nTime = as.numeric(args[6])
#=========================================#


#testing
# rm(list=ls())
# i=1
# j = sprintf("%04d", i)
# yearsBeforeAB=2
# metric = "thickness"
# analysname = "analysnum1matchF0procLONGpredAB0singleT0reproduce1"
# nTime = 2
# subset.size=20; jj = 1; totaltests=1000
# start = (jj*subset.size)-subset.size+1
# end = (start + subset.size)-1
# print(paste("running for inputs", start, "-", end))



print("-----------------")
print("loading packages")
print("-----------------")
library("dplyr"); library("data.table"); library("broom")



for (i in start:end) {
  
  iii = sprintf("%04d", i)
  
  b = "/cluster/projects/p274/projects/p040-ad_change/Berkeley"
  ufile = paste0("/U_brainslopes_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, ".csv")
  input = paste0("resample_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, "-", iii, ".csv")
  
  resdir = paste0(b, "/results/cortmaps_resample/", analysname, "/", metric, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime)
  U.all.tmp = fread(paste0(resdir, ufile))
  ROIs.all = names(U.all.tmp)[1:364]
  ifile = file.path(resdir, input)
  
  
  print("-----------------")
  print(paste("loading", ifile))
  print("-----------------")
  sample = fread(ifile)
  
  
  corthickrois = ROIs.all[grepl("thickness", ROIs.all, ignore.case = T)]
  rois = corthickrois
  
  
  for (ii in 1:length(rois)) {
    
    if (ii == 1) { 
      modlist = modlist_centcor = list() 
    }
    print(ii)
    roi = rois[ii]
    
    sample$brainvarSlope = NULL
    sample = left_join(sample,
              U.all.tmp %>% dplyr::select(subject_id,
                                          all_of(roi)) %>% rename(brainvarSlope = roi)
    )
    modlist[[ii]] = tidy(lm(scale(brainvarSlope) ~ scale(convgroup) + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = sample))
    modlist_centcor[[ii]] = tidy(lm(scale(brainvarSlope) ~ scale(convgroup) + meanCentiloids + subject_sex + cohort + meanAge + nTimepoints + intervalFirstlast, data = sample))
    
  }
  
  
  #write maps
  cortmap = data.frame(rois, bind_rows(modlist) %>% filter(term == "scale(convgroup)")) %>% mutate(
    sig = ifelse(p.value <.05, 1, 0),
    N_conv = unname(table(sample$convgroup)[2]),
    N_never = unname(table(sample$convgroup)[1]),
    yearsBeforeAB = yearsBeforeAB,
    ifile = ifile,
    i = i
  )
  
  cortmap_centcor = data.frame(rois, bind_rows(modlist_centcor) %>% filter(term == "scale(convgroup)")) %>% mutate(
    sig = ifelse(p.value <.05, 1, 0),
    N_conv = unname(table(sample$convgroup)[2]),
    N_never = unname(table(sample$convgroup)[1]),
    yearsBeforeAB = yearsBeforeAB,
    ifile = ifile,
    i = i
  )
  
  
  print("-----------------")
  print("writing cortmaps")
  print("-----------------")
  
  map1 = paste0(resdir, "/results_cortmap_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, "-", iii, ".csv")
  map2 = paste0(resdir, "/results_cortmapcentcor_yearsBeforeAB", yearsBeforeAB,"nTime", nTime, "-", iii, ".csv")
  
  fwrite(cortmap, map1)
  fwrite(cortmap_centcor, map2)

}
quit()
