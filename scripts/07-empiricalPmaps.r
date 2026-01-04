#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: get empirical p-value maps using null models from regional wild bootstrap resampling
#========================================================================================#


#================ INPUTS =================#
args = commandArgs(TRUE)
yearsBeforeAB = as.numeric(args[1])
analysname = as.character(args[2])
metric = as.character(args[3])
nTime = as.numeric(args[4])
# j = sprintf("%04d", k)
#=========================================#

library("data.table"); library("dplyr"); library("ggplot2"); library("MetBrewer")
#testing
# rm(list=ls())
# i=1
# j = sprintf("%04d", i)
# yearsBeforeAB=1
# metric = "thickness"
# analysname = "analysnum1matchF0procLONGpredAB0singleT0reproduce1"
# nTime = 2



samplerois = c("lh_bankssts_thickness.aparcnative71",
               "lh_caudalanteriorcingulate_thickness.aparcnative71",
               "lh_caudalmiddlefrontal_thickness.aparcnative71",
               "lh_cuneus_thickness.aparcnative71",
               "lh_entorhinal_thickness.aparcnative71",
               "lh_fusiform_thickness.aparcnative71",
               "lh_inferiorparietal_thickness.aparcnative71",
               "lh_inferiortemporal_thickness.aparcnative71",
               "lh_isthmuscingulate_thickness.aparcnative71",
               "lh_lateraloccipital_thickness.aparcnative71",
               "lh_lateralorbitofrontal_thickness.aparcnative71",
               "lh_lingual_thickness.aparcnative71",
               "lh_medialorbitofrontal_thickness.aparcnative71",
               "lh_middletemporal_thickness.aparcnative71",
               "lh_parahippocampal_thickness.aparcnative71",
               "lh_paracentral_thickness.aparcnative71",
               "lh_parsopercularis_thickness.aparcnative71",
               "lh_parsorbitalis_thickness.aparcnative71",
               "lh_parstriangularis_thickness.aparcnative71",
               "lh_pericalcarine_thickness.aparcnative71",
               "lh_postcentral_thickness.aparcnative71",
               "lh_posteriorcingulate_thickness.aparcnative71",
               "lh_precentral_thickness.aparcnative71",
               "lh_precuneus_thickness.aparcnative71",
               "lh_rostralanteriorcingulate_thickness.aparcnative71",
               "lh_rostralmiddlefrontal_thickness.aparcnative71",
               "lh_superiorfrontal_thickness.aparcnative71",
               "lh_superiorparietal_thickness.aparcnative71",
               "lh_superiortemporal_thickness.aparcnative71",
               "lh_supramarginal_thickness.aparcnative71",
               "lh_frontalpole_thickness.aparcnative71",
               "lh_temporalpole_thickness.aparcnative71",
               "lh_transversetemporal_thickness.aparcnative71",
               "lh_insula_thickness.aparcnative71",
               "lh_MeanThickness_thickness.aparcnative71",
               "rh_bankssts_thickness.aparcnative71",
               "rh_caudalanteriorcingulate_thickness.aparcnative71",
               "rh_caudalmiddlefrontal_thickness.aparcnative71",
               "rh_cuneus_thickness.aparcnative71",
               "rh_entorhinal_thickness.aparcnative71",
               "rh_fusiform_thickness.aparcnative71",
               "rh_inferiorparietal_thickness.aparcnative71",
               "rh_inferiortemporal_thickness.aparcnative71",
               "rh_isthmuscingulate_thickness.aparcnative71",
               "rh_lateraloccipital_thickness.aparcnative71",
               "rh_lateralorbitofrontal_thickness.aparcnative71",
               "rh_lingual_thickness.aparcnative71",
               "rh_medialorbitofrontal_thickness.aparcnative71",
               "rh_middletemporal_thickness.aparcnative71",
               "rh_parahippocampal_thickness.aparcnative71",
               "rh_paracentral_thickness.aparcnative71",
               "rh_parsopercularis_thickness.aparcnative71",
               "rh_parsorbitalis_thickness.aparcnative71",
               "rh_parstriangularis_thickness.aparcnative71",
               "rh_pericalcarine_thickness.aparcnative71",
               "rh_postcentral_thickness.aparcnative71",
               "rh_posteriorcingulate_thickness.aparcnative71",
               "rh_precentral_thickness.aparcnative71",
               "rh_precuneus_thickness.aparcnative71",
               "rh_rostralanteriorcingulate_thickness.aparcnative71",
               "rh_rostralmiddlefrontal_thickness.aparcnative71",
               "rh_superiorfrontal_thickness.aparcnative71",
               "rh_superiorparietal_thickness.aparcnative71",
               "rh_superiortemporal_thickness.aparcnative71",
               "rh_supramarginal_thickness.aparcnative71",
               "rh_frontalpole_thickness.aparcnative71",
               "rh_temporalpole_thickness.aparcnative71",
               "rh_transversetemporal_thickness.aparcnative71",
               "rh_insula_thickness.aparcnative71",
               "rh_MeanThickness_thickness.aparcnative71")
# "pcaL",
# "pcaR")


# load in null distributions and calculate empirical p-value
# compute empirical p-values --------
for (jj in 1:length(samplerois)) {
  
  if (jj == 1) {
    tmpnullmodlist = tmpnullmodlist_centcor = list()
    empirical_p = empirical_t = c()
    empirical_pcorr = empirical_tcorr = c()
    empirical_dens = empirical_denslog = list()
    
    empirical_p_centcor = empirical_t_centcor = c()
    empirical_pcorr_centcor = empirical_tcorr_centcor = c()
    empirical_dens_centcor = empirical_denslog_centcor = list()
  }
  
  # jj = 27
  print(jj)
  sampleroi = samplerois[jj]
  
  odir = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/results/cortmaps_resample/", analysname, "/", metric, "_nullmod/yearsBeforeAB", yearsBeforeAB, "nTime", nTime,
                "/", sampleroi)
  ifile1 = file.path(odir, paste0("nullresult_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-",sampleroi, ".csv"))
  ifile2 = file.path(odir, paste0("nullresultcentcor_yearsBeforeAB", yearsBeforeAB, "nTime", nTime, "-",sampleroi, ".csv"))
  nulldist = fread(ifile1)
  nulldist_centcor = fread(ifile2)
  
  #compute empirical p
  tmpnullmodlist[[jj]] = nulldist
  tmpnullmodlist_centcor[[jj]] = nulldist_centcor
  tmpdf = nulldist
  tmpdf2 = nulldist_centcor
  
  
  #load in cortmaps for observed p values
  idir = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/results/cortmaps/", analysname, "/", metric)
  map1 = paste0(metric, "_cortmapAB", yearsBeforeAB, "_nTime", nTime, "_", analysname, ".csv")
  map2 = paste0(metric, "_cortmapAB", yearsBeforeAB, "_nTime", nTime, "_", analysname, "_centcor.csv")
  
  cortmap = fread(file.path(idir, map1))
  cortmap_centcor = fread(file.path(idir, map2))
  
  observed_p = unname(unlist(cortmap[jj, "p.value"]))
  observed_t = unname(unlist(cortmap[jj, "statistic"]))
  
  observed_p_centcor = unname(unlist(cortmap_centcor[jj, "p.value"]))
  observed_t_centcor = unname(unlist(cortmap_centcor[jj, "statistic"]))
  
  ggplot(tmpdf) + geom_density(aes(x = -log10(p.value))) +
    geom_vline(xintercept = unlist(-log10(cortmap[jj, "p.value"])), col="red")
  
  pal = met.brewer("Derain", n=8)
  
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
  
  empirical_dens[[jj]] = ggplot(tmpdf) + 
    geom_vline(xintercept = unlist(cortmap[jj, "p.value"]), col="#f8e025", size =1) +
    geom_density(aes(x = p.value), col="black", fill = "NA", alpha = 0.8, size=0.5) +
    # geom_vline(xintercept = (cortmap[jj, "p.value"]), col="#EFC86E", size =1) + #f8e025
    labs(x = "p-value (null distribution)") +
    theme_classic() + mytheme
  
  empirical_denslog[[jj]] = ggplot(tmpdf) + 
    geom_vline(xintercept = unlist(observed_t), col="#f8e025", size =1) +
    geom_density(aes(x = statistic), col="black", fill = "NA", alpha = 0.8, size=0.5) +
    labs(x = "t (null distribution)") +
    theme_classic() + mytheme
  
  empirical_denslog_centcor[[jj]] = ggplot(tmpdf2) + 
    geom_vline(xintercept = unlist(observed_t_centcor), col="#f8e025", size =1) +
    geom_density(aes(x = statistic), col="black", fill = "NA", alpha = 0.8, size=0.5) +
    labs(x = "t (null distribution)") +
    theme_classic() + mytheme
  
  
  #compute empirical p
  (empirical_p[jj] <- (sum(abs(nulldist$p.value) <= abs(observed_p))) / (length(nulldist$p.value)))
  
  #now equal to
  mean(as.vector(abs(nulldist$statistic)) >= abs(observed_t))
  
  (empirical_t[jj] <- (sum(abs(nulldist$statistic) >= abs(observed_t))) / (length(nulldist$statistic)))
  (empirical_tcorr[jj] <- (sum(abs(nulldist$statistic) >= abs(observed_t)) + 1) / (length(nulldist$statistic) + 1))
  
  (empirical_pcorr[jj] <- (sum(abs(nulldist$p.value) <= abs(observed_p)) + 1) / (length(nulldist$p.value) + 1))
  (empirical_tcorr[jj] <- (sum(abs(nulldist$statistic) >= abs(observed_t)) + 1) / (length(nulldist$statistic) + 1))
  
  (empirical_p_centcor[jj] <- (sum(abs(nulldist_centcor$p.value) <= abs(observed_p_centcor))) / (length(nulldist_centcor$p.value)))
  (empirical_t_centcor[jj] <- (sum(abs(nulldist_centcor$statistic) >= abs(observed_t_centcor))) / (length(nulldist_centcor$p.value)))
  
  (empirical_pcorr_centcor[jj] <- (sum(abs(nulldist_centcor$p.value) <= abs(observed_p_centcor)) + 1) / (length(nulldist_centcor$p.value) + 1))
  (empirical_tcorr_centcor[jj] <- (sum(abs(nulldist_centcor$statistic) >= abs(observed_t_centcor)) + 1) / (length(nulldist_centcor$statistic) + 1))
  
  
  savefigs=1
  if (savefigs) {
    plotdir = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/plots/empirical_densplots/", analysname, "/yearsBeforeAB", yearsBeforeAB, "nTime", nTime)
    if (!dir.exists(plotdir)) { 
      system(paste("mkdir -p", plotdir))
    }
    sigflag = ifelse(observed_p < .05, 1, 0)
    esigflag = ifelse(empirical_p[jj] < .05, 1, 0)
    sigflag_centcor = ifelse(observed_p_centcor < .05, 1, 0)
    esigflag_centcor = ifelse(empirical_p_centcor[jj] < .05, 1, 0)
    filename1 = paste0(plotdir, "/pEmpirical-sig", sigflag,"-", "esig", esigflag, "-", samplerois[jj],".pdf")
    filename2 = paste0(plotdir, "/pEmpiricalcentcor-sig", sigflag_centcor,"-", "esig", esigflag_centcor, "-", samplerois[jj],".pdf")
    
    print(paste("saving", filename1))
    ggsave(plot = empirical_denslog[[jj]] + theme(axis.ticks = element_line()) + ggtitle(samplerois[[jj]]),
           filename = filename1,
           width=9, height=6, units="cm", dpi=600,
           device = cairo_pdf)
    
    ggsave(plot = empirical_denslog_centcor[[jj]] + theme(axis.ticks = element_line()) + ggtitle(samplerois[[jj]]),
           filename = filename2,
           width=9, height=6, units="cm", dpi=600,
           device = cairo_pdf)
  }
  
}


empirical_pmap1 = data.frame(samplerois, p.value = empirical_p, 
                             statistic = empirical_t, 
                             empirical_pcorr,
                             empirical_tcorr
)
empirical_pmap2 = data.frame(samplerois, p.value = empirical_p_centcor, 
                             statistic = empirical_t_centcor, 
                             empirical_pcorr_centcor,
                             empirical_tcorr_centcor
)

writeMaps = 1
if (writeMaps) {
  cortmapdir = paste0("/cluster/projects/p274/projects/p040-ad_change/Berkeley/results/cortmaps/", analysname, "/", metric, "_null")
  if (!dir.exists(cortmapdir)) {
    dir.create(cortmapdir)
  }
  print("writing empirical map")
  emp_map1 = paste0(cortmapdir, "/empirical_", metric, "_cortmapAB", yearsBeforeAB, "_nTime", nTime, "_", analysname, ".csv")
  print(paste("writing", emp_map1))
  write.table(empirical_pmap1,
              emp_map1, 
              row.names = F, col.names = T, quote = F)
  emp_map2 = paste0(cortmapdir, "/empiricalcentcor_", metric, "_cortmapAB", yearsBeforeAB, "_nTime", nTime, "_", analysname, ".csv")
  print(paste("writing", emp_map2))
  write.table(empirical_pmap2,
              emp_map2, 
              row.names = F, col.names = T, quote = F)
}

