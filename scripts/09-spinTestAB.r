#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose:  run spatial correlation between effect size maps at each cutoff
#           and 3 AB deposition maps
#========================================================================================#


rm(list=ls())


#---load packages
loadPackages = function() {
  packages = c("here", "dplyr","magrittr","corrplot","ggrepel","data.table","stringr","MetBrewer")
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  print(sapply(packages, require, character.only = T))
  print(sapply(packages, function(p) as.character(packageVersion(p))))
}
loadPackages()
here()



mainlist = list()
p_mainlist = list()
p_mainlistright = list()
df_means_mainlist = list()


for (amyloidmap in 1:3) {
  
  
  pb = txtProgressBar(min=1, max=6, style=3)
  df_means = list()
  spinlist = list()
  
  for (analysisnum in 1:6) {
    print(paste("###########", analysisnum, "###########"))
    setTxtProgressBar(pb,analysisnum)
    if (analysisnum == 1) {
      # main (all MRI)
      analysname = "analysnum1matchF0procLONGpredAB0singleT0reproduce1"
      plottext = "F0 predAB0 T0"
    }
    if (analysisnum == 2) {
      # sensitivity 1 (all MRI + matchfollow)
      analysname = "analysnum2matchF1procLONGpredAB0singleT0reproduce1" 
      plottext = "F1 predAB0 T0"
    }
    if (analysisnum == 3) {
      # sensitivity 2 (1.5T)
      analysname = "analysnum3matchF0procLONGpredAB0singleT1reproduce1"
      plottext = "F0 predAB0 T1"
    }
    if (analysisnum == 4) {
      # sensitivity 3 (1.5T + matchfollow)
      analysname = "analysnum4matchF1procLONGpredAB0singleT1reproduce1"
      plottext = "F1 predAB0 T1"
    }
    if (analysisnum == 5) {
      # sensitivity 4 (1.5T + matchfollow T2 slopes)
      analysname = "analysnum4matchF1procLONGpredAB0singleT1reproduce1"
      plottext = "F1 predAB0 T1 (slope)"
    }
    if (analysisnum == 6) {
      # sensitivity 4 (1.5T + matchfollow T2 slopes + level)
      analysname = "analysnum4matchF1procLONGpredAB0singleT1reproduce1"
      plottext = "F1 predAB0 T1 (intcor)"
    }
    
    
    
    runit=1
    if (runit) {
      base = "/cluster/projects/p274/projects/p040-ad_change/Berkeley"
      b = file.path(base, paste0("results/cortmaps/", analysname, "/thickness"))
      if (analysisnum == 5) {
        mapAB1 = fread(file.path(paste0(b, "/thickness_cortmapAB1_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB2 = fread(file.path(paste0(b, "/thickness_cortmapAB2_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB3 = fread(file.path(paste0(b, "/thickness_cortmapAB3_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB4 = fread(file.path(paste0(b, "/thickness_cortmapAB4_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB5 = fread(file.path(paste0(b, "/thickness_cortmapAB5_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB6 = fread(file.path(paste0(b, "/thickness_cortmapAB6_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB7 = fread(file.path(paste0(b, "/thickness_cortmapAB7_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB8 = fread(file.path(paste0(b, "/thickness_cortmapAB8_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB9 = fread(file.path(paste0(b, "/thickness_cortmapAB9_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB10 = fread(file.path(paste0(b, "/thickness_cortmapAB10_nTime2_", analysname, "_Yslopoys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        
        mapAB1_c = fread(file.path(paste0(b, "/thickness_cortmapAB1_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB2_c = fread(file.path(paste0(b, "/thickness_cortmapAB2_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB3_c = fread(file.path(paste0(b, "/thickness_cortmapAB3_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB4_c = fread(file.path(paste0(b, "/thickness_cortmapAB4_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB5_c = fread(file.path(paste0(b, "/thickness_cortmapAB5_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB6_c = fread(file.path(paste0(b, "/thickness_cortmapAB6_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB7_c = fread(file.path(paste0(b, "/thickness_cortmapAB7_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB8_c = fread(file.path(paste0(b, "/thickness_cortmapAB8_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB9_c = fread(file.path(paste0(b, "/thickness_cortmapAB9_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB10_c = fread(file.path(paste0(b, "/thickness_cortmapAB10_nTime2_", analysname, "_Yslopoyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        
      } else if (analysisnum == 6) {
        mapAB1 = fread(file.path(paste0(b, "//thickness_cortmapAB1_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB2 = fread(file.path(paste0(b, "//thickness_cortmapAB2_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB3 = fread(file.path(paste0(b, "/thickness_cortmapAB3_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB4 = fread(file.path(paste0(b, "/thickness_cortmapAB4_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB5 = fread(file.path(paste0(b, "/thickness_cortmapAB5_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB6 = fread(file.path(paste0(b, "/thickness_cortmapAB6_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB7 = fread(file.path(paste0(b, "/thickness_cortmapAB7_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB8 = fread(file.path(paste0(b, "/thickness_cortmapAB8_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB9 = fread(file.path(paste0(b, "/thickness_cortmapAB9_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB10 = fread(file.path(paste0(b, "/thickness_cortmapAB10_nTime2_", analysname, "_Yslointooys.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        
        mapAB1_c = fread(file.path(paste0(b, "/thickness_cortmapAB1_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB2_c = fread(file.path(paste0(b, "/thickness_cortmapAB2_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB3_c = fread(file.path(paste0(b, "/thickness_cortmapAB3_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB4_c = fread(file.path(paste0(b, "/thickness_cortmapAB4_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB5_c = fread(file.path(paste0(b, "/thickness_cortmapAB5_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB6_c = fread(file.path(paste0(b, "/thickness_cortmapAB6_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB7_c = fread(file.path(paste0(b, "/thickness_cortmapAB7_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB8_c = fread(file.path(paste0(b, "/thickness_cortmapAB8_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB9_c = fread(file.path(paste0(b, "/thickness_cortmapAB9_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB10_c = fread(file.path(paste0(b, "/thickness_cortmapAB10_nTime2_", analysname, "_Yslointooyscentcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        
        
      } else {
        mapAB1 = fread(file.path(paste0(b, "//thickness_cortmapAB1_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB2 = fread(file.path(paste0(b, "//thickness_cortmapAB2_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB3 = fread(file.path(paste0(b, "/thickness_cortmapAB3_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB4 = fread(file.path(paste0(b, "/thickness_cortmapAB4_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB5 = fread(file.path(paste0(b, "/thickness_cortmapAB5_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB6 = fread(file.path(paste0(b, "/thickness_cortmapAB6_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB7 = fread(file.path(paste0(b, "/thickness_cortmapAB7_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB8 = fread(file.path(paste0(b, "/thickness_cortmapAB8_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB9 = fread(file.path(paste0(b, "/thickness_cortmapAB9_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB10 = fread(file.path(paste0(b, "/thickness_cortmapAB10_nTime2_", analysname, ".csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        
        mapAB1_c = fread(file.path(paste0(b, "/thickness_cortmapAB1_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB2_c = fread(file.path(paste0(b, "/thickness_cortmapAB2_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB3_c = fread(file.path(paste0(b, "/thickness_cortmapAB3_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB4_c = fread(file.path(paste0(b, "/thickness_cortmapAB4_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB5_c = fread(file.path(paste0(b, "/thickness_cortmapAB5_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB6_c = fread(file.path(paste0(b, "/thickness_cortmapAB6_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB7_c = fread(file.path(paste0(b, "/thickness_cortmapAB7_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB8_c = fread(file.path(paste0(b, "/thickness_cortmapAB8_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB9_c = fread(file.path(paste0(b, "/thickness_cortmapAB9_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        mapAB10_c = fread(file.path(paste0(b, "/thickness_cortmapAB10_nTime2_", analysname, "_centcor.csv"))) %>% filter(!rois %in% c("pcaL", "pcaR"),!str_detect(rois, "Mean"))
        
      }
      
      mapsAB = data.frame(
        rois = mapAB1$rois,
        mapAB1 = mapAB1$statistic,
        mapAB2 = mapAB2$statistic,
        mapAB3 = mapAB3$statistic,
        mapAB4 = mapAB4$statistic,
        mapAB5 = mapAB5$statistic,
        mapAB6 = mapAB6$statistic,
        mapAB7 = mapAB7$statistic,
        mapAB8 = mapAB8$statistic,
        mapAB9 = mapAB9$statistic,
        mapAB10 = mapAB10$statistic
      )
      mapsAB_centcor = data.frame(
        rois = mapAB1$rois,
        mapAB1_c = mapAB1_c$statistic,
        mapAB2_c = mapAB2_c$statistic,
        mapAB3_c = mapAB3_c$statistic,
        mapAB4_c = mapAB4_c$statistic,
        mapAB5_c = mapAB5_c$statistic,
        mapAB6_c = mapAB6_c$statistic,
        mapAB7_c = mapAB7_c$statistic,
        mapAB8_c = mapAB8_c$statistic,
        mapAB9_c = mapAB9_c$statistic,
        mapAB10_c = mapAB10_c$statistic
      )
      
      
      # compute correlation matrix across loaded effect size maps (NB! non-independent)
      cor_mat = cor(mapsAB[,2:ncol(mapsAB)], use = "pairwise.complete.obs")
      cor_mat_centcor = cor(mapsAB_centcor[,2:ncol(mapsAB)], use = "pairwise.complete.obs")
      
      # plot it
      corrplot::corrplot(cor_mat, method = "circle", addCoef.col = "black", number.cex = 0.7)
      corrplot::corrplot(cor_mat_centcor, method = "circle", addCoef.col = "black", number.cex = 0.7)
      
      
      path_normmaps = file.path(base, "results/normMaps")
      readNormMaps = function(map_name, statmap = F) {
        df_map = data.table::fread(file.path(path_normmaps, map_name))
        df_map$rois = df_map$plotrois
        
        # stat to test against
        if (statmap == F) {
          df_map$istat = df_map$SUVR
        } else {
          df_map$istat = df_map$d
        }
        return(df_map)
      }
      
      
      # av SUVR across all ABpos timepoints of diagnosed cases in ADNI -----
      if (amyloidmap == 1) {
        abmap = readNormMaps("amyloid_map_ABdiagnosed_suvrmean.csv")
      }
      # group diff map in your sample -----
      if (amyloidmap == 2) {
        abmap = readNormMaps("amyloid_map_convVnonconv_suvrdiff_d.csv", statmap = T)
      }
      # group diff map between all ABpos timepoints of diagnosed cases in ADNI and nonconverters in your sample -----
      if (amyloidmap == 3) {
        abmap = readNormMaps("amyloid_map_ABdiagnosedVnonconv_suvrdiff_d.csv", statmap = T)
      }
      
      
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
      
      
      spinTest = function(xmap, ymap, color_by = xmap) {
        
        # testing
        # xmap = abmap
        # ymap = mapAB1
        
        
        wd = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/scripts/rotate_parcellation"
        
        
        # # Extract hemisphere coordinates
        source(file.path(wd, "R", "rotate.parcellation.R"))
        source(file.path(wd, "R", "perm.sphere.p.R"))
        coord.l <- read.table(file.path(wd, "sphere_left_coordinates.txt")) %>% as.matrix()
        coord.r <- read.table(file.path(wd, "sphere_right_coordinates.txt")) %>% as.matrix()
        
        
        # set up perm
        perm.id = rotate.parcellation(coord.l, coord.r, 1000)
        # #save(perm.id, file = "perm.id")
        
        
        #ensure ordered the same
        xmap$label = gsub("_thickness.aparcnative71", "", xmap$rois)
        ymap$label = gsub("_thickness.aparcnative71", "", ymap$rois)
        xmap$label %in% ymap$label
        identical(xmap$label, ymap$label)
        ymap <- ymap[match(xmap$label, ymap$label), ]
        identical(xmap$label, ymap$label)
        
        
        # NB! input stat got mapped on loading amyloid maps
        x = xmap$istat
        y = ymap$d
        cor.test(x, y)
        plot(x, y)
        
        
        # Spin test
        sigp = perm.sphere.p(x,y,perm.id,corr.type='pearson')
        
        
        # Spin result p-val
        sigp$p_value
        
        
        # Extract data
        empirical_correlation <- unname(unlist(sigp$empirical_correlation))
        permuted_correlations <- unname(unlist(sigp$permuted_correlations))
        (p_empirical_twosided <- mean(abs(permuted_correlations) >= abs(empirical_correlation)))
        (p_empirical_onesided <- mean(permuted_correlations >= empirical_correlation))
        
        
        (sum(abs(sigp$permuted_correlation) >= abs(empirical_correlation))) / (length(sigp$permuted_correlation))
        
        
        # Plot histogram of permuted correlations
        # hist(permuted_correlations, breaks = 30, col = "lightblue",
        #      main = "Permutation Test Results-Thickness",
        #      xlab = "Correlation", ylab = "Frequency")
        
        # # Add observed correlation
        # abline(v = empirical_correlation, col = "firebrick", lwd = 2)
        # text(x = empirical_correlation, y = max(table(cut(permuted_correlations, 30))),
        #      labels = paste("Observed:", round(empirical_correlation, 2)),
        #      pos = 4, col = "firebrick")
        
        
        region_names <- xmap$label
        
        
        # Create abbreviations for region names
        region_abbrev <- c(
          "lh_bankssts" = "L_bSTS",
          "lh_caudalanteriorcingulate" = "L_cACC",
          "lh_caudalmiddlefrontal" = "L_cMFG",
          "lh_cuneus" = "L_CUN",
          "lh_entorhinal" = "L_ENT",
          "lh_frontalpole" = "L_FP",
          "lh_fusiform" = "L_FUS",
          "lh_inferiorparietal" = "L_iPAR",
          "lh_inferiortemporal" = "L_iTEMP",
          "lh_insula" = "L_INS",
          "lh_isthmuscingulate" = "L_iCC",
          "lh_lateraloccipital" = "L_lOCC",
          "lh_lateralorbitofrontal" = "L_lOFC",
          "lh_lingual" = "L_LING",
          "lh_medialorbitofrontal" = "L_mOFC",
          "lh_middletemporal" = "L_mTEMP",
          "lh_paracentral" = "L_PARC",
          "lh_parahippocampal" = "L_PHC",
          "lh_parsopercularis" = "L_pOPER",
          "lh_parsorbitalis" = "L_pORB",
          "lh_parstriangularis" = "L_pTRI",
          "lh_pericalcarine" = "L_PCAL",
          "lh_postcentral" = "L_POSTC",
          "lh_posteriorcingulate" = "L_pCC",
          "lh_precentral" = "L_PREC",
          "lh_precuneus" = "L_PCUN",
          "lh_rostralanteriorcingulate" = "L_rACC",
          "lh_rostralmiddlefrontal" = "L_rMFG",
          "lh_superiorfrontal" = "L_sFG",
          "lh_superiorparietal" = "L_sPAR",
          "lh_superiortemporal" = "L_sTEMP",
          "lh_supramarginal" = "L_SMAR",
          "lh_temporalpole" = "L_TP",
          "lh_transversetemporal" = "L_tTEMP",
          "rh_bankssts" = "R_bSTS",
          "rh_caudalanteriorcingulate" = "R_cACC",
          "rh_caudalmiddlefrontal" = "R_cMFG",
          "rh_cuneus" = "R_CUN",
          "rh_entorhinal" = "R_ENT",
          "rh_frontalpole" = "R_FP",
          "rh_fusiform" = "R_FUS",
          "rh_inferiorparietal" = "R_iPAR",
          "rh_inferiortemporal" = "R_iTEMP",
          "rh_insula" = "R_INS",
          "rh_isthmuscingulate" = "R_iCC",
          "rh_lateraloccipital" = "R_lOCC",
          "rh_lateralorbitofrontal" = "R_lOFC",
          "rh_lingual" = "R_LING",
          "rh_medialorbitofrontal" = "R_mOFC",
          "rh_middletemporal" = "R_mTEMP",
          "rh_paracentral" = "R_PARC",
          "rh_parahippocampal" = "R_PHC",
          "rh_parsopercularis" = "R_pOPER",
          "rh_parsorbitalis" = "R_pORB",
          "rh_parstriangularis" = "R_pTRI",
          "rh_pericalcarine" = "R_PCAL",
          "rh_postcentral" = "R_POSTC",
          "rh_posteriorcingulate" = "R_pCC",
          "rh_precentral" = "R_PREC",
          "rh_precuneus" = "R_PCUN",
          "rh_rostralanteriorcingulate" = "R_rACC",
          "rh_rostralmiddlefrontal" = "R_rMFG",
          "rh_superiorfrontal" = "R_sFG",
          "rh_superiorparietal" = "R_sPAR",
          "rh_superiortemporal" = "R_sTEMP",
          "rh_supramarginal" = "R_SMAR",
          "rh_temporalpole" = "R_TP",
          "rh_transversetemporal" = "R_tTEMP"
        )
        
        if (color_by == "x") {
          sig_regions <- xmap$label[xmap$p.value < .05]
        } else {
          sig_regions <- ymap$label[ymap$p.value < .05]
        }
        
        df <- data.frame(
          x,
          y,
          Region = xmap$label,
          significant = xmap$label %in% sig_regions
        )
        
        
        (p_spatcor_labelled <- ggplot(df, aes(x = x, y = y)) +
            geom_point(data = subset(df, !significant), color = "black", size = 2) +
            geom_point(data = subset(df, significant), color = "cadetblue4", size = 3) +
            geom_smooth(method = "lm",
                        color = "firebrick",
                        se = FALSE) + 
            geom_text_repel(aes(label = region_abbrev[Region]),
                            size = 3,
                            max.overlaps = Inf) +
            
            theme_minimal() +
            theme(
              axis.title = element_text(size = 25),
              axis.text = element_text(size = 20),
              plot.title = element_text(size = 25),
              axis.line = element_line(color = "black"),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black")
            )
        )
        
        (p_spatcor_colored <- ggplot(df, aes(x = x, y = y)) +
            geom_point(data = subset(df, !significant), color = "black", size = 2) +
            geom_point(data = subset(df, significant), color = "cadetblue4", size = 3) +
            geom_smooth(method = "lm",
                        color = "firebrick",
                        se = FALSE) + 
            theme_minimal() +
            theme(
              axis.title = element_text(size = 25),
              axis.text = element_text(size = 20),
              plot.title = element_text(size = 25),
              axis.line = element_line(color = "black"),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black")
            )
        )
        (p_spatcor_bland <- ggplot(df, aes(x = x, y = y)) +
            geom_point(color = "cadetblue4", size = 3) +
            geom_smooth(method = "lm",
                        color = "firebrick",
                        se = FALSE) + 
            theme_minimal() +
            theme(
              axis.title = element_text(size = 25),
              axis.text = element_text(size = 20),
              plot.title = element_text(size = 25),
              axis.line = element_line(color = "black"),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black")
            )
        )
        
        
        (p_spatcor_newbland <- ggplot(df, aes(x = y, y = x)) +
            geom_point(color = "black", size = 3) +
            geom_smooth(method = "lm",
                        color = "#DC4D35",
                        # color = "#A26296",
                        
                        se = FALSE, size = 2) +
            theme_classic() + mytheme +
            theme(
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank()
            ) +
            labs(y = "SUVR", x = "Effect (d)")
          
        )
        print(p_spatcor_newbland)
        
        # ggsave(plot = p_spatcor_newbland,
        #        filename = paste0(base, "/paper2/figs_yearsBeforeAB/p_spatcor_newbland_AB8_abmap1_T2slopes.png"),
        #        width=9, #12,
        #        height=9, #11,
        #        units="cm", dpi=600)
        
        
        return = list(
          spinP = sigp$p_value,
          empirical_correlation = empirical_correlation,
          permuted_correlations = permuted_correlations,
          p_empirical_twosided = p_empirical_twosided,
          p_empirical_onesided = p_empirical_onesided,
          df = df,
          p_spatcor_labelled = p_spatcor_labelled,
          p_spatcor_colored = p_spatcor_colored,
          p_spatcor_bland = p_spatcor_bland
        )
      }
      
      
      
      # effect size maps spin-tested against ABmap
      print("-----------------1-----------------")
      spinabmap_1 = spinTest(abmap, mapAB1, color_by = "y")
      print("-----------------2-----------------")
      spinabmap_2 = spinTest(abmap, mapAB2, color_by = "y")
      print("-----------------3-----------------")
      spinabmap_3 = spinTest(abmap, mapAB3, color_by = "y")
      print("-----------------4-----------------")
      spinabmap_4 = spinTest(abmap, mapAB4, color_by = "y")
      print("-----------------5-----------------")
      spinabmap_5 = spinTest(abmap, mapAB5, color_by = "y")
      print("-----------------6-----------------")
      spinabmap_6 = spinTest(abmap, mapAB6, color_by = "y")
      print("-----------------7-----------------")
      spinabmap_7 = spinTest(abmap, mapAB7, color_by = "y")
      print("-----------------8-----------------")
      spinabmap_8 = spinTest(abmap, mapAB8, color_by = "y")
      print("-----------------9-----------------")
      spinabmap_9 = spinTest(abmap, mapAB9, color_by = "y")
      print("-----------------10-----------------")
      spinabmap_10 = spinTest(abmap, mapAB10, color_by = "y")
      
      
      make_spin_df <- function(prefix = "spinabmap_", yearsBeforeAB = 1:10) {
        n <- length(yearsBeforeAB)
        data.frame(
          yearsBeforeAB = yearsBeforeAB,
          empirical_correlation = sapply(seq_len(n), function(i) get(paste0(prefix, i))$empirical_correlation),
          spinP = sapply(seq_len(n), function(i) get(paste0(prefix, i))$spinP)
        ) %>%
          dplyr::mutate(spinsignificant = spinP < 0.05)
      }
      
      
      gg_spinab = make_spin_df(yearsBeforeAB = 1:10)
      df_means[[analysisnum]] = data.frame(mean = mean(gg_spinab$empirical_correlation))
      
      
      gg_spinab$spinaalpha = as.factor(ifelse(gg_spinab$spinsignificant == T, 0, 1))
      gg_spinab$spinaalpha2 = as.factor(ifelse(gg_spinab$spinsignificant == T, 1, 0))
      gg_spinab$spinaalpha2[gg_spinab$spinaalpha2 == 0] = NA
      
      
      spinlist[[analysisnum]] = gg_spinab %>% mutate(analysisname = plottext)
      pal = rev(met.brewer("Signac", n = 10, type = "continuous"))[c(1, 2, 3, 4, 7, 8)]
      
      
      if (analysisnum == 1) {
        
        (p_spinmain = ggplot(gg_spinab, aes(y = empirical_correlation, x = yearsBeforeAB)) +
           geom_line(col=pal[1], size=1, alpha = 0.5) +
           scale_alpha_discrete(c(NA, 1)) +
           geom_hline(yintercept = 0, col = "lightgrey", linetype = 2) +
           geom_point(data=gg_spinab[!is.na(gg_spinab$spinaalpha2), ], col="black", size=3, shape = 21, fill = pal[analysisnum], stroke = 0.75) +
           scale_x_continuous(breaks = seq(1, 10)) +
           theme_classic() + mytheme + 
           annotate("text", x = 7, y = 1, label = plottext, color = pal[analysisnum], size = 5, hjust = 0) +
           geom_segment(data = df_means[[analysisnum]],
                        aes(x = 10.2, xend = 11.3, y = mean, yend = mean),
                        color = pal[analysisnum],
                        inherit.aes = FALSE,
                        size = 0.5,
                        linetype = 2)
        )
        
        (p_spinmain_right = ggplot(gg_spinab, aes(y = empirical_correlation, x = yearsBeforeAB)) +
            geom_line(col=pal[1], size=1, alpha = 0.5) +
            scale_alpha_discrete(c(NA, 1)) +
            scale_y_continuous(position = "right", limits = c(-0.3, 1), breaks = seq(-0.2, 1.0, by=0.2)) +
            geom_hline(yintercept = 0, col = "lightgrey", linetype = 2) +
            geom_point(data=gg_spinab[!is.na(gg_spinab$spinaalpha2), ], col="black", size=3, shape = 21, fill = pal[analysisnum], stroke = 0.75) +
            scale_x_continuous(breaks = seq(1, 10)) +
            theme_classic() + mytheme + 
            annotate("text", x = 7, y = 1, label = plottext, color = pal[analysisnum], size = 5, hjust = 0) +
            geom_segment(data = df_means[[analysisnum]],
                         aes(x = -0.1, xend = 0.8, y = mean, yend = mean),
                         color = pal[analysisnum],
                         inherit.aes = FALSE,
                         size = 0.5,
                         linetype = 2)
        )
        
      } else {
        
        (p_spinmain = p_spinmain + 
           geom_line(data = gg_spinab, aes(y = empirical_correlation, x = yearsBeforeAB), col=pal[analysisnum], size=1,  alpha = 0.5) +
           scale_alpha_discrete(c(0, 1)) +
           geom_point(data=gg_spinab[!is.na(gg_spinab$spinaalpha2), ], col="black", size=3, shape = 21, fill = pal[analysisnum], stroke = 0.75) +
           annotate("text", x = 7, y = 1 - (analysisnum - 1) * 0.05, label = plottext, color = pal[analysisnum], size = 5, hjust = 0) +
           geom_segment(data = df_means[[analysisnum]],
                        aes(x = 10.2, xend = 11.3, y = mean, yend = mean),
                        color = pal[analysisnum],
                        inherit.aes = FALSE,
                        size = 0.5,
                        linetype = 2))
        
        (p_spinmain_right = p_spinmain + 
            geom_line(data = gg_spinab, aes(y = empirical_correlation, x = yearsBeforeAB), col=pal[analysisnum], size=1,  alpha = 0.5) +
            scale_alpha_discrete(c(0, 1)) +
            scale_y_continuous(position = "right", limits = c(-0.3, 1), breaks = seq(-0.2, 1.0, by=0.2)) +
            geom_point(data=gg_spinab[!is.na(gg_spinab$spinaalpha2), ], col="black", size=3, shape = 21, fill = pal[analysisnum], stroke = 0.75) +
            annotate("text", x = 7, y = 1 - (analysisnum - 1) * 0.05, label = plottext, color = pal[analysisnum], size = 5, hjust = 0) +
            geom_segment(data = df_means[[analysisnum]],
                         aes(x = -0.1, xend = 0.8, y = mean, yend = mean),
                         color = pal[analysisnum],
                         inherit.aes = FALSE,
                         size = 0.5,
                         linetype = 2))
        
        
      }
      
    }
  }
  
  mainlist[[amyloidmap]] = bind_rows(spinlist, .id = "analysismap")
  p_mainlist[[amyloidmap]] = p_spinmain
  p_mainlistright[[amyloidmap]] = p_spinmain_right
  df_means_mainlist[[amyloidmap]] = bind_rows(df_means) %>% mutate(analysisname = plottext)
  
  
  # save('mainlist','p_mainlist','p_mainlistright', 'df_means_mainlist',
  #      file = paste0(base, "/spatCor_abmap", amyloidmap, ".Rda")
  # )
  
  
  # (p_spinmain = p_spinmain +
  #   scale_y_continuous(position = "right", limits = c(-0.3, 1), breaks = seq(-0.2, 1.0, by=0.2))
  # )
  
  
  p_mainlist[[amyloidmap]] = p_mainlist[[amyloidmap]] + scale_y_continuous(position = "left", limits = c(-0.3, 1), breaks = seq(-0.2, 1.0, by=0.2))
  
  
  ggsave(plot = p_mainlist[[amyloidmap]] + labs(y = "r") + 
           theme(axis.text = element_text(size = 12),
                 axis.title.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12)
           ),
         filename = paste0(base, "/paper2/figs_yearsBeforeAB/p_spatcorABeffects_slope_abmap", amyloidmap, "D_left.png"),
         width=10, #12,
         height=8, #11,
         units="cm", dpi=600)
  
  
  # ggsave(plot = p_mainlistright[[amyloidmap]] + labs(y = "r") + theme(axis.text = element_text(size = 12),
  #                                                       axis.title.x = element_text(size = 12),
  #                                                       axis.title.y = element_text(size = 12)),
  #        filename = paste0(base, "/paper2/figs_yearsBeforeAB/p_spatcorABeffects_slope_abmap", amyloidmap, "_right.png"),
  #        width=10, #12, 
  #        height=8, #11, 
  #        units="cm", dpi=600)
  
}
