#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Reproduce Extended Data Fig. 4 (Sensitivity analyses correcting for additional covariates)
#          Script is fully executable
#========================================================================================#

rm(list=ls())
library(stringr); library(tidyverse); library(magrittr); library(here); library(patchwork); library(ggpubr)


#--------------------------
base = here()
resdir0 = here(base, "sourceData/extendedData4/moreCovariates")
modelfiles01 <- list.files(path = resdir0, pattern = "thickness_cortmapAB.*_MAINrevALLCOV\\.csv", full.names = TRUE) # main models with reduced sample to match models ran with more covariates
modelfiles02 <- list.files(path = resdir0, pattern = "thickness_cortmapAB.*_revALLCOV\\.csv", full.names = TRUE) # models with more covariates


print("loading all model outputs")
df_models01 = read_delim(modelfiles01, id = "src", delim = " ") |>
  mutate(src = basename(src)) 

df_models02 = read_delim(modelfiles02, id = "src", delim = " ") |>
  mutate(src = basename(src)) 


# check terms in models
unique(df_models01$term)
unique(df_models02$term)


extractYearsBeforeAB = function(sourcefile) {
  tmp = str_sub(sourcefile, 20, 21)
  tmp = as.numeric(gsub("_", "", tmp))  
  return(tmp)
}
df_models01$yearsBeforeAB = extractYearsBeforeAB(df_models01$src)
df_models02$yearsBeforeAB = extractYearsBeforeAB(df_models02$src)


# original effect in model with same N
df_MAINeffectorig = df_models01 %>% filter(term == "convgroup")


# new effect in more adjusted model
df_MAINeffect = df_models02 %>% filter(term == "convgroup")
df_APOEeffect = df_models02 %>% filter(term == "carrier")
df_EDUeffect = df_models02 %>% filter(term == "scale(edu)")
df_WMHeffect = df_models02 %>% filter(term == "scale(wmhSlope)")


cor(df_MAINeffectorig$statistic, df_MAINeffect$statistic)
# plot(df_MAINeffectorig$statistic, df_MAINeffect$statistic)


df_compare = data.frame(df_MAINeffectorig,
                        df_MAINeffect %>% select(estimate, statistic, p.value, d) %>% rename_with(~ paste0("mod02_", .))
)


mytheme = theme(plot.background = element_rect(fill = "white"),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                title = element_text(size=12),
                text = element_text(color = "black", size = 12, family="Nimbus Sans Narrow"),
                plot.title = element_text(hjust = 0.5),
                axis.title.y = element_text(color = "black", size = 12, vjust =-1),
                axis.title.x = element_text(color = "black", size = 12, vjust = -2),
                axis.text = element_text(color = "black", size = 12),
                legend.key.size = unit(1,"cm"))


pal = MetBrewer::MetPalettes$Signac[[1]]


plotDots = function(imap, x, y, title, xtitle, ytitle) {
  
  correl = round(cor(imap[x], imap[y])[[1]], 3)
  
  (
    ggcor = ggplot(imap, aes_string(x = x, y = y)) +
      geom_point(aes(col = as.factor(yearsBeforeAB))) +
      scale_color_manual(values = pal[1:10]) +
      scale_fill_manual(values = pal[1:10]) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "blue") +
      
      geom_smooth(
        method = "lm",
        col = "black",
        se = F
      ) +
      theme_classic() + mytheme +
      coord_fixed() +
      theme(legend.position = "right",
            axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
            axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0))
      ) +
      labs(title = title,
           color = "Years Before obs. Aβ+",
           x = xtitle,
           y = ytitle)
  )
  
  return(list(ggcor = ggcor,
              correl = correl
  )
  )
}


# change in signifcance
# num sig in initial model
sum(df_compare$p.value<.05)
sum(df_compare$mod02_p.value<.05)

#125 of 150 associations remained significant
sum(df_compare$mod02_p.value< .05 & df_compare$p.value<.05)

#how many became sig with new model - 6
sum(df_compare$mod02_p.value< .05 & df_compare$p.value>.05)

# 25 lost signifcance
sum(df_compare$mod02_p.value> .05 & df_compare$p.value<.05)

# which are these
df_compare$roiyears = paste0(df_compare$rois, "_", df_compare$yearsBeforeAB)
df_compare$roi = gsub("thickness.aparcnative71", "", df_compare$rois)
df_compare$roi = stringi::stri_reverse(
  substring(stringi::stri_reverse(df_compare$roi), 2)
)
df_compare$hemi = ifelse(substr(df_compare$roi, 1, 1) == "l", "L", "R")
df_compare$roi = paste0(df_compare$roi, " ", df_compare$hemi)
df_compare$roi = gsub("lh_", "", df_compare$roi)
df_compare$roi = gsub("rh_", "", df_compare$roi)

tmplost = df_compare[which(df_compare$mod02_p.value> .05 & df_compare$p.value<.05),]
tmpgain = df_compare[which(df_compare$mod02_p.value< .05 & df_compare$p.value>.05),]

table(tmplost$roi)
table(tmpgain$roi)

tmplosttab = as.data.frame(table(tmplost$roi))
tmpgaintab = as.data.frame(table(tmpgain$roi))


# plot bar chart
(p_lost = ggplot(tmplosttab, aes(x = Freq, y = Var1)) +
    geom_bar(stat = "identity", fill = "#1E4E79") +
    coord_cartesian(xlim = c(1, 10)) +
    scale_x_continuous(breaks = 1:10, labels = 1:10) +
    
    theme_classic() + mytheme +
    # coord_fixed() +
    theme(legend.position = "right",
          axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
          axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0))
    ) +
    labs(title = paste0("Lost significance (", nrow(tmplost), ")"),
         color = "Years Before obs. Aβ+",
         x = "Freq",
         y = NULL) + theme(axis.line.y = element_blank(), axis.ticks = element_blank())
)


(p_gain = ggplot(tmpgaintab, aes(x = Freq, y = Var1)) +
    geom_bar(stat = "identity", fill = "#A52A2A") +
    
    theme_classic() + mytheme +
    coord_cartesian(xlim = c(1, 10)) +
    scale_x_continuous(breaks = 1:10, labels = 1:10) +
    
    theme(legend.position = "right",
          axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
          axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0))
    ) +
    labs(title = paste0("Gained significance (", nrow(tmpgain), ")"),
         color = "Years Before obs. Aβ+",
         x = "Freq",
         y = NULL) + theme(axis.line.y = element_blank(), axis.ticks = element_blank())
)

p_lost = p_lost + theme(axis.title.x = element_blank(), axis.text = element_text(size = 14), axis.title = element_text(size = 20), plot.title = element_text(size = 20))
p_gain = p_gain + theme(axis.title.x = element_blank(), axis.text = element_text(size = 14), axis.title = element_text(size = 20), plot.title = element_text(size = 20))

combined_lostgain <-  p_lost / p_gain
combined_lostgain = cowplot::plot_grid(p_lost, p_gain, rel_heights = c(1, 0.5), nrow = 2, align = T)


# ggsave(plot = combined_lostgain,
#        filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_main_v_review_lostgain_largertext_corectfullscale.pdf",
#        width=13, height=20, units="cm", dpi=600, device=cairo_pdf)


# plot correlation of effects sizes
p_correl1 = plotDots(df_compare, x = "mod02_statistic", y = "statistic", title = "T statistic", xtitle = "+ APOE carrier + edu + WMH slope", ytitle = "main model")
p_correl2 = plotDots(df_compare, x = "mod02_d", y = "d", title = "Cohens d", xtitle = "+ APOE carrier + edu + WMH slope", ytitle = "main model")
p_correl1$ggcor = p_correl1$ggcor + annotate("text", x = -2, y = 3, label = paste("R =", p_correl1$correl), color = "black", hjust = 0)
p_correl2$ggcor = p_correl2$ggcor + annotate("text", x = -0.2, y = .35, label = paste("R =", p_correl2$correl), color = "black", hjust = 0)  


# plot(-log10(df_compare$mod02_p.value), -log10(df_compare$p.value))
p_both = ggarrange(
  p_correl1$ggcor + theme(axis.title.x = element_text(size = 15),
                          axis.title.y = element_text(size = 15)),
  p_correl2$ggcor + theme(axis.title.x = element_text(size = 15),
                          axis.title.y = element_text(size = 15)),
  common.legend = 2,
  align = "hv"
)
p_both

(p_logp = 
    ggplot(df_compare, aes(x = -log10(mod02_p.value), y = -log10(p.value))) +
    geom_point(aes(col = as.factor(yearsBeforeAB))) +
    scale_color_manual(values = pal[1:10]) +
    scale_fill_manual(values = pal[1:10]) +
    geom_smooth(
      method = "lm",
      col = "black",
      se = F
    ) +
    theme_classic() + mytheme +
    coord_fixed() +
    theme(legend.position = "right",
          axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
          axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0))
    ) +
    labs(title = "-log10(p)",
         color = "Years Before obs. Aβ+",
         x = "+ APOE carrier + edu + WMH slope",
         y = "main model") + 
    geom_hline(yintercept = 1.3, linetype = 2, col = "lightgrey") + geom_vline(xintercept = 1.3, linetype = 2, col = "lightgrey")
)


p_bothp = ggarrange(
  p_logp + theme(axis.title.x = element_text(size = 15),
                 axis.title.y = element_text(size = 15)),
  p_logp + theme(axis.title.x = element_text(size = 15),
                 axis.title.y = element_text(size = 15)),
  common.legend = 2,
  align = "hv"
)


savefigs = 0
if (savefigs) {
  ggsave(plot = p_both,
  filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/p_main_v_review.pdf",
  width=20, height=20, units="cm", dpi=600, device=cairo_pdf)
}


# apply FDR correction per effect map
correctEffectMap = function(imap) {
  
  omap = imap
  omap$FDRind = 0
  omap$pInd = 0
  omap$pInd[omap$p.value<.05] = 1
  omap$alphaP = 0.3
  omap$alphaP[omap$pInd == 1] = 1
  
  omap$ci = omap$std.error*1.96
  omap$ci_lwr = omap$estimate-omap$ci
  omap$ci_upr = omap$estimate+omap$ci
  
  
  # FDR correct per map (as in main analysis)
  for (ii in 1:10) {
    if (ii==1) {
      OMAP=c()
    }
    tmp = omap[omap$yearsBeforeAB == ii,]
    tmpbh = sgof::BH(c(tmp$p.value))
    fdthresh = max(tmpbh$data[tmpbh$Adjusted.pvalues < .05])
    tmp$FDRind = ifelse(tmp$p.value <= fdthresh, 1, 0)
    OMAP %<>% rbind(., tmp)
  }
  
  return(OMAP)
}


df_APOEeffect = correctEffectMap(df_APOEeffect)
df_EDUeffect = correctEffectMap(df_EDUeffect)
df_MAINeffect = correctEffectMap(df_MAINeffect)
df_WMHeffect = correctEffectMap(df_WMHeffect)


df_APOEeffect$p.value[df_APOEeffect$pInd == 1]
df_APOEeffect[df_APOEeffect$pInd == 1,]



common_max = max(
  c(
    df_APOEeffect$ci_upr,
    df_EDUeffect$ci_upr,
    df_MAINeffect$ci_upr,
    df_WMHeffect$ci_upr
  )
) + 0.1

common_min = min(
  c(
    df_APOEeffect$ci_lwr,
    df_EDUeffect$ci_lwr,
    df_MAINeffect$ci_lwr,
    df_WMHeffect$ci_lwr
  )
) - 0.1


#original
df_MAINeffectorig$yearsBeforeAB = factor(df_MAINeffectorig$yearsBeforeAB, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
df_MAINeffectorig %<>% arrange(yearsBeforeAB)
df_MAINeffectorig = correctEffectMap(df_MAINeffectorig)



effectPlot = function(dat, effectString, yString, alphaFac = "fdr") {
  
  # gwas style p-value plot
  (punivar_gwasp = dat %>%
     
     ggplot(.) +
     scale_x_discrete() +
     scale_alpha_discrete(range = c(0, 1)) +
     scale_x_discrete(expand = c(0.05, 0.05)) + #stopo cutting plot at top
     scale_color_manual(values = pal) +
     geom_point(aes(x=seq(1,nrow(dat)), y=-log10(p.value), col = factor(yearsBeforeAB)),alpha = dat$alphaP) +
     
     geom_hline(yintercept = 1.3) +
     coord_cartesian(ylim = c(0, 5)) +
     
     labs(title = effectString, x = NULL, y="-log10(p)") +
     theme_classic() + mytheme +
     
     theme(text = element_text(color = "black", size = 18, family="Nimbus Sans Narrow"),
           plot.title = element_text(hjust = 0.5),
           axis.ticks = element_line(),
           axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
           axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,20,0)),
           axis.text = element_text(color = "black", size = 18)) +
     theme(axis.text.x = element_text(size = 8, angle=90, hjust = 1),
           # axis.line.x.top = element_line())
     ) +
     
     labs(y = "-log10(p)")
  )
  
  if (alphaFac == "fdr") {
    punivar_gwasp = punivar_gwasp + 
      geom_point(aes(x=seq(1,nrow(dat)), y=-log10(p.value),alpha=factor(FDRind)),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75)
  }
  if (alphaFac == "p") {
    punivar_gwasp = punivar_gwasp + 
      geom_point(aes(x=seq(1,nrow(dat)), y=-log10(p.value),alpha=factor(pInd)),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75)
  }
  
  
  # gwas style effect size plot
  (punivar_gwasest = dat %>%
      
      ggplot(.) +
      scale_x_discrete(position = "top") +
      geom_hline(yintercept = 0, linetype=2, size=1) +
      scale_alpha_discrete(range = c(0, 1)) +
      scale_x_discrete(expand = c(0.05, 0.05)) + #stopo cutting plot at top
      geom_pointrange(aes(alpha=factor(pInd), x=seq(1,nrow(dat)), y=estimate, ymin=estimate-ci,ymax=estimate+ci, col = factor(yearsBeforeAB)),position=position_dodge(width=0.3)) +
      
      scale_color_manual(values = pal) +
      coord_cartesian(ylim = c(common_min, common_max)) +
      labs(title = effectString, x = NULL, y=expression(beta)) +
      
      theme_classic() + mytheme +
      
      theme(text = element_text(color = "black", size = 18, family="Nimbus Sans Narrow"),
            plot.title = element_text(hjust = 0.5),
            axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
            axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,20,0)),
            axis.text = element_text(color = "black", size = 18)) +
      theme(axis.text.x = element_text(size = 8, angle=90, hjust = 1),
      )
    
    
  )
  
  if (alphaFac == "fdr") {
    punivar_gwasest = punivar_gwasest + 
      geom_point(aes(x=seq(1,nrow(dat)), y=estimate,alpha=factor(FDRind)),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75)
  }
  if (alphaFac == "p") {
    punivar_gwasest = punivar_gwasest + 
      geom_point(aes(x=seq(1,nrow(dat)), y=estimate,alpha=factor(pInd)),col = "black", position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75)
  }
  
  
  return(list(p_partial = punivar_gwasp + theme(legend.position = "none"),
              p_beta = punivar_gwasest + theme(legend.position = "none")))
}


p_APOE = effectPlot(df_APOEeffect, "APOE carriership", "-log10(p)")
p_EDU = effectPlot(df_EDUeffect, "Education", "-log10(p)", alphaFac = "fdr")
p_WMH = effectPlot(df_WMHeffect, "White Matter Hypointensity change", "-log10(p)")
p_MAINPARTIAL = effectPlot(df_MAINeffect, "Converters v Aβ- (re-adjusted)", "-log10(p)", alphaFac = "none")
p_MAINORIG = effectPlot(df_MAINeffectorig, "Converters v Aβ- (main)", "-log10(p)", alphaFac = "none")



nSigOrig = sum(df_MAINeffectorig$p.value < .05)
nSigNew = sum(df_MAINeffect$p.value < .05)


ggarrange(
  p_MAINORIG$p_partial + theme(axis.title.x = element_text(size = 15),
                               axis.title.y = element_text(size = 15)) + annotate("text", x = 10, y = 5, label = paste("N sig =", nSigOrig), color = "black", hjust = 0),
  p_MAINPARTIAL$p_partial + theme(axis.title.x = element_text(size = 15),
                                  axis.title.y = element_text(size = 15)) + annotate("text", x = 10, y = 5, label = paste("N sig =", nSigNew), color = "black", hjust = 0)
)


(p_effects = ggarrange(
  p_APOE$p_partial + theme(axis.title.x = element_text(size = 15),
                           axis.title.y = element_text(size = 15)),
  p_EDU$p_partial + theme(axis.title.x = element_text(size = 15),
                          axis.title.y = element_text(size = 15)),
  p_WMH$p_partial + theme(axis.title.x = element_text(size = 15),
                          axis.title.y = element_text(size = 15)),
  common.legend = 2,
  align = "hv",
  ncol = 5
)
)


combinePlots = function(p_obj){
  punivar_gwasp = p_obj$p_partial
  punivar_gwasest = p_obj$p_beta
  punivar_gwasest + scale_x_discrete(position = "top") + theme(axis.line.x = element_line())
  
  punivar1 = punivar_gwasest + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
  punivar2 = punivar_gwasp + theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
  punivar1 = punivar1 + theme(legend.position = "none") + theme(axis.line = element_line()) 
  punivar2 = punivar2 + theme(legend.position = "none") + theme(axis.line = element_line()) 
  # library(patchwork)
  combined_plot <- punivar2 / punivar1
  
  return(combined_plot)
}


p_MAINORIG_combined = combinePlots(p_MAINORIG)
p_MAINPARTIAL_combined = combinePlots(p_MAINPARTIAL)
p_APOE_combined = combinePlots(p_APOE)
p_WMH_combined = combinePlots(p_WMH)
p_EDU_combined = combinePlots(p_EDU)


# plots in Extended Data 4
p_MAINORIG_combined
p_MAINPARTIAL_combined
p_APOE_combined
p_EDU_combined
p_WMH_combined

