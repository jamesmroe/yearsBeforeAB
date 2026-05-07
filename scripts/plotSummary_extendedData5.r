#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Reproduce Extended Data Fig. 5 (Summary of the number of reported significant effects and their directions across main and supplementary analyses (raw differences uncorrected for Aβ levels))
#          Script is fully executable
#========================================================================================#

# summary effect plots across all analyses
rm(list=ls())
library(stringr); library(tidyverse); library(magrittr); library(here); library(patchwork); library(ggpubr); library(data.table); library(dplyr); library(ggplot2)
base = here()
df_powered_slope = fread(file.path(base, "sourceData/extendedData5/df_powered_slope.csv"))
df_powered_level = fread(file.path(base, "sourceData/extendedData5/df_powered_level.csv"))
df_reduced_slope = fread(file.path(base, "sourceData/extendedData5/df_reduced_slope.csv"))
df_reduced_level = fread(file.path(base, "sourceData/extendedData5/df_reduced_level.csv"))
df_predicted_slope = fread(file.path(base, "sourceData/extendedData5/df_predicted_slope.csv"))
df_predicted_level = fread(file.path(base, "sourceData/extendedData5/df_predicted_level.csv"))


# ROIS from anterior to posterior
dk_ap = c(
  # Global means
  "lh_MeanThickness","rh_MeanThickness",
  
  # Frontal (most anterior)
  "lh_frontalpole","rh_frontalpole",
  "lh_parsorbitalis","rh_parsorbitalis",
  "lh_parstriangularis","rh_parstriangularis",
  "lh_parsopercularis","rh_parsopercularis",
  "lh_lateralorbitofrontal","rh_lateralorbitofrontal",
  "lh_medialorbitofrontal","rh_medialorbitofrontal",
  "lh_rostralmiddlefrontal","rh_rostralmiddlefrontal",
  "lh_caudalmiddlefrontal","rh_caudalmiddlefrontal",
  "lh_superiorfrontal","rh_superiorfrontal",
  "lh_precentral","rh_precentral",
  "lh_paracentral","rh_paracentral",
  
  # Cingulate / insula (anterior → posterior)
  "lh_rostralanteriorcingulate","rh_rostralanteriorcingulate",
  "lh_caudalanteriorcingulate","rh_caudalanteriorcingulate",
  "lh_insula","rh_insula",
  
  # Temporal (anterior → posterior)
  "lh_temporalpole","rh_temporalpole",
  "lh_entorhinal","rh_entorhinal",
  "lh_parahippocampal","rh_parahippocampal",
  "lh_fusiform","rh_fusiform",
  "lh_superiortemporal","rh_superiortemporal",
  "lh_middletemporal","rh_middletemporal",
  "lh_inferiortemporal","rh_inferiortemporal",
  "lh_transversetemporal","rh_transversetemporal",
  "lh_bankssts","rh_bankssts",
  
  # Parietal / medial posterior
  "lh_postcentral","rh_postcentral",
  "lh_supramarginal","rh_supramarginal",
  "lh_inferiorparietal","rh_inferiorparietal",
  "lh_superiorparietal","rh_superiorparietal",
  "lh_precuneus","rh_precuneus",
  "lh_posteriorcingulate","rh_posteriorcingulate",
  "lh_isthmuscingulate","rh_isthmuscingulate",
  
  # Occipital / posterior-most
  "lh_lingual","rh_lingual",
  "lh_cuneus","rh_cuneus",
  "lh_pericalcarine","rh_pericalcarine",
  "lh_lateraloccipital","rh_lateraloccipital"
)
length(dk_ap)



# example FDR corrected positive effects
df_powered_slope$rois[df_powered_slope$sigdirection=="positive_fdrsig" & !is.na(df_powered_slope$sigdirection)]
df_powered_slope$p.value[df_powered_slope$sigdirection=="positive_fdrsig" & !is.na(df_powered_slope$sigdirection)]


df_powered_slopesum = df_powered_slope %>%
  filter(!is.na(sigdirection)) %>%
  group_by(roi, sigdirection) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(roi = factor(roi, levels = dk_ap)) %>%
  arrange(roi)

df_powered_levelsum = df_powered_level %>%
  filter(!is.na(sigdirection)) %>%
  group_by(roi, sigdirection) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(roi = factor(roi, levels = dk_ap)) %>%
  arrange(roi)

df_predicted_slopesum = df_predicted_slope %>%
  filter(!is.na(sigdirection)) %>%
  group_by(roi, sigdirection) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(roi = factor(roi, levels = dk_ap)) %>%
  arrange(roi)

df_predicted_levelsum = df_predicted_level %>%
  filter(!is.na(sigdirection)) %>%
  group_by(roi, sigdirection) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(roi = factor(roi, levels = dk_ap)) %>%
  arrange(roi)

df_reduced_slopesum = df_reduced_slope %>%
  filter(!is.na(sigdirection)) %>%
  group_by(roi, sigdirection) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(roi = factor(roi, levels = dk_ap)) %>%
  arrange(roi)

df_reduced_levelsum = df_reduced_level %>%
  filter(!is.na(sigdirection)) %>%
  group_by(roi, sigdirection) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(roi = factor(roi, levels = dk_ap)) %>%
  arrange(roi)


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




makeBarPlot = function(df) { 
  # Prepare data for plotting
  df_plot = df %>%
    mutate(
      # ensure anatomical ordering
      roi = factor(roi, levels = dk_ap),
      
      # determine effect direction
      direction = case_when(
        grepl("^positive", sigdirection) ~ "positive",
        grepl("^negative", sigdirection) ~ "negative"
      ),
      
      # assign signed counts for diverging bars
      n_signed = ifelse(direction == "positive", n, -n)
    )
  
  missing_rois = setdiff(dk_ap, as.character(df_plot$roi))
  
  df_nosig = data.frame(
    roi = factor(missing_rois, levels = dk_ap),
    sigdirection = "nosig",
    n = 0,
    direction = "nosig",
    n_signed = 0,
    stringsAsFactors = FALSE
  )
  df_plot = rbind(df_plot, df_nosig)
  df_plot$roi = factor(df_plot$roi, levels = dk_ap)
  df_plot = df_plot[order(df_plot$roi), ]
  
  length(unique(df_plot$roi))
  
  # symmetric y-axis limits
  ylim = max(abs(df_plot$n_signed), na.rm = TRUE)
  
  length(unique(df_plot$roi))
  
  # Diverging stacked bar plot
  (p_barcount = ggplot(df_plot, aes(x = roi, y = n_signed, fill = sigdirection)) +
      geom_col(width = 0.8) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      scale_fill_manual(
        values = c(
          "positive_fdrsig"   = "#3a0f0b",
          "positive_sigonly" = "#e7b1aa",
          "negative_sigonly" = "#aad6e7",
          "negative_fdrsig"  = "#0f2a35"
        ),
        breaks = c(
          "positive_fdrsig",
          "positive_sigonly",
          "negative_sigonly",
          "negative_fdrsig"
        )
        
      ) +
      scale_y_continuous(
        name = "Number of significant effects",
        # limits = c(-15, ylim),
        # expand = expansion(mult = c(0.05, 0.05))
      ) +
      theme_classic() + mytheme +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0))
      )
  )
  
  return(list(p_barcount = p_barcount))
  
}

p_slope = makeBarPlot(df_powered_slopesum)
p_level = makeBarPlot(df_powered_levelsum)

p_reduced_slope = makeBarPlot(df_reduced_slopesum)
p_reduced_level = makeBarPlot(df_reduced_levelsum)

p_predicted_slope = makeBarPlot(df_predicted_slopesum)
p_predicted_level = makeBarPlot(df_predicted_levelsum)



# slope analyses
p_slope$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(-6, 40),
    breaks = seq(0, 40, by = 10)
  )

p_reduced_slope$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(-3, 46),
    breaks = seq(0, 46, by = 10)
  )

p_predicted_slope$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(-2, 10),
    breaks = seq(0, 10, by = 2)
  )


# level analyses
p_level$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(-6, 40),
    breaks = seq(0, 40, by = 10)
  )

p_reduced_level$p_barcount + theme(legend.position = "none")  +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(-3, 46),
    breaks = seq(0, 40, by = 10)
  )

p_predicted_level$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(-2, 10),
    breaks = seq(0, 10, by = 2)
  )


saveplots=0
if (saveplots) {
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/pX_sigeffects_slopepowered.pdf",
         plot = p_slope$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
           scale_y_continuous(
             limits = c(-6, 40),
             breaks = seq(0, 40, by = 10)
           ),
         width = 36,
         height = 18,
         units = "cm",
         device = cairo_pdf
  )
  
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/pX_sigeffects_slopereducedextranodup.pdf",
         plot = p_reduced_slope$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
           scale_y_continuous(
             limits = c(-3, 46),
             breaks = seq(0, 46, by = 10)
           ),
         width = 36,
         height = 18,
         units = "cm",
         device = cairo_pdf
  )
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/pX_sigeffects_slopepredicted.pdf",
         plot = p_predicted_slope$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
           scale_y_continuous(
             limits = c(-2, 10),
             breaks = seq(0, 10, by = 2)
           ),
         width = 36,
         height = 18,
         units = "cm",
         device = cairo_pdf
  )
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/pX_sigeffects_levelpowered.pdf",
         plot = p_level$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
           scale_y_continuous(
             limits = c(-6, 40),
             breaks = seq(0, 40, by = 10)
           ),
         width = 36,
         height = 18,
         units = "cm",
         device = cairo_pdf
  )
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/pX_sigeffects_levelreducedextranodup.pdf",
         plot = p_reduced_level$p_barcount + theme(legend.position = "none")  +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
           scale_y_continuous(
             limits = c(-3, 46),
             breaks = seq(0, 40, by = 10)
           ),
         width = 36,
         height = 18,
         units = "cm",
         device = cairo_pdf
  )
  
  ggsave(filename = "/cluster/projects/p274/projects/p040-ad_change/Berkeley/paper2/figs_yearsBeforeAB/pX_sigeffects_levelpredicted.pdf",
         plot = p_predicted_level$p_barcount + theme(legend.position = "none") +scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
           scale_y_continuous(
             limits = c(-2, 10),
             breaks = seq(0, 10, by = 2)
           ),
         width = 36,
         height = 18,
         units = "cm",
         device = cairo_pdf
  )
  
}

