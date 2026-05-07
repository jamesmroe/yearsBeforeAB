# write all source data to file

rm(list=ls())

# Load packages ----
library("here")
library("magrittr")
library("tidyverse")
library("sgof")
library("openxlsx")


# results dir ----
resdir = here("sourceData")

loadResults = function(files) {
  df_models = lapply(files, function(f) {
    read_table(f) |> mutate(src = basename(f))
  }) |> bind_rows()
  
  # Clean and derive variables ----
  df_models$term.1 = NULL
  df_models = df_models |>
    mutate(
      yearsBeforeAB = str_split_fixed(src, "_", n = 3)[, 2] |>
        str_extract("(?<=AB)\\d+") |>
        as.integer(),
      centcor  = grepl("centcor", src),
      Yintonly = grepl("Yintonly", src),
      type = case_when(
        centcor & Yintonly  ~ "level_centcor",
        centcor & !Yintonly ~ "slope_centcor",
        !centcor & Yintonly ~ "level",
        TRUE                ~ "slope"
      ),
      ABtype = case_when(
        grepl("predAB1", src) ~ "predicted",
        grepl("predAB0", src) ~ "observed",
        TRUE                  ~ NA_character_
      )
    ) |>
    arrange(
      yearsBeforeAB,
      factor(type, levels = c("slope", "slope_centcor", "level", "level_centcor"))
    )
  return(df_models)
}


# fig 1 --------------------------------------------
files = list.files(path = paste0(resdir, "/fig1/cortmaps"), pattern = "thickness_cortmap*", full.names = TRUE)
df_fig1 = loadResults(files)


# fig 2 --------------------------------------------
files = list.files(path = paste0(resdir, "/fig2/cortmaps"), pattern = "thickness_cortmap*", full.names = TRUE)
df_fig2 = loadResults(files)


# fig 3  --------------------------------------------
# Load and label AB maps ----
map1 = read_table(file.path(resdir, "fig3/normMaps/amyloid_map_conv_suvrMean.csv")) |> mutate(type = "amyloid_map_conv_suvrMean.csv")
map2 = read_table(file.path(resdir, "fig3/normMaps/amyloid_map_ABdiagnosed_suvrmean.csv")) |> mutate(type = "amyloid_map_ABdiagnosed_suvrmean")
map3 = read_table(file.path(resdir, "fig3/normMaps/amyloid_map_convVnonconv_suvrdiff_d.csv")) |> mutate(type = "amyloid_map_convVnonconv_suvrdiff_d")
map4 = read_table(file.path(resdir, "fig3/normMaps/amyloid_map_ABdiagnosedVnonconv_suvrdiff_d.csv")) |> mutate(type = "amyloid_map_ABdiagnosedVnonconv_suvrdiff_d")
map1$SUVR_Z_within = NULL; map1$SUVR_Z_across = NULL
map2$SUVR_Z_within = NULL; map2$SUVR_Z_across = NULL


# Load spin test results ----
spinfiles = list.files(path = file.path(resdir, "fig3/spinTestAB"), pattern = "*.csv", full.names = TRUE)
df_spin = lapply(spinfiles, function(f) {
  read_delim(f, delim = " ", quote = "\"", show_col_types = FALSE) |> mutate(src = basename(f))
}) |> bind_rows() |>
  mutate(
    analysis = str_extract(src, "(?<=abmap\\d_)(Sensitivity\\d+|main)"),
  )


# fig 4 --------------------------------------------
files = list.files(path = paste0(resdir, "/fig4/cortmaps"), pattern = "thickness_cortmap*", full.names = TRUE)
df_fig4 = loadResults(files)


# fig 5 --------------------------------------------
df_models_fig5_1 = lapply(paste0(resdir, "/fig5/predTraj.csv"), function(f) {
  read_table(f) |> mutate(src = basename(f))
}) |> bind_rows()

df_models_fig5_2 = lapply(paste0(resdir, "/fig5/cortmapPredTraj.csv"), function(f) {
  read_table(f) |> mutate(src = basename(f))
}) |> bind_rows()

df_models_fig5_3 = lapply(paste0(resdir, "/fig5/rankABpredTraj.csv"), function(f) {
  read_table(f) |> mutate(src = basename(f))
}) |> bind_rows()

df_models_fig5_4 = lapply(paste0(resdir, "/fig5/df_rankABthick.csv"), function(f) {
  read_table(f) |> mutate(src = basename(f))
}) |> bind_rows()



# Extended Data Figure 1 --------------------------------------------
files = list.files(path = paste0(resdir, "/extendedDataFig1/cortmaps"), pattern = "thickness_cortmap*", full.names = TRUE)
df_extfig1 = loadResults(files)


# Extended Data Figure 2 --------------------------------------------
files = list.files(path = paste0(resdir, "/extendedDataFig2/cortmaps"), pattern = "thickness_cortmap*", full.names = TRUE)
df_extfig2 = loadResults(files)


# Extended Data Figure 3a-d --------------------------------------------
files = list.files(path = paste0(resdir, "/extendedDataFig3/cortmaps/matchF0"), pattern = "thickness_cortmap*", full.names = TRUE)
df_extfig3_left = loadResults(files)


# Extended Data Figure 3a-d --------------------------------------------
files = list.files(path = paste0(resdir, "/extendedDataFig3/cortmaps/matchF1"), pattern = "thickness_cortmap*", full.names = TRUE)
df_extfig3_right = loadResults(files)

# Extended Data Figure 4 --------------------------------------------
modelfiles01 <- list.files(path = paste0(resdir, "/extendedDataFig4/moreCovariates"), pattern = "thickness_cortmapAB.*_MAINrevALLCOV\\.csv", full.names = TRUE) # main models with reduced sample to match models ran with more covariates
modelfiles02 <- list.files(path = paste0(resdir, "/extendedDataFig4/moreCovariates"), pattern = "thickness_cortmapAB.*_revALLCOV\\.csv", full.names = TRUE) # models with more covariates

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
df_extfig4Orig = df_models01 %>% mutate(modelType = "mainOriginal")
df_extfig4More = df_models02 %>% mutate(modelType = "moreCovariates")


# Extended Data Figure 5 --------------------------------------------
df_powered_slope = fread(file.path(resdir, "extendedDataFig5/df_powered_slope.csv")) %>% mutate(analysisType = "slope effects (full sample analyses)") %>% select(analysisType, everything())
df_powered_level = fread(file.path(resdir, "extendedDataFig5/df_powered_level.csv")) %>% mutate(analysisType = "level effects (full sample analyses)") %>% select(analysisType, everything())
df_reduced_slope = fread(file.path(resdir, "extendedDataFig5/df_reduced_slope.csv")) %>% mutate(analysisType = "slope effects (reduced sample analyses)") %>% select(analysisType, everything())
df_reduced_level = fread(file.path(resdir, "extendedDataFig5/df_reduced_level.csv")) %>% mutate(analysisType = "level effects (reduced sample analyses)") %>% select(analysisType, everything())
df_predicted_slope = fread(file.path(resdir, "extendedDataFig5/df_predicted_slope.csv")) %>% mutate(analysisType = "slope effects (predicted AB)") %>% select(analysisType, everything())
df_predicted_level = fread(file.path(resdir, "extendedDataFig5/df_predicted_level.csv")) %>% mutate(analysisType = "level effects (predicted AB)") %>% select(analysisType, everything())
df_extfig5 = rbind(df_powered_slope,
      df_powered_level,
      df_reduced_slope,
      df_reduced_level,
      df_predicted_slope,
      df_predicted_level)


# write to excel source file ----
makeExcel = 1
if (makeExcel) {
  
  wb <- createWorkbook()
  header_style <- createStyle(textDecoration = "bold")
  
  # Add sheets
  sheets = c("Fig_1e-h",
             "Fig_2b-e",
             "Fig_3a",
             "Fig_3b",
             "Fig_3c",
             "Fig_3d",
             "Fig_3e-f",
             "Fig_4b-i",
             "Fig_5_predTraj",
             "Fig_5c",
             "Fig_5e",
             "Fig_5f-h",
             "ExtFig1",
             "ExtFig2",
             "ExtFig3a-d",
             "ExtFig3e-h",
             "ExtFig4_original",
             "ExtFig4_moreCov",
             "ExtFig5a-f"
  )
  data_to_write = c("df_fig1",
                    "df_fig2",
                    "map1",
                    "map2",
                    "map3",
                    "map4",
                    "df_spin",
                    "df_fig4",
                    "df_models_fig5_1",
                    "df_models_fig5_2",
                    "df_models_fig5_3",
                    "df_models_fig5_4",
                    "df_extfig1",
                    "df_extfig2",
                    "df_extfig3_left",
                    "df_extfig3_right",
                    "df_extfig4Orig",
                    "df_extfig4More",
                    "df_extfig5"
                    
  )
  length(sheets) == length(data_to_write)
  
  for (sheetnum in 1:length(sheets)) {
    addWorksheet(wb, sheets[sheetnum])
  }
  
  header_style <- createStyle(textDecoration = "bold")
  
  for (sheetnum in 1:length(sheets)) {
    writeData(
      wb,
      sheet = sheets[sheetnum],
      get(data_to_write[sheetnum])
    )
  }
  
  
  for (sheetnum in 1:length(sheets)) {
    ncols = ncol(get(data_to_write[sheetnum]))
    addStyle(wb, sheet = sheets[sheetnum], style = header_style, rows = 1, 
             cols = 1:ncols,
             gridExpand = TRUE)
  }
  saveWorkbook(wb, file = file.path(here(), "/sourceData.xlsx"), overwrite = TRUE)
  file.path(here(), "/sourceData.xlsx")
}
