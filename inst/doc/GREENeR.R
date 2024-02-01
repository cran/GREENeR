## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
# load the GREENeR package
library(GREENeR)

## -----------------------------------------------------------------------------
# load the topological sequence of catchments and complementary info
data(catch_data_TN)
head(catch_data_TN)

## -----------------------------------------------------------------------------
# load the sources of TN for each year and catchment
data(annual_data_TN)
head(annual_data_TN)

# load the sources of TP for each year and catchment
data(annual_data_TP)
head(annual_data_TP)

## -----------------------------------------------------------------------------
# load the geographical information of the basin (require for some functionalities)
data(sh_file)
head(sh_file)

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the barplot for the Lay TN and TP Scenarios
#  input_plot(annual_data_TN, sh_file, "Lay", "B")
#  input_plot(annual_data_TP, sh_file, "Lay", "B")

## ---- echo=FALSE, out.width='100%', fig.cap="Barplot for the Lay TN scenario."----
knitr::include_graphics('figures/barplot_LayTN.png')

## ---- echo=FALSE, out.width='100%', fig.cap="Barplot for the Lay TP scenario."----
knitr::include_graphics('figures/barplot_LayTP.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the density plots for the Lay TN and TP Scenarios
#  input_plot(annual_data_TN, "Lay", "D",coef_SD=1, sh_file)
#  input_plot(annual_data_TP, "Lay", "D",coef_SD=1, sh_file)

## ---- echo=FALSE, out.width='100%', fig.cap="Density plot for the Lay TN scenario."----
knitr::include_graphics('figures/densityplot_LayTN.png')

## ---- echo=FALSE, out.width='100%', fig.cap="Density plot for the Lay TP scenario."----
knitr::include_graphics('figures/densityplot_LayTP.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the time serie plot type 1 (areas)
#  input_Tserie(catch_data_TN, annual_data_TN, "Lay", "gr1")
#  input_Tserie(catch_data_TP, annual_data_TP, "Lay", "gr1")

## ---- echo=FALSE, out.width="45%", out.height="25%",fig.cap="Time series plot type 1 (areas) for TN and TP.",fig.show='hold',fig.align='center'----
knitr::include_graphics(c("figures/timeserieTN_type1.png","figures/timeserieTP_type1.png"))

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the time serie plot type 2 (lines)
#  input_Tserie(catch_data_TN, annual_data_TN,  "Lay", "gr2")
#  input_Tserie(catch_data_TP, annual_data_TP,  "Lay", "gr2")

## ---- echo=FALSE, out.width="45%", out.height="25%",fig.cap="Time series plot type 2 (lines) for TN and TP.",fig.show='hold',fig.align='center'----
knitr::include_graphics(c("figures/timeserieTN_type2.png","figures/timeserieTP_type2.png"))

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # the time serie plot type 3 (areas by year and km2)
#  input_Tserie_area(catch_data_TN, annual_data_TN, sh_file, "Lay", "gr1")
#  input_Tserie_area(catch_data_TP, annual_data_TP, sh_file, "Lay", "gr1")

## ---- echo=FALSE, out.width="45%", out.height="25%",fig.cap="Time series plot type 3 (areas by year and km2) for TN and TP.",fig.show='hold',fig.align='center'----
knitr::include_graphics(c("figures/timeserieTN_type3.png","figures/timeserieTP_type3.png"))

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the time serie plot type 4 (by km2 and Shreve)
#  input_Tserie_area(catch_data_TN, annual_data_TN, sh_file, "Lay", "gr2")
#  input_Tserie_area(catch_data_TP, annual_data_TP, sh_file, "Lay", "gr2")

## ---- echo=FALSE, out.width="45%", out.height="25%",fig.cap="Time series plot type 4, thatcompares the levels of nutrient inputs in three zones of the basin.",fig.show='hold',fig.align='center'----
knitr::include_graphics(c("figures/timeserieTN_type4.png","figures/timeserieTP_type4.png"))

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the Input Load Map by source type 1 (kt/year)
#  input_maps(catch_data_TN, annual_data_TN, sh_file, "gr1", legend_position = 3)

## ---- echo=FALSE, out.width='100%', fig.cap="Input map by source type 1 for TN scenario"----
#![Input maps.](figures/input_maps.png){width=75%} #300px for HTML
knitr::include_graphics('figures/input_mapsTN_type1.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the Input Load Map by source type 1 (kt/year)
#  input_maps(catch_data_TP, annual_data_TP, sh_file, "gr1", legend_position = 3)

## ---- echo=FALSE, out.width='100%', fig.cap="Input map by source type 1 for TP scenario"----
#![Input maps.](figures/input_maps.png){width=75%} #300px for HTML
knitr::include_graphics('figures/input_mapsTP_type1.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the Input Load Map by source type 2 (kt/y/km2)
#  input_maps(catch_data_TN, annual_data_TN, sh_file, "gr2", legend_position = 3)

## ---- echo=FALSE, out.width='100%', fig.cap="Input map by source type 2 for TN scenario"----
knitr::include_graphics('figures/input_mapsTN_type2.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # the Input Load Map by source type 2 (kt/y/km2)
#  input_maps(catch_data_TP, annual_data_TP, sh_file, "gr2", legend_position = 3)

## ---- echo=FALSE, out.width='100%', fig.cap="Input map by source type 2 for TP scenario"----
knitr::include_graphics('figures/input_mapsTP_type2.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  LakeRetent_plot(catch_data_TN)

## ---- echo=FALSE, out.width='100%', fig.cap="Lake retention values distribution for TN scenario"----
knitr::include_graphics('figures/LakeRetent_plot.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  TN_ref_values <- annual_data_TN[!is.na(annual_data_TN$YearlyMass),]
#  nrow(TN_ref_values)
#  table(TN_ref_values$HydroID,TN_ref_values$YearValue)
#  
#  references_plot(annual_data_TN)

## ---- echo=FALSE, out.width='100%', fig.cap="Observed references distribution for TN scenario"----
knitr::include_graphics('figures/references.png')

## ---- message=FALSE-----------------------------------------------------------
# Parameter for the calibration of the model
# the lower limits for all params (alpha_P, alpha_L, sd_coef)
low <- c(10, 0.000, 0.1) 
# the upper limits for all params (alpha_P, alpha_L, sd_coef)  
upp <- c(50, 0.08, 0.9)       

# number of iterations
n_iter <- 200 

# years in which the model should be executed
years <- c(2003:2009)

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # execution of the calibration
#  DF_calib_TN <- calib_green(catch_data_TN, annual_data_TN, n_iter, low, upp, years)
#  DF_calib <- data.frame(DF_calib)

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Generating the box plots
#  rateBS <-  5  #  percent of parameters selected from the whole calibration set
#  calib_boxplot(DF_calib_TN, rateBS)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/calibration_boxplotTN.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Generating the dot plots
#  Gof_mes <- "NSE"
#  calib_dot(DF_calib_TN, Gof_mes)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/calibration_dotTN.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Generating the scatter plots
#  Gof_mes <- "NSE"
#  scatter_plot(DF_calib_TN, Gof_mes)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/calibration_scatterTN.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  Gof_mes <- "NSE" # according the NSE goodness of fit metric
#  NSE_bestParams <- select_params(DF_calib_TN, Gof_mes)
#  NSE_bestParams
#  Param_NSE2 <- as.numeric(NSE_bestParams[2:4])

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  Gof_mes <- "rNSE"
#  rNSE_bestParams <- select_params(DF_calib_TN,Gof_mes)
#  rNSE_bestParams
#  Param_rNSE2 <- as.numeric(rNSE_bestParams[2:4])

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # annual scatter plot comparing observed vs modeled loads by year
#  label_plot <- "NSE best params for TN in the Lay"
#  Themx <- max(annual_data_TN$YearlyMass, na.rm = TRUE)
#  simobs_annual_plot(catch_data_TN, annual_data_TN, Param_NSE2[1], Param_NSE2[2], Param_NSE2[3], years, label_plot, Themx)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/annual_scatter_vs_modeled_loadsTN.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  setPlabels <- c("NSE", "rNSE")
#  label_plot <- "Comparing two sets of parameters for Lay TN"
#  compare_calib(catch_data_TN, annual_data_TN, Param_NSE2[1], Param_NSE2[2], Param_NSE2[3], Param_rNSE2[1], Param_rNSE2[2], Param_rNSE2[3], years, label_plot, setPlabels)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/comparing_twosets_paramsTN.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Computing the nutrient balance
#  years <- c(2003:2009)
#  Nut_bal_TN <- region_nut_balance(catch_data_TN, annual_data_TN, Param_NSE2[1], Param_NSE2[2], Param_NSE2[3], years, atm_coeff = 0.38)
#  # Plot the sankey plot with the result of the balance
#  sank <- N4_sankey(Nut_bal_TN)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/n4sankey.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  years <- c(2003:2009)
#  basin_sa <- green_shares(catch_data_TN, annual_data_TN, Param_NSE2[1], Param_NSE2[2], Param_NSE2[3], years)

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # The title of the plot
#  plot_title <- "Time series Load Output for the Lay Basin"
#  outFile <- tempfile(pattern = "OutletLoad", fileext = ".csv")
#  nutrient_tserie(basin_sa, plot_title, "gr1", outFile)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/timeserie_loadoutputLay_type1.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Output Load Basin average time series (lines)
#  nutrient_tserie(basin_sa, plot_title, "gr2")

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/timeserie_loadoutputLay_type2.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Total Load in the Basin Outlet by source apportionment time series (lines)
#  nutrient_tserie_darea(basin_sa, sh_file, "Lay Basin", "gr1")

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/timeserie_totalloadLay.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  style <- "log10"
#  nutrient_maps(basin_sa, sh_file, "gr1", style, legend_position = 3)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/nutrient_maps_outputLoadsType1.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Basin Output Total Load  Maps
#  nutrient_maps(basin_sa, sh_file, style, "gr2", legend_position = 3)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/nutrient_maps_outputLoadsType2.png')

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  # Basin Output Specific Load by km2 Maps
#  nutrient_maps(basin_sa, sh_file, style, "gr3", legend_position = 3)

## ---- echo=FALSE, out.width='100%'--------------------------------------------
knitr::include_graphics('figures/nutrient_maps_outputLoadsType3.png')

