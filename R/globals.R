load("./data/annual_data_TN.rda")
load("./data/catch_data_TN.rda")
load("./data/annual_data_TP.rda")
load("./data/catch_data_TP.rda")
load("./data/sh_file.rda")
globalVariables(c("AA", "SumDS", "SumPSCmp1", "CatchToRiver", "BB", "Year",
                  "CatchLoad", "ForestFraction", "HydroID", "InvNrmRain",
                  "LakeFrRet", "NrmLengthKm", "Shreve", "To_catch", "YearValue",
                  "YearlyNMass", ":=", ".", "%dopar%", "alpha_p", "alpha_l",
                  "sd_C", "num_cores", "n_iter", "task", "launch_green",
                  "ObsLoad", "shrLevel", "value", "DF_LoadAvg_Type_Shr_Year",
                  "DF_LoadAvg_Type_Year", "Total", "variable", "scale_barRef",
                  "scale_barTextS", "scale_barRefs", "scale_barTexts",
                  "DrainAreaS", "Atm", "Atm2Diff", "Atm2DiffAgri",
                  "Atm2LandRet", "Bg", "CatchLakeRet", "CatchRivRet",
                  "DiffAgri2LandRet", "Diffuse", "Fix", "LandRet", "Man", "Min",
                  "Ps", "SD2Diff", "SD2LandRet", "Sd", "Soil", "YearlyMass",
                  "AllDesting", "AllOrig", "BasinID",  "DiffuseAgri",
                  "NextDownID", "RemaiHydroId", "WholeHydroID", "PredictLoad"))
