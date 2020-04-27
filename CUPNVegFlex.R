# CHECK FOR NA.RM RESULTS FOR SUM, MEAN, ETC.
# In shapefiles--Warning in bind_rows_(x, .id) :
#binding character and factor vector, coercing into character vector
# When i tried to read the temp-out qaqc trees, $TreeCounts$ECDF_plot
# Error: NULL value passed as symbol address

### Load libraries -----
packages <- c("plyr", "magrittr", "reshape2", "plotrix", "cowplot", "lubridate", "tidyverse", "scales", "stringr", "rgdal", "here")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})

### CALLED FUNCTIONS -----
FuncReformatMC <- function (Dat_df) {
  # Function to fill in missing MC data
  #
  # Args:
  #   Dat_df:  A data frame with the raw module-corner data
  # Returns:
  #   A data frame with the corrected module-corner data
  Dat_df %<>%
    dplyr::mutate(M1C2 = ifelse(!is.na(M1C2), M1C2, ifelse(is.na(M1C4), 0, 1))) %>%
    dplyr::mutate(M1C4 = ifelse(!is.na(M1C4), M1C4, ifelse(M1C2==0, 0, 1))) %>%
    dplyr::mutate(M2C2 = ifelse(!is.na(M2C2), M2C2, ifelse(is.na(M2C3), 0, 1))) %>%
    dplyr::mutate(M2C3 = ifelse(!is.na(M2C3), M2C3, ifelse(M2C2==0, 0, 1))) %>%
    dplyr::mutate(M3C2 = ifelse(!is.na(M3C2), M3C2, ifelse(is.na(M3C4), 0, 1))) %>%
    dplyr::mutate(M3C4 = ifelse(!is.na(M3C4), M3C4, ifelse(M3C2==0, 0, 1))) %>%
    dplyr::mutate(M4C2 = ifelse(!is.na(M4C2), M4C2, ifelse(is.na(M4C3), 0, 1))) %>%
    dplyr::mutate(M4C3 = ifelse(!is.na(M4C3), M4C3, ifelse(M4C2==0, 0, 1)))
  return(Dat_df)
}

FuncFormatVeg <- function () {
  # Function to format vegetation data for EDA.
  #
  # Returns:
  #   A formatted data frame of the vegetation data
  VegTemp_df <- read_csv(here::here("Data_in", "VegDat.csv"))
  TempCG_df <- read_csv(here::here("Data_in", "FinalCG.csv"))
  
  VegTemp_df$CoverPercent <- sapply(as.character(VegTemp_df$CoverClass), switch,
           `1` = sqrt(0.01*0.1),
           `2` = sqrt(0.1*1),
           `3` = sqrt(1*2),
           `4` = sqrt(2*5),
           `5` = sqrt(5*10),
           `6` = sqrt(10*25),
           `7` = sqrt(25*50),
           `8` = sqrt(50*75),
           `9` = sqrt(75*95),
           `10` = sqrt(95*100)) # calculate geometric mean for each cover class interval (replaces the arithmetic mean that is exported from Access)
  
  PlotCG_df <- TempCG_df %>% # for each plot, gives the short CG
    dplyr::select(EventID=event_name_date_calc, CGShort=lastCommGroup) %>%
    dplyr::mutate(Plot = substr(EventID, start=1, stop=7)) %>% # 
    dplyr::select(Plot, CGShort) %>%
    distinct()
  if (length(DupPlots <- PlotCG_df$Plot[duplicated(PlotCG_df$Plot)]) > 0) stop("Multiple community groups provided for these plots: ", DupPlots) # >>>>>>>>>>>>>>>>>>>>>>>> CHECK THIS ERROR MESSAGE
  
  VegTemp2_df <- VegTemp_df %>%
    dplyr::rename(Year = EventYear, M1C2=PresMod1Corn1, M1C4=PresMod1Corn2, M2C2=PresMod2Corn1, M2C3=PresMod2Corn2, M3C2=PresMod3Corn1, M3C4=PresMod3Corn2, M4C2=PresMod4Corn1, M4C3=PresMod4Corn2) %>%
    dplyr::mutate(GRTS3 = formatC(GRTSOrder, width=3, zero.print=TRUE, flag="0"), # create PlotGRTS with leading zero's to GRTS numbers to ensure 3 digits
           Plot = paste0(SubPark, GRTS3),  # create unique ParkGRTS
           Cycle = cut(Year, breaks=seq(2010, max(Year, na.rm = TRUE) + 4, by = 5)), # create Cycle classification
           Q = ifelse(QAQCPlot, "Q", ""),
           EventID = paste0(Plot, Q, " - ", Year),
           Year = as.integer(Year),
           QAQCPlot = as.logical(QAQCPlot),
           CoverPercent = c(0.05, 0.55, 1.50, 3.5, 7.50, 17.50, 37.50, 62.50, 85.00, 97.50)[CoverClass], # midpoint of cover class ranges--some of values in Access database are incorrect
           PlantName = gsub(" ", "_", PlantName, fixed=TRUE),
           PlantGenus = gsub('_.*$','', PlantName),
           PlantConfirmed = regexpr("_", PlantName) > 0) %>% # plant identified to species level
    left_join(PlotCG_df, by = "Plot")
  
  VegTemp2_df$BotanistName[VegTemp2_df$BotanistName=="_"] <- NA #  replace BotanistName dash cells with NA
  VegTemp2_df$DifficultID <- apply(VegTemp2_df["PlantName"], 1, FUN=function(x) any(startsWith(x, c("Carex", "Poaceae", "Carya", "Galium", "Sanicula")))) # TRUE if the plant name starts with any of the vector elements
  VegTemp2_df$DifficultID[VegTemp2_df$PlantName=="Galium_circaezans"] <- FALSE
  
  missing_plots_main <- unique(PlotCG_df$Plot)[!unique(PlotCG_df$Plot) %in% unique(VegTemp2_df$Plot)] # <<<<<<< THROW ERROR MESSAGE--THESE PLOTS ARE IN PlotCG, but not in the VegStats data
  missing_PlotCG <- unique(VegTemp2_df$Plot)[!unique(VegTemp2_df$Plot) %in% unique(PlotCG_df$Plot)] # <<<<<<< THROW ERROR MESSAGE--THESE PLOTS ARE IN VegStats data, but not in PlotCG
  VegTemp2_df$CGShort[VegTemp2_df$Plot=="SHIL010"] <- "G159" # <<<<<<<<<<<<<<< ADD MANUALLY--CUPN NEEDS TO FIX THIS
  VegTemp2_df$CGShort[VegTemp2_df$Plot=="CUGA001"] <- "G165" # <<<<<<<<<<<<<<< ADD MANUALLY--CUPN NEEDS TO FIX THIS

  VegDat_df <- FuncReformatMC(Dat_df = VegTemp2_df) %>%
    dplyr::select(Cycle, Year, QAQCPlot, Park, SubPark, Plot, EventID, OptimalSize = optimalSize, PlantName, PlantGenus, PlantConfirmed, DifficultID, ExoticPlant, HighPriorityExotic, CoverPercent, CoverClass, M1C2, M1C4, M2C2, M2C3, M3C2, M3C4, M4C2, M4C3, BotanistName, Latitude, Longitude, CGShort, CommunityGroup) %>%
    arrange(Park, SubPark, Year, Plot, QAQCPlot)
  
  sel_options <- list(
    Park = unique(VegDat_df$Park),
    SubPark = unique(VegDat_df$SubPark),
    Cycle = unique(VegDat_df$Cycle))
  saveRDS(sel_options, "Temp_out/sel_options.RDS")
    
  # write.csv(VegDat_df, "Temp_out/VegDat_cleaned.csv", row.names = FALSE)
  write_csv(VegDat_df, here::here("Temp_out", "VegDat_cleaned.csv"))
  rm(VegTemp_df, VegTemp2_df)
  return(VegDat_df)
}

FuncMapBounds <- function () {
  CUPNmap <- readOGR(here::here("Map_files", "CUPN_shapefiles", "CUPN.shp"))
  CUPNmap <- spTransform(CUPNmap, CRS("+init=epsg:4326")) # make sure shapefiles are in same projection as google maps
  CUPNmap@data$id <- rownames(CUPNmap@data)
  CUPNmapbounds_df <- CUPNmap %>% broom::tidy() %>%  # convert to data frame with boundary information
    left_join(CUPNmap@data[,c("SubUnitLab", "id")], by = "id")
  saveRDS(CUPNmapbounds_df, "Map_files/CUPNmapbounds.RDS")
}
 
FuncMapRpt <- function(Dat_df, CUPNmapbounds, cyc, metric) {
  SPExoticsMap_list <- list()
  for (SP in unique(Dat_df$SubPark)) {
    incProgress(1/10, detail = SP)
    CUPNmapbounds_SP <- subset(CUPNmapbounds, SubUnitLab==SP)
    
    MapDat_df <- FuncMaps(
      Dat_df = Dat_df,
      p = SP,
      cyc = cyc,
      metric = metric)
    
    map_center <- c(mean(c(min(CUPNmapbounds_SP$long), max(CUPNmapbounds_SP$long))), mean(c(min(CUPNmapbounds_SP$lat), max(CUPNmapbounds_SP$lat))))
    map_zoom <- min(MaxZoom(lonrange = range(CUPNmapbounds_SP$long), latrange = range(CUPNmapbounds_SP$lat)))
    
    basemap <- leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
      addProviderTiles("Esri.WorldImagery") %>%
      setView(lng = map_center[1], lat = map_center[2], zoom = map_zoom)
    pieces <- unique(unlist(subset(CUPNmapbounds_SP, select = "piece")))
    for (p in pieces) {
      d <- subset(CUPNmapbounds_SP, piece==p)
      basemap %<>% 
        addPolylines(lat = d$lat, lng = d$long, weight = 2, color = "black", opacity = 1)
    }
    
    SPExoticsMap_list[[SP]] <- basemap %>%
      addMinicharts(
        lng = MapDat_df$Longitude, lat = MapDat_df$Latitude,
        type = "bar",
        width = 25,
        height = 25,
        chartdata = as.matrix(MapDat_df[, c("Native", "Exotic")]),
        colorPalette = c("purple", "white"),
        opacity = 1
      )
  } # end SubParks
  return(SPExoticsMap_list)
}

FuncOHSanderPlot <- function(dat_df, comm_type) {
  # Function to generate plots of Sander metric
  #
  # Args:
  #   Dat_df:  The final Sander data
  # Returns:
  #   Plots
  ggplot(subset(dat_df, CGtype==comm_type), aes(x=reorder(SubPark,`50%`),  y=`50%`)) +
    geom_point() +
    geom_linerange(aes(ymin=`5%`, ymax=`95%`)) +
    geom_hline(yintercept=1, linetype="dashed", color="blue") +
    ylim(-0.1, 2.1) +
    labs(y="Count of oak-hickory saplings per 10m2", x="SubPark", title=paste0("Oak-hickory sapling densities in ", comm_type, " communities")) +
    theme_bw(base_size = 10)
}

FuncOHWilcox <- function(Dat_df, metric, CycLevs = CycLevs) {
  # Function to generate histogram of QAQC for oak-hickory regeneration metrics
  #
  # Args:
  #   Dat_df:  The final OH data
  #   metric:  "PartialOSI", "log10ImpScore"
  #   CycLevs:  Vector of survey cycles, e.g., "(2010,2015]" "(2015,2020]"
  # Returns:
  #   Wilcoxon signed-rank test table by burn status. Use only plots with partial OSI calculated for both cycles. Therefore, the count of resurveys in this table is less than actual number resurveyed (b/c sometimes a plot that was resurveyed had OSI=NaN for one of the surveys)
  names(Dat_df)[names(Dat_df)==metric] <- "Metric"
  Datwide_df <- Dat_df %>%
    filter(QAQCPlot == FALSE) %>%
    dplyr::select(SubPark, Plot, CGShort, Cycle, Metric) %>%
    spread(key = Cycle, value = Metric)
  Datwide_df$delta <- unlist(Datwide_df[max(CycLevs)] - Datwide_df[min(CycLevs)])
  
  Unburned = setdiff(unique(Dat_df$Plot), unique(Dat_df$Plot[Dat_df$PriorBurn==TRUE]))
  
  RegenBurned_list <- vector("list", length = length(unique(Dat_df$SubPark)))
  for (s in unique(Dat_df$SubPark)) {
    SubDat_df <- Datwide_df %>%
      filter(SubPark==s) %>% # for counting sample size per survey, keep ones that were not resurveyed
      dplyr::mutate(Unburned = Plot %in% Unburned)
    unburned_N <- nrow(subset(SubDat_df, !is.na(delta) & Unburned==TRUE))
    # Wilcox test only uses plots that were resurveyed. Calculate only if at least 6 unburned resurveyed plots.
    if(unburned_N  >= 6) { # ifelse doesn't work for some reason
      unburned_stat <- wilcox.test( 
        unlist(subset(SubDat_df, !is.na(delta) & Unburned==TRUE, select = max(CycLevs))),
        unlist(subset(SubDat_df, !is.na(delta) & Unburned==TRUE, select = min(CycLevs))),
        paired=TRUE, exact=TRUE, conf.int = TRUE, conf.level=0.9)
      unburned_CI90 <- round(as.numeric(unburned_stat$conf.int), 2)
      unburned_p <- round(as.numeric(unburned_stat$p.value), 2)
    } else {
      unburned_CI90 <- c(NA, NA)
      unburned_p <- NA
    }
    # Specify in notes that it's a two-sided 90% CI. 'dat2' is first because desired test statistic is 'dat2 - dat1'
    
    burned_N <- nrow(subset(SubDat_df, !is.na(delta) & Unburned==FALSE))
    if(burned_N  >= 6) {
      burned_stat <- wilcox.test( 
        unlist(subset(SubDat_df, !is.na(delta) & Unburned==FALSE, select = max(CycLevs))),
        unlist(subset(SubDat_df, !is.na(delta) & Unburned==FALSE, select = min(CycLevs))),
        paired=TRUE, exact=TRUE, conf.int = TRUE, conf.level=0.9)
      burned_CI90 <- round(as.numeric(burned_stat$conf.int), 2)
      burned_p <- round(as.numeric(burned_stat$p.value), 2)
    } else {
      burned_CI90 <- c(NA, NA)
      burned_p <- NA
    }
    
    RegenBurned_list[[s]] <- c(
      s,
      sum(!is.na(SubDat_df$delta)),
      unburned_N,
      unburned_CI90,
      unburned_p,
      burned_N,
      burned_CI90,
      burned_p)
  }
  BurnWilcox_df <- data.frame(do.call("rbind", RegenBurned_list))
  
  for(x in 2:ncol(BurnWilcox_df)) {
    BurnWilcox_df[,x] <- as.numeric(as.character(BurnWilcox_df[,x]))
  }
  names(BurnWilcox_df) <- c("SubPark", "N_resurveys", "UNBURNED_N", "Low90CI", "High90CI", "P_value", "BURNED_N", "Low90CI", "High90CI", "P_value")
  return(BurnWilcox_df)
}

# FuncOHPlotHist <- function(Dat_df, metric) {
#   # Function to generate histogram of QAQC for oak-hickory regeneration metrics
#   #
#   # Args:
#   #   Dat_df:  The final OH data
#   #   metric:  "PartialOSI", "log10ImpScore"
#   # Returns:
#   #   Histogram of OH metric QAQC
#   RelOH_df <- Dat_df %>%
#     filter(PlotCycle %in% unique(Dat_df$PlotCycle[Dat_df$QAQCPlot==TRUE]))
#   RelOHWide_df <- spread(data = RelOH_df[, c("Plot", "Cycle", "QAQCPlot", metric)], key = QAQCPlot, value = metric)    
#   names(RelOHWide_df) <- c("Plot", "Cycle", "Q1", "Q2") # <<<<<<<<<<< MAY NEED TO GENERALIZE TO >2 SURVEY CYCLES
#   RelOHWide_df$Diff <- RelOHWide_df$Q1 - RelOHWide_df$Q2 # only for oak-hickory plots <<<<<<<<<<<<<< THIS MAY DIFFER FROM BEFORE
#   
#   OHHist_plot <- ggplot(RelOHWide_df, aes(x=Diff)) +
#     geom_histogram(binwidth = 0.1, center=0, colour="black") +
#     geom_vline(xintercept=0, color="red") +
#     labs(x=paste0("Difference in Oak-Hickory ", switch(metric, "PartialOSI" = "partial OSI", "log10ImpScore" = "log10(Importance)")), y="# of QAQC surveys", subtitle=paste0("Across ALL Community Groups, 90% CI: ", round(quantile(RelOHWide_df$Diff, probs=c(.05,.95)), 2)[1], " to ", round(quantile(RelOHWide_df$Diff, probs=c(.05,.95)), 2)[2])) +
#     theme_bw(base_size = 10)
#   return(OHHist_plot)
# }
# 
# FuncOHBoxplot <- function(Dat_df, metric) {
#   # Function to generate boxplots of oak-hickory regeneration metric by Subpark-Burn
#   #
#   # Args:
#   #   Dat_df:  The final OH data
#   #   metric:  "PartialOSI", "log10ImpScore"
#   # Returns:
#   #   Boxplots by SubPark-Burn
#   names(Dat_df)[names(Dat_df)==metric] <- "Metric"
#   OH_Boxplots <- ggplot(data=subset(Dat_df, QAQCPlot==FALSE), aes(y=Metric, x=SubPark, fill=PriorBurn)) +
#      geom_boxplot() +
#      scale_fill_manual(values=c("gray", "red")) +
#      labs(title=paste0(metric," by SubPark and prescribed burn history"), subtitle="(Oak-hickory communities only; Excludes QAQC data)", y=metric) +
#      theme_bw(base_size = 10)
#   return(OH_Boxplots)
# }
# 
# FuncOHBurnTiming <- function(Dat_df, metric) {
#   # Function to generate scatterplots of oak-hickory regeneration metric by time-since-burn
#   #
#   # Args:
#   #   Dat_df:  The final OH data
#   #   metric:  "PartialOSI", "log10ImpScore"
#   # Returns:
#   #   Scatterplots showing relationship between time-since-burn and regeneration metric, in Oak-Hickory communities
#   names(Dat_df)[names(Dat_df)==metric] <- "Metric"
#   TimeSinceBurn_df <- Dat_df %>%
#     filter(QAQCPlot==FALSE & !is.na(MostRecentBurn) & !is.na(SurveyDate)) %>%
#     dplyr::mutate(YrsSinceBurn = as.numeric((SurveyDate - MostRecentBurn)/365))
#   max_x <- ceiling(max(TimeSinceBurn_df$YrsSinceBurn))
#   BurnTiming_plot <- ggplot(data=TimeSinceBurn_df, aes(x = YrsSinceBurn, y = Metric)) +
#     geom_point() +
#     labs(x="Years since most recent prescribed burn", subtitle=paste0("Relationship between ", metric, " and time-since-prescribed burn, by Community Group"), y=metric) +
#     scale_x_continuous(labels=seq(0, max_x, by = 2), breaks=seq(0, max_x, by = 2), limits=c(0, max_x)) +
#     geom_vline(xintercept=seq(1, max_x, by = 1), color="red", linetype="dotted") +
#     geom_vline(xintercept = 0, color="red") +
#     theme_bw(base_size = 10) +
#     facet_wrap(~CGShort)
#   return(BurnTiming_plot)
# }

FuncHerbsCover <- function (Dat_df) {
  # Function to summarize herb match types. Use only for plot-level analysis, and all species or hicov (not rolled to genus) 
  # 
  # Args:
  #   Dat_df:  A data frame of the cleaned herb data
  # Returns:
  #   A list object that can be printed via an .RMD file
  HerbsCC_df <- Dat_df
  HerbsCC_df$CC[is.na(HerbsCC_df$CC)] <- 0 # Zero for now
  HerbsCC_df <- HerbsCC_df %>%
    dplyr::select(PlotID, QAQC, PlantName, CC) %>%
    spread(key=QAQC, value=CC)
  HerbsCC_df <- HerbsCC_df[complete.cases(HerbsCC_df),]
  colnames(HerbsCC_df)[3:4] <- c("Nsurv", "Qsurv")
  HerbsCC_df$Nsurv[HerbsCC_df$Nsurv==0] <- NA # convert back to NA
  HerbsCC_df$Qsurv[HerbsCC_df$Qsurv==0] <- NA
  HerbsCC_df$AbsDiff <- abs(HerbsCC_df$Qsurv - HerbsCC_df$Nsurv)
  HerbsCC_df <- as.data.frame(HerbsCC_df)
  return(HerbsCC_df)
}
    
FuncHerbsMC <- function (Dat_df) {
  # Function to fill in module presence-absence for herb data
  # 
  # Args:
  #   Dat_df:  A data frame with the formatted module-corner data
  # Returns:
  #   A data frame with presence-absence data at the module level
  #
  Dat2_df <- Dat_df %>%
    rowwise() %>%
    dplyr::mutate(InPlot=sum(M1C2, M1C4, M2C2, M2C3, M3C2, M3C4, M4C2, M4C3)) %>%
    filter(InPlot > 0) %>% # only include the ones that were detected in at least one module-corner (sometimes a record has a cover class but wasn't detected in any module-corner)
    rowwise() %>%
    dplyr::mutate(InPlot = InPlot > 0, # this is a given (should be TRUE always because already filtered data for this)
           M1=sum(M1C2, M1C4) > 0, # presence-absence at the module level
           M2=sum(M2C2, M2C3) > 0,
           M3=sum(M3C2, M3C4) > 0,
           M4=sum(M4C2, M4C3) > 0) %>%
    dplyr::mutate_at(vars(names(Dat_df[,grep("M.C.", colnames(Dat_df))])),funs(.>1)) %>% # converts to presence-absence at the 10m2 plot level
    ungroup() %>%
    dplyr::mutate(PlotID=as.factor(PlotID)) %>%
    dplyr::select(PlotID, QAQC, PlantName, InPlot, everything())
  return(Dat2_df)
}

FuncHerbsMatch <- function (Dat_df, group_vec) {
  # Function to summarize herb match types 
  # 
  # Args:
  #   Dat_df:  A match data frame
  #   group_vec:  Vector of variables to group by for summary, e.g., c("PlotID", "10SqM")
  # Returns:
  #   A list object that can be printed via an .RMD file
  
  HM_df <- Dat_df
  HM_df[, c("Nsurv", "Qsurv")][is.na(HM_df[, c("Nsurv", "Qsurv")])] <- FALSE
  HM_df$MatchTT <- HM_df$Nsurv + HM_df$Qsurv == 2
  HM_df$MatchType <- paste(HM_df$Nsurv, HM_df$Qsurv, sep="-")
  
  HS_df <- HM_df %>%
    group_by(.dots = group_vec) %>%
    dplyr::summarize(
      Nonly = sum(MatchType=="TRUE-FALSE"),
      Qonly = sum(MatchType=="FALSE-TRUE"),
      NsurvSpp = sum(Nsurv),
      QsurvSpp = sum(Qsurv),
      PSturn = round(((Nonly + Qonly)/(NsurvSpp + QsurvSpp)) * 100),
      RichDiff = QsurvSpp - NsurvSpp,
      PercRichDiff = round((abs(RichDiff)/mean(c(NsurvSpp, QsurvSpp))) * 100), 
      UniqueSpp = sum(Nsurv | Qsurv), # any species detected at least once 
      MatchTrue = sum(MatchTT),
      MatchByMean = round(MatchTrue/mean(c(NsurvSpp, QsurvSpp)) * 100),
      MatchByUnique = round((MatchTrue/UniqueSpp) * 100)) %>%
    dplyr::mutate_if(is.numeric, as.integer) %>%
    ungroup() %>%
    as.data.frame()
  
  Temp_x <- abs(HS_df$RichDiff)
  HS_tab <- data.frame(table( factor(Temp_x, levels = min(Temp_x):max(Temp_x))))
  HS_tab$CDF <- round(cumsum(HS_tab$Freq)/sum(HS_tab$Freq)*100, 1)
  names(HS_tab) <- c("Difference in Counts", "Frequency", "CDF (%)")
  
  return_list <- list(HM_df, HS_df, HS_tab)
  return(return_list)
  }

FuncColMat <- function (QAQC_df, categ_vec, plot_title) {
  # Function to generate a color matrix of cover category changes. Numbers in cells are frequencies. 
  # 
  # Args:
  #   QAQC_df:  A data frame with two columns of data (called Value1, Value2), one is QAQC category and the other is non-QAQC category. Missing entries should be NA.
  #   categ_vec:  Vector of categories (groups), sorted
  # Returns:
  #   ggplot
  #
  QAQC_df[is.na(QAQC_df)] <- "XXX"
  categ_vec2 <- c(categ_vec, "XXX")
  QAQC_df[, "Value1"] <- factor(QAQC_df[, "Value1"], levels=categ_vec2)
  QAQC_df[, "Value2"] <- factor(QAQC_df[, "Value2"], levels=categ_vec2)
  
  QAQC_mat <- as.matrix(QAQC_df)
  QAQC_mat <- t(apply(QAQC_mat,1,sort, decreasing = TRUE)) # For upper triangular, within each row, the smaller number is first
  QAQC_squaremat <- table(factor(QAQC_mat[,1], levels=categ_vec2), factor(QAQC_mat[,2], categ_vec2))
  QAQC_final_df <- data.frame(QAQC_squaremat) 
  QAQC_final_df <- QAQC_final_df[as.integer(QAQC_final_df$Var2) <= as.integer(QAQC_final_df$Var1), ] # get rid of below-diagonal elements
  
  plot_ColMat <- 
    ggplot(QAQC_final_df, aes(x = Var1, y = Var2, fill = log(Freq))) +
    geom_tile(color = "black") +
    geom_text(aes(label = Freq), na.rm = TRUE) +
    scale_fill_gradient(low = "yellow", high = "dodgerblue3", na.value = "white") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(categ_vec2)) +
    labs(title = paste0(plot_title, " (N = ", sum(QAQC_squaremat), ")")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          panel.grid = element_blank())
  
  return(plot_ColMat)
}

FuncPlotTreesECDF <- function (dat, plot_vals) {
  # Function to generate a cumulative distribution curve. 
  # 
  # Args:
  #   dat:  A QAQC data frame.
  #   plot_vals:  Additional information for plot axes labels, etc.
  # Returns:
  #   png of plot
  plot.new()
  dev.control(displaylist="enable")
  
  plot.ecdf(abs(dat[["Metric"]]), col.01line="white", main = paste0("CDF for Difference in ", plot_vals$metric_nam), ylab="Proportion of paired observations", xlab=paste0(plot_vals$lab, " in ", plot_vals$metric_nam), cex.lab=1, cex.main=1)
  
  ecdf_plot <- recordPlot()
  invisible(dev.off())
  return(ecdf_plot)
}

FuncPlotHerbsRichUnmatched <- function (HM_df, MinDetect, ScaleLab) {
  # Function to generate a bar plot of proportion of species detections missed in QAQC surveys.
  # 
  # Args:
  #   HM_df:  Herb match data frame.
  #   MinDetect:  Minimum number of detections for a species to be included
  #   ScaleLab:  Character string indicating the scale of evaluation (options are "100SqM", "10SqM", "1SqM")
  # Returns:
  #   Bar plot
  #
  HFspecies <- which(table(HM_df$PlantName) >= MinDetect) # species detected on at least X plot/10SqM/1SqM-surveys
  HFmissed <- HM_df %>%
    filter(PlantName %in% names(HFspecies)) %>%
    group_by(PlantName) %>%
    dplyr::summarize(Count = n(),
              AbsMissed = Count - sum(MatchTT),
              PropMissed = as.integer(round((sum(!MatchTT)/Count) * 100))) %>%
    dplyr::mutate(PlantName = paste0(PlantName, " (", Count, ")"))
  
  BarMissed_plot <- ggplot() +
    geom_bar(data = subset(HFmissed, PropMissed > 50), aes(x=reorder(PlantName, -PropMissed), y = PropMissed), width=.5, stat="identity") +
    theme_bw(base_size=10) +
    ylim(0, 100) +
    geom_hline(yintercept=50, color="red") +
    labs(y="Proportion of Detections Unmatched (%)", x="Species", title= paste0("Proportion of ", ScaleLab, "Survey Detections Unmatched"), subtitle=paste0("(Species with at least ", MinDetect, " detections and >50% unmatched)")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    theme(plot.title = element_text(size=10), legend.title = element_blank())
  
  return(BarMissed_plot)
}

FuncPlotHerbsRichDiff_Hist <- function (HS_df, ScaleLab) {
  # Histogram showing difference between standard and QAQC herb richness counts.
  # 
  # Args:
  #   HS_df:  Herb species data frame.
  #   ScaleLab:  Character string indicating the scale of evaluation (options are "100SqM", "10SqM", "1SqM")
  # Returns:
  #   Histogram
  #
  HistRichDiff_plot <- 
    ggplot(HS_df, aes(x=RichDiff)) +
    geom_bar() +
    geom_vline(xintercept=0, col="red") +
    labs(subtitle="(Positive difference means QAQC richness is higher)", x=paste0("Difference in Herb ", ScaleLab, " Richness"),  y=paste0("Number of ", ScaleLab, "-Surveys")) +
    theme_bw(base_size = 12)
  
  return(HistRichDiff_plot)
}

FuncPlotHerbsRichDiff_Scatter <- function (HS_df, ScaleLab) {
  # Scatterplot showing difference between standard and QAQC herb richness counts.
  # 
  # Args:
  #   HS_df:  Herb species data frame.
  #   ScaleLab:  Character string indicating the scale of evaluation (options are "100SqM", "10SqM", "1SqM")
  # Returns:
  #   HScatterplot
  #
  axesmax <- max(HS_df$NsurvSpp, HS_df$QsurvSpp, na.rm = TRUE) + 1
  countmax <- HS_df %>%
    dplyr::count(NsurvSpp, QsurvSpp) %>%
    dplyr::top_n(n, n=1) %>%
    dplyr::select(n)
  countmax <- as.integer(unlist(countmax))
  
  PointsRichDiff_plot <-
    ggplot(HS_df, aes(x=NsurvSpp, y=QsurvSpp)) +
    geom_count() +
    scale_size_continuous(breaks = seq(1, max(5, countmax), by = 2), limits = c(1, max(5, countmax)), range = c(1, 5), name = "# obs") +
    geom_abline(slope=1, col="red") +
    xlim(0, axesmax) +
    ylim(0, axesmax) +
    coord_fixed() + # preserve the 1:1 aspect
    labs(x=paste0("Herb ", ScaleLab, " Richness in Standard Survey"), y=paste0("Herb ", ScaleLab, " Richness in QAQC Survey"), subtitle = "(Point size is scaled to the number of\noverlapping observations)") +
    theme_bw(base_size = 12)
  
  return(PointsRichDiff_plot)
}

FuncPlotHerbsECDF <- function (HS_df, ScaleLab) {
  # Function to generate a cumulative distribution curve. 
  # 
  # Args:
  #   HS_df:  Herb species data frame.
  #   ScaleLab:  Character string indicating the scale of evaluation (options are "100SqM", "10SqM", "1SqM")
  #   plot_vals:  Additional information for plot axes labels, etc.
  
  # Returns:
  #   png of plot
  plot.new()
  dev.control(displaylist="enable")
  
  plot.ecdf(abs(HS_df[["RichDiff"]]), col.01line="white", main = paste0("CDF for Difference in Herb ", ScaleLab, " Richness"), ylab="Proportion of paired observations", xlab="Absolute difference in herb richness", cex.lab=1, cex.main=1)
  
  ecdf_plot <- recordPlot()
  invisible(dev.off())
  return(ecdf_plot)
}

FuncBarHist <- function (dat, plot_vals, alt_xlab = NA) {
  p <- ggplot(data = dat, aes(x = Metric)) +
    scale_x_continuous(limits = c(-1,1)*plot_vals$range_max, breaks = seq(-1*plot_vals$range_max, plot_vals$range_max, by = plot_vals$range_int)) +
    theme_bw(base_size = 12)
  
  if (plot_vals$type == "discrete") { # if count data, return a bar plot
    p + geom_bar(fill = "white", color = "black", na.rm = TRUE) +
      geom_vline(xintercept = 0, col = "red") +
      labs(x = ifelse(is.na(alt_xlab), "Difference between surveys, q - n", alt_xlab), y = "Count")
  } else {
    p + geom_histogram(fill = "white", color = "black", binwidth = plot_vals$range_int, na.rm = TRUE) +
      geom_vline(xintercept = 0, col = "red") +
      labs(x = ifelse(is.na(alt_xlab), "% difference between surveys, (q - n)/mean", alt_xlab), y = "Count")
  }
}

FuncQAQCTreesFormat <- function(NDat_df, QDat_df, cname, NANA_mismatch) {
  # Function to format tree data for QAQC analysis
  # 
  # Args:
  #   NDat_df: The data file with non-QAQC tree data. Required columns are 'TreeID', 'QTreeID', and the variable of interest
  #   QDat_df: The data file with QAQC tree data. Required columns are 'TreeID', 'NTreeID', and the variable of interest
  #   cname:  Name of the column containing the variable of interest
  #   NANA_mismatch: TRUE if NA/NA should be considered a mismatch; FALSE if NA/NA should be reported and then omitted from analysis.
  #
  # Returns:
  #   Data frame for analysis
  #   ID for cases where NA for both surveys
  TempDat_df <- merge(subset(NDat_df, select=c("QTreeID", cname)), subset(QDat_df, select=c("TreeID", cname)), by.x="QTreeID", by.y="TreeID") # merge the data by Tree ID
  names(TempDat_df) <- c("QAQCtreeID", "Value1", "Value2")
  
  if(NANA_mismatch) {
    warn_missing <- NA
  } else {
    # ERROR CHECKING
    warn_missing <- TempDat_df$QAQCtreeID[is.na(TempDat_df$Value1)&is.na(TempDat_df$Value2)]
    TempDat_df <- TempDat_df[!TempDat_df$QAQCtreeID %in% warn_missing,] # remove cases where NA for both surveys
  }
  
  TempDat_df$Match <- TempDat_df$Value1==TempDat_df$Value2
  TempDat_df$Match[is.na(TempDat_df$Match)] <- FALSE # If at least one of the values was NA (not recorded), then Match = FALSE
  # return_list <- list(TempDat_df = TempDat_df, warn_missing = warn_missing) # it seems that none of the variables have warnings and when they do, it's a lot. If NA/NA is NOT a required field and therefore is NOT a mismatch, do I even need to add a warning? 
  return(TempDat_df)
}

FuncPlotTreesValDiff_Hist <- function (RT_df, metric_nam) {
  # Histogram showing difference between standard and QAQC herb richness counts.
  # 
  # Args:
  #   RT_df:  Tree QAQC data frame
  #   metric_nam:  Name of metric being evaluated
  # Returns:
  #   Histogram
  #
  RT_df$ValDiff <- RT_df$Qsurv - RT_df$Nsurv
  HistValDiff_plot <- 
    ggplot(RT_df, aes(x=ValDiff)) +
    geom_bar() +
    scale_x_continuous(breaks= pretty_breaks()) +
    geom_vline(xintercept=0, col="red") +
    labs(subtitle="(Positive difference means QAQC value was higher)", x=paste0("Difference in ", metric_nam),  y="Number of paired observations") +
    theme_bw(base_size = 12)
  
  return(HistValDiff_plot)
}

FuncPlotTreesValDiff_Scatter <- function (RT_df, metric_nam) {
  # Scatterplot showing difference between standard and QAQC herb richness counts.
  # 
  # Args:
  #   RT_df:  Tree QAQC data frame
  #   metric_nam:  Name of metric being evaluated
  # Returns:
  #   Scatterplot
  #
  axesmax <- max(RT_df$Nsurv, RT_df$Qsurv, na.rm = TRUE) + 1
  countmax <- RT_df %>%
    dplyr::count(Nsurv, Qsurv) %>%
    dplyr::top_n(n, n=1) %>%
    dplyr::select(n)
  countmax <- as.integer(unlist(countmax))
  
  PointsValDiff_plot <-
    ggplot(RT_df, aes(x=Nsurv, y=Qsurv)) +
    geom_count() +
    scale_size_continuous(breaks = seq(1, max(7, countmax), by = round(max(7, countmax)/7)), limits = c(1, max(7, countmax)), range = c(1, 7), name = "# obs") +
    geom_abline(slope=1, col="red") +
    xlim(0, axesmax) +
    ylim(0, axesmax) +
    coord_fixed() + # preserve the 1:1 aspect
    labs(x=paste0(metric_nam, " in Standard Survey"), y=paste0(metric_nam, " in QAQC Survey"), subtitle = "(Point size is scaled to the number of\noverlapping observations)") +
    theme_bw(base_size = 12)
  
  return(PointsValDiff_plot)
}

FuncPSA <- function(categ_vec, DatWide_df, c1, c2, missing_val=TRUE) {
  # Function to calculate proportion of specific agreement for categorical data (includes binary & nominal). For ordered categorical data, use FuncOrderedPSA. PSA measures the probability, if one observer detects something (positive), that the other observer will detect it also. Calculate 90% CI only if there are at least 5 detections (fairly arbitrary threshold).
  # 
  # Args:
  #   categ_vec:  Vector of categories (groups), sorted
  #   DatWide_df, c1, c2:  Dataframe with the columns c1 and c2 that should be compared for a match. The values in these columns are the recorded categories. Missing entries should be NA.
  #   missing_val: TRUE if there is a possibility of NAs in data
  #
  # Returns:
  #   Data frame summary
  #
  if(missing_val) {
    DatWide_df[, c1][is.na(DatWide_df[, c1])] <- "XXX"
    DatWide_df[, c2][is.na(DatWide_df[, c2])] <- "XXX"
    categ_vec <- c(categ_vec, "XXX")
  }
  FuncPSA_df <- function(Xcateg_vec, XDatWide_df, Xc1, Xc2) {
    PSA_list <- list()
    PSA_list <- lapply(Xcateg_vec, function(x) {
      num <- 2*nrow(XDatWide_df[XDatWide_df[, Xc1]==x & XDatWide_df[, Xc2]==x,])
      denom <- sum(XDatWide_df[, Xc1]==x) + sum(XDatWide_df[, Xc2]==x)
      return(list(Grp = x, Num = num, Denom = denom, PSA = num/denom))
    })
    Out_df <- data.frame(do.call(rbind, PSA_list))
    Out_df$Grp <- as.factor(unlist(Out_df$Grp))
    Out_df$Num <- as.integer(Out_df$Num)
    Out_df$Denom <- as.integer(Out_df$Denom)
    Out_df$PSA <- round(as.numeric(Out_df$PSA)*100) 
    return(Out_df)
  }
  PSA_df <- FuncPSA_df(categ_vec, DatWide_df, c1, c2) # This is the PSA by group
  
  # Now bootstrap 90% confidence intervals
  Prop_tab <- prop.table(table(DatWide_df[, c1], DatWide_df[, c2])) # This is the matrix of proportions
  PropMelt_tab <- setNames(melt(Prop_tab), c("Grp1", "Grp2", "Prob"))
  PropMelt_tab <- PropMelt_tab[PropMelt_tab$Prob>0,]
  Boot_list <- vector("list", length = 1000)
  for(b in 1:1000) {
    PropSample_tab <- PropMelt_tab[sample(1:nrow(PropMelt_tab), size=nrow(DatWide_df), replace=TRUE, prob=PropMelt_tab$Prob), c("Grp1", "Grp2")] # A bootstrapped sample
    BootTemp_df <- FuncPSA_df(categ_vec, PropSample_tab, "Grp1", "Grp2")
    BootTemp_df$PSA[is.na(BootTemp_df$PSA)] <- 0
    Boot_list[[b]] <- BootTemp_df$PSA
  }
  Boot_df <- data.frame(do.call(cbind, Boot_list))
  BootCI_df <- as.data.frame(t(apply(Boot_df, MARGIN=1, FUN=function(x) quantile(as.numeric(x), probs=c(0.05, 0.95)))))
  PSA_df <- cbind(PSA_df, BootCI_df)
  PSA_df$`5%`[PSA_df$Denom < 5] <- NA # Don't calculate if less than 5 positive detections of the species
  PSA_df$`95%`[PSA_df$Denom < 5] <- NA
  PSA_df <- PSA_df[order(PSA_df$Grp),] 
  return(PSA_df)
}

FuncPlotPSA <- function(PSA_df, xlabel) {
  # Function to create PSA plot
  # 
  # Args:
  #   PSA_df:  Dataframe output from FuncPSA
  #   xlabel: Text for x-axis label (category)
  #
  # Returns:
  #   PSA point plot with bootstrapped 90% CI's
  PlotPSA <- ggplot(PSA_df[!is.na(PSA_df$`5%`),], aes(x=Grp,  y=PSA)) +
    geom_point(na.rm = TRUE) +
    geom_point(data = PSA_df[is.na(PSA_df$`5%`),], aes(x=Grp, y=PSA), shape=1, na.rm = TRUE) + 
    geom_linerange(aes(ymin=`5%`, ymax=`95%`), na.rm = TRUE) +
    ylim(0, 100) +
    scale_x_discrete(labels=paste0(PSA_df$Grp,"\n(",PSA_df$Denom,")")) +
    labs(y="Proportion of specific agreement (%)", x=xlabel, subtitle="(If <5 positive detections, point is open circle without a 90% CI)") +
    theme_bw(base_size = 12)
  return(PlotPSA)
}

FuncPlotThreshMetric <- function(dat_df, by_cegl, by_cegl_park = NULL, y_nam, plot_labs = list(x = NULL, y = NULL, subtitle = NULL), combo, add95) {
  # Plots of data with threshold-colored background. One metric, many subparks.
  #
  # Args:
  #   dat_df:  A data frame with the raw data
  #   by_cegl: TRUE if data are summarized by cegl. In this case, facet by SubunitCode (otherwise x= SubunitCode)
  #   by_cegl_park: If by_cegl is TRUE, then plot for only one Park (possibly with multiple SubParks)
  #   y_nam: Column name for y-axis variable
  #   plot_labs: list with x, y, subtitle for plots
  #   combo: dataframe with limits information
  #   add95: add 95% CI
  #   
  # Returns:
  #   Plots
  # 
  num_subparks <- NA
  names(dat_df)[names(dat_df)==y_nam] <- "y"
  
  if(y_nam == "avg") {
    names(dat_df)[names(dat_df)=="CI.avg_q5"] <- "low95"
    names(dat_df)[names(dat_df)=="CI.avg_q95"] <- "high95"
    names(dat_df)[names(dat_df)=="CI.avg_q16"] <- "low68"
    names(dat_df)[names(dat_df)=="CI.avg_q84"] <- "high68"
  }
  if(y_nam == "medn") {
    names(dat_df)[names(dat_df)=="CI.medn_q5"] <- "low95"
    names(dat_df)[names(dat_df)=="CI.medn_q95"] <- "high95"
    names(dat_df)[names(dat_df)=="CI.medn_q16"] <- "low68"
    names(dat_df)[names(dat_df)=="CI.medn_q84"] <- "high68"
  }
  
  if(by_cegl) {
    summary_df <- dat_df %>%
      dplyr::select(SubPark = SubunitCode, CommGroup = bin, StartYr = panel_init_year, N = countplot, y, low68, high68, low95, high95) %>%
      arrange(SubPark, CommGroup, StartYr)
    
    dat_df$binwrap = str_wrap(dat_df$bin, width = 15)
    names(dat_df)[names(dat_df)=="binwrap"] <- "x"
    num_subparks <- length(unique(dat_df$SubunitCode)) # if plot by cegl, how many subparks?
  } else {
    summary_df <- dat_df %>%
      dplyr::select(SubPark = SubunitCode, StartYr = panel_init_year, N = countplot, y, low68, high68, low95, high95) %>%
      arrange(SubPark, StartYr)
    names(dat_df)[names(dat_df)=="SubunitCode"] <- "x"
  }
  names(summary_df)[names(summary_df) == "y"] <- ifelse(y_nam=="medn", "median", "mean")
  summary_df$N <- as.integer(summary_df$N)

  axes_max = max(dat_df$high95, dat_df$y, na.rm=TRUE)
  combo$max[combo$max==99999] <- max(axes_max, 1.2*max(combo$max[combo$max!=99999]), na.rm = TRUE)

  p <- ggplot(data = dat_df, aes(x=x, y=y, fill = panel_init_year)) + 
    geom_errorbar(aes_string(ymin = "low68", ymax = "high68"), width = 0, size = 0) + # just to set the stage
    {if(by_cegl & num_subparks>1) {
      facet_grid(.~SubunitCode)
    }}
  
  # Add colored background
  if ("good_lims" %in% rownames(combo)) { 
    p <- p +
      geom_rect(aes(xmin = 0.5,
                    xmax = length(unique(dat_df$x))+0.5,
                    ymin = combo["good_lims", "min"], 
                    ymax = combo["good_lims", "max"]),
                fill="#33FF33", alpha=.3)}
  
  if ("caution_lims" %in% rownames(combo)) { 
    p <- p +
      geom_rect(aes(xmin = 0.5,
                    xmax = length(unique(dat_df$x))+0.5,
                    ymin = combo["caution_lims", "min"], 
                    ymax = combo["caution_lims", "max"]),
                fill="#FFFF00", alpha=.3)}
  if ("sigconcern_lims" %in% rownames(combo)) { 
    p <- p +
      geom_rect(aes(xmin = 0.5,
                    xmax = length(unique(dat_df$x))+0.5,
                    ymin = combo["sigconcern_lims", "min"], 
                    ymax = combo["sigconcern_lims", "max"]),
                fill="#FF3333", alpha=.3)}
  
  p <- p + 
  {if(add95)
  {geom_errorbar(aes_string(ymin = "low95", ymax = "high95"), width = 0, size = 0.25, linetype = "dashed", position = position_dodge(width = 0.5))}} +
    geom_errorbar(aes_string(ymin = "low68", ymax = "high68"), width = 0, size = 0.5, position = position_dodge(width = 0.5)) +
    geom_point(size = 4, shape = 21, position = position_dodge(width = 0.5)) +
    labs(x = plot_labs$x, y = plot_labs$y, subtitle = plot_labs$subtitle) +
    scale_fill_grey() +
    theme_bw(base_size = 11) +
    theme(legend.title = element_blank(),
          legend.position = "top",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return_list <- list(summary_df = summary_df, plot = p)
  return(return_list)
}

FuncPlotThreshPark <- function(dat_df, y_nam, plot_title, combo, add95) {
  # Plots of data with threshold-colored background. One park, many metrics.
  #
  # Args:
  #   dat_df:  A data frame with the raw data for that subpark, all metrics
  #   y_nam: Column name for y-axis variable
  #   plot_title: title for plot
  #   combo: dataframe with limits information for that subpark, all metrics
  #   add95: add 95% CI
  #   
  # Returns:
  #   List of point plots, one for each subpark in the park.
  # 
  names(dat_df)[names(dat_df)==y_nam] <- "y"

  if(y_nam == "avg") {
    names(dat_df)[names(dat_df)=="CI.avg_q5"] <- "low95"
    names(dat_df)[names(dat_df)=="CI.avg_q95"] <- "high95"
    names(dat_df)[names(dat_df)=="CI.avg_q16"] <- "low68"
    names(dat_df)[names(dat_df)=="CI.avg_q84"] <- "high68"
  }
  if(y_nam == "medn") {
    names(dat_df)[names(dat_df)=="CI.medn_q5"] <- "low95"
    names(dat_df)[names(dat_df)=="CI.medn_q95"] <- "high95"
    names(dat_df)[names(dat_df)=="CI.medn_q16"] <- "low68"
    names(dat_df)[names(dat_df)=="CI.medn_q84"] <- "high68"
  }
  
  # Summary table for output
  summary_df <- dat_df %>%
    dplyr::select(SubPark = SubunitCode, Parameter = param, StartYr = panel_init_year, N = countplot, y, low68, high68, low95, high95) %>%
    arrange(SubPark, Parameter, StartYr)
  names(summary_df)[names(summary_df) == "y"] <- ifelse(y_nam=="medn", "median", "mean")
  summary_df$N <- as.integer(summary_df$N)
  
  names(dat_df)[names(dat_df)=="SubunitCode"] <- "x"
  
  axes_max_df <- dat_df %>%
    group_by(param) %>%
    dplyr::summarize(axes_max = max(high95, dat_df$y, na.rm=TRUE))
  
  for(p in combo$param) {
    combo$max[combo$param==p & combo$max==99999] <- max(axes_max_df$axes_max[axes_max_df$param==p], 1.2*max(combo$max[combo$param==p & combo$max!=99999]), na.rm = TRUE)
  }
  
  plots_list <- list()
  for(i in unique(dat_df$x)) { # separate plot per subpark
    p <- ggplot(data = subset(dat_df, x==i), aes(x=x, y=y, fill = panel_init_year)) +
      geom_errorbar(aes_string(ymin = "low68", ymax = "high68"), width = 0, size = 0) + # just to set the stage
      facet_wrap(.~ y_short_label, scales = "free_y", nrow = 1, , labeller = labeller(y_short_label = label_wrap_gen(17)))
    
    # Add colored background
    if ("good_lims" %in% unique(combo$lim)) {
      temp_dat <- dat_df %>%
        left_join(subset(combo, lim=="good_lims", select = c(param, min, max)), by = "param") %>%
        select(y, y_short_label, x, min, max) %>%
        group_by(y_short_label, x, min, max) %>%
        dplyr::summarize(y = max(y)) %>%
        distinct()
      temp_dat <- temp_dat[complete.cases(temp_dat),]
      xmax_len <- length(unique(temp_dat$x))+0.5
      p <- p +
        geom_rect(
          data = temp_dat,
          aes(xmin = 0.5, xmax = xmax_len, ymin = min, ymax = max), fill="#33FF33", alpha = 0.9)}
    
    if ("caution_lims" %in% combo$lim) {
      temp_dat <- dat_df %>%
        left_join(subset(combo, lim=="caution_lims", select = c(param, min, max)), by = "param") %>%
        group_by(y_short_label, x, min, max) %>%
        dplyr::summarize(y = max(y)) %>%
        distinct()
      temp_dat <- temp_dat[complete.cases(temp_dat),]
      p <- p +
        geom_rect(
          data = temp_dat,
          aes(xmin = 0.5, xmax = length(unique(temp_dat$x))+0.5, ymin = min, ymax = max), fill="#FFFF00", alpha = 0.9)}
    
    if ("sigconcern_lims" %in% combo$lim) {
      temp_dat <- dat_df %>%
        left_join(subset(combo, lim=="sigconcern_lims", select = c(param, min, max)), by = "param") %>%
        group_by(y_short_label, x, min, max) %>%
        dplyr::summarize(y = max(y)) %>%
        distinct()
      temp_dat <- temp_dat[complete.cases(temp_dat),]
      p <- p +
        geom_rect(
          data = temp_dat,
          aes(xmin = 0.5, xmax = length(unique(temp_dat$x))+0.5, ymin = min, ymax = max), fill="#FF3333", alpha=0.9)}
    p <- p +
    {if(add95)
    {geom_errorbar(aes_string(ymin = "low95", ymax = "high95"), width = 0, size = 0.25, linetype = "dashed", position = position_dodge(width = 0.5))}} +
      geom_errorbar(aes_string(ymin = "low68", ymax = "high68"), width = 0, size = 0.5, position = position_dodge(width = 0.5)) +
      geom_point(size = 4, shape = 21, position = position_dodge(width = 0.5)) +
      labs(subtitle = i) +
      scale_fill_grey() +
      theme_bw(base_size = 11) +
      theme(legend.position = "top",
            legend.title = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"),
            plot.title = element_text(vjust = -9),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.text.x = element_text(angle = 90, hjust = 0),
            plot.background = element_rect(color = "black"))
    plots_list[[i]] <- p
  }
  title <- ggdraw() + draw_label(plot_title, fontface = "bold")
  if(length(plots_list) > 1) {
    x <- plot_grid(plotlist = plots_list, nrow = 1)
    p <- plot_grid(title, x, ncol = 1, rel_heights = c(0.1, 1))
  } else {
    x <- plot_grid(title, plots_list[[1]], ncol = 1, rel_heights = c(0.1, 1))
    p <- plot_grid(x, NULL, nrow = 1)
  }
  return_list <- list(summary_df = summary_df, plot = p)
  return(return_list)
}

FuncThreshSub <- function (cegl_bin, rawdata, y_label, y_short_label, field_name, cutoffs, cutoffs_one_is_sigconcern, cutoffs_high_is_good) {
  
  BootFunc <- function(dat, N) {
    if (N > 1) {
      boot_median <- boot_mean <- list()
      for (r in 1:1000) {
        boot_median[[r]] <- median(sample(x=dat, size=length(dat), replace=TRUE))
        boot_mean[[r]] <- mean(sample(x=dat, size=length(dat), replace=TRUE)) # mean for one bootstrapped sample
      } # end of bootstrap replicates
      b_medianCI <- quantile(unlist(boot_median), probs=c(0.025, 0.16, 0.84, 0.975))
      names(b_medianCI) <- c("medn_q5", "medn_q16", "medn_q84", "medn_q95")
      b_meanCI <- quantile(unlist(boot_mean), probs=c(0.025, 0.16, 0.84, 0.975)) # <<<<<<<< CHECK ON THIS--EQUIVALENT TO 68% AND 95% AROUND THE MEAN <<<<<<<<<
      names(b_meanCI) <- c("avg_q5", "avg_q16", "avg_q84", "avg_q95")
      boot <- c(b_medianCI, b_meanCI)
    } else {
      boot <- rep(NA, 8)}
    
    return(boot)
  }
  
  #standardize case (Access is not case sensitive, R is)
  rawdata$SubunitCode<-as.factor(toupper(rawdata$SubunitCode))
  rawdata$UnitCode<-as.factor(toupper(rawdata$UnitCode))
  #set the attribute as "y" to graph throughout
  rawdata$y<-rawdata[[field_name]]
  
  # if there is CEGL data, join to bins, otherwise don't
  have_cegl <- FALSE
  if ("CEGL" %in% names(rawdata)) { # if the column even exists
    if (length(rawdata$CEGL) != 0 ) { # and if there are data in it
      have_cegl <- TRUE
      ## join to CEGL / bin
      rawdata_cegl<-left_join(rawdata,cegl_bin,by="CEGL")
      
      #overwrite old with new CEGL-bin added 
      rawdata<-rawdata_cegl
    } # rawdata has CEGL data
  } # rawdata has CEGL column
  
  # remove NA and exceedingly high values - this could cause count issues later if null data exported
  rawdata<-subset(rawdata,y<9000000) %>%
    mutate_if(is.factor, as.character)
  
  ##summarize the data for some functions
  template <- rawdata %>%
    select(UnitCode,SubunitCode,panel_init_year) %>%
    distinct()
  subunit_summ <- list()
  for (i in 1:nrow(template)) {
    subdat <- subset(rawdata, UnitCode==template$UnitCode[i] & SubunitCode==template$SubunitCode[i] & panel_init_year==template$panel_init_year[i])$y
    subunit_summ[[i]] <- c(
      unlist(template[i,]),
      countplot= length(subdat),
      avg=mean(subdat),
      medn=median(subdat), 
      CI=unlist(BootFunc(subdat, N= length(subdat))))
  }
  subunit_summ_df <- as_tibble(do.call("rbind", subunit_summ))
  subunit_summ_df[, 4:14] <- sapply(subunit_summ_df[, 4:14],as.numeric)
  
  # CEGL summaries
  if(have_cegl) {
    template_cegl <- rawdata %>%
      select(UnitCode,SubunitCode,bin,panel_init_year) %>%
      distinct()
    template_cegl <- template_cegl[complete.cases(template_cegl),]
    cegl_summ <- list()
    for (i in 1:nrow(template_cegl)) {
      subdat <- subset(rawdata, UnitCode==template_cegl$UnitCode[i] & SubunitCode==template_cegl$SubunitCode[i] & bin==template_cegl$bin[i] & panel_init_year==template_cegl$panel_init_year[i])$y
      cegl_summ[[i]] <- c(
        unlist(template_cegl[i,]),
        countplot= length(subdat),
        avg=mean(subdat),
        medn=median(subdat), 
        CI=unlist(BootFunc(subdat, length(subdat))))
    }
    cegl_summ_df <- as_tibble(do.call("rbind", cegl_summ))
    cegl_summ_df[, 5:15] <- sapply(cegl_summ_df[, 5:15],as.numeric)
  } else {
    cegl_summ_df <- NA
  }
  
  if (length(cutoffs)==1) {
    # if there is only one cutoff given, cutoffs_one_is_sigconcern should be TRUE if the cutoff is between
    #    caution and significant concern, otherwise, it's between good and caution
    if (cutoffs_one_is_sigconcern==TRUE) {
      #cutoff between significant concern and caution
      good_lims<-NULL
      if (cutoffs_high_is_good == TRUE) {
        sigconcern_lims<-c(min(rawdata$y,0),cutoffs[1]) # high is good
        caution_lims<-c(cutoffs[1],99999) # high is good
      } else {
        sigconcern_lims<-c(cutoffs[1],99999) #high is bad
        caution_lims<-c(0,cutoffs[1]) # high is good
      }
    } else {  
      #cutoff between caution and good
      sigconcern_lims<-NULL
      if (cutoffs_high_is_good == TRUE) {
        caution_lims<-c(min(rawdata$y,0),cutoffs[1])
        good_lims<-c(cutoffs[1],99999) #high is good
      } else {
        #high is bad
        good_lims<-c(min(rawdata$y,0),cutoffs[1])
        caution_lims<-c(cutoffs[1],99999) #high is good
      }
    }  #which "one" is this
  } else { #two cutoffs
    
    if (cutoffs[1]<cutoffs[2]) {
      ##small is good
      good_lims<-c(min(rawdata$y,0),cutoffs[1])
      caution_lims<-cutoffs
      if (99999 < cutoffs[2]) {
        sigconcern_lims<-NULL ##c(cutoffs[2],cutoffs[2])
      } else {
        # sigconcern_lims<-c(cutoffs[2],max(ylim_upper,ymax,99999))
        sigconcern_lims<-c(cutoffs[2],99999)
      }
    } else {
      ##small is bad
      sigconcern_lims<-c(min(rawdata$y,0),cutoffs[2])
      caution_lims<-c(cutoffs[2],cutoffs[1])
      
      if (99999 < cutoffs[1]) {
        good_lims<-NULL # c(cutoffs[1],cutoffs[1])
      } else {
        good_lims<-c(cutoffs[1],99999)
      }
    } }
  
  ##combine labels into combo object
  combo <- NULL
  combo_colors <- NULL
  combo_text <-NULL
  if (length(good_lims)!=0) { 
    combo <- rbind(combo,good_lims)
    combo_colors <- c("#33FF33")
    combo_text <-c("Good")
    
  }
  if (length(caution_lims)!=0) {
    combo <- rbind(combo,caution_lims)
    combo_colors <- c(combo_colors,c("#FFFF00"))
    combo_text <- c(combo_text ,c("Caution"))
    
  }
  if (length(sigconcern_lims)!=0) {
    combo <- rbind(combo,sigconcern_lims)
    combo_colors <- c(combo_colors,c("#FF3333"))
    combo_text <- c(combo_text ,c("Significant\nConcern"))
    
  }
  
  combo<-as.data.frame(combo)
  
  names(combo)<-c("min","max")
  combo$color <- combo_colors  
  combo$text <-combo_text
  
  return_list <- list(
    subunit_summ_df = subunit_summ_df,
    cegl_summ_df = cegl_summ_df,
    combo = combo,
    y_label = y_label,
    y_short_label = tools::toTitleCase(y_short_label),
    cutoffs_one_is_sigconcern = cutoffs_one_is_sigconcern,
    cutoffs_high_is_good = cutoffs_high_is_good)
  return(return_list)
  
}

### MAIN FUNCTIONS -----
FuncThreshMain <- function() {
  thresh_list <- list()
  cegl_bin <- read_csv(here::here("Data_in", "cegl_x_bin.csv")) %>%
    select(CEGL, bin = gp1)
  
  # One chart per park-metric, x-axis is year ----
  #cwd (coarse woody debris)
  cwd<-read_csv(here::here("Data_in", "qry_stat_cwd_byPanelYear.csv"))
  thresh_list$cwd <- FuncThreshSub(cegl_bin, cwd,"cwdPercentOfLiveVol","Coarse Woody Debris","cwdPercentOfLiveVol", c(15,5), cutoffs_one_is_sigconcern=NULL, cutoffs_high_is_good=NULL)
  
  #snag LARGE DENSITY
  snag<-read_csv(here::here("Data_in", "qry_snagDensity_event_panel.csv"))
  thresh_list$snag <- FuncThreshSub(cegl_bin, snag,"plot_large_Snag_per_ha 25cm","Large snag density","plot_large25_Snag_per_ha", c(8,3),cutoffs_one_is_sigconcern=NULL, cutoffs_high_is_good=NULL)
  
  #HPE species count
  hpe2<-read_csv(here::here("Data_in", "qry_event_HPExotic_count_byPanelYear.csv"))
  thresh_list$hpe2 <- FuncThreshSub(cegl_bin, hpe2,"HPexotic_species_present_per_plot_v2","Exotic species #","HPexotic_species_present",c(0.5,3.5),cutoffs_one_is_sigconcern=NULL, cutoffs_high_is_good=NULL)
  
  #tree mortality and growth rate
  growth_mort<-read_csv(here::here("Data_in", "qry_tree_eventGrowthMort_panel.csv"))
  growth_mort$UnitCode<-growth_mort$Unit_Code
  
  thresh_list$tree_growth_mort <- FuncThreshSub(cegl_bin, growth_mort,"Class_2-4_Annual_Tree_Mortality_Percent","Tree Mortality","AnnualEventMortalityPercent",c(1.6),cutoffs_one_is_sigconcern=FALSE, cutoffs_high_is_good=FALSE)
  
  growth_mort$relGrowthRateAsPercentOfRegion <- 100*growth_mort$EventAvgRelBAGrowth/ mean(growth_mort$EventAvgRelBAGrowth)
  thresh_list$tree_growth_rate <- FuncThreshSub(cegl_bin, growth_mort,"Class_2-4 Tree Growth Rate as Percent of Region Average","Tree Growth Rate","relGrowthRateAsPercentOfRegion", c(60),cutoffs_one_is_sigconcern=FALSE, cutoffs_high_is_good=TRUE)
  saveRDS(thresh_list, "Temp_out/ThreshList.RDS")
}

FuncMapDatMain <- function () {
  # Function to run map summaries on vegetation data
  MapOut_list <- list()
  
  VegMap_df <- FuncFormatVeg()
  CycLevs <- sort(unique(VegMap_df$Cycle), decreasing = FALSE)
  TempLevs <- as.numeric(substring(CycLevs, 7, 10))
  CycLabels <- paste0("[", TempLevs - 4, ",", TempLevs, "]")
  
  RC_df <- VegMap_df %>%
    dplyr::mutate(Type = ifelse(ExoticPlant, "Exotic", "Native")) %>%
    filter(QAQCPlot==FALSE) %>% # non-QAQC only <<<<<<<
    dplyr::select(Park, SubPark, Longitude, Latitude, Plot, Cycle, Year, Type, HighPriorityExotic, PlantName, PlantConfirmed, CoverPercent) 
  
  # Top 3 cover exotics per SubPark-Cycle
  Exotics_tab <- RC_df %>% # Identify the exotics found on the most plots
  dplyr::filter(Type=="Exotic") %>%
    group_by(SubPark, Cycle, PlantName) %>%
    dplyr::summarise(Freq = n(),
                     Cov = sum(CoverPercent))
  HiCov3 <- Exotics_tab %>% 
    group_by(SubPark, Cycle) %>%
    top_n(n = 3, wt = Cov) %>%
    dplyr::mutate(rank = rank(desc(Cov), ties.method = "first")) %>%
    dplyr::select(-Freq, -Cov) 
  toprank <- max(HiCov3$rank)
  
  HiCov3 %<>% 
    spread(key = rank, value = PlantName) 
  HiCov3$Top3 <- apply( HiCov3[ , 3:ncol(HiCov3)] , 1 , paste , collapse = ", " )
  HiExoticCov3_tab <- subset(HiCov3, select=c(SubPark, Cycle, Top3))
  HiExoticCov3_tab$Top3 <- gsub("_", " ", HiExoticCov3_tab$Top3)
  HiExoticCov3_tab$Top3 <- gsub(", NA", "", HiExoticCov3_tab$Top3) 
      
  # Top 3 frequency exotics per SubPark-Cycle
  HiFreq3 <- Exotics_tab %>% 
    group_by(SubPark, Cycle) %>%
    top_n(n = 3, wt = Freq) %>%
    dplyr::mutate(rank = rank(desc(Freq), ties.method = "first")) %>%
    dplyr::select(-Freq, -Cov) 
  toprank <- max(HiFreq3$rank)
  
  HiFreq3 %<>% 
    spread(key = rank, value = PlantName) 
  HiFreq3$Top3 <- apply( HiFreq3[ , 3:ncol(HiFreq3)] , 1 , paste , collapse = ", " )
  HiExoticFreq3_tab <- subset(HiFreq3, select=c(SubPark, Cycle, Top3))
  HiExoticFreq3_tab$Top3 <- gsub("_", " ", HiExoticFreq3_tab$Top3)
  HiExoticFreq3_tab$Top3 <- gsub(", NA", "", HiExoticFreq3_tab$Top3) 
  
  # Exotic cover
  RC_Template_df <- RC_df %>%
    dplyr::select(Park, SubPark, Longitude, Latitude, Plot, Cycle, Year, CoverPercent) %>%
    group_by(Park, SubPark, Longitude, Latitude, Plot, Cycle, Year) %>%
    dplyr::summarize(TotCover = sum(CoverPercent))
    
    RC_Template2_df <- rbind(cbind(RC_Template_df, Type=rep("Native", nrow(RC_Template_df))), cbind(RC_Template_df, Type=rep("Exotic", nrow(RC_Template_df))))
  
  RC_NativeExotic_df <- RC_df %>%
    dplyr::select(Plot, Cycle, Year, Type, PlantConfirmed, CoverPercent) %>%
    group_by(Plot, Cycle, Year, Type) %>%
    dplyr::mutate(SpeciesRichness = sum(PlantConfirmed),
           PercCover = sum(CoverPercent)) %>%
    ungroup() %>%
    dplyr::select(-CoverPercent, -PlantConfirmed) %>%
    distinct() %>%
    right_join(RC_Template2_df, by=c("Plot", "Type", "Cycle", "Year"))
  RC_NativeExotic_df[, c("SpeciesRichness", "PercCover")][is.na(RC_NativeExotic_df[, c("SpeciesRichness", "PercCover")])] <- 0
  RC_NativeExotic_df$RelCover = round((RC_NativeExotic_df$PercCover/RC_NativeExotic_df$TotCover)*100, 2)
  
  Exotics_popup <- RC_df %>%
    filter(Type == "Exotic") %>%
    dplyr::select(SubPark, Plot, Cycle, Year, PlantName, CoverPercent) %>% 
    left_join(subset(RC_NativeExotic_df, Type=="Exotic", select = c(SubPark, Plot, Cycle, Year, TotCover)), by = c("SubPark", "Plot", "Cycle", "Year")) %>%
    mutate(RelCover = round((CoverPercent/TotCover)*100, 1)) %>%
    dplyr::arrange(desc(CoverPercent)) %>%
    group_by(SubPark, Plot, Cycle, Year) %>% 
    dplyr::mutate(SumExoticCov = sum(RelCover),
      Poptext = paste0(Plot, " (", Year, ")", "<br/>", paste0(SumExoticCov, "% EXOTIC COVER"), "<br>-----<br>", paste0(PlantName, " (", RelCover, "%)", collapse = "<br>"))) %>% 
    dplyr::select(SubPark, Plot, Cycle, Year, Poptext)
  Exotics_popup$Poptext <- gsub("_", " ", Exotics_popup$Poptext)
  Exotics_popup$Poptext <- gsub("(0%)", "trace", Exotics_popup$Poptext)
  
  ExoticCov_df <- RC_NativeExotic_df %>%
    filter(Type=="Exotic") %>%
    left_join(Exotics_popup, by = c("SubPark", "Plot", "Cycle", "Year")) %>%
    mutate(Poptext = ifelse(is.na(Poptext), paste0(Plot, " (", Year, ")", "<br/>", "0% EXOTIC COVER"), Poptext)) %>%
    dplyr::select(SubPark, Cycle, Plot, Longitude, Latitude, RelCover, Poptext) %>%
    dplyr::rename(Metric = RelCover)
  
  # Exotic richness
  ExoticRich_df <- RC_NativeExotic_df %>%
    filter(Type=="Exotic") %>%
    dplyr::mutate(Poptext = paste0(Plot, " (", Year, "):", "<br/>", SpeciesRichness, " exotic herb species")) %>%
    dplyr::select(SubPark, Cycle, Plot, Longitude, Latitude, SpeciesRichness, Poptext) %>%
    dplyr::rename(Metric = SpeciesRichness)
  
  # Native richness
  NativeRich_df <- RC_NativeExotic_df %>%
    filter(Type=="Native") %>%
    dplyr::mutate(Poptext = paste0(Plot, " (", Year, "):", "<br/>", SpeciesRichness, " native herb species")) %>%
    dplyr::select(SubPark, Cycle, Plot, Longitude, Latitude, SpeciesRichness, Poptext) %>%
    dplyr::rename(Metric = SpeciesRichness)
  
  # Cover for high priority exotics
  RC_HPtemp_df <- VegMap_df %>%
    filter(QAQCPlot==FALSE & HighPriorityExotic) %>% # non-QAQC only <<<<<<<
    dplyr::select(Plot, Cycle, Year, PlantName, CoverPercent)

  RC_TemplateHP_df <- merge(RC_Template_df, unique(RC_HPtemp_df$PlantName)) %>%
    dplyr::mutate(y = as.character(y))
  names(RC_TemplateHP_df)[names(RC_TemplateHP_df)=="y"] <- "PlantName"
  
  RC_HP_df <- RC_HPtemp_df %>%
    right_join(RC_TemplateHP_df, by=c("Plot", "PlantName", "Cycle", "Year"))
  RC_HP_df[, c("CoverPercent")][is.na(RC_HP_df[, c("CoverPercent")])] <- 0
  
  RC_HP_df %<>%
    dplyr::mutate(RelCover = round((CoverPercent/TotCover)*100, 2),
                  SpeciesRichness = as.numeric(RelCover > 0),
                  Poptext = paste0(Plot, " (", Year, "):", "<br/>", round(RelCover, 1), "% cover")) %>%
    dplyr::select(SubPark, Cycle, Plot, Longitude, Latitude, PlantName, RelCover, Poptext) %>%
    dplyr::rename(Metric = RelCover)
  RC_HP_df$PlantName <- gsub("_", " ", RC_HP_df$PlantName)
  
  MapOut_list <- list(
    HiExoticCov3 = HiExoticCov3_tab,
    HiExoticFreq3 = HiExoticFreq3_tab,
    ExoticRich = ExoticRich_df,
    ExoticCov = ExoticCov_df, 
    HPExoticCov = RC_HP_df,
    NativeRich = NativeRich_df)
  saveRDS(MapOut_list, "Temp_out/MapDatOut.RDS")
}

FuncQAQCSitesMain <- function() {
  # Function to run site QAQC analyses.
  QAQCSitesOut_list <- list()
  
  # Stand height -----
  SiteTemp_df <- read_csv(here::here("Data_in", "SiteMetrics.csv"))
  
  Site_df <- SiteTemp_df %>%
    dplyr::rename(Plot = location_name, Year = event_year, StandHeight = Stand_Height) %>%
    dplyr::select(Plot, QAQC, Year, StandHeight) %>%
    spread(key = QAQC, value = StandHeight) %>%
    dplyr::rename(Nsurv = `-1`, Qsurv = `0`) %>%
    dplyr::mutate(PlotYear = paste0(Plot, "_", Year)) # In case a single plot has QAQC survey in multiple survey cycles
  Site_df <- Site_df[complete.cases(Site_df), ]
  
  Site_df %<>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Diff = Qsurv - Nsurv, # negative value means Qsurv smaller
           Metric = round((Diff/Mean)*100),
           LT5perc = abs(Metric) < 5,
           LT10perc = abs(Metric) < 10)
  
  QAQCSitesOut_list[["StandHeight"]][["FinalDat_df"]] <- Site_df
  QAQCSitesOut_list[["StandHeight"]][["plot_vals"]] <- list(type = "continuous", lab = "% difference", v1 = 5, v2 = 10, range_max = 5*ceiling(max(abs(Site_df$Metric))/5) + 5, range_int = 5)
  
  StandHeightTemp_tab <- round(ecdf(abs(Site_df$Metric))(seq(0, as.numeric(QAQCSitesOut_list[["StandHeight"]][["plot_vals"]]["range_max"]), by = 5))*100)
  StandHeightTemp_tab <- data.frame(cbind(Diff = paste0("< = ", seq(0, as.numeric(QAQCSitesOut_list[["StandHeight"]][["plot_vals"]]["range_max"]), by = 5), "%"), CumPercStands = StandHeightTemp_tab))
  QAQCSitesOut_list[["StandHeight"]][["Diff_tab"]] <- StandHeightTemp_tab
  
  # Coarse woody debris counts by Plot -----
  CWD_df <- read_csv(here::here("Data_in", "CWD.csv"))
  # TRANSECT INFORMATION (X OR Y) MISSING FOR A LOOK007 2012 NON-QAQC RECORD, SO FOR THE TRANSECT-SPECIFIC ANALYSIS IT WAS COUNTED AS A MISMATCH <<<<<<<<
  CWD_df %<>%
    dplyr::select(Plot = location_name, Year = event_year, QAQC, Transect, Diameter, Decay = Decay_Class) %>%
    dplyr::mutate(PlotYear = paste0(Plot, "_", Year)) %>%
    dplyr::select(-Plot, -Year)
  
  TempCWD_df <- expand.grid(PlotYear = unique(CWD_df$PlotYear), QAQC = unique(CWD_df$QAQC), Transect = c("X", "Y"), stringsAsFactors = FALSE)
  CWD_df %<>%
    full_join(TempCWD_df, by=c("PlotYear", "QAQC", "Transect"))
  
  CWD_ByPlot_df <- CWD_df %>%
    group_by(PlotYear, QAQC) %>%
    dplyr::summarize(N = n()) %>%
    spread(key = QAQC, value = N) %>%
    dplyr::rename(Qsurv = `-1`, Nsurv = `0`) 
  
  CWD_ByPlot_df %<>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Metric = Qsurv - Nsurv,
           MatchBy1 = abs(Metric) <= 1,
           MatchBy2 = abs(Metric) <= 2)
  
  QAQCSitesOut_list[["CWDPlotCounts"]][["CWD_orig_df"]] <- CWD_df
  QAQCSitesOut_list[["CWDPlotCounts"]][["FinalDat_df"]] <- as.data.frame(CWD_ByPlot_df)
  QAQCSitesOut_list[["CWDPlotCounts"]][["plot_vals"]] <- list(type = "discrete", lab = "Absolute difference", v1 = 1, v2 = 2, range_max = (max(abs(CWD_ByPlot_df$Metric))) + 1, range_int = 1)
  QAQCSitesOut_list[["CWDPlotCounts"]][["Diff_tab"]] <- 
    table(abs(CWD_ByPlot_df$Metric))
  
  # Coarse woody debris counts by Transect -----
  CWD_ByTransect_df <- CWD_df %>%
    group_by(PlotYear, QAQC, Transect) %>%
    dplyr::summarize(N = n()) %>%
    spread(key = QAQC, value = N) %>%
    dplyr::rename(Qsurv = `-1`, Nsurv = `0`) %>%
    filter(!is.na(Transect)) # if there is no transect name, remove it because we don't know what transect it should be assigned to
  
  CWD_ByTransect_df %<>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Metric = Qsurv - Nsurv,
           MatchBy1 = abs(Metric) <= 1,
           MatchBy2 = abs(Metric) <= 2)
  
  QAQCSitesOut_list[["CWDTransectCounts"]][["FinalDat_df"]] <- as.data.frame(CWD_ByTransect_df)
  QAQCSitesOut_list[["CWDTransectCounts"]][["plot_vals"]] <- list(type = "discrete", lab = "Absolute difference", v1 = 1, v2 = 2, range_max = (max(abs(CWD_ByTransect_df$Metric))) + 1, range_int = 1)
  QAQCSitesOut_list[["CWDTransectCounts"]][["Diff_tab"]] <- 
    table(abs(CWD_ByTransect_df$Metric))
  
  
  # COMMUNITY GROUP (OWN DASHBOARD TAB) -----
  # REQUIRED FIELD; NA/NA IS A MISMATCH
  CG <- read_csv(here::here("Data_in", "Community.csv"))
  CG$QAQC <- as.logical(CG$QAQC) # Change to T/F
  CG %<>%
    dplyr::rename(Plot = location_name, Year=event_year, CC = Community_Code, CG  = CommGroup) %>% # use standard variable names
    dplyr::mutate(PlotYear = paste0(Plot, "-", Year)) %>% # In case a single plot has QAQC survey in multiple survey cycles
    dplyr::select(PlotYear, QAQC, CG, CC, Fit)
  
  # ERROR CHECKING
  noCGmatch <- names(table(CG$PlotYear)[table(CG$PlotYear)!=2])
  if(length(noCGmatch)>0) cat("\n********** BEGIN ERROR CHECKING **********\nThese survey events do not have exactly two entries (QAQC & non-QAQC) for evaluating Community Group match:\n\n ", noCGmatch, "\n\n********** END ERROR CHECKING **********\n")
  
  # Format the data
  CG <- droplevels(CG)
  all_CG <- sort(unique(CG$CG)) # vector of CG codes
  CCWide_df <- dcast(CG, PlotYear ~ QAQC, value.var="CC")
  names(CCWide_df)[2:3] <- c("CC1", "CC2")
  CCWide_df$CCMatch <- CCWide_df$CC1==CCWide_df$CC2
  CCWide_df$CCMatch[is.na(CCWide_df$CCMatch)] <- FALSE # If at least one of the CC values was NA (not recorded), then Match = FALSE
  
  CGWide_df <- dcast(CG, PlotYear ~ QAQC, value.var="CG")
  names(CGWide_df)[2:3] <- c("CG1", "CG2")
  CGWide_df$CGMatch <- CGWide_df$CG1==CGWide_df$CG2
  CGWide_df$CGMatch[is.na(CGWide_df$CGMatch)] <- FALSE # If one of the CG values was NA (not recorded), then Match = FALSE
  
  FitWide_df <- dcast(CG, PlotYear ~ QAQC, value.var="Fit")
  names(FitWide_df)[2:3] <- c("Fit1", "Fit2")
  FitWide_df[,2:3] <- FitWide_df[,2:3]==1 # TRUE for Fit=1, FALSE for Fit != 1
  FitWide_df[,2:3][is.na(FitWide_df[,2:3])] <- FALSE # If one of the Fit values was NA (not recorded), then FALSE for Fit (i.e., means Fit is not 1)
  FitWide_df$FitMatch <- ifelse(FitWide_df$Fit1 & FitWide_df$Fit2, "11", ifelse(!FitWide_df$Fit1 & !FitWide_df$Fit2, "++", "1+/+1")) # FitMatch of '11' means both were Fit=1; '++' means both were Fit>1; otherwise '1+/+1'
  
  # Summary table
  GroupSummary <- merge(CGWide_df[, c("PlotYear", "CGMatch")], CCWide_df[, c("PlotYear", "CCMatch")], by="PlotYear", all=TRUE)
  GroupSummary <- merge(GroupSummary, FitWide_df[, c("PlotYear", "FitMatch")], by="PlotYear", all=TRUE)
  GroupSummary$FitMatch <- factor(GroupSummary$FitMatch, levels=c("11", "1+/+1", "++"))
  
  # Proportion of specific agreement by CG
  # For each CG, P(CG) is (2*agreement)/(number of cases in which at least one called it)
  CG_PSA_df <- FuncPSA(categ_vec = all_CG, DatWide_df = CGWide_df, c1="CG1", c2="CG2", missing_val = TRUE)
  
  QAQCSitesOut_list[["CG_list"]] <- list(
    all_CG = all_CG, # vector of CG codes
    GroupSummary = GroupSummary, # for summary tables
    CGWide_df = CGWide_df, # for color matrix
    CG_PSA_df = CG_PSA_df)
  
  # SITE DISTURBANCE (OWN DASHBOARD TAB) -----
  DistCol <- sort(c("Disturbance_Open_Grown","Disturbance_Homesite","Disturbance_Fire_Recent","Disturbance_Fire_Historic","Disturbance_Deer_Browse","Disturbance_Logging","Disturbance_Plowing","Disturbance_Pasture","Disturbance_Erosion"))  # The names of disturbance columns      
  
  # Extract the relevant columns and format
  SiteDisturbSummary_list <- SiteDisturbPSA_list <- list()
  for (d in DistCol) {
    SD_df <- SiteTemp_df %>%
      dplyr::select(Plot = location_name, QAQC, Year = event_year, d) %>%
      spread(key = QAQC, value = d) %>%
      dplyr::rename("Value1" = `-1`, "Value2" = `0`)
    
    # ERROR CHECKING <<<<<<<<<<<<<<<<<<<<
    NoDistMatch <- SD_df[is.na(SD_df$Value1)|is.na(SD_df$Value2), c("Plot", "Year")]
    if(nrow(NoDistMatch)>0) cat("\n********** BEGIN ERROR CHECKING **********\nThese site-years are missing an entry for the disturbance. These will automatically be counted as a mismatch (even if both groups had NA) for the disturbance ", d, "\n\n********** END ERROR CHECKING **********\n")
    SD_df$MatchType <- ifelse(SD_df$Value1 & SD_df$Value2, "11", ifelse(!SD_df$Value1 & !SD_df$Value2, "00", "10/01")) # 1 = TRUE, 0 = FALSE
    SD_df$MatchType[is.na(SD_df$Value1)|is.na(SD_df$Value2)] <- "10/01" # if any of the entries is NA, count as mismatch
    TempSDMatchCount <- SD_df %>% dplyr::count(MatchType)
    SiteDisturbSummary_list[[d]] <- data.frame(SD = d, TempSDMatchCount)
    SD_df <- as.data.frame(SD_df)
    SiteDisturbPSA_list[[d]] <- FuncPSA(categ_vec = c(TRUE, FALSE), DatWide_df = SD_df, c1 = "Value1", c2 = "Value2", missing_val = TRUE)
  }
  
  # Summarize lists in data frames
  SiteDisturbSummary_df <- do.call("rbind", SiteDisturbSummary_list)
  SiteDisturbSummary_df <- spread(SiteDisturbSummary_df, key=MatchType, value=n)
  SiteDisturbSummary_df[is.na(SiteDisturbSummary_df)] <- 0
  SiteDisturbSummary_df$SD <- as.character(SiteDisturbSummary_df$SD)
  
  SiteDisturbPSA_df <- do.call("rbind", SiteDisturbPSA_list)
  SiteDisturbPSA_df %<>% 
    rownames_to_column(var="SD") %>%
    filter(!is.nan(PSA) & Grp %in% c("TRUE","XXX")) %>%
    dplyr::select(SD, PSA, Denom, `5%`, `95%`) %>%
    dplyr::mutate(SD=str_replace_all(SD, ".[1-9]", ""))
  
  SiteDisturbPSA_df$PSA[is.na(SiteDisturbPSA_df$PSA)] <- 0
  
  SiteDisturbPSASummary_df <- data.frame(SD=DistCol, stringsAsFactors = FALSE) %>%
    left_join(SiteDisturbSummary_df, by="SD") %>%
    full_join(SiteDisturbPSA_df, by="SD") %>%
    dplyr::rename("FALSE-FALSE"=`00`, "TRUE-FALSE"=`10/01`, "TRUE-TRUE"=`11`) %>%
    dplyr::mutate(SD=str_replace_all(SD, "Disturbance_", ""))
  
  SiteDisturbPSASummary_df$SD <- str_replace_all(SiteDisturbPSASummary_df$SD, "Disturbance_", "")
  SiteDisturbPSASummary_df$SD <- factor(SiteDisturbPSASummary_df$SD)
  SiteDisturbPSASummary_df <- droplevels(SiteDisturbPSASummary_df)
  
  # This is the summary PSA plot
  SiteDisturbPSA_plot <- ggplot(SiteDisturbPSASummary_df, aes(x=SD,  y=PSA)) +
    geom_point(na.rm = TRUE) +
    geom_linerange(aes(ymin=`5%`, ymax=`95%`), na.rm = TRUE) +
    ylim(0, 100) +
    scale_x_discrete(labels=paste0(SiteDisturbPSASummary_df$SD, "\n(", SiteDisturbPSASummary_df$Denom, ")")) +
    labs(y="Proportion of specific agreement (%)", x="Site Disturbance", subtitle="(90% CI included for N > 5)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = .6))
  
  SiteDisturbPSASummary_df$Denom <- NULL
  
  QAQCSitesOut_list[["SiteDisturb_list"]] <- list(
    SiteDisturbPSASummary_df = SiteDisturbPSASummary_df,
    SiteDisturbPSA_plot = SiteDisturbPSA_plot)
  
  saveRDS(QAQCSitesOut_list, "Temp_out/QAQCSitesOut.RDS")
}

FuncQAQCTreesMain <- function() {
  # Function to run tree QAQC analyses.
  TreesN_df <- read_csv(here::here("Data_in", "TreesNorm.csv"))
  TreesQ_df <- read_csv(here::here("Data_in", "TreesQAQC.csv"))

  # Data formatting -----
  names(TreesN_df) <- gsub("Condition_Damage_Code_", "Damage", names(TreesN_df))
  names(TreesQ_df) <- gsub("Condition_Damage_Code_", "Damage", names(TreesQ_df))
  
  QAQCTreesDat_list <- list()
  
  TreesN_df %<>% dplyr::mutate(
    Foliage1=paste0(Foliage_Condition_1,".",Foliage_Percent_Class_1),
    Foliage2=paste0(Foliage_Condition_2,".",Foliage_Percent_Class_2),
    Foliage3=paste0(Foliage_Condition_3,".",Foliage_Percent_Class_3),
    Foliage4=paste0(Foliage_Condition_4,".",Foliage_Percent_Class_4),
    Foliage5=paste0(Foliage_Condition_5,".",Foliage_Percent_Class_5))
  TreesQ_df %<>% dplyr::mutate(
    Foliage1=paste0(Foliage_Condition_1,".",Foliage_Percent_Class_1),
    Foliage2=paste0(Foliage_Condition_2,".",Foliage_Percent_Class_2),
    Foliage3=paste0(Foliage_Condition_3,".",Foliage_Percent_Class_3),
    Foliage4=paste0(Foliage_Condition_4,".",Foliage_Percent_Class_4),
    Foliage5=paste0(Foliage_Condition_5,".",Foliage_Percent_Class_5))

  TreesN_df %<>% dplyr::select(TreeID = Tree_ID, DBH, Status = Tree_Status, Crown = Crown_Class, Vigor, Dieback, Damage1:Damage5, Foliage1:Foliage5, SnagHeight = Snag_Height, SnagDecay = Dead_Tree_Decay_Class, Arborist = ArboristName, QTreeID = qaqcPlotTree_ID)
  TreesQ_df %<>% dplyr::select(TreeID = Tree_ID, DBH, Status = Tree_Status, Crown = Crown_Class, Vigor, Dieback, Damage1:Damage5, Foliage1:Foliage5, SnagHeight = Snag_Height, SnagDecay = Dead_Tree_Decay_Class, Arborist = ArboristName, NTreeID = normalPlotTree_ID)
  
  # Warnings -----
  # QAQC - Trees
  WarnQAQCTrees_list <- sapply(c("NoQAQC", "NoNormal"), function(x) NULL) # <<<<<<< FILL IN MORE
  WarnQAQCTrees_list$NoQAQC <- setdiff(TreesQ_df$NTreeID, TreesN_df$TreeID) # >>>>>> ERROR message to use: These individual trees do not have ANY corresponding non-QAQC data for evaluating match in any metric. They will be deleted.
  QAQCTreesDat_list$TreesQ <- TreesQ_df[!TreesQ_df$NTreeID %in% WarnQAQCTrees_list$NoQAQC, ]; QAQCTreesDat_list$TreesQ <- droplevels(QAQCTreesDat_list$TreesQ)
  
  WarnQAQCTrees_list$NoNormal <- setdiff(TreesN_df$QTreeID, TreesQ_df$TreeID) # >>>>>> ERROR message to use: These individual trees do not have ANY corresponding QAQC data for evaluating match in any metric. They will be deleted.
  QAQCTreesDat_list$TreesN <- TreesN_df[!TreesN_df$QTreeID %in% WarnQAQCTrees_list$NoNormal, ]; QAQCTreesDat_list$TreesN <- droplevels(QAQCTreesDat_list$TreesN)
  
  # These are the QAQC trees recorded as Status=1 in both surveys
  Stat1Trees_df <- subset(QAQCTreesDat_list$TreesN, Status==1, select=c("QTreeID")) %>%
    inner_join(subset(QAQCTreesDat_list$TreesQ, Status==1, select=c("TreeID")), by=c("QTreeID"="TreeID")) # these are the trees with Status=1 in both surveys
  QAQCTreesDat_list$TreesNSub1 <- subset(QAQCTreesDat_list$TreesN, QTreeID %in% Stat1Trees_df$QTreeID) 
  QAQCTreesDat_list$TreesQSub1 <- subset(QAQCTreesDat_list$TreesQ, TreeID %in% Stat1Trees_df$QTreeID) 
  saveRDS(QAQCTreesDat_list, "Temp_out/QAQCTreesDat.RDS")
  
  QAQCTreesOut_list <- list()
  
  # START WITH THE CATEGORICAL TREE DATA (OWN DASHBOARD TAB) -----
  # Tree status -----
  QAQCTreesOut_list[["Status"]][["FinalDat_df"]] <- FuncQAQCTreesFormat(NDat_df = QAQCTreesDat_list$TreesN, QDat_df = QAQCTreesDat_list$TreesQ, cname = "Status", NANA_mismatch = TRUE) 
  QAQCTreesOut_list[["Status"]][["categ_vec"]] <- 0:4
  QAQCTreesOut_list[["Status"]][["PSA_df"]] <- FuncPSA(categ_vec = 0:4, DatWide_df = QAQCTreesOut_list[["Status"]][["FinalDat_df"]], c1 = "Value1", c2 = "Value2", missing_val = TRUE)
  
  # Tree crown -----
  QAQCTreesOut_list[["Crown"]][["FinalDat_df"]] <- FuncQAQCTreesFormat(NDat_df = QAQCTreesDat_list$TreesNSub1, QDat_df = QAQCTreesDat_list$TreesQSub1, cname="Crown", NANA_mismatch = TRUE) 
  QAQCTreesOut_list[["Crown"]][["categ_vec"]] <- 1:5
  QAQCTreesOut_list[["Crown"]][["PSA_df"]] <- FuncPSA(categ_vec = 1:5, DatWide_df = QAQCTreesOut_list[["Crown"]][["FinalDat_df"]], c1 = "Value1", c2 = "Value2", missing_val = TRUE)
  
  # Tree vigor -----
  QAQCTreesOut_list[["Vigor"]][["FinalDat_df"]] <- FuncQAQCTreesFormat(NDat_df = QAQCTreesDat_list$TreesNSub1, QDat_df = QAQCTreesDat_list$TreesQSub1, cname="Vigor", NANA_mismatch = TRUE) 
  QAQCTreesOut_list[["Vigor"]][["categ_vec"]] <- 1:5
  QAQCTreesOut_list[["Vigor"]][["PSA_df"]] <- FuncPSA(categ_vec = 1:5, DatWide_df = QAQCTreesOut_list[["Vigor"]][["FinalDat_df"]], c1 = "Value1", c2 = "Value2", missing_val = TRUE)
  
  # Tree dieback -----
  QAQCTreesOut_list[["Dieback"]][["FinalDat_df"]] <- FuncQAQCTreesFormat(NDat_df = QAQCTreesDat_list$TreesNSub1, QDat_df = QAQCTreesDat_list$TreesQSub1, cname = "Dieback", NANA_mismatch = FALSE) # This one is false for NANA_mismatch because Dieback is not a required field 
  QAQCTreesOut_list[["Dieback"]][["categ_vec"]] <- 1:5
  QAQCTreesOut_list[["Dieback"]][["PSA_df"]] <- FuncPSA(categ_vec = 1:5, DatWide_df = QAQCTreesOut_list[["Dieback"]][["FinalDat_df"]], c1 = "Value1", c2 = "Value2", missing_val = TRUE)
  
  # Snag decay -----
  # These are the QAQC trees recorded as Status=2 in both surveys
  Stat2Trees_df <- subset(QAQCTreesDat_list$TreesN, Status==2, select = c("QTreeID")) %>%
    inner_join(subset(QAQCTreesDat_list$TreesQ, Status==2, select = c("TreeID")), by = c("QTreeID"="TreeID")) # these are the trees with Status=2 in both surveys
  TreesNSub2 <- subset(QAQCTreesDat_list$TreesN, QTreeID %in% Stat2Trees_df$QTreeID) 
  TreesQSub2 <- subset(QAQCTreesDat_list$TreesQ, TreeID %in% Stat2Trees_df$QTreeID) 
  
  QAQCTreesOut_list[["SnagDecay"]][["FinalDat_df"]] <- FuncQAQCTreesFormat(NDat_df = TreesNSub2, QDat_df = TreesQSub2, cname = "SnagDecay", NANA_mismatch = TRUE) 
  QAQCTreesOut_list[["SnagDecay"]][["categ_vec"]] <- 1:5
  QAQCTreesOut_list[["SnagDecay"]][["PSA_df"]] <- FuncPSA(categ_vec = 1:5, DatWide_df = QAQCTreesOut_list[["SnagDecay"]][["FinalDat_df"]], c1 = "Value1", c2 = "Value2", missing_val = TRUE)
  
  # NOW THE COUNT/CONTINUOUS TREE DATA (OWN DASHBOARD TAB) -----
  
  # Tree counts -----
  # Used data in 'all_trees', regardless of status
  
  TreeCounts_df <- read_csv(here::here("Data_in", "AllTrees.csv"))
  # >>>>>>>> ERROR: Warning: 1 parsing failure.
  #row                 col           expected actual                   file
  #1453 Foliage_Condition_3 1/0/T/F/TRUE/FALSE      N 'Data_in/AllTrees.csv'
  TreeCounts_df %<>%
    dplyr::rename(Plot = location_name, Year = event_year, QAQC = qaqc, TreeID=Tree_Number) %>%
    dplyr::mutate(PlotYear = paste0(Plot, "_", Year)) %>% # In case a single plot has QAQC survey in multiple survey cycles
    dplyr::select(PlotYear, QAQC, Module, TreeID) %>%
    group_by(PlotYear, QAQC, Module) %>%
    dplyr::summarize(N = n()) %>%
    ungroup() %>%
    spread(key = QAQC, value = N) %>%
    dplyr::rename(Qsurv = `-1`, Nsurv = `0`)
  TempCounts_df <- expand.grid(PlotYear = unique(TreeCounts_df$PlotYear), Module=1:4, stringsAsFactors = FALSE)
  TreeCounts_df %<>%
    right_join(TempCounts_df, by=c("PlotYear", "Module"))
  TreeCounts_df[is.na(TreeCounts_df)] <- 0
  
  TreeCounts_df %<>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Metric = Qsurv - Nsurv, # This metric is the difference in counts
           MatchBy1 = abs(Metric) <= 1,
           MatchBy2 = abs(Metric) <= 2)
  
  x <- abs(TreeCounts_df$Metric)
  TreeCounts_tab <- data.frame(table( factor(x, levels = min(x):max(x))))
  TreeCounts_tab$CDF <- round(cumsum(TreeCounts_tab$Freq)/sum(TreeCounts_tab$Freq)*100, 1)
  names(TreeCounts_tab) <- c("Difference in counts", "Frequency", "CDF (%)")
  comment(TreeCounts_tab) = paste0("MAX. COUNT IN A SURVEY WAS ", max(TreeCounts_df$Qsurv, TreeCounts_df$Nsurv))
  attributes(TreeCounts_tab)$comment2 = paste0("% OF 0-0 PLOTS WAS ", round((nrow(subset(TreeCounts_df, Qsurv==0 & Nsurv==0))/nrow(TreeCounts_df))*100), "%")
  
  QAQCTreesOut_list[["TreeCounts"]][["FinalDat_df"]] <- as.data.frame(TreeCounts_df)
  QAQCTreesOut_list[["TreeCounts"]][["Diff_tab"]] <- TreeCounts_tab # count of trees in each match category
  QAQCTreesOut_list[["TreeCounts"]][["plot_vals"]] <- list(metric_nam = "Count of Trees", type = "discrete", lab = "Absolute difference", v1 = 0, v2 = 1, range_max = max(abs(TreeCounts_df$Metric)) + 1, range_int = 1)
  QAQCTreesOut_list[["TreeCounts"]][["ECDF_plot"]] <- FuncPlotTreesECDF(dat = QAQCTreesOut_list[["TreeCounts"]][["FinalDat_df"]], plot_vals = QAQCTreesOut_list[["TreeCounts"]][["plot_vals"]])
  
  # Seedling / sapling counts ----
  SeedSapCounts_df <- read_csv(here::here("Data_in", "SeedlingsSaplings.csv"))
  SeedSapCounts_df %<>%
    dplyr::mutate(MC = paste0("M", Module, "C", Corner),
           PlotYear = paste(location_name, event_year, sep="_")) %>%
    dplyr::select(PlotYear, QAQC, MC, Seed5_15=Seedling_Ht_5_15, Seed15_30=Seedling_Ht_15_30, Seed30_50=Seedling_Ht_30_50, Seed50_137=Seedling_Ht_50_137, Sap0_1=Sapling_DBH_0_1, Sap1_2.5=Sapling_DBH_1_2half, Sap2.5_5=Sapling_DBH_2half_5, Sap5_10=Sapling_DBH_5_10)
  TempSeedSapCounts <- expand.grid(PlotYear = unique(SeedSapCounts_df$PlotYear), QAQC = unique(SeedSapCounts_df$QAQC), MC = unique(SeedSapCounts_df$MC), stringsAsFactors = FALSE)
  SeedSapCounts_df %<>%
    right_join(TempSeedSapCounts, by=c("PlotYear", "QAQC", "MC"))
  SeedSapCounts_df[is.na(SeedSapCounts_df)] <- 0
  SeedSapCounts_df %<>%
    group_by(PlotYear, QAQC, MC) %>%
    dplyr::summarize_at(.vars=c("Seed5_15", "Seed15_30",  "Seed30_50", "Seed50_137", "Sap0_1", "Sap1_2.5", "Sap2.5_5", "Sap5_10"), funs(n=sum)) %>%
    ungroup()
  names(SeedSapCounts_df) <- c("PlotYear", "QAQC", "MC", "Seed5_15", "Seed15_30", "Seed30_50", "Seed50_137", "Sap0_1", "Sap1_2.5", "Sap2.5_5", "Sap5_10")
  
  for (i in names(SeedSapCounts_df)[4:11]) {
    SeedSapSub_df <- SeedSapCounts_df %>%
      dplyr::select(PlotYear, QAQC, MC, i) %>%
      spread(key=QAQC, value=i) %>%
      dplyr::rename(Qsurv=`-1`, Nsurv=`0`) %>%
      dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
             Metric = Qsurv - Nsurv)
    
    x <- abs(SeedSapSub_df$Metric)
    SeedSapSub_tab <- data.frame(table( factor(x, levels = min(x):max(x))))
    SeedSapSub_tab$CDF <- round(cumsum(SeedSapSub_tab$Freq)/sum(SeedSapSub_tab$Freq)*100, 1)
    names(SeedSapSub_tab) <- c("Difference in counts", "Frequency", "CDF (%)")
    comment(SeedSapSub_tab) = paste0("MAX. COUNT IN A SURVEY WAS ", max(SeedSapSub_df$Qsurv, SeedSapSub_df$Nsurv))
    attributes(SeedSapSub_tab)$comment2 = paste0("% OF 0-0 PLOTS WAS ", round((nrow(subset(SeedSapSub_df, Qsurv==0 & Nsurv==0))/nrow(SeedSapSub_df))*100), "%")

    QAQCTreesOut_list[[i]][["FinalDat_df"]] <- as.data.frame(SeedSapSub_df)
    QAQCTreesOut_list[[i]][["Diff_tab"]] <- SeedSapSub_tab
    QAQCTreesOut_list[[i]][["plot_vals"]] <- list(metric_nam = paste0("Count of ", i), type = "discrete", lab = "Absolute difference", v1 = 1, v2 = 3, range_max = max(abs(SeedSapSub_df$Metric)) + 1, range_int = 1)
    QAQCTreesOut_list[[i]][["ECDF_plot"]] <- FuncPlotTreesECDF(dat = QAQCTreesOut_list[[i]][["FinalDat_df"]], plot_vals = QAQCTreesOut_list[[i]][["plot_vals"]])
    }
  
    # Canopy cover by Plot-Year -----
    Canopy_df <- read_csv(here::here("Data_in", "Densiometer.csv"))
    Canopy_df %<>%
      dplyr::rename(Plot = location_name, Year = event_year, Dxn=Direction, Densio=Spherical_Densiometer_Tally) %>%
      dplyr::mutate(PlotYear = paste0(Plot, "_", Year)) %>% # In case a single plot has QAQC survey in multiple survey cycles
      dplyr::select(PlotYear, QAQC, Module, Dxn, Densio)
    
    Canopy_ByPlot_df <- Canopy_df %>%
      dplyr::select(-Dxn) %>%
      group_by(PlotYear, QAQC) %>%
      dplyr::summarize(AvgDensio = mean(Densio)) %>%
      spread(key = QAQC, value = AvgDensio) %>%
      dplyr::rename(Qsurv = `-1`, Nsurv =`0`)
    
    Canopy_ByPlot_df %<>%
      dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
             Metric = Qsurv - Nsurv,
             MatchBy1 = abs(Metric) <= 1,
             MatchBy3 = abs(Metric) <= 3)
    
    QAQCTreesOut_list[["CanopyPlotCounts"]][["FinalDat_df"]] <- as.data.frame(Canopy_ByPlot_df)
    QAQCTreesOut_list[["CanopyPlotCounts"]][["plot_vals"]] <- list(metric_nam = "Plot-Level Densiometer Count", type = "continuous", lab = "Absolute difference", v1 = 1, v2 = 3, range_max = 3*ceiling(max(abs(Canopy_ByPlot_df$Metric))/3) + 3, range_int = 3)
    QAQCTreesOut_list[["CanopyPlotCounts"]][["Diff_tab"]] <- table(cut(abs(Canopy_ByPlot_df$Metric), breaks = c(0, 1, seq(3, as.numeric(QAQCTreesOut_list[["CanopyPlotCounts"]][["plot_vals"]]["range_max"]), by = 3))))
    
  # Canopy cover by Module-Year -----
  Canopy_ByMod_df <- Canopy_df %>%
    dplyr::select(-Dxn) %>%
    group_by(PlotYear, QAQC, Module) %>%
    dplyr::summarize(AvgDensio = mean(Densio)) %>%
    spread(key = QAQC, value = AvgDensio) %>%
    dplyr::rename(Qsurv = `-1`, Nsurv =`0`)
  
  Canopy_ByMod_df %<>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Metric = Qsurv - Nsurv,
           MatchBy1 = abs(Metric) <= 1,
           MatchBy3 = abs(Metric) <= 3)
  
  QAQCTreesOut_list[["CanopyModCounts"]][["FinalDat_df"]] <- as.data.frame(Canopy_ByMod_df)
  QAQCTreesOut_list[["CanopyModCounts"]][["plot_vals"]] <- list(metric_nam = "Module-Level Densiometer Count", type = "continuous", lab = "Absolute difference", v1 = 1, v2 = 3, range_max = 3*ceiling(max(abs(Canopy_ByMod_df$Metric))/3) + 3, range_int = 3)
  QAQCTreesOut_list[["CanopyModCounts"]][["Diff_tab"]] <- table(cut(abs(Canopy_ByMod_df$Metric), breaks = c(0, 1, seq(3, as.numeric(QAQCTreesOut_list[["CanopyModCounts"]][["plot_vals"]]["range_max"]), by = 3))))

  # Tree DBH -----
  # REQUIRED FIELD for live (status 1) and standing dead (status 2) trees
  # Collected on trees with DBH > 10
  # Only use trees for which a DBH was recorded in both surveys
  # NOTE: MACA008-2016, tree#16 (JUNIVIR) with status=2 but no DBH recorded. This tree was deleted from analysis, leaving 526 trees for analysis. <<<<<<<<<<<<<<<< ADD TO WARNINGS?
  
  # These are the QAQC trees recorded as Status=1 OR 2 in both surveys
  Stat12Trees_df <- subset(QAQCTreesDat_list$TreesN, Status==1 | Status==2, select=c("QTreeID")) %>% # these are the trees with Status=1 or 2 in both surveys
    inner_join(subset(QAQCTreesDat_list$TreesQ, Status==1 | Status==2, select=c("TreeID")), by=c("QTreeID"="TreeID"))
  TreesNSub12 <- subset(QAQCTreesDat_list$TreesN, QTreeID %in% Stat12Trees_df$QTreeID) 
  TreesQSub12 <- subset(QAQCTreesDat_list$TreesQ, TreeID %in% Stat12Trees_df$QTreeID) 
  
  # Summary table
  DBH_df <- FuncQAQCTreesFormat(NDat_df = TreesNSub12, QDat_df = TreesQSub12, cname = "DBH", NANA_mismatch = FALSE) 
  DBH_df <-DBH_df[complete.cases(DBH_df),] # remove the one tree that was missing DBH in non-QAQC survey <<<<<<<<<<<<<<< ADD TO WARNINGS?
  
  DBH_df %<>%
    dplyr::select(-QAQCtreeID) %>%
    dplyr::rename(Nsurv = Value1, Qsurv = Value2) %>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Metric = Qsurv - Nsurv, # negative value means Qsurv smaller
           MatchBy1 = abs(Metric) <= 1,
           MatchBy2 = abs(Metric) <= 2)
  
  x <- abs(DBH_df$Metric)
  DBH_tab <- data.frame(table( factor(x, levels = min(x):max(x))))
  DBH_tab$CDF <- round(cumsum(DBH_tab$Freq)/sum(DBH_tab$Freq)*100, 1)
  names(DBH_tab) <- c("Difference in DBH", "Frequency", "CDF (%)")
  
  QAQCTreesOut_list[["DBH"]][["FinalDat_df"]] <- DBH_df
  QAQCTreesOut_list[["DBH"]][["Diff_tab"]] <- DBH_tab
  QAQCTreesOut_list[["DBH"]][["plot_vals"]] <- list(metric_nam = "Tree DBH", type = "discrete", lab = "Absolute difference", v1 = 0, v2 = 1, range_max = max(abs(DBH_df$Metric)) + 1, range_int = 1)
  QAQCTreesOut_list[["DBH"]][["ECDF_plot"]] <- FuncPlotTreesECDF(dat = QAQCTreesOut_list[["DBH"]][["FinalDat_df"]], plot_vals = QAQCTreesOut_list[["DBH"]][["plot_vals"]])
  
  # Snag height -----
  # REQUIRED FIELD for standing dead (status 2) trees
  # Only used trees for which a DBH was recorded in both surveys
  # NOTE: Two trees classified as Status =2 did not have a snag height recorded: (1) MACA012-2017, tree#11 (JUNIVIR), (2) KIMO011-2012, tree#19 (OXYDARB). These were omitted from analysis, leaving 62 dead trees for analysis. <<<<<<<<<<<<<<<<<<<<<<<< ADD TO WARNINGS?
  
  # Summary table
  SnagHeight_df <- FuncQAQCTreesFormat(NDat_df = TreesNSub2, QDat_df = TreesQSub2, cname = "SnagHeight", NANA_mismatch = FALSE) 
  SnagHeight_df <- SnagHeight_df[complete.cases(SnagHeight_df),] # remove the two trees that were missing DBH in QAQC survey <<<<<<< ADD TO WARNINGS?
  SnagHeight_df %<>%
    dplyr::select(-QAQCtreeID) %>%
    dplyr::rename(Nsurv = Value1, Qsurv = Value2) %>%
    dplyr::mutate(Mean = (Nsurv + Qsurv)/2,
           Diff = Qsurv - Nsurv, # negative value means Qsurv smaller
           Metric = round((Diff/Mean)*100),
           LT10perc = abs(Metric) < 10)
  
  QAQCTreesOut_list[["SnagHeight"]][["FinalDat_df"]] <- SnagHeight_df
  QAQCTreesOut_list[["SnagHeight"]][["plot_vals"]] <- list(metric_nam = "Snag Height", type = "continuous", lab = "% difference", v1 = 10, v2 = 20, range_max = 10*ceiling(max(abs(SnagHeight_df$Metric))/10) + 10, range_int = 10)
  
  SnagHeightTemp_tab <- round(ecdf(abs(SnagHeight_df$Metric))(seq(0, as.numeric(QAQCTreesOut_list[["SnagHeight"]][["plot_vals"]]["range_max"]), by = 10))*100)
  SnagHeightTemp_tab <- data.frame(cbind(Diff = paste0("< = ", seq(0, as.numeric(QAQCTreesOut_list[["SnagHeight"]][["plot_vals"]]["range_max"]), by = 10), "%"), CumPercSnags = SnagHeightTemp_tab))
  QAQCTreesOut_list[["SnagHeight"]][["Diff_tab"]] <- SnagHeightTemp_tab
  
   # TREE CONDITION DAMAGE CODE (OWN DASHBOARD TAB) -----
  # Extract the relevant columns and format
  TreesNCond_df <- QAQCTreesDat_list$TreesNSub1 %>%
    dplyr::select(TreeID, QTreeID, Arborist, Damage1:Damage5) %>%
    gather(key = Temp, value = DamageCode, Damage1:Damage5, na.rm = TRUE) %>% # convert to long
    dplyr::select(-Temp) %>%
    dplyr::mutate(DamageCode = toupper(DamageCode)) %>% # change any lowercase DC to uppercase
    spread(key = DamageCode, value = DamageCode)
  
  TreesNCond_df <- subset(QAQCTreesDat_list$TreesNSub1, select=c("TreeID", "QTreeID", "Arborist")) %>% # include even trees with zero damage codes
    left_join(TreesNCond_df, by = c("TreeID", "QTreeID", "Arborist"))
  TreesNCond_df[, 4:ncol(TreesNCond_df)] <- !is.na(TreesNCond_df[, 4:ncol(TreesNCond_df)]) # For each unique tree, TRUE if that DC was reported    
  
  TreesQCond_df <- QAQCTreesDat_list$TreesQSub1 %>%
    dplyr::select(TreeID, NTreeID, Arborist, Damage1:Damage5) %>%
    gather(key = Temp, value = DamageCode, Damage1:Damage5, na.rm = TRUE) %>%
    dplyr::select(-Temp) %>%
    dplyr::mutate(DamageCode = toupper(DamageCode)) %>% # change any lowercase DC to uppercase
    spread(key = DamageCode, value = DamageCode)
  
  TreesQCond_df <- subset(QAQCTreesDat_list$TreesQSub1, select=c("TreeID", "NTreeID", "Arborist")) %>% # include even trees with zero damage codes
    left_join(TreesQCond_df, by = c("TreeID", "NTreeID", "Arborist"))
  TreesQCond_df[, 4:ncol(TreesQCond_df)] <- !is.na(TreesQCond_df[, 4:ncol(TreesQCond_df)]) # For each unique tree, TRUE if that DC was reported    
  
  all_DC <- sort(unique(c(colnames(TreesNCond_df)[4:ncol(TreesNCond_df)],colnames(TreesQCond_df)[4:ncol(TreesQCond_df)]))) # these are all the possible damage codes
  all_arborists <- sort(unique(c(TreesNCond_df$Arborist, TreesQCond_df$Arborist)))
  
  # Summarize for each Arborist, count of survey-trees
  TreeCondArbN_df <- data.frame(Arborist=all_arborists) %>%
    rowwise() %>%
    dplyr::mutate(NumTrees = sum(TreesNCond_df$Arborist==Arborist) + sum(TreesQCond_df$Arborist==Arborist)) 
  
  # For each DC, evaluate QAQC
  TreeCondSummary_list <- TreeCondPSA_list <- TreeCondArb_list <- list()
  for (x in all_DC) {
    # Extract the relevant non-QAQC data
    if(x %in% colnames(TreesNCond_df)) { # if the DC exists in non-QAQC data
      CondNDat_df <- subset(TreesNCond_df, select = c("TreeID", "QTreeID", "Arborist", x)) # then extract that column along with tree ID information
    } else {
      CondNDat_df <- add_column(subset(TreesNCond_df, select=c("TreeID", "QTreeID", "Arborist")), Temp = FALSE) # otherwise, create the data frame with FALSE for all entries of that DC
      names(CondNDat_df)[names(CondNDat_df)=="Temp"] <- x
    }
    
    # Extract the relevant QAQC data
    if(x %in% colnames(TreesQCond_df)) { # if the DC exists in QAQC data
      CondQDat_df <- subset(TreesQCond_df, select = c("TreeID", "NTreeID", "Arborist", x)) # then extract that column along with tree ID information
    } else {
      CondQDat_df <- add_column(subset(TreesQCond_df, select=c("TreeID", "NTreeID", "Arborist")), Temp = FALSE)
      names(CondQDat_df)[names(CondQDat_df)=="Temp"] <- x
    }
    
    # Merge the QAQC and non-QAQC data
    TreeCond_df <- CondNDat_df %>%
      left_join(CondQDat_df, by=c("QTreeID" = "TreeID")) %>%
      dplyr::select(-TreeID, -QTreeID, -NTreeID)
    TreeCond_df[is.na(TreeCond_df)] <- FALSE # now it's two columns of logicals, where TRUE means the disturbance was recorded. One column is for QAQC, the other is non-QAQC.
    names(TreeCond_df) <- c("Arborist1", "Value1", "Arborist2", "Value2")
    
    # Tree summary
    TreeCond_df$MatchType <- ifelse(TreeCond_df$Value1 & TreeCond_df$Value2, "11", ifelse(!TreeCond_df$Value1 & !TreeCond_df$Value2, "00", "10/01")) # 1 = TRUE, 0 = FALSE
    TempMatchCount <- TreeCond_df %>% dplyr::count(MatchType)
    TreeCondSummary_list[[x]] <- data.frame(DC = x, TempMatchCount)
    TreeCond_df <- as.data.frame(TreeCond_df)
    TreeCondPSA_list[[x]] <- FuncPSA(categ_vec = c(TRUE, FALSE), DatWide_df = TreeCond_df, c1 = "Value1", c2 = "Value2", missing_val = FALSE)
    
    # Summary by Arborist-DC
    TreeCondArb_list[[x]] <- data.frame(Arborist = all_arborists) %>%
      rowwise() %>%
      dplyr::mutate(DC = x,
             ArbTrue = sum(TreeCond_df$Arborist1==Arborist & TreeCond_df$Value1) + sum(TreeCond_df$Arborist2==Arborist & TreeCond_df$Value2),
             ArbDenom = nrow(TreeCond_df[TreeCond_df$Arborist1==Arborist & (TreeCond_df$Value1|TreeCond_df$Value2),]) + nrow(TreeCond_df[TreeCond_df$Arborist2==Arborist & (TreeCond_df$Value1|TreeCond_df$Value2),]))
  } # end of damage code loop
  
  # Summarize lists in data frames
  TreeCondSummary_df <- do.call("rbind", TreeCondSummary_list)
  TreeCondSummary_df <- spread(TreeCondSummary_df, key=MatchType, value=n)
  TreeCondSummary_df$DC <- as.character(TreeCondSummary_df$DC)
  TreeCondSummary_df[is.na(TreeCondSummary_df)] <- 0
  
  TreeCondPSA_df <- do.call("rbind", TreeCondPSA_list)
  TreeCondPSA_df <- TreeCondPSA_df %>% 
    rownames_to_column(var="DC") %>%
    filter(Grp==TRUE) %>%
    dplyr::select(DC, PSA, Denom, `5%`, `95%`) %>%
    dplyr::mutate(DC=str_replace_all(DC, ".1", ""))
  TreeCondPSA_df$PSA[is.na(TreeCondPSA_df$PSA)] <- 0
  TreeCondPSASummary_df <- data.frame(DC=unique(all_DC), stringsAsFactors = FALSE) %>%
    left_join(TreeCondSummary_df, by="DC") %>%
    full_join(TreeCondPSA_df, by="DC") %>%
    dplyr::rename("FALSE-FALSE"=`00`, "TRUE-FALSE"=`10/01`, "TRUE-TRUE"=`11`)
  TreeCondPSASummary_df$DC <- factor(TreeCondPSASummary_df$DC, levels=all_DC)
  TreeCondPSASummary_df <- droplevels(TreeCondPSASummary_df) # <<<<<< SAVE RDS
  
  # This is the summary PSA plot
  TreeCondPSA_plot <- ggplot(TreeCondPSASummary_df, aes(x=DC,  y=PSA)) +
    geom_point(na.rm = TRUE) +
    geom_linerange(aes(ymin=`5%`, ymax=`95%`), na.rm = TRUE) +
    ylim(0, 100) +
    scale_x_discrete(labels=paste0(TreeCondPSASummary_df$DC,"\n(",TreeCondPSASummary_df$Denom,")")) +
    labs(y="Proportion of specific agreement (%)", x="Tree Condition", subtitle="(90% CI included for N > 5)") +
    theme_bw(base_size = 10)
  TreeCondPSASummary_df$Denom <- NULL
  
  TreeCondArb_df <- do.call("rbind", TreeCondArb_list)
  TreeCondArb_df$DC <- factor(TreeCondArb_df$DC, levels=all_DC)
  TreeCondArb_df$percentTRUE <- as.integer((TreeCondArb_df$ArbTrue/TreeCondArb_df$ArbDenom)*100)
  
  QAQCTreesOut_list[["TreeCond_list"]] <- list(
    all_DC = all_DC, # vector of damage codes
    TreeCondArbN_df = TreeCondArbN_df, # count of survey-trees per Arborist
    TreeCondPSASummary_df = TreeCondPSASummary_df, # PSA summary data
    TreeCondPSA_plot = TreeCondPSA_plot, # the summary PSA plot
    TreeCondArb_df = TreeCondArb_df) # for 'Proportion of Specific True' by 'Arborist' plots 
  
  # TREE FOLIAGE CONDITION (OWN DASHBOARD TAB) -----
  # ONLY USE TREES WITH STATUS=1 IN BOTH SURVEYS. BECAUSE THIS IS NOT A REQUIRED FIELD, "NA/NA" IS TABULATED SEPARATELY FROM MATCH AND MISMATCH
  
  # Extract the relevant columns and format
  NTreeFol_df <- QAQCTreesDat_list$TreesNSub1 %>%
    dplyr::select(QTreeID, Foliage1:Foliage5) %>%
    gather(key = Temp, value = FCcombined, Foliage1:Foliage5, na.rm = TRUE) %>% # convert to long
    dplyr::select(-Temp) %>%
    filter(FCcombined!="NA.NA") %>%
    dplyr::mutate(FC = gsub("\\..*", "", FCcombined), 
           FClevel = gsub("^*.\\.", "", FCcombined)) %>% # split apart the foliage code and foliage category
    dplyr::select(-FCcombined)
  
  QTreeFol_df <- QAQCTreesDat_list$TreesQSub1 %>%
    dplyr::select(TreeID, Foliage1:Foliage5) %>%
    gather(key = Temp, value = FCcombined, Foliage1:Foliage5, na.rm = TRUE) %>% # convert to long
    dplyr::select(-Temp) %>%
    filter(FCcombined!="NA.NA") %>%
    dplyr::mutate(FC = gsub("\\..*", "", FCcombined), 
           FClevel = gsub("^*.\\.", "", FCcombined)) %>% # split apart the foliage code and foliage category
    dplyr::select(-FCcombined)
  
  all_FC <- sort(unique(c(NTreeFol_df$FC, QTreeFol_df$FC)))
  all_arborists <- sort(unique(c(QAQCTreesDat_list$TreesNSub1$Arborist, QAQCTreesDat_list$TreesQSub1$Arborist)))
  
  FCSummary_list <- FC_PSA_list <- FC_PSA_plotlist <- list()
  
  for (f in all_FC) {
    nFC_df <- subset(QAQCTreesDat_list$TreesNSub1, select = c("QTreeID", "Arborist")) %>% # make sure all relevant trees are included
      left_join(subset(NTreeFol_df, FC==f), by = "QTreeID") %>%
      dplyr::mutate(FC = f) 
    qFC_df <- subset(QAQCTreesDat_list$TreesQSub1, select = c("TreeID", "Arborist")) %>% # make sure all relevant trees are included
      left_join(subset(QTreeFol_df, FC==f), by="TreeID") %>%
      dplyr::mutate(FC = f) %>%
      dplyr::rename("QTreeID" = "TreeID")
    FC_df <- nFC_df %>%
      left_join(qFC_df, by=c("QTreeID", "FC")) %>%
      dplyr::select(-QTreeID, -FC)
    names(FC_df) <-c("Arborist1", "Value1", "Arborist2", "Value2")
    FC_df$Value1[FC_df$Value1 == "NA"] <- NA
    FC_df$Value1 <- as.integer(FC_df$Value1)
    FC_df$Value2[FC_df$Value2 == "NA"] <- NA
    FC_df$Value2 <- as.integer(FC_df$Value2)
    
    FC_df$MatchType <- ifelse(is.na(FC_df$Value1) & is.na(FC_df$Value2), "NA/NA", ifelse(FC_df$Value1==FC_df$Value2, "MATCH", "MISMATCH"))
    FC_df$MatchType[is.na(FC_df$MatchType)] <- "MISMATCH"
    Temp <- FC_df %>% dplyr::count(MatchType)
    FCSummary_list[[f]] <- data.frame(FC = f, Temp)
    FC_df <- as.data.frame(FC_df)
    
    FC_PSA_list[[f]] <- FuncPSA(categ_vec = 1:5, DatWide_df = FC_df, c1="Value1", c2="Value2", missing_val = TRUE)
    FC_PSA_plotlist[[f]] <- FuncPlotPSA(PSA_df = FC_PSA_list[[f]], xlabel = paste0("Foliage Condition: ", f))
  }
  
  FCSummary_tab <- do.call("rbind", FCSummary_list)
  FCSummary_tab <- spread(FCSummary_tab, key = MatchType, value = n)
  FCSummary_tab[is.na(FCSummary_tab)] <- 0
  
  QAQCTreesOut_list[["FoliageCond_list"]] <- list(
    all_FC = all_FC, # vector of foliage condition codes
    FCSummary_tab = FCSummary_tab, # summary table
    FC_PSA_list = FC_PSA_list, 
    FC_PSA_plotlist = FC_PSA_plotlist)
  
  saveRDS(QAQCTreesOut_list, "Temp_out/QAQCTreesOut.RDS")
}

FuncQAQCHerbsMain <- function () {
  # Function to run herb QAQC analyses.
  HerbsTemp_df <- read_csv(here::here("Data_in", "Herbs.csv")) %>%
    filter(Species_Original!="Carya sp.") # <<<<<<<<<<<<<<<<<<<<<< THIS IS A TEMPORARY FIX SO I CAN CONTINUE WORKING ON DATABASE--THERE ARE DUPLICATE RECORD PROBLEMS WITH THIS SPECIES
  herbs_to_combine <- c("VIOL", "VITI", "CELT", "CARY", "CARE", "SANI", "AGRI", "ULMU", "SMIL", "DICH", "QUER", "POTE", "ELEP", "RUBU", "DESM", "GALI") # species within these herb genera are difficult to ID, so are analyzed at genus level for some analyses
  
  # Herb species richness, data formatting -----
  HerbsTemp_df %<>% # rename columns
    dplyr::select(Plot=location_name, QAQC, Year=event_year, PlantName=Species_Original, CC=Cover_Class, M1C2=Presence_Column_1, M1C4=Presence_Column_2, M2C2=Presence_Column_3, M2C3=Presence_Column_4, M3C2=Presence_Column_5, M3C4=Presence_Column_6, M4C2=Presence_Column_7, M4C3=Presence_Column_8) %>%
    dplyr::mutate(QAQC= QAQC==-1, PlotID= paste0(Plot, "_", Year), # '-1' means TRUE for QAQC
           RowID = as.numeric(row.names(.))) %>% # add the row number so records can be easily referenced back to original data
    dplyr::select(RowID, PlotID, everything())
  
  # Identify which records to include in base analysis -----
  HerbsTemp_df[grep("sp.", HerbsTemp_df$PlantName, invert = FALSE), "BaseAnalysis"] <- FALSE # these were defined only to genus. These affect QAQC results.
  HerbsTemp_df$BaseAnalysis[is.na(HerbsTemp_df$PlantName)] <- FALSE # records where plant name not provided
  HerbsTemp_df$BaseAnalysis[is.na(HerbsTemp_df$BaseAnalysis)] <- TRUE # these will be included in the base analysis
  
  # Identify herbs with high cover in a base analysis (evaluated at species level) -----
  HerbsHiCov_df <- HerbsTemp_df %>% # At least one of the surveys called it CC 2 or above for that plot.
    dplyr::select(PlotID, PlantName, CC) %>%
    filter(CC > 1) %>%
    dplyr::select(-CC) %>%
    distinct() %>%
    dplyr::mutate(BaseHiCov = TRUE)
  
  HerbsTemp2_df <- HerbsTemp_df %>%
    left_join(HerbsHiCov_df, by=c("PlotID", "PlantName")) %>%
    dplyr::mutate(BaseHiCov = replace_na(BaseHiCov, FALSE))
  HerbsTemp2_df$BaseHiCov[HerbsTemp2_df$BaseAnalysis == FALSE] <- NA # NA if the record shouldn't be included in baseline analysis anyway
  
  # Herb species that should be rolled up to genus level for some analyses -----
  HerbsTemp2_df$Genus <- substr(HerbsTemp2_df$PlantName, start=1, stop=4) # The genus is the first four letters in the plant name
  HerbsTemp2_df$Genus <- toupper(HerbsTemp2_df$Genus)
  HerbsTemp2_df$PlantNameCOMB <- HerbsTemp2_df$PlantName
  HerbsTemp2_df$PlantNameCOMB[HerbsTemp2_df$Genus %in% herbs_to_combine] <- HerbsTemp2_df$Genus[HerbsTemp2_df$Genus %in% herbs_to_combine]
  
  HerbsTemp2_df[grep("sp.", HerbsTemp2_df$PlantNameCOMB, invert = FALSE), "RolledAnalysis"] <- FALSE # these will NOT be included in the rolled analysis
  HerbsTemp2_df$RolledAnalysis[is.na(HerbsTemp2_df$PlantNameCOMB)] <- FALSE # records where plant name not provided
  HerbsTemp2_df$RolledAnalysis[is.na(HerbsTemp2_df$RolledAnalysis)] <- TRUE # these will be included in the rolled analysis
  HerbsTemp2_df$PlantNameCOMB[HerbsTemp2_df$RolledAnalysis == FALSE] <- NA
  HerbsTemp2_df <- dplyr::select(HerbsTemp2_df, -Plot, -Year, -Genus)
  
  # Fill-in the module-corner data for baseline analysis -----
  HerbsTemp3_df <- FuncReformatMC(Dat_df = subset(HerbsTemp2_df))
  HerbsCleanBase_df <- FuncHerbsMC(Dat_df = HerbsTemp3_df)
  
  # Fill-in the module-corner data for rolled analysis -----
  HerbsRolledTemp_df <- HerbsTemp2_df %>%
    filter(RolledAnalysis == TRUE) %>%
    dplyr::select(-PlantName) %>% 
    group_by(PlotID, QAQC, PlantNameCOMB) %>% # keep the maximum MC value (1,2,3) for species rolled up by genus
    dplyr::rename(PlantName = PlantNameCOMB) %>%
    dplyr::summarize(M1C2 = max(M1C2, na.rm=TRUE),
              M1C4 = max(M1C4, na.rm=TRUE),
              M2C2 = max(M2C2, na.rm=TRUE),
              M2C3 = max(M2C3, na.rm=TRUE),
              M3C2 = max(M3C2, na.rm=TRUE),
              M3C4 = max(M3C4, na.rm=TRUE),
              M4C2 = max(M4C2, na.rm=TRUE),
              M4C3 = max(M4C3, na.rm=TRUE)) %>%
    ungroup()
  
  HerbsRolledTemp_df[HerbsRolledTemp_df=="-Inf"] <- NA
  HerbsRolledTemp2_df <- FuncReformatMC(Dat_df = HerbsRolledTemp_df)
  HerbsCleanRolled_df <- FuncHerbsMC(Dat_df = HerbsRolledTemp2_df)
  QAQCHerbsDat_list <- list(
    HerbsBase = HerbsCleanBase_df, 
    HerbsHiCov = subset(HerbsCleanBase_df, BaseHiCov==TRUE),
    HerbsRolled = HerbsCleanRolled_df)
  saveRDS(QAQCHerbsDat_list, "Temp_out/QAQCHerbsDat.RDS")
  
  # Summarize the herb QAQC data -----
  # Initialize lists
  HM100SqM_list <- HS100SqM_list <- HS100SqMTable_list <- HM10SqM_list <- HS10SqM_list <- HS10SqMTable_list <- HM1SqM_list <- HS1SqM_list <- HS1SqMTable_list <- list()
  
  for (x in c("HerbsBase", "HerbsHiCov", "HerbsRolled")) {
    HerbsSubdat_df <- QAQCHerbsDat_list[[x]]
    
    # 100SqM-level summary
    Herbs_100SqM <- HerbsSubdat_df %>%
      dplyr::select(PlotID, QAQC, PlantName, InPlot) %>%
      spread(key=QAQC, value=InPlot) # these are all TRUE, so will have to then replace NA's with FALSE
    colnames(Herbs_100SqM)[3:4] <- c("Nsurv", "Qsurv")
    
    Temp100SqM_list <- FuncHerbsMatch(Dat_df = Herbs_100SqM, group_vec = "PlotID")
    
    HM100SqM_list[[x]] <- Temp100SqM_list[[1]] # this one is used to generate the herb match report in .Rmd
    HS100SqM_list[[x]] <- Temp100SqM_list[[2]] # this one is used to plot herb richness histograms in .Rmd
    HS100SqMTable_list[[x]] <- Temp100SqM_list[[3]]
    
    # 10SqM-level summary
    Herbs10SqMTemp_df <- HerbsSubdat_df[, c("PlotID", "QAQC", "PlantName", "M1", "M2", "M3", "M4")]
    Herbs10SqM_df <- melt(data=Herbs10SqMTemp_df, measure.vars=c("M1", "M2", "M3", "M4"), variable.name="10SqM") %>%
      spread(key=QAQC, value=value) 
    colnames(Herbs10SqM_df)[4:5] <- c("Nsurv", "Qsurv")
    
    Temp10SqM_list <- FuncHerbsMatch(Dat_df = Herbs10SqM_df, group_vec = c("PlotID", "`10SqM`"))
    
    HM10SqM_list[[x]] <- Temp10SqM_list[[1]] # this one is used to generate the herb match report in .Rmd
    HS10SqM_list[[x]] <- Temp10SqM_list[[2]] # this one is used to plot herb richness histograms in .Rmd
    HS10SqMTable_list[[x]] <- Temp10SqM_list[[3]]
    
    # 1SqM-level summary
    Herbs1SqMTemp_df <- HerbsSubdat_df[, c("PlotID", "QAQC", "PlantName", "M1C2", "M1C4", "M2C2", "M2C3", "M3C2", "M3C4", "M4C2", "M4C3")]
    Herbs1SqM_df <- melt(data=Herbs1SqMTemp_df, measure.vars=c("M1C2", "M1C4", "M2C2", "M2C3", "M3C2", "M3C4", "M4C2", "M4C3"), variable.name="1SqM") %>%
      spread(key=QAQC, value=value) 
    colnames(Herbs1SqM_df)[4:5] <- c("Nsurv", "Qsurv")
    
    Temp1SqM_list <- FuncHerbsMatch(Dat_df = Herbs1SqM_df, group_vec = c("PlotID", "`1SqM`"))
    
    HM1SqM_list[[x]] <- Temp1SqM_list[[1]] # this one is used to generate the herb match report in .Rmd
    HS1SqM_list[[x]] <- Temp1SqM_list[[2]] # this one is used to plot herb richness histograms in .Rmd
    HS1SqMTable_list[[x]] <- Temp1SqM_list[[3]]
  }
  
  QAQCHerbsOut_list <- list(
    HM100SqM_list = HM100SqM_list, # for matched species reports
    HS100SqM_list = HS100SqM_list, # for histograms
    HS100SqMTable_list = HS100SqMTable_list,
    HM10SqM_list = HM10SqM_list,
    HS10SqM_list = HS10SqM_list,
    HS10SqMTable_list = HS10SqMTable_list,
    HM1SqM_list = HM1SqM_list,
    HS1SqM_list = HS1SqM_list,
    HS1SqMTable_list = HS1SqMTable_list)
  saveRDS(QAQCHerbsOut_list, "Temp_out/QAQCHerbsOut.RDS")
}

FuncOakHickMain <- function () {
  # Function to run oak-hickory regeneration analyses on Oak-Hickory plots only
  VegRegen_df <- read_csv(here::here("Temp_out", "VegDat_cleaned.csv")) %>% # <<<<<<<<<<<< GENERATE AN ERROR IF THIS IS NOT AVAILABLE --ACTUALLY NEED TO MAKE SURE ALSO THAT IT'S UPDATED
    dplyr::select(SubPark, Cycle, Year, Plot, Latitude, Longitude, QAQCPlot, CGShort, EventID) %>%
    distinct()
  VegRegen_df$Cycle <- factor(VegRegen_df$Cycle)
  
  CycLevs <- sort(levels(VegRegen_df$Cycle), decreasing = FALSE)
  TempLevs <- as.numeric(substring(CycLevs, 7, 10)) # grab the end-year for each survey cycle
  CycLabels <- paste0("[", TempLevs - 4, ",", TempLevs, "]") 
  
  PlotCG_link <- unique(VegRegen_df[c("Plot", "CGShort")])
  
  # If there are survey events in the 'RegenSS_df' and 'RegenTree_df' data that are not in the 'main' data, need to fill in missing plot information
  FuncMissingRegen <- function(Dat_df, PlotCG_link, TreeDat = FALSE) {
    MissingEvents <- unique(Dat_df$EventID[is.na(Dat_df$Plot)])
    for (event in MissingEvents) {
      which_row <- which(Dat_df$EventID==event)
      Dat_df[which_row, "SubPark"] <- substr(event, start=1, stop=4)
      Dat_df[which_row, "Plot"] <- temp_plot <- substr(event, start=1, stop=7)
      Dat_df[which_row, "CGShort"] <- ifelse(temp_plot %in% PlotCG_link$Plot, PlotCG_link$CGShort[PlotCG_link$Plot==temp_plot], NA) # <<<<<<<<<<<<< PUTS NA FOR CGSHORT IF IT DOESN'T EXIST IN PLOTCG_LINK, WHICH IS CREATED FROM MAIN VEGETATION DATA--OTHERWISE ERRORS B/C CUGA010 IS MISSING
      Dat_df[which_row, "Year"] <- temp_yr <- as.integer(gsub('.*- ','', event))
      Dat_df[which_row, "Cycle"] <- cut(temp_yr, breaks=seq(2010, max(Dat_df$Year, na.rm = TRUE) + 4, by = 5))
      Dat_df[which_row, "QAQCPlot"] <- FALSE # <<<<<<<<<<<<< ASSUMING THIS IS TRUE...NEED TO GENERALIZE
      if(TreeDat) Dat_df[which_row, "TreeID"] <- gsub("NA", temp_plot, Dat_df[which_row, "TreeID"])
    }
    AddToMain <- Dat_df %>%
      filter(EventID %in% MissingEvents) %>%
      dplyr::select(SubPark, Cycle, Year, Plot, QAQCPlot, CGShort, EventID) %>%
      distinct() # add these rows of data to the main dataset
    return_list <- list(Dat_df, AddToMain)
    return(return_list)
  }
  
  # Format seedling/sapling data -----
  RegenSS_df <- read_csv(here::here("Data_in", "OakHickSS.csv")) %>%
    dplyr::rename(EventID = event_name_date_calc) %>%
    dplyr::mutate(MC = paste0("M", Module, "C", Corner)) %>%
    left_join(VegRegen_df, by = "EventID") %>%
    dplyr::select(SubPark, EventID, Plot, CGShort, Year, Cycle, QAQCPlot, Module, MC, PlantName=Plant_Name, Seed5_50=Seedling_Ht_5_50, Seed5_15=Seedling_Ht_5_15, Seed15_30=Seedling_Ht_15_30, Seed30_50=Seedling_Ht_30_50, Seed50_137=Seedling_Ht_50_137, Sap0_1=Sapling_DBH_0_1, Sap1_2.5=Sapling_DBH_1_2half, Sap2.5_5=Sapling_DBH_2half_5, Sap5_10=Sapling_DBH_5_10) %>%
    dplyr::mutate(OHSpecies = gsub(' .*$','', PlantName) %in% c("Quercus", "Carya"))
  wrong_MC <- RegenSS_df[!RegenSS_df$MC %in% c("M1C2", "M1C4", "M2C2", "M2C3", "M3C2", "M3C4", "M4C2", "M4C3"), c("EventID", "MC")] # >>>>>>>>>>>>>>>>>>>>> THROW ERROR--THESE PLOT SURVEYS HAVE INCORRECT MODULE-CORNER DESIGNATIONS IN THE SEEDLING/SAPLING DATA (WON'T AFFECT CALCS BUT THEY SHOULD CORRECT IT)
  RegenSS_df <- RegenSS_df[!RegenSS_df$EventID == "CARL001 - 2016", ] # >>>>>> TWO CORNERS WERE NOT SURVEYED, SO DELETE THIS SURVEY EVENT -- NEED TO DISCUSS WITH CUPN
  RegenSS_df$MC[RegenSS_df$MC=="M2C4"] <- "M2C3" # <<<<<<<<<<<< SEEMS THEY HAVE TWO WRONG CORNERS. THERE ARE ALSO OTHERS THAT ARE INCORRECT.
  
  MissingTemp <- FuncMissingRegen(Dat_df = RegenSS_df, PlotCG_link = PlotCG_link, TreeDat = FALSE)
  RegenSS_df <- MissingTemp[[1]]
  RegenSS_df[is.na(RegenSS_df)] <- 0
  VegRegen_df <- rbind.fill(VegRegen_df, MissingTemp[[2]])
  
  # Format tree DBH data -----
  RegenTree_df <- read_csv(here::here("Data_in", "OakHickTree.csv")) %>% # tree dbh data, status=1 means alive, status=2 means dead. Some data import issue with "Snag_Height" but we are not using those data anyway in this analysis <<<<<<<<<<<<<< THROWS A WARNING
    dplyr::rename(EventID = event_name_date_calc) %>%
    left_join(VegRegen_df, by = "EventID") %>%
    # filter(Plot %in% oak_hick_plots) %>% # only keep data from Oak-Hickory plots
    dplyr::mutate(TreeID = paste0(Plot, "_", Tree_Number)) %>%
    dplyr::select(EventID, SubPark, Plot, CGShort, Year, Cycle, QAQCPlot, Module, TreeID, PlantName = Plant_Name, DBH, Status = Tree_Status) %>%
    dplyr::mutate(OHSpecies = gsub(' .*$','', PlantName) %in% c("Quercus", "Carya"))
  RegenTree_df$Status[RegenTree_df$Status %in% c(0,4)] <- 2 # Status of 0 and 4 also means tree is naturally dead
  RegenTree_df <- RegenTree_df[!RegenTree_df$EventID == "CARL001 - 2016", ] # >>>>>> TWO CORNERS WERE NOT SURVEYED, SO DELETE THIS SURVEY EVENT
  MissingTemp <- FuncMissingRegen(Dat_df = RegenTree_df, PlotCG_link = PlotCG_link, TreeDat = TRUE)
  RegenTree_df <- MissingTemp[[1]]
  VegRegen_df <- rbind.fill(VegRegen_df, MissingTemp[[2]])
  
  oak_hick_plots <- unique(VegRegen_df$Plot[VegRegen_df$CGShort %in% c("G015", "G165", "G601", "G159", "G180", "G166", "G012")]) # vector of plots that are oak-hickory.
  ruderal_plots <- unique(VegRegen_df$Plot[VegRegen_df$CGShort %in% c("G030", "G031", "G552", "G553", "G557", "G583")])
  
  # Calculate Sander threshold -----
  # Threshold value is 1 per 10m2 box
  # Make sure no more than eight MC's per plot
  MC_count <- RegenSS_df %>%
    dplyr::select(EventID, MC) %>%
    distinct() %>%
    dplyr::count(EventID)
  if(nrow(over_MC <- MC_count[MC_count$n > 8 | is.na(MC_count$n),]) > 0) stop(paste0("These survey events have more than 8 module-corners: \n", over_MC)) # <<<<<<<<<< ERROR-CHECK
  
  SanderSap_df <- RegenSS_df %>%
    dplyr::select(SubPark, EventID, Plot, CGShort, Year, Cycle, QAQCPlot, Module, MC, PlantName, OHSpecies, Sap0_1, Sap1_2.5, Sap2.5_5) %>% #only keep the counts from relevant sapling sizes. These are all counted in the 10m2 plots.
    dplyr::mutate(CGtype=ifelse(Plot %in% oak_hick_plots, "OakHick", ifelse(Plot %in% ruderal_plots, "Ruderal", "Other"))) %>%
    filter(CGtype!="Other")
  SanderSap_df <- droplevels(SanderSap_df)
  
  # All communities... count only the oak-hickory species
  SanderSapAgg_df <- SanderSap_df %>%
    filter(OHSpecies==TRUE) %>% # only work with oak-hickory species
    dplyr::mutate(SapSum = Sap0_1 + Sap1_2.5 + Sap2.5_5) %>%
    group_by(EventID, SubPark, CGtype, CGShort, Plot, Year, Cycle, QAQCPlot, MC) %>% # aggregate by plot-corner
    dplyr::summarise(SapCount = sum(SapSum)) %>%
    group_by(EventID) %>%
    dplyr::mutate(MCid = 1:n()) %>%
    dplyr::select(-MC)
  SanderTemplate_df <- SanderSap_df %>%
    select(EventID, SubPark, CGtype, CGShort, Plot, Year, Cycle, QAQCPlot) %>%
    distinct()
  SanderTemplate_df <- merge(SanderTemplate_df, as.data.frame(1:8))
  names(SanderTemplate_df)[names(SanderTemplate_df)=="1:8"] <- "MCid"
  
  SanderFinal_df <- SanderSapAgg_df %>%
    right_join(SanderTemplate_df, by=c("EventID", "SubPark", "CGtype", "CGShort", "Plot", "Year", "Cycle", "QAQCPlot", "MCid"))
  SanderFinal_df$SapCount[is.na(SanderFinal_df$SapCount)] <- 0
  
  # Hierarchical bootstrap to obtain confidence intervals for each Park in Cycle 1 (not doing Cycle 2 because not all plots surveyed in Cycle 2 yet). QAQC data are included in the bootstrapping. <<<<<<<<<<<<<<<
  
  SanderSapCycle1 <- SanderFinal_df %>%
    ungroup() %>%
    filter(Cycle=="(2010,2015]") %>%
    dplyr::select(SubPark, CGtype, CGShort, Plot, SapCount)
  SanderSapCycle1 <- droplevels(SanderSapCycle1)
  
  Sander_Nplots <- SanderSapCycle1 %>%
    group_by(CGtype, SubPark) %>%
    dplyr::summarise(Nplots=n_distinct(Plot))
  Sander_Nplots <- as.data.frame(Sander_Nplots)
  
  BootFunc <- function(dat, reps=2000) {
    dat <- droplevels(dat)
    boot_mean <- list()
    for (r in 1:reps) {
      boot_counts <- list()
      boot_plots <- sample(x=unique(dat$Plot), size=length(unique(dat$Plot)), replace=TRUE) # bootstrapped plots from that SubPark
      for (i in 1:length(boot_plots)) {
        boot_counts[[i]] <- sample(x=dat$SapCount[dat$Plot==boot_plots[i]], size=8, replace=TRUE) # the bootstrapped counts--8 per bootstrapped plot
      }
      boot_mean[[r]] <- mean(unlist(boot_counts)) # mean box count of saplings for one bootstrapped sample of the SubPark
    } # end of bootstrap replicates
    return_list <- list(
      unlist(boot_mean),
      quantile(unlist(boot_mean), probs=c(0.05, 0.5, 0.95))
    )
    return(return_list) # returns the vector of bootstrapped means and also a vector of the median and 90% confidence interval limits
  }
  
  SanderBootMeans <- SanderBootQuants <- list()
  for(t in c("OakHick", "Ruderal")) { # separately for OakHick & Ruderal
    for (p in unique(SanderSapCycle1$SubPark)) { # separately for each SubPark
      dat <- subset(SanderSapCycle1, CGtype==t & SubPark==p)
      if(nrow(dat)>0) {
        temp <- BootFunc(dat)
        SanderBootMeans[[paste0(t, "_", p)]] <- temp[[1]]
        SanderBootQuants[[paste0(t, "_", p)]] <- c(
          Type=t,
          SubPark=p,
          temp[[2]])
      }
    }
  }
  
  SanderBootSummary_df <- do.call("rbind", SanderBootQuants)
  SanderBootSummary_df <- as.data.frame(SanderBootSummary_df, stringsAsFactors = FALSE)
  rownames(SanderBootSummary_df) <- NULL
  for(x in 3:ncol(SanderBootSummary_df)) {
    SanderBootSummary_df[,x] <- round(as.numeric(as.character(SanderBootSummary_df[,x])), 2)
  }
  names(SanderBootSummary_df)[names(SanderBootSummary_df)=="Type"] <- "CGtype"
  attr(SanderBootSummary_df$CGtype, "names") <- attr(SanderBootSummary_df$SubPark, "names") <- NULL
  
  # Sander summaries
  SanderSummary_df <- Sander_Nplots %>%
    full_join(SanderBootSummary_df, by=c("CGtype", "SubPark")) %>%
    arrange(CGtype, `50%`)
  
  # Now only working with oak-hickory plots -----
  OHDat_df <- VegRegen_df %>%
    filter(Plot %in% oak_hick_plots)
  
  # Calculate OSI -----
  # These are the oak-hickory plots non-QAQC surveyed more than once
  OHResurveySP_df <- OHDat_df %>%
    filter(QAQCPlot==FALSE) %>% # non-QAQC Oak-Hickory plots only
    group_by(SubPark, Plot) %>%
    dplyr::summarize(NumSurveys=n()) %>%
    group_by(SubPark) %>%
    dplyr::summarize(NumOHResurveys = sum(NumSurveys > 1, na.rm=TRUE)) # <<<<<<<<<<< ORIGINALLY THIS WAS EQUAL 2, NOW I JUST SAY > 1--DOES THAT WORK? CHANGE TEXT
  
  OHPlotN_df <- OHDat_df %>% # Count of Oak-Hickory plots and plot resurveys by SubPark
    dplyr::select(SubPark, Plot, Cycle, QAQCPlot) %>%
    distinct() %>%
    group_by(SubPark) %>%
    dplyr::summarize(NumOHPlots = n_distinct(Plot),
                     NumOHQAQC = sum(QAQCPlot==TRUE, na.rm=TRUE)) %>%
    full_join(OHResurveySP_df, by="SubPark") %>%
    dplyr::select(SubPark, NumOHPlots, NumOHResurveys, NumOHQAQC)
  
  OHresurvey_N <- table(subset(OHDat_df, QAQCPlot==FALSE)$Plot)
  OHresurvey_plots <- names(OHresurvey_N)[OHresurvey_N > 1] # these are the Oak-Hickory plots with non-QAQC data for more than one survey cycle
  
  # Calculate the additional sapling counts in Oak-Hickory plots
  OHSapAdd_df <- RegenTree_df %>%
    filter(Plot %in% oak_hick_plots & DBH >=10 & DBH <= 15 & Status==1 & QAQCPlot == FALSE) %>% # has a "live" status and a DBH <= 15cm
    group_by(EventID, SubPark, CGShort, Plot, Cycle, OHSpecies, Module) %>%
    dplyr::summarize(Trees10_15=n())
  OHSpecies <- c(TRUE, FALSE)
  Module <- 1:4
  OHSapTemplate_df <- merge(data.frame(unique(OHSapAdd_df[c("EventID", "SubPark", "CGShort", "Plot", "Cycle")]), stringsAsFactors = FALSE), as.data.frame(Module))
  OHSapTemplate_df <- merge(data.frame(unique(OHSapTemplate_df[c("EventID", "SubPark", "CGShort", "Plot", "Cycle", "Module")]), stringsAsFactors = FALSE), as.data.frame(OHSpecies))
  OHSapAdd_df <- left_join(OHSapTemplate_df, OHSapAdd_df, by=c("EventID", "SubPark", "CGShort", "Plot", "Cycle", "OHSpecies", "Module"))
  OHSapAdd_df[is.na(OHSapAdd_df)] <- 0
  
  # Calculate seedling and sapling ratios for OSI -----
  # For each plot there should be no more than 8 unique MC. 
  OH_SSTemp_df <- RegenSS_df %>%
    dplyr::filter(Plot %in% oak_hick_plots) %>%
    dplyr::mutate(CUPNseed = dplyr::select(., Seed5_50:Seed50_137) %>% rowSums,
                  Sap0_2.5 = dplyr::select(., Sap0_1:Sap1_2.5) %>% rowSums,
                  Sap2.5_10 = dplyr::select(., Sap2.5_5:Sap5_10) %>% rowSums) %>%
    filter(QAQCPlot == FALSE) %>%
    dplyr::select(EventID, SubPark, CGShort, Plot, Cycle, MC, OHSpecies, CUPNseed, Sap0_2.5, Sap2.5_10) %>%
    group_by(EventID, SubPark, CGShort, Plot, Cycle, MC, OHSpecies) %>%
    dplyr::summarize(CUPNseed = sum(CUPNseed, na.rm=TRUE), 
                     Sap0_2.5 = sum(Sap0_2.5, na.rm=TRUE), # these are seedlings, by Woodall's definition (but recorded as saplings by CUPN)
                     Sap2.5_10 = sum(Sap2.5_10, na.rm=TRUE)) %>% # these are a subset of the saplings by Woodall's definition (but need to add small trees)
    ungroup()
  MC <- unique(OH_SSTemp_df$MC)
  OH_SSTemplate_df <- merge(data.frame(unique(OH_SSTemp_df[c("EventID", "SubPark", "CGShort", "Plot", "Cycle")]), stringsAsFactors = FALSE), as.data.frame(MC, stringsAsFactors = FALSE))
  OH_SSTemplate_df <- merge(OH_SSTemplate_df, as.data.frame(OHSpecies)) # make sure for every survey event, there is a record for both OakHick and non-OakHick species
  OH_SSTemp2_df <- left_join(OH_SSTemplate_df, OH_SSTemp_df, by=c("EventID", "SubPark", "CGShort", "Plot", "Cycle", "MC", "OHSpecies"))
  OH_SSTemp2_df[is.na(OH_SSTemp2_df)] <- 0
  
  # Bootstrap counts for OSI
  OSIboot_out <- vector("list", length = 1000)
  for (B in 1:1000) { # 1000 bootstrap samples
    
    OSIboot_samp <- OH_SSTemp2_df %>%
      group_by(Plot, Cycle, OHSpecies) %>%
      sample_n(size = 8, replace = TRUE)
    
    # Bootstrap the SapAdd counts
    OHSapAddSamp_df <- OHSapAdd_df %>%
      group_by(Plot, Cycle, OHSpecies) %>%
      sample_n(size = 4, replace = TRUE) %>%
      group_by(EventID, SubPark, CGShort, Plot, Cycle, OHSpecies) %>%
      dplyr::summarize(SapAdd = sum(Trees10_15, na.rm=TRUE))
    
    OSIboot_temp <- OSIboot_samp %>%
      group_by(EventID, SubPark, CGShort, Plot, Cycle, OHSpecies) %>%
      dplyr::summarize(CUPNseed = sum(CUPNseed, na.rm=TRUE), 
                       Sap0_2.5 = sum(Sap0_2.5, na.rm=TRUE),
                       Sap2.5_10 = sum(Sap2.5_10, na.rm=TRUE)) %>%
      ungroup() %>%
      left_join(OHSapAddSamp_df, by=c("EventID", "SubPark", "CGShort", "Plot", "Cycle", "OHSpecies")) # add the tree data
    
    OSIboot_temp$SapAdd[is.na( OSIboot_temp$SapAdd)] <- 0
    
    OSIboot_temp %<>%
      dplyr::mutate(SeedOSI = round((CUPNseed * 50) + (Sap0_2.5 * 5)), # for CUPNseed, it's 8m2 seedling counts extrapolated to 400m2; for sap0_2.5, it's 80m2 counts extrapolated to 400m2
                    SapOSI = round((Sap2.5_10 * 5) + SapAdd)) # for sap2.5_10, it's 80m2 extrapolated to 400m2; Trees10_15 is counts for the 400m2 plot
    OSIboot_temp$OHSpecies[OSIboot_temp$OHSpecies] <- "OH"
    OSIboot_temp$OHSpecies[OSIboot_temp$OHSpecies!="OH"] <- "notOH"
    OSIboot_temp$OHSpecies <- factor(OSIboot_temp$OHSpecies)
    
    OHSeed_df <- dcast(OSIboot_temp, EventID + SubPark + CGShort + Plot + Cycle ~ OHSpecies, value.var=c("SeedOSI"))
    OHSeed_df$SeedRatio <- OHSeed_df$OH/(OHSeed_df$OH+OHSeed_df$notOH)
    OHSeed_df$SeedRatio[is.na(OHSeed_df$SeedRatio)] <- 0 # If the denominator is zero (so undefined ratio), count it as 0
    
    OHSap_df <- dcast(OSIboot_temp, EventID + SubPark + CGShort + Plot + Cycle ~ OHSpecies, value.var=c("SapOSI"))
    OHSap_df$SapRatio <- OHSap_df$OH/(OHSap_df$OH+OHSap_df$notOH)
    OHSap_df$SapRatio[is.na(OHSap_df$SapRatio)] <- 0 # If the denominator is zero (so undefined ratio), count it as 0
    
    OSIboot_out[[B]] <- full_join(OHSeed_df, OHSap_df, by=c("EventID", "SubPark", "CGShort", "Plot", "Cycle")) %>%
      dplyr::select(EventID, SubPark, CGShort, Plot, Cycle, SeedRatio, SapRatio) %>%
      dplyr::mutate(PartialOSI = SeedRatio + SapRatio)
  } # end 1000 bootstraps
  
  # Calculate bootstrap 95% CI
  OSIboot_df <- do.call("rbind", OSIboot_out) %>%
    group_by(EventID, SubPark, CGShort, Plot, Cycle) %>%
    dplyr::summarize(MeanOSI = mean(PartialOSI, na.rm=TRUE),
                     Low95OSI = quantile(PartialOSI, 0.025),
                     High95OSI = quantile(PartialOSI, 0.975))
  
  # Calculate seedling and sapling categories for Importance -----
  # For each plot there should be no more than 8 unique MC. 
  # DON'T CALCULATE IMPORTANCE FOR PLOTS THAT COUNTED 'Seed5_50' rather than the smaller categories
  
  # Calculate the Oak-Hickory tree counts in Oak-Hickory plots
  OHTreeAdd_df <- RegenTree_df %>%
    filter(Plot %in% oak_hick_plots & Status==1 & QAQCPlot == FALSE & OHSpecies == TRUE) %>% # has a "live" status and a DBH <= 15cm
    select(-Year, -QAQCPlot, -OHSpecies, -TreeID, -PlantName, -DBH, -Status) %>%
    group_by(EventID, SubPark, CGShort, Plot, Cycle, Module) %>%
    tally()
  
  OHTreeTemplate_df <- merge(data.frame(unique(OHTreeAdd_df[c("EventID", "SubPark", "CGShort", "Plot", "Cycle")]), stringsAsFactors = FALSE), as.data.frame(Module))
  OHTreeAdd_df <- left_join(OHTreeTemplate_df, OHTreeAdd_df, by=c("EventID", "SubPark", "CGShort", "Plot", "Cycle", "Module")) %>%
    select(-Module)
  OHTreeAdd_df[is.na(OHTreeAdd_df)] <- 0
  
  # # Don't calculate Importance for earlier surveys with difference seedling categories <<<<<<<<<<<<<
  # SkipImport_df <- RegenSS_df %>%
  #   filter(QAQCPlot==FALSE & Cycle=="(2010,2015]") %>%
  #   dplyr::mutate(CUPNseed_subcat = select(.,"Seed5_15:Seed30_50") %>% rowSums) %>%
  #   group_by(Plot, Year, Cycle) %>%
  #   dplyr::summarize(Sum_Seed5_50 = sum(Seed5_50),
  #                    Sum_CUPNseed_subcat = sum(CUPNseed_subcat)) %>%
  #   filter(Sum_CUPNseed_subcat==0) %>% # these are the ones should not calculate Importance for
  #   mutate(CalcImport = FALSE) %>%
  #   select(-Sum_Seed5_50, -Sum_CUPNseed_subcat)
  
  OHImport_SSTemp_df <- RegenSS_df %>%
    filter(Plot %in% oak_hick_plots & OHSpecies & QAQCPlot==FALSE) %>%
    # left_join(SkipImport_df, by = c("Plot", "Year", "Cycle")) %>%
    # filter(is.na(CalcImport)) %>%
    # select(-Seed5_50, -CalcImport, -Year, -QAQCPlot, -Module, -PlantName, -OHSpecies) %>% # These are the ones to calculate
    select(-Seed5_50, -Year, -QAQCPlot, -Module, -PlantName, -OHSpecies) %>%
    dplyr::mutate(Import_1pt = dplyr::select(., Seed5_15:Seed15_30) %>% rowSums,
                  Import_50pt = dplyr::select(., Sap0_1:Sap5_10 ) %>% rowSums) %>%
    dplyr::rename(Import_2pt = Seed30_50,
                  Import_10pt = Seed50_137) %>%
    select(-Seed5_15, -Seed15_30, -Sap0_1, -Sap1_2.5, -Sap2.5_5, -Sap5_10) %>%
    group_by(EventID, SubPark, CGShort, Plot, Cycle, MC) %>%
    dplyr::summarize(Import_1pt = sum(Import_1pt, na.rm=TRUE),
                     Import_2pt = sum(Import_2pt, na.rm=TRUE),
                     Import_10pt = sum(Import_10pt, na.rm=TRUE),
                     Import_50pt = sum(Import_50pt, na.rm=TRUE))
  
  OHImport_SSTemplate_df <- merge(data.frame(unique(OHImport_SSTemp_df[c("EventID", "SubPark", "CGShort", "Plot", "Cycle")]), stringsAsFactors = FALSE), as.data.frame(MC, stringsAsFactors = FALSE))
  OHImport_SSTemp2_df <- left_join(OHImport_SSTemplate_df, OHImport_SSTemp_df, by=c("EventID", "SubPark", "CGShort", "Plot", "Cycle", "MC"))
  OHImport_SSTemp2_df[is.na(OHImport_SSTemp2_df)] <- 0
  
  # Bootstrap counts for Importance
  Importboot_out <- vector("list", length = 1000)
  for (B in 1:1000) { # 1000 bootstrap samples
    
    Importboot_samp <- OHImport_SSTemp2_df %>%
      group_by(Plot, Cycle) %>%
      sample_n(size = 8, replace = TRUE)
    
    Importboot_temp <- Importboot_samp %>%
      group_by(EventID, SubPark, CGShort, Plot, Cycle) %>%
      dplyr::summarize(Import_1ptScore = sum(Import_1pt, na.rm=TRUE)*(10000/8), 
                       Import_2ptScore = sum(Import_2pt, na.rm=TRUE)*2*(10000/8),
                       Import_10ptScore = sum(Import_10pt, na.rm=TRUE)*10*(10000/8),
                       Import_50ptScore = sum(Import_50pt, na.rm=TRUE)*50*(10000/80)) %>%
      ungroup()
    
    # Bootstrap the tree counts
    Importboot_treesamp <- OHTreeAdd_df %>%
      group_by(Plot, Cycle) %>%
      sample_n(size = 4, replace = TRUE)
    
    Importboot_treetemp <- Importboot_treesamp %>%
      group_by(EventID, SubPark, CGShort, Plot, Cycle) %>%
      dplyr::summarize(ImportTrees_50ptScore = sum(n, na.rm=TRUE)*50*(10000/400)) # trees get score of 50, then multiplied out to a hectare
    
    Importboot2_temp <- Importboot_temp %>% left_join(Importboot_treetemp, by = c("EventID", "SubPark", "CGShort", "Plot", "Cycle"))
    Importboot2_temp[is.na(Importboot2_temp)] <- 0
    Importboot2_temp$Import_50ptScore <- Importboot2_temp$Import_50ptScore + Importboot2_temp$ImportTrees_50ptScore
    Importboot2_temp$ImportTrees_50ptScore <- NULL
    Importboot_out[[B]] <- Importboot2_temp %>% 
      mutate(logImportance = log(Import_1ptScore + Import_2ptScore + Import_10ptScore + Import_50ptScore + 1)) %>%
      select(-Import_1ptScore, -Import_2ptScore, -Import_10ptScore, -Import_50ptScore)
  } # end 1000 bootstraps
  
  # Calculate bootstrap 95% CI
  Importboot_df <- do.call("rbind", Importboot_out) %>%
    group_by(EventID, SubPark, CGShort, Plot, Cycle) %>%
    dplyr::summarize(
      MeanImport = mean(logImportance, na.rm=TRUE),
      Low95Import = quantile(logImportance, 0.025),
      High95Import = quantile(logImportance, 0.975))
  
  OHtemp_df <- full_join(OSIboot_df, Importboot_df, by = c("EventID", "SubPark", "CGShort", "Plot", "Cycle"))
  OHtemp_df %<>% mutate_if(is.numeric, round, 2)
  
  # Add in specific survey dates -----
  Burn_df <- read_csv(here::here("Data_in", "BurnPlots.csv"))
  PlotDates <- gather(Burn_df[, c("Plot", "Date_Plot_Established", "Date_Plot_Revisited")], key = SurveyType, value = SurveyDate, -Plot) %>%
    dplyr::select(-SurveyType) %>%
    dplyr::mutate(SurveyDate = mdy(SurveyDate),
                  Year = year(SurveyDate), 
                  QAQCPlot = as.logical(FALSE)) %>%
    filter(!is.na(Year) & Plot %in% oak_hick_plots) %>%
    distinct()
  
  OHDat_df2 <- OHDat_df %>% # add the specific survey dates to these plots
    full_join(PlotDates, by=c("Plot", "Year", "QAQCPlot"))
  OHDat_df2[is.na(OHDat_df2$Plot) | is.na(OHDat_df2$Year) | is.na(OHDat_df2$QAQCPlot),] # <<<<<<<<<<<<< ERROR CHECK--data frame of entries in Burn Data that weren't in OH Veg Data
  
  OH_df <- OHtemp_df %>%
    left_join(OHDat_df2, by=c("SubPark", "EventID", "Plot", "CGShort", "Cycle"))
  MissingYears <- which(is.na(OH_df$Year)) # <<<<<<<<<<<< ERROR CHECK FOR PLOT-SURVEYS MISSING IN THE MAIN DATA
  MissingYears
  OH_df[MissingYears, "Year"] <- substr(OH_df[MissingYears, "EventID"], nchar(OH_df[MissingYears, "EventID"])-3, nchar(OH_df[MissingYears, "EventID"]))   
  OH_df$Year <- as.integer(OH_df$Year)
  
  # Add in the burn dates -----
  Burns <- Burn_df %>%
    dplyr::select(Plot, BurnDate = Fire_Start_Date) %>%
    filter(!is.na(BurnDate)) %>%
    distinct() %>%
    dplyr::mutate(BurnDate = mdy(BurnDate))
  
  OH_df$MostRecentBurn <- OH_df$PriorBurnCycle <- OH_df$PriorBurn <- NA
  class(OH_df$MostRecentBurn) <- "Date"
  for (i in 1:nrow(OH_df)) { 
    if(OH_df[i, "QAQCPlot"]==FALSE) {
      AllBurns_df <- subset(Burns, Plot==as.character(OH_df[i, "Plot"]), select=c("BurnDate")) # These are all the prescribed burn dates for that plot
      AllBurns_df$Use <- ifelse(is.na(AllBurns_df$BurnDate), NA, AllBurns_df$BurnDate < OH_df[i,]$SurveyDate) # Logical indicating if burn occurred prior to survey
      OH_df[i, "PriorBurn"] <- sum(AllBurns_df$Use, na.rm=TRUE) > 0
      RecentBurn <- ifelse(OH_df[i, "PriorBurn"], ymd(max(AllBurns_df$BurnDate[AllBurns_df$Use==TRUE], na.rm=TRUE)), NA)
      if(!is.na(RecentBurn)) class(RecentBurn) <- "Date"
      OH_df[i, "MostRecentBurn"] <- RecentBurn
      if(!is.na(RecentBurn)) {
        CycleNumber <- as.numeric(OH_df[i, "Cycle"])
        if(CycleNumber==1) OH_df[i, "PriorBurnCycle"] <- TRUE # If this is Cycle1, then it's counted in PriorBurnCycle if there is ANY recorded prior burn, no matter how long ago...
        if(CycleNumber > 1) {
          PriorSurvDate <- unique(OH_df$SurveyDate[OH_df$Plot==as.character(OH_df[i, "Plot"]) & OH_df$Cycle==levels(OH_df$Cycle)[CycleNumber-1] & OH_df$QAQCPlot==FALSE]) # If this is post-Cycle 1, then it's counted in PriorBurnCycle if it occurred after the prior survey date and before this cycle's survey date
          if(length(PriorSurvDate)>0) { # if there is a prior survey date
            OH_df[i, "PriorBurnCycle"] <- RecentBurn > PriorSurvDate
          } else {
            OH_df[i, "PriorBurnCycle"] <- RecentBurn > ymd(OH_df[i,]$SurveyDate-1825) # this is equivalent to five years; takes care of plots that didn't have Cycle 1 surveys but had Cycle 2 surveys
          }
        }
      }
    }
  }
  
  # Generate output data -----
  OH_df$Cycle <- as.character(OH_df$Cycle)
  OHfinal_template <- expand.grid(Plot = unique(OH_df$Plot), Cycle = unique(OH_df$Cycle), stringsAsFactors = FALSE) # make sure all cycles for each plot are represented
  OHfinal_df <- left_join(OHfinal_template, OH_df, by = c("Plot", "Cycle"))
  OHfinal_df$PriorBurnCycle[is.na(OHfinal_df$PriorBurnCycle)] <- FALSE
  OHfinal_df$CycLabels <- paste0("[", as.numeric(substring(OHfinal_df$Cycle, 7, 10)) - 4, ",", as.numeric(substring(OHfinal_df$Cycle, 7, 10)), "]")
  
  saveRDS(OHfinal_df, "Temp_out/OakHickDat.RDS")
  OakHickOut_list <- list(
    OH_CycLevs = CycLevs,
    OH_CycLabels = CycLabels,
    OH_Sander_df = SanderSummary_df,
    OHfinal_df = OHfinal_df)
  saveRDS(OakHickOut_list, "Temp_out/OakHickOut.RDS")
}