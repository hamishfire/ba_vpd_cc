# Hamish Clarke 2021.03.23 hamishc@uow.edu.au

# This script runs a generalised linear model of the probability of 
#  fire as function of vapour pressure deficit (VPD)

# load libraries
library(data.table)
library(tidyverse)
library(rgdal)
library(Hmisc)

set.seed(42)

# specify decimal function
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# calculate vpd function
calc_vpd <- function(temp, dewpt) {
  # From temperature we can compute saturation vapour pressure (per Tetens equation)
  es <- 0.6108*exp(17.27*temp/(temp + 237.3))  # Temp in DegC, es in kPa
  
  # From dewpoint temperature and temperature we can calculate RH, 
  # per http://www.bom.gov.au/climate/averages/climatology/relhum/calc-rh.pdf 
  rh <- 100*exp(1.8096 + (17.2694*dewpt)/(237.3+dewpt)) /  
    exp(1.8096 + (17.2694*temp)/(237.3+temp)) # dewpt temp in DegC, temp in DegC, rh in %
  
  # From saturation vapour pressure and relative humidity we can compute vpd
  vpd <- es*(1-rh/100) # es in kPa, rh in %, vpd in kPa
  return(vpd)
}

### SET DIRECTORIES, CONSTANTS 
base_diri <- "PATH/TO/BASE/DIRECTORY/" # base
burn_diri <- "PATH/TO/BURN/DIR/" # burn data
abs_diri <- "PATH/TO/ABS/DIR/" # absence data
diro <- "PATH/TO/DIRO/" # out directory

biome_diri <- "PATH/TO/BIOME/" # Terrestrial biomes

winstr <- sprintf("Win%.02d",1:21) # 21 windows; See MODIS C6 Burned Area User Guide 1.2
yrnum <- 2003:2020                 # 18 years

biome_num <- c(1:14,98:99)
biome_name <- c("Tropical and Subtropical Moist Broadleaf Forests",
                "Tropical and Subtropical Dry Broadleaf Forests",
                "Tropical and Subtropical Coniferous Forests",
                "Temperate Broadleaf and Mixed Forests",
                "Temperate Coniferous Forests",
                "Boreal Forests/Taiga",
                "Tropical and subtropical grasslands, savannas, and shrublands",
                "Temperate Grasslands, Savannas, and Shrublands",
                "Flooded Grasslands and Savannas",
                "Montane Grasslands and Shrublands",
                "Tundra",
                "Mediterranean Forests, Woodlands, and Scrub",
                "Deserts and Xeric Shrublands",
                "Mangroves",
                "Lakes",
                "Rock and Ice")

win_name <- c("Alaska","Canada","USA (Conterminous)","Central America","South America (North)",
              "South America (Central)","South America (South)","Europe","West and North Africa",
              "Central and North Africa","East Africa and Arabian Peninsula","Southern Africa (North)",
              "Southern Africa (South)","Madagascar","Russia and Central Asia 1",
              "Russia and Central Asia 2","Russia (Kamachatka)","South Asia","South East Asia",
              "Australia","New Zealand","Azores","Cape Verde Island","Hawaii")

##### GET BIOME 
# from WWF Terrestrial Ecoregions of the World

big_list <- list() # put all windows, biomes in here

# Loop through windows
for (wctr in 1:length(winstr)) {

  print(paste0("Now processing ", winstr[wctr]))
  
    # load all files into single big dataframe
    pres_data <- 
      list.files(paste0(burn_diri,winstr[wctr]), pattern="*.csv", full.names = TRUE) %>% 
      map_df(~fread(.))
    
    # remove incomplete BIOME cases
    pres_data <- pres_data[complete.cases(pres_data$biome),]
    
    pres_data$date <- as.Date(pres_data$date)
    
    # remove BA post Feb 29 2020
    pres_data <- pres_data[pres_data$date<"2020-03-01",]

    # adjust units for vpd equation - celsius required
    pres_data$era5_tmax <- pres_data$era5_tmax-273.15
    pres_data$era5_dmin <- pres_data$era5_dmin-273.15
    
    # load absence data
    abs_data <- 
      list.files(paste0(abs_diri,winstr[wctr]), pattern="*.csv", full.names = TRUE) %>% 
      map_df(~fread(.))
    
    # remove incomplete BIOME cases
    abs_data <- abs_data[complete.cases(abs_data$biome),]
    
    # remove points burnt in last ~5 years
    abs_data <- abs_data[is.na(abs_data$tsf),]
    
    abs_data$date <- as.Date(abs_data$date)
    
    # adjust units for vpd equation - celsius required
    abs_data$era5_tmax <- abs_data$era5_tmax-273.15
    abs_data$era5_dmin <- abs_data$era5_dmin-273.15
    
    # compute vpd
    pres_data$vpd <- calc_vpd(pres_data$era5_tmax, pres_data$era5_dmin)
    abs_data$vpd <- calc_vpd(abs_data$era5_tmax, abs_data$era5_dmin)
    
    # divide into biomes
    ubiome_pres <- unique(pres_data$biome)
    ubiome_abs <- unique(abs_data$biome)
    ubiome <- intersect(ubiome_pres, ubiome_abs)
    ubiome <- sort(ubiome)
    
    # Print exceptions
    if (length(setdiff(ubiome_abs,ubiome_pres))!=0) {
      print(paste0("Biome ", setdiff(ubiome_abs,ubiome_pres)," omitted because missing from presence data"))  
    }
    
    if (length(setdiff(ubiome_pres,ubiome_abs))!=0) {
      print(paste0("Biome ", setdiff(ubiome_pres,ubiome_abs)," omitted because missing from absence data"))  
    }
    rm(list=c("ubiome_pres","ubiome_abs"))
    

  win_list <- list() # all biome results for each window here
  
  # loop through biomes
  for (bctr in 1:length(ubiome)) { 

    print(paste0("Now processing biome ", ubiome[bctr]))

    logitmod_results.df<-data.frame(matrix(nrow=1,ncol=16))
    colnames(logitmod_results.df)<-c("window","biome",
                                     "nFires", "nPresence", "nAbsence",
                                     "logitmodCoef1", "logitmodCoef2",
                                     "VPD_Pfire05","VPD_Pfire05_upr","VPD_Pfire05_lwr",
                                     "Accuracy", "DevExpl",
                                     "confmat1_1", "confmat1_2","confmat2_1", "conf_mat2_2")


    # remove NA 
    pres_vpd <- pres_vpd[complete.cases(pres_vpd$vpd),]
    abs_vpd_all <- abs_vpd_all[complete.cases(abs_vpd_all$vpd),]

    # skip if sample size too small
    if (nrow(pres_vpd)<100) {
      print("Less than 100 samples; skipping...")
      next
    }
    
    # equalise sample size
    if (nrow(abs_vpd_all)<nrow(pres_vpd)) {
      print("WARNING: Fewer absence points than presence points")
      print(paste0("  Ratio = ",specify_decimal(nrow(abs_vpd_all)/nrow(pres_vpd),2),":1"))
      abs_vpd <- abs_vpd_all
    } else {
      abs_vpd <- abs_vpd_all[sample(nrow(abs_vpd_all), size=nrow(pres_vpd), replace = F),]
    }
    rm(abs_vpd_all)      
    
    mdata <- rbind(pres_vpd,abs_vpd)

    # logistic model
    logit_mod <- glm(formula = fire ~ vpd , data = mdata , family = binomial )
    
    summary(logit_mod)
    
    # uncertainty per https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
    ilink <- family(logit_mod)$linkinv
  
    ## add fit and se.fit on the **link** scale
    ndata <- setNames(as_tibble(predict(logit_mod, mdata, se.fit = TRUE)[1:2]),
                                       c('fit_link','se_link'))
    ## create the interval and backtransform
    ndata <- mutate(ndata,
                    fit_resp  = ilink(fit_link),
                    right_upr = ilink(fit_link + (2 * se_link)),
                    right_lwr = ilink(fit_link - (2 * se_link)))
  
    conf_mat_pred <- table( mdata$fire,ndata$fit_resp >0.5)
    conf_mat_upr <- table( mdata$fire,ndata$right_upr>0.5)
    conf_mat_lwr <- table( mdata$fire,ndata$right_lwr>0.5)
    
    # working on upper and lower bounds...
    se <- sqrt(diag(vcov(logit_mod)))
    
    # get coefficients and accuracy
    logitmod_results.df[1] <- winstr[wctr]
    logitmod_results.df[2] <- ubiome[bctr]
    logitmod_results.df[3] <- nrow(mdata)
    logitmod_results.df[4] <- nrow(pres_vpd)
    logitmod_results.df[5] <- nrow(abs_vpd)
    logitmod_results.df[6] <- coefficients(logit_mod)[1]
    logitmod_results.df[7] <- coefficients(logit_mod)[2]
    logitmod_results.df[8] <- -(coefficients(logit_mod)[1]/coefficients(logit_mod)[2])
    logitmod_results.df[9] <- -(coefficients(logit_mod)[1] + se[1]) / (coefficients(logit_mod)[2] - se[2])
    
    logitmod_results.df[10] <- -(coefficients(logit_mod)[1]-se[1]) / (coefficients(logit_mod)[2] +se[2])
    logitmod_results.df[11] <- (conf_mat_pred[1,1]+conf_mat_pred[2,2])/sum(conf_mat_pred)
    dev_expl <- 100*(logit_mod$null.deviance-logit_mod$deviance)/logit_mod$null.deviance 
    print(dev_expl)
    logitmod_results.df[12] <- dev_expl
    rm(dev_expl)
    logitmod_results.df[13] <- conf_mat_pred[1,1]
    logitmod_results.df[14] <- conf_mat_pred[1,2]
    logitmod_results.df[15] <- conf_mat_pred[2,1]
    logitmod_results.df[16] <- conf_mat_pred[2,2]
  
    VPD_Pfire05 <- logitmod_results.df[8]  # VPD_Pfire05 from the model fit
    VPD_Pfire05.UPR <- logitmod_results.df[9]  # UPR 95% CI VPD_Pfire05 from the 100 x model fit
    VPD_Pfire05.LWR <- logitmod_results.df[10]  # LWR 95% CI VPD_Pfire05 from the 100 x model fit
    

    rm(list=c("pres_vpd","abs_vpd","mdata"))  
    rm(list=c("conf_mat_lwr","conf_mat_upr","conf_mat_pred","poly.yy","poly.xx","se",
              "VPD_Pfire05","VPD_Pfire05.LWR","VPD_Pfire05.UPR",
              "logit_mod","ndata"))
    
    win_list[[bctr]] <- logitmod_results.df
    rm(logitmod_results.df)
  } # end biome loop 
    
    win_df <- do.call(rbind,win_list) # join biome DFs into window DF
    rm(win_list)
    
    big_list[[wctr]] <- win_df 
    rm(list=c("pres_data","abs_data"))
} # end window loop

big_df <- do.call(rbind, big_list) # join window DFs into one big DF
rm(big_list)
filo_ext <- "output.csv"
diro_csv <- diro

filo_csv <-paste0(diro_csv,filo_ext)
write.table(big_df, file=filo_csv, sep=",", row.names=FALSE, col.names=TRUE)
