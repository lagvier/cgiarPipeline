#accumulative_temp
#accumulative_rh
#accumulative_ppt
#accumulative_dpv (posible)
# library(nasapower)
# library(pacman)

# p_load(nasapower,readr,tidyverse,plotly)
#
# # head(result$data$pheno)
# # unique(result$data$pheno[,c("Env","Year")])
# #
# # # Upload pheno data -------------------------------------------------------------
# # bioflow <- read_csv("D:/OneDrive - CGIAR/2024/HACKATHON FEB 2024/bioflow.csv")
# # n<-names(bioflow %>% dplyr::select(-all_of(trait)))
# # #unique(bioflow$Field_Location)
#
# # Parameters --------------------------------------------------------------
# trait<-c("Yield_Mg_ha")
# #f_location<-c("NYH3")
# #y<-c(2021)
# LAT<- -27.48
# LONG<-151.81
#
# obj<-data()
# df<-obj$data$pheno
#
# date_planted<-min(df$Date_Planted)
# date_harvest<-max(df$Date_Harvested)
#
# temporal<-"hourly" #daily,monthly

# Funtion for download Climate from NASA  ---------------------------------
nasaPowerExtraction <- function(LAT,LONG,date_planted,date_harvest){

# Climate data ------------------------------------------------------------
    WTH<-nasapower::get_power(community = "ag", #c("ag", "re", "sb")
                            lonlat = c(LONG,LAT), #Decimal degrees
                            pars = c("RH2M","T2M","PRECTOTCORR"),# "T2M_MAX","T2M_MIN"
                            dates = c(date_planted, date_harvest),#YYYY-MM-DD
                            temporal_api = temporal) %>%
          dplyr::mutate(datetime=ISOdate(YEAR, MO, DY,HR),
                date=as.Date(datetime))

# descriptive  --------------------------------------------------------
  descriptive<-summary(dplyr::select(WTH,RH2M,T2M,PRECTOTCORR))

# Outputs -----------------------------------------------------------------
  output<- list(WTH = WTH,
                descriptive = descriptive)
  return(output)
}









# d<-np(LAT,LONG,date_planted,date_harvest) #test example
#
# d$TS_RH
#
#
#
#
# # Filter by... ------------------------------------------------------------
# bioflow <- bioflow %>%
#   dplyr::select(all_of(n),trait) %>%
#   dplyr::filter(Year %in% y) %>%
#   #dplyr::filter(Field_Location %in% f_location) %>%
#   dplyr::mutate("pl_date.{trait}":=Date_Planted,
#                 "ev_date.{trait}":=Date_Harvested,
#                 "days.{trait}"   :=Date_Harvested-Date_Planted)
# #table(bioflow$days.Yield_Mg_ha)
# arg<-unique(bioflow$days.Yield_Mg_ha)
# # max(bioflow$Date_Harvested)
# # min(bioflow$Date_Planted)
#
# # Intern parameters -------------------------------------------------------
# date1<-min(bioflow$Date_Planted)
# date2<-max(bioflow$Date_Harvested)
#
# # Upload climate data  ----------------------------------------------------
# WTH<-nasapower::get_power(community = "ag", #c("ag", "re", "sb")
#                       lonlat = c(LONG,LAT), #Decimal degrees
#                       pars = c("RH2M","T2M","PRECTOTCORR"),# "T2M", "PRECTOTCORR","T2M_MAX","T2M_MIN"),
#                       dates = c(date1, date2),#YYYY-MM-DD
#                       temporal_api = "hourly")
# wth<-WTH %>%
#   dplyr::mutate(date=ISOdate(YEAR, MO, DY,HR),
#          date=as.Date(date),"days.{trait}":=date-date1) %>%
#   dplyr::filter(days.Yield_Mg_ha>=min(bioflow$days.Yield_Mg_ha))
#
# #View(wth)
# #names(wth)
#
#
#
# # b_clim <- wth%>%
# #   #dplyr::filter(date<=date2) %>%
# #   mutate(days.Yield_Mg_ha=as.factor(days.Yield_Mg_ha)) %>%
# #   ggplot(aes(x=days.Yield_Mg_ha,y=RH2M,colour=days.Yield_Mg_ha))+
# #   geom_boxplot()+
# #   # geom_label_repel(
# #   #   mapping = aes(label = ifelse(best, Lines, NA)),
# #   #   size = 2,
# #   #   max.overlaps = 50,
# #   #   alpha = 0.9, fill = "yellow")+
# #   #facet_wrap(~TRT,nrow = 2)+
# #   stat_summary(
# #     fun = median,
# #     geom = 'line',
# #     aes(group = days.Yield_Mg_ha , colour = days.Yield_Mg_ha),
# #     position = position_dodge(width = 0.9) #this has to be added
# #   )+xlab(" ")+
# #   ylab("Relative Humidity (%)")+
# #   labs(title="",colour="")+
# #   theme(legend.position = "none",axis.text.x=element_text(angle = 60))
# #
# # b_clim





