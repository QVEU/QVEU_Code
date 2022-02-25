#NEON_optimizationplot.R
#Patrick Dolan - 12/1/21
#Purpose: plotting optimization data from Neon optimizaton protocol. 

library(ggplot2)
library(data.table)

optimTable<-fread("/Volumes/dolanpt$/My Documents/Research_Admin_QVEU/Lab_Documents/Neon_Optimization.csv")#point to your optimization table (see file in sharepoint folder)

ggplot(optimTable)+geom_tile(aes(pulse_width,voltage,fill=transfection_eff))+facet_grid(~pulse_number)+scale_fill_viridis_c()+theme_light()+geom_text(aes(pulse_width,voltage+50,label=sample))