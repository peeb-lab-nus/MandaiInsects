## PACKAGES ============
library(readxl)
library(ggplot2)
library(plyr)

## DIRECTORIES ============
main.dir <- "~/Dropbox/Projects/Special Projects/Mandai"
fig.dir <- file.path(main.dir, "figures")

## Checking trap metadata
x <- read_excel(file.path(main.dir, "MandaiAllSeq_Trapwise_JP_17Aug2021_Ver3_edited.xlsx"),
           sheet = "Sample List")


x <- read_excel(file.path(main.dir, "MandaiAllSeq_Trapwise_JP_17Aug2021_Ver3_edited.xlsx"),
                sheet = "All Seq")

# Create a subset with only malaise trap metadata
x2 <- subset(x, type == "MT")

# Plot trapping frequency for each sampling locality
trap_freq_plot <- ggplot() + geom_point(aes(y = TrapID, x = date), data = x2)
ggsave(filename = file.path(fig.dir, "trap_freq.pdf"), 
       trap_freq_plot,
       width = 6, 
       height = 6)

# Calculate temporal coverage of trapping
trapping_duration <- ddply(.data = x2,
                           .variables = .(TrapID),
                           .fun = summarise, 
                           sampling_duration = as.numeric(diff(range(date))),
                           n_samples = length(date),
                           sampling_dens = (n_samples * 100 )/ 52   )
trapping_duration

trapping_duration <- ddply(.data = x2,
                           .variables = .(TrapID),
                           .fun = function(x){
                             sampling_duration = as.numeric(diff(range(x$date)))
                             n_samples = length(x$date)
                             sampling_dens = (n_samples * 100 )/ 52  
                             return(data.frame(sampling_duration,
                                               n_samples,
                                               sampling_dens))
                           })

## 
otu <- read_excel(file.path(main.dir, "MandaiAllSeq_Trapwise_JP_17Aug2021_Ver3_edited.xlsx"),
                  sheet = "All Seq")
otu2 <- subset(otu, `Trap Type` == "MT")
all(unique(otu2$TrapID_type_date) %in% x2$TrapID_type_date)

#~/Dropbox/Projects/Special Projects/Mandai/

## Spatial stuff
library(terra)
trees <- vect(x = file.path(main.dir, "TreesMPD_202001.shp"))

ggplot(data = as.data.frame(trees)) + geom_point(aes(y = Lat, x = Long))
range(trees$Long)

