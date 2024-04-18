

## PACKAGES ============
library(readxl)
library(ggplot2)
library(plyr)

## DIRECTORIES ============
main.dir <- "~/Dropbox/Projects/Special Projects/Mandai"
fig.dir <- file.path(main.dir, "figures")

## Checking trap metadata
otu <- read_excel(file.path(main.dir, "MandaiAllSeq_Trapwise_JP_17Aug2021_Ver3_edited.xlsx"),
                  sheet = "All Seq")

names(otu)
ddply(.data = otu, 
      .variables = .(trapid, `Coll. Date`),
      .fun = summarise,
      otu_div = length(unique(`3%id`)),
      abundance = length(`3%id`))


## Expanding
library(reshape2)
head(otu)
x <- acast(data = otu, formula = trapid ~ `3%id`, fun.aggregate = length,
           value.var = "3%id")
x_pa <- vegan::decostand(x, method = "pa")
rowSums(x_pa)
colSums(x_pa)

plot(vegan::specaccum(x, method = "rarefaction"))
?specaccum

View(x_pa)

subset(otu, trapid == "MIS-L07" & `3%id` == "ESCO-SG007-A")


View(x)
