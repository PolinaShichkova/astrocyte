library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

# Load FPKM data with supplementary info
fpkm_dt = fread("/Users/polina/Documents/EPFL/TA_GB/genomics_bioinf_project/01_agingAstrocyte/fpkm_mmc3.txt")
colnames(fpkm_dt)
# if you use data from GSE, check which columns do you need. I used data from supplementary of the paper. 
#GSE and Supp xlsx are almost the same, but there is some additional info in Supp xlsx
fpkm_dt <- fpkm_dt[,-41] # drop duplicated column Gene name
colnames(fpkm_dt)[5:40] # columns with the main data


sapply(fpkm_dt, mode) # check data types in columns


# consider to add pseudocounts to avoid NaN problem with log and to decrease the effect of noise of low-expressed genes. 
# I use pseudcount of 0.1 here

fpkm_dt_s = fpkm_dt %>%
  select( colnames(fpkm_dt)[5:40] ) + 0.1

fpkm_dt_s = cbind(fpkm_dt_s, fpkm_dt$`Gene name`)

fpkm_dt_sf = fpkm_dt_s[apply(fpkm_dt_s[,c(1:36)]>30.1,1,any),] # filter low-expressed
fpkm_dt_L = log(fpkm_dt_sf[,c(1:36)]) # log transform

fpkm_dt_L = cbind(fpkm_dt_L, fpkm_dt_sf$V2) #fpkm_dt_L = cbind(fpkm_dt_L, fpkm_dt_s3$`Gene name`)


##

# check normalization

# I would like to compare MC samples in both young and old ages ("4mo MC1"  "4mo MC2"  "4mo MC3"  "2yo MC 1" "2yo MC 2" "2yo MC 3")

library(GGally)

fpkmData <-fpkm_dt[,5:40]

pdf("output/pairsAstrocyte_MC_12may2019.pdf") 
ggpairs(fpkmData[,c(1:3,16:18)])  # from library(GGally)
#pairs(fpkmData[,c(1:3,16:18)], pch = 19) # similar fun from base R
dev.off() 

pdf("output/pairsAstrocyte_yHTHoCB_12may2019.pdf") 
ggpairs(fpkmData[,c(10:12,25:27)])
dev.off() 

fpkmData <-fpkm_dt_L[,1:36] # log transformed data

pdf("output/pairsAstrocyte_logMC_12may2019.pdf") 
ggpairs(fpkmData[,c(1:3,16:18)])
dev.off() 

pdf("output/pairsAstrocyte_logyHTHoCB_12may2019.pdf") 
ggpairs(fpkmData[,c(10:12,25:27)])
dev.off() 


# another way to plot the data, with * significance (you will need to explain it in your reports if you choose to plot this way) 
library(PerformanceAnalytics)

fpkmData <-fpkm_dt[,5:40]

pdf("output/pairsAstrocyte_perfAn_MC_12may2019.pdf") 
chart.Correlation(fpkmData[,c(1:3,16:18)])  # from library(GGally)
#pairs(fpkmData[,c(1:3,16:18)], pch = 19) # similar fun from base R
dev.off() 

pdf("output/pairsAstrocyte_perfAn_yHTHoCB_12may2019.pdf") 
chart.Correlation(fpkmData[,c(10:12,25:27)])
dev.off() 

fpkmData <-fpkm_dt_L[,1:36] # log transformed data

pdf("output/pairsAstrocyte_perfAn_logMC_12may2019.pdf") 
chart.Correlation(fpkmData[,c(1:3,16:18)])
dev.off() 

pdf("output/pairsAstrocyte_perfAn_logyHTHoCB_12may2019.pdf") 
chart.Correlation(fpkmData[,c(10:12,25:27)])
dev.off() 


## ! What can you say from the plots above?

######

# reproduce fig 2b

library(gplots)   # contains the heatmap.2 package

fpkmData <-fpkm_dt_L[,1:36] # log transformed data
scaled <- scale(fpkmData, center = TRUE, scale = TRUE) # for z scores
tr = t(scaled) # transpose

# heatmap with ward.D2

pdf("output/2b_heatmap.pdf") 
heatmap.2(tr, 
          cexRow=0.5, cexCol=0.95, 
          scale="none", 
          trace="none",
          col = bluered(100),
          hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off() 

##############

# figure 7c (part), simplified demo

fpkm_dt_s = fpkm_dt %>%
  select( colnames(fpkm_dt)[5:40] ) + 0.1

fpkm_dt_s = cbind(fpkm_dt_s, fpkm_dt$`Gene name`)
colnames(fpkm_dt_s)[37] = "Gene name"

fpkm_dt_s = cbind(fpkm_dt_s, fpkm_dt$`Annotation/Divergence`) # find gene names in this column
colnames(fpkm_dt_s)[38] = "Annotation/Divergence"

fpkm_dt_sf = fpkm_dt_s[apply(fpkm_dt_s[,c(1:36)]>0.1,1,any),]  # here changed threshold to find all genes I'm interested in

#find proper gene in fpkm_dt_L$Annotation
Thbs1 = fpkm_dt_sf[fpkm_dt_sf$`Annotation/Divergence` %like% "Thbs1", ]
Thbs2 = fpkm_dt_sf[fpkm_dt_sf$`Annotation/Divergence` %like% "Thbs2", ]
Thbs4 = fpkm_dt_sf[fpkm_dt_sf$`Annotation/Divergence` %like% "Thbs4", ]
Sparcl1 = fpkm_dt_sf[fpkm_dt_sf$`Annotation/Divergence` %like% "Sparcl1", ]

synapseFormation = rbind(Thbs1,Thbs2,Thbs4,Sparcl1)
colnames(synapseFormation)

rowMedians(synapseFormation[c("4mo MC1","4mo MC2","4mo MC3","4mo SC 1","4mo SC 2","4mo SC 3","4mo VC1","4mo VC2","4mo VC3" )])

synapseFormation = synapseFormation %>% 
  mutate(mean_ctx = rowMeans(select(., c("4mo MC1","4mo MC2","4mo MC3","4mo SC 1","4mo SC 2","4mo SC 3","4mo VC1","4mo VC2","4mo VC3" ))),
         mean_cb = rowMeans(select(., c("4mo CB1", "4mo CB2","4mo CB3"  ))),
         mean_hth = rowMeans(select(., c("4mo HTH1","4mo HTH2","4mo HTH3" ))))


synapseFormation$hth_ctx = log2(synapseFormation$mean_hth/synapseFormation$mean_ctx)
synapseFormation$cb_ctx = log2(synapseFormation$mean_cb/synapseFormation$mean_ctx)

#####

library(formattable)
library(DT)

# create 19 breaks and 20 rgb color values ranging from blue to red
brks <- quantile(synapseFormation[,c('hth_ctx', 'cb_ctx')], probs = seq(.05, .95, .05), na.rm = TRUE)

clrs_red <- round(seq(255, 0, length.out = length(brks)/2), 0) %>%
{paste0("rgb(255,", ., ",", ., ")")}

clrs_blue <- round(seq(0, 255, length.out = length(brks)/2 + 1), 0) %>%
{paste0("rgb(", ., ",", ., ", 255 )")}

clrs = c(clrs_blue,clrs_red[2:length(clrs_red)])

datatable(synapseFormation[,c('hth_ctx', 'cb_ctx')]) %>% formatStyle(names(synapseFormation[,c('hth_ctx', 'cb_ctx')]), backgroundColor = styleInterval(brks, clrs))








