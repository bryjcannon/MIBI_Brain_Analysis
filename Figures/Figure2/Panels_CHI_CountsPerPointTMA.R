library(EBImage)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(viridisLite)
library(scales)
library(circlize)

allTifsData = data.frame()

outputFolder = ("/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/RResults/HealthyTMA/ADPanel/StructuralMarkers/")
dir.create(outputFolder, showWarnings= F, recursive = T)

##TMA Healthy 1024 HIPCA1 (pt9,nov2018),Midbrain(pt2,feb2019), Striatum(pt5,feb2019), Locus Coeruleus(pt8,nov2018), Medulla(pt1 (pt5,nov2018)), Cerebellum (pt6,nov2018)
#setwd("/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/HiResScan_TMA/TMADataFromNov&Feb/AllADPanelTMAHealthyHiRes/no_fftnoiseolder")
#dataSize=1024
#pointList=c(9,2,5,8,1,6)
# N = length(pointList)
# path = list.files(paste0("Point",pointList[1],"/TIFs/"))
# markers = gsub("\\.tif", "", basename(path))

##TMA Healthy ADpanel 512 HIPDG (pt21-24,may2019), Striatum(pt41-44,may2019), Midbrain(pt5-8,may2019), Locus Coeruleus(p45-48,may2019), Medulla(pt54-57,may2019), Cerebellum (pt64-67,may2019)
setwd("/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/AllADPanelTMAData/NoAuNoNoiseFFT_new/no_fftnoise_2X")
dataSize=512
 pointList=c(21,22,23,24,41,42,43,44,5,6,7,8,45,46,47,48,54,55,56,57,64,65,66,67) #All Healthy Brain areas
# pointList= c(13,14,15,16,25,26,27,28,1,2,3,4,9,10,11,12,21,22,23,24,5,6,7,8,41,42,43,44,45,46,47,48,54,55,56,57,64,65,66,67)#ALL H and AD and LDB
# pointList= c(64,65,66,67,5,6,7,8,21,22,23,24,13,14,15,16,25,26,27,28,1,2,3,4,9,10,11,12)# H(Cerebellum,SN HIP) and AD (HIP MTG) and LDB (SN EC)
# pointList= c(5,6,7,8,21,22,23,24,13,14,15,16,25,26,27,28,1,2,3,4,9,10,11,12)# H(SN HIP) and AD (HIP MTG) and LDB (SN EC)
# pointList= c(21,22,23,24,13,14,15,16)# H(HIP) and AD (HIP)
# pointList= c(13,14,15,16,25,26,27,28)# AD(HIP) and AD (MTG)
# pointList= c(21,22,23,24,13,14,15,16,25,26,27,28)# H(HIP), AD(HIP) and AD (MTG)
# pointList= c(5,6,7,8,21,22,23,24)# H(SN HIP)
# pointList= c(13,14,15,16,25,26,27,28)#  AD (HIP MTG)
# pointList= c(1,2,3,4,9,10,11,12)# LDB (SN EC)
N = length(pointList)
path = list.files(paste0("Point",pointList[1],"/TIFs/"))
markers = gsub("\\.tif", "", basename(path))

##TMA Healthy PDPanel 512 HIPDG (pt21-24,june2019),Midbrain(pt9-12,june2019), Striatum(pt17-20,june2019), Locus Coeruleus(p5-8,june2019), Medulla(pt1-4,may2019), Cerebellum (pt63-66,may2019)
#setwd("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/AllPDPanelTMAData/no_background_Au_pASyn2x_PD") # Option 1
setwd("/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/AllPDPanelTMAData/CleanedNoNoiseNoBGHSNLBDSNLBDEC/no_fftnoisebleadremoval") #Option 2
dataSize=512
# pointList=c(21,22,23,24,9,10,11,12,17,18,19,20,5,6,7,8,1,2,3,4,63,64,65,66)#All Healthy Brain areas
# # # #pointList=c(9,10,11,12,21,22,23,24,33,34,35,36,29,30,31,32,37,38,39,40,41,42,43,44)# H(SN HIP) and AD (HIP MTG) and LDB (SN EC) 
# # #pointList=c(37,38,39,40,41,42,43,44)# LDB (SN EC) 
# # #pointList=c(9,10,11,12,37,38,39,40)# H(SN) and LDB (SN) 
 pointList=c(9,10,11,12,37,38,39,40,41,42,43,44)# H(SN) and LDB (SN)  and LBD (EC)
# pointList=c(9,10,11,12,21,22,23,24)# H(SN) H(Hip)
# #pointList=c(33,34,35,36,29,30,31,32)# AD (HIP MTG)  
# #pointList=c(37,38,39,40,41,42,43,44)# LDB(SN) LBD(EC)
N = length(pointList)
path = list.files(paste0("Point",pointList[1],"/TIFs/"))
markers = gsub("\\.tif", "", basename(path))

##TMA AD PD Disease vs Healthy 1024 HIPCA1_AD(pt11,nov2018), MTG_AD(pt7,feb2019), EC_LBD(pt2=pt3,aug2019), midbrain_LBD(pt1=pt4,aug2019),HIPCA1_H(pt9,feb2019),Midbrain_H(pt2,feb2019), Striatum_H(pt5,feb2019), Locus Coeruleus_H(pt8,nov2018), Medulla_H(pt1 (pt5,nov2018)), Cerebellum_H(pt6,nov2018)
# setwd("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/HiResScan_TMA/ADPDHealthyTMA")
# dataSize=1024
# pointList=c(1,2,3,4,5,6,7,8,9,11)
# path = list.files("Point11/TIFs/")
# N = length(pointList)
# path = list.files(paste0("Point",pointList[1],"/TIFs/"))
# markers = unique(intersect(gsub("\\.tif", "", basename(path)),gsub("\\.tif", "", basename(pathPD))))

##TMA AD Disease ADPanel 512 HIPCA1-AD (pt13-16),MTG_AD(pt25-28), Midbrain_LBD(pt1-4), EC_LBD(pt9-12)
# setwd("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/NoAuNoNoiseFFT_new/no_fftnoise")
# dataSize=512
# #pointList=c(13,14,15,16,25,26,27,28,1,2,3,4,9,10,11,12) #ALL DISEASE
# #pointList=c(13,14,15,16,25,26,27,28) #AD DISEASE
# pointList=c(1,2,3,4,9,10,11,12) #PD DISEASE
# N = length(pointList)
# path = list.files(paste0("Point",pointList[1],"/TIFs/"))
# markers = gsub("\\.tif", "", basename(path))

##TMA PD Disease PDPanel 512 HIPCA1-AD (pt33-36),MTG_AD(pt29-32), Midbrain_LBD(pt37-40), EC_LBD(pt41-44)
# setwd("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/AllPDPanelTMAData/no_background_Au_pASyn2x_PD")
# dataSize=512
# #pointList=c(33,34,35,36,29,30,31,32,37,38,39,40,41,42,43,44) #ALL DISEASE
# #pointList=c(33,34,35,36,29,30,31,32) #AD DISEASE
# pointList=c(37,38,39,40,41,42,43,44) #PD DISEASE
# N = length(pointList)
# path = list.files(paste0("Point",pointList[1],"/TIFs/"))
# markers = gsub("\\.tif", "", basename(path))

tifData <- data.frame(matrix(NA, N, length(markers)))
colnames(tifData) = markers
rownames(tifData) = paste0("Point",pointList )
index=1
for (PointNumber in pointList){
  path = paste0("Point", PointNumber, "/TIFs")
  
  allTifs =  list.files(path = path, pattern = paste0(".*tif$"), full.names = T)
  
# Import original TIFs to be used for training and get total signal count per channel
  for (tifName in allTifs) {
    if (!(gsub("\\.tif", "", basename(tifName)) %in% markers))
      next
    print(paste0("Processing: ", tifName))
    tif_image <- readImage(tifName)
    
# Normalise signal intensites before summing by:
    
# 1. Unit vector calculation for narmalisation across different runs (Batch effect). need to change the way data is imported. wont work in this format
    # for (i in 1:N){
    #    r = sqrt(sum(plot_data[i,] * plot_data[i,], na.rm=T))
    #   plot_data[i,] = plot_data[i,]/r
    #  }
    
# 2.1  99th Percentile narmalisation for normalisation individual marker expression
    fcsquantile <- quantile(tif_image, 0.99)
      
# 2.2 If all values are zero, set quantile to 1 to avoid Undefined values
    fcsquantile[fcsquantile == 0] <- 1
    # between zero and max signal/99th percentile
    tif_image <- t(apply(tif_image, 1, function(x) {
      x / fcsquantile
    })) 
    
# creats table of pixel intensity and counts
    countTable = data.frame(table(tif_image))
    if (length(countTable$tif_image) == 1) # when zeros only in tif image sets sum to zero
      sumCounts = 0
    else {
      countTable$Counts = round(as.numeric(as.character(countTable$tif_image)) / as.numeric(as.character(countTable[2,1])))
      countTable$tif_image = NULL
      
      #have single values
      countTable = data.frame(countTable)
      # colnames(countTable) = countTable[2,]
      # countTable = countTable[1,]
      sumCounts = sum(apply(countTable, 1, function(x) {
        x[1] * x[2]
      }))
    
      ##to calculate total mean tif intensity
      #sumCounts = sumCounts/(dataSize*dataSize) 
      
      ##to calculate log transform total mean tif intensity
      sumCounts = log(sumCounts+1)/(dataSize*dataSize)
      
      ##to calculate mean non-zero tif intensity
      #sumCounts = sumCounts/(sum(countTable[2:nrow(countTable),"Freq"])) 
     }
    
    tifData[index, gsub("\\.tif", "", basename(tifName))] = sumCounts
  }
  index=index+1
}
    
tifData$TIFs = NULL

##For Single AD or PD Panel All Good Markers Healthy & Diseased Brain Region Heatmaps
plot_data = tifData

## To report All Good Markers
#exclude = "Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|PanApoE2E3E4" 

##To report Diease AD Markers
#exclude = "Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|PanApoE2E3E4|CD47|Calbindin|Calretinin|Presenilin1NTF|Presenilin1|CD31|CD33Lyo|CD45|CD56Lyo|CD105|GFAP|Iba1|MAG|MAP2|MBP|MCT1|PanGAD6567|Parvalbumin|PSD95|Synaptophysin|TH|VGAT|VGLUT1|VGLUT2|TotalTau|HistoneH3Lyo" 

## To report Diease PD Markers
#exclude = "Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|PanApoE2E3E4|CD31|CD33Lyo|CD45|CD47|CD56Lyo|CD105|DAT|GFAP|Iba1|MAG|MAP2|MBP|MCT1|PanGAD6567|Parvalbumin|PSD95|Synaptophysin|VGAT|VGLUT1|VGLUT2|TotalTau|HistoneH3Lyo|PARK7DJ1|ParkinPARK2|PINK1PARK6|LRRK2Parkin8"

## To report Common PD Basic
#exclude = "Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|PanApoE2E3E4|ApoE4|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|LRRK2Parkin8|PolyubiK63|PolyubiK48|MFN2|pTDP43|PanASyn|PHF1Tau|PanAmyloidbeta1724|Amyloidbeta142|Amyloidbeta140|8OHGuano" 

## To report Common AD Basic
exclude = "Ca40|C12|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|PanApoE2E3E4|ApoE4|PolyubiK63|PolyubiK48|MFN2|pTDP43|PHF1Tau|PanAmyloidbeta1724|Amyloidbeta142|Amyloidbeta140|8OHGuano"

exclude = grep(exclude, colnames(plot_data))
plot_data = plot_data[,-exclude]
allTifsData = rbind(allTifsData, plot_data)

##For Heatmap from two different runs
# plot_data = tifData
# plot_data$SampleType = "TMA_HealthyAD" # need to change when joining runs Control1, Control2
# exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|PanApoE2E3E4"
# exclude = grep(exclude, colnames(plot_data))
# plot_data = plot_data[,-exclude]
# rownames(plot_data) = paste0(rownames(plot_data), "_" , plot_data$SampleType)
# rownames(plot_data) = gsub("1_PD", "_PD", rownames(plot_data))
# allTifsData = rbind(allTifsData, plot_data)
# plot_data$SampleType = NULL

##For Joint AD and PD Panel Healthy Heatmaps
# plot_data = tifData
# plot_data$SampleType = "PD"
# exclude = "C12|Ca40|Na23|Si28|Ta181|totallon|empty113|empty139|Au197|Background|Reelin|SERT|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|PanApoE2E3E4|MFN2|pTDP43|TREM2|EEA1|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|LRRK2Parkin8|DAT|PanAsyn|PHF1Tau|Amyloidbeta142|Amyloidbeta140|Calbindin|Calretinin|PanASyn|Parvalbumin"
# exclude = grep(exclude, colnames(plot_data))
# plot_data = plot_data[,-exclude]
# allTifsData = rbind(allTifsData, plot_data)

##For Single AD or PD Panel Common Marker Healthy Brain Region Heatmaps
# plot_data = tifData
# exclude = "C12|Ca40|Na23|Si28|Ta181|totallon|empty113|empty139|Au197|Background|Reelin|SERT|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|PanApoE2E3E4|ApoE4|MFN2|pTDP43|TREM2|EEA1|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|LRRK2Parkin8|DAT|PanAsyn|PHF1Tau|PanAmyloidbeta1724|Amyloidbeta142|Amyloidbeta140|PanASyn|8OHGuano"
# exclude = grep(exclude, colnames(plot_data))
# plot_data = plot_data[,-exclude]
# allTifsData = rbind(allTifsData, plot_data)

##For ScatterPlot of AD and PD markers
#exclude = "C12|Ca40|Na23|Si28|Ta181|totallon|empty113|empty139|Au197|Background|Reelin|SERT|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|PanApoE2E3E4|MFN2|pTDP43|TREM2|EEA1|MBP|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|LRRK2Parkin8|DAT|PanAsyn"
#exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|PanApoE2E3E4|CD31|CD33Lyo|CD45|CD47|CD56Lyo|CD105|GFAP|Iba1|MAG|MAP2|MBP|MCT1|PSD95|Synaptophysin|TotalTau|HistoneH3Lyo|DAT|PanGAD6567|Parvalbumin|VGAT|VGLUT1|VGLUT2|PARK7DJ1|ParkinPARK2|PINK1PARK6|LRRK2Parkin8|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|MFN2|TREM2|EEA1|Calbindin|Calretinin|Parvalbumin"
#exclude = grep(exclude, colnames(tifData))

##PD Markers Only
#exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|PanApoE2E3E4|CD31|CD33Lyo|CD45|CD47|CD56Lyo|CD105|DAT|GFAP|Iba1|MAG|MAP2|MBP|MCT1|PanGAD6567|Parvalbumin|PSD95|Synaptophysin|VGAT|VGLUT1|VGLUT2|TotalTau|HistoneH3Lyo|PARK7DJ1|ParkinPARK2|PINK1PARK6|Calbindin|Calretinin"
#exclude = grep(exclude, colnames(tifData))

##AD Markers Only
 # exclude = "C12|Ca40|Na23|Si28|Ta181|totallon|empty113|empty139|Au197|Background|Reelin|SERT|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|PanApoE2E3E4|MFN2|TREM2|EEA1Calbindin|Calretinin|Parvalbumin"
 # exclude = grep(exclude, colnames(tifData))

#plot_data = tifData[,-exclude]
boxplotMarkers = c("PanASyn", "pASyn","PanAmyloidbeta1724", "PHF1Tau")
#plot_data = tifData[, c("PanAmyloidbeta1724", "PHF1Tau","ApoE4","8OHGuano")]#AD Disease
plot_data = tifData[, boxplotMarkers] #PD Disease ,"ApoE4","8OHGuano"
plot_data$SampleType = "Healthy"#___________________________________________________________________________________________change
#plot_data$SampleType = "AD"#___________________________________________________________________________________________change
#plot_data$SampleType = "PD"#___________________________________________________________________________________________change
allTifsData = rbind(allTifsData, plot_data)
#_________________________________________________________________________________________________________________________________________________Run till here for each sample H AD PD

#tifData =  tifData[, c("PanAmyloidbeta1724", "PHF1Tau","ApoE4","8OHGuano")]#AD Disease
#tifData =  tifData[, c("PanASyn", "pASyn","8OHGuano","PanAmyloidbeta1724", "PHF1Tau","ApoE4")]#PD Disease
#plot_data =  plot_data[, c("PanASyn", "pASyn")]

tifData = allTifsData
scatterPlotdata =tifData
scatterPlotdata$SampleType = allTifsData$SampleType

##for having same sample numbers in ScatterboxPlot
# healthyData = subset(scatterPlotdata, SampleType == "Healthy")
# scatterPlotdata = subset(scatterPlotdata, SampleType != "Healthy")
# scatterPlotdata = rbind(scatterPlotdata,healthyData[sample(1:nrow(healthyData), (nrow(scatterPlotdata)/2)), ])
##for having same sample numbers in ScatterboxPlot

scatterPlotdata = melt(scatterPlotdata, ids = c("SampleType"))
scatterPlotdata$SampleType = as.factor(scatterPlotdata$SampleType)
scatterPlotdata$SampleType = ordered(scatterPlotdata$SampleType, levels =c("Healthy", "AD", "PD"))
ggplot(scatterPlotdata,aes(x=SampleType, y=value)) +
  geom_boxplot(outlier.alpha = 0) +
  #geom_violin() +
  geom_point() +
  facet_wrap (~variable, scale = "free_y") + 
  scale_y_continuous(trans=log10_trans(),
          breaks = trans_breaks("log10", function(x) 10^x),
          labels = trans_format("log10", math_format(10^.x)))+
          #limits = c(0.00001,0.00006)) +
  #ylim(c(0,0.00005)) +
  #ylim(c(0.00003,0.00004)) +
  theme_bw(15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename= paste0(outputFolder,"/ScatterBoxPlot_ADPanel.pdf"), width = 6, height = 6)

#Paired t-test
#Test normality of distribution
# Shapiro-Wilk normality test for the differences and paired t=test
statTestResults = data.frame(matrix(NA,length(boxplotMarkers), 4))
colnames(statTestResults) = c("Normality","Healthy_AD","Healthy_PD","AD_PD")
rownames(statTestResults) = boxplotMarkers
for (marker in boxplotMarkers){
  res = shapiro.test(allTifsData[,marker])
  statTestResults [marker,1] = res$p.value
  # Compute t-test for Healthy_AD
  test_data = subset(allTifsData, SampleType %in% c("Healthy", "AD"))
  res = t.test(as.formula(paste0(marker,"~SampleType")), test_data, paired = TRUE)
  statTestResults [marker,2] = res$p.value
  # Compute t-test for Healthy_PD
  test_data = subset(allTifsData, SampleType %in% c("Healthy", "PD"))
  res = t.test(as.formula(paste0(marker,"~SampleType")), test_data, paired = TRUE)
  statTestResults [marker,3] = res$p.value
  # Compute t-test for AD_PD
  test_data = subset(allTifsData, SampleType %in% c("AD", "PD"))
  res = t.test(as.formula(paste0(marker,"~SampleType")), test_data, paired = TRUE)
  statTestResults [marker,4] = res$p.value
}

##EXCLUSIONS OF MARKERS LIST
##For Healthy using AD and PD Panel (512 image)
#exclude = "C12|Ca40|Na23|Si28|Ta181|totallon|empty113|empty139|Au197|Background|Reelin|SERT|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|PanApoE2E3E4|MFN2|pTDP43|TREM2|EEA1|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|LRRK2Parkin8|DAT|PanAsyn|PHF1Tau|Amyloidbeta142|Amyloidbeta140"
##For Healthy Brain Tissue using AD Panel (512 image)
#exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|pTDP43|PolyubiK63|PolyubiK48|PanApoE2E3E4|MFN2|PHF1Tau|Amyloidbeta142|Amyloidbeta140|Presenilin1NTF"
#exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|PanApoE2E3E4"
##For Healthy Brain Tissue using PD Panel (512 image)
#exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|DAT|EEA1|LRRK2Parkin8|MFN2|PanApoE2E3E4|PanAsyn|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|SERT|PolyubiK63|pTDP43|PHF1Tau|Amyloidbeta142"
#exclude = "C12|Ca40|Na23|Si28|Ta181|Au197|totallon|empty113|empty139|Background|Reelin|SERT|EEA1|PanApoE2E3E4"
##For AD and PD Panel (1024 image)
#exclude = exclude = "C12|Ca40|Na23|Si28|Ta181|totallon|empty113|empty139|Au197|Background|Reelin|SERT|SorLA|Presenilin1NTF|Presenilin1|PolyubiK63|PolyubiK48|ApoE4|PanApoE2E3E4|MFN2|pTDP43|TREM2|EEA1|MBP|PARK7DJ1|ParkinPARK2|PINK1PARK6|pASyn|LRRK2Parkin8|DAT|PanAsyn|PHF1Tau|Amyloidbeta142|Amyloidbeta140|HistoneH3Lyo"

#plot histograms of all channels all points
plot_dataHist = melt(plot_data[,colnames(plot_data)])#all channels all points
colnames(plot_dataHist) = c("Markers","LogCellIntensities")
ggplot(plot_dataHist, aes(x=LogCellIntensities)) +
  geom_histogram(bins = 30)+
  facet_wrap(Markers~., scales="free")
ggsave(filename = paste0(outputFolder,"/MarkersHistograms.pdf"), width = 20, height = 20)


# To plot dendrogram & heatmaps
set.seed(240214)
hc = hclust(dist(plot_data), "ward.D2")
pdf(paste0(outputFolder,"/dendrogram_AD.pdf"), width = 10, height = 4, onefile = FALSE)
ggdendrogram(hc,rotate = FALSE, size = 2)
dev.off()

callback=function(x,data){
  hc=hclust(dist(data), "ward.D2")
  dend=reorder(as.dendrogram(hc),wts=data)
  as.hclust(dend,"ward.D2")
}
set.seed(240214)

## Scaled by column (able to compare Points). "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), and "viridis" (or "D", the default option).
## Cluster_cols = F to not cluster the brain regions

# UnScaled heatmap. 
# pdf(paste0(outputFolder,"/Heatmap_MeanIntensityUnScaled_AD.pdf"), width = 7, height = 7, onefile = FALSE)
# p=pheatmap((plot_data), cluster_rows = T, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "99th Percentile Normalised \nLog transform total mean tif intensities \nUnscaled scaled",
#            col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
dev.off()
pdf(paste0(outputFolder, "/HeatMap_AllMarkers_VIRIDIS2.pdf"), width = 7, height = 7, onefile = F)
pheatmap((plot_data), cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "99th Percentile Normalised \nLog transform total mean tif intensities \nUnscaled scaled", 
         color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
dev.off()

# RowScaled  (able to compare expression of markers within one point)
scaled_data = apply(plot_data, 1, function(x) (x - min(x,na.rm=T))/diff(range(x,na.rm=T))) #puts all values 0-1
#scaled_data = scaled_data[p$tree_row$order,p$tree_col$order] #applies clustering of Basic Heatmap

pdf(paste0(outputFolder,"/Heatmap_MeanIntensityRowScale_AD.pdf"), width = 7, height = 7, onefile = FALSE)
p=pheatmap((scaled_data), cluster_rows = T, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "99th Percentile Normalised \nLog transform total mean tif intensities \nPoint/Region scaled",
           color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
          #col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

# ColumnScaled (able to compare expression of markers across different points).
scaled_data = apply(plot_data, 2, function(x) (x - min(x,na.rm=T))/diff(range(x,na.rm=T)))#puts all values 0-1
#scaled_data = scaled_data[p$tree_row$order,p$tree_col$order] #applies clustering of Basic Heatmap

pdf(paste0(outputFolder,"/Heatmap_MeanIntensityColumnScale_AD.pdf"), width = 7, height = 7, onefile = FALSE)

p=pheatmap((scaled_data), cluster_rows = T, cluster_cols = T,  clustering_callback = callback, border_color = FALSE, scale = "none",main = "99th Percentile Normalised \nLog transform total mean tif intensities \nMarker scaled",
           color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
           #col = colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

# Z-Score ColumnScaled (able to compare expression of markers across different points).
scaled_data = apply(plot_data, 2, function(x) (x - min(x,na.rm=T))/diff(range(x,na.rm=T)))#puts all values 0-1
#scaled_data = scaled_data[p$tree_row$order,p$tree_col$order] #applies clustering of Basic Heatmap

pdf(paste0(outputFolder,"/Heatmap_zscore_ColumnScale_AD.pdf"), width = 7, height = 7, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000"))(100) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap(t(scaled_data), cluster_rows = T, cluster_cols = T,  clustering_callback = callback, border_color = FALSE, scale = "row",main = "Z-scores \nper FOV \n",
           col = colfunc, cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

# Transposed RowScaled (able to compare expression of markers within one point)
scaled_data = apply(plot_data, 1, function(x) (x - min(x,na.rm=T))/diff(range(x,na.rm=T))) #puts all values 0-1
#scaled_data = scaled_data[p$tree_row$order,p$tree_col$order] #applies clustering of Basic Heatmap

pdf(paste0(outputFolder,"/Heatmap_MeanIntensity_TransposedRowScale_AD.pdf"), width = 7, height = 7, onefile = FALSE)
p=pheatmap(t(scaled_data), cluster_rows = T, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "99th Percentile Normalised \nLog transform total mean tif intensities \nPoint/Region scaled",
           color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
           #col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

# Transposed ColumnScaled (able to compare expression of markers across different points).
scaled_data = apply(plot_data, 2, function(x) (x - min(x,na.rm=T))/diff(range(x,na.rm=T)))#puts all values 0-1
#scaled_data = scaled_data[p$tree_row$order,p$tree_col$order] #applies clustering of Basic Heatmap

pdf(paste0(outputFolder,"/Heatmap_MeanIntensity_TransposedColumnScale_AD.pdf"), width = 7, height = 7, onefile = FALSE)

p=pheatmap(t(scaled_data), cluster_rows = T, cluster_cols = T,  clustering_callback = callback, border_color = FALSE, scale = "none", main = "99th Percentile Normalised \nLog transform total mean tif intensities \nMarker scaled",
           color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
           #col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

write.csv(tifData, paste0(outputFolder, "/CountsPerPixel.csv"))


## Extraheatmap settings
#pdf(paste0(outputFolder,"/PLACENEWNAME.pdf"), width = 7, height = 14, onefile = FALSE)
#breakList = c(0.0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc = colorRamp2(breakList, c("darkblue", "white", "darkgreen"))	
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808"))
#p=pheatmap((plot_data), cluster_rows = T, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "Intensity summary, row scaled",
           #color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
           #col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
          #col=colfunc(length(breakList), cellwidth = 10, cellheight = 10)
#dev.off()
#dev.off()






