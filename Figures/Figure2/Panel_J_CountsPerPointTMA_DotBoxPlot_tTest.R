library(EBImage)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(viridisLite)
library(scales)

allTifsData = data.frame()

outputFolder = ("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/RResults/BoxScatterplotADScaled")
dir.create(file.path(outputFolder), showWarnings = FALSE, recursive = T)

##TMA Healthy ADpanel 512 HIPDG (pt21-24,may2019), Striatum(pt41-44,may2019), Midbrain(pt5-8,may2019), Locus Coeruleus(p45-48,may2019), Medulla(pt54-57,may2019), Cerebellum (pt64-67,may2019)
setwd("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/AllADPanelTMAData/NoAuNoNoiseFFT_new/no_fftnoise_2X")
dataSize=512
pointList= c(5,6,7,8,21,22,23,24)# H(SN HIP)
pointList= c(13,14,15,16,25,26,27,28)#  AD (HIP MTG)
pointList= c(1,2,3,4,9,10,11,12)# LDB (SN EC)
N = length(pointList)
path = list.files(paste0("Point",pointList[1],"/TIFs/"))
markers = gsub("\\.tif", "", basename(path))

##TMA Healthy PDPanel 512 HIPDG (pt21-24,june2019),Midbrain(pt9-12,june2019), Striatum(pt17-20,june2019), Locus Coeruleus(p5-8,june2019), Medulla(pt1-4,may2019), Cerebellum (pt63-66,may2019)
setwd("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/AllPDPanelTMAData/CleanedNoNoiseNoBGHSNLBDSNLBDEC/no_fftnoisebleadremoval")
dataSize=512
pointList=c(9,10,11,12,21,22,23,24)# H(SN) H(Hip)
pointList=c(33,34,35,36,29,30,31,32)# AD (HIP MTG)  
pointList=c(37,38,39,40,41,42,43,44)# LDB(SN) LBD(EC)
N = length(pointList)
path = list.files(paste0("Point",pointList[1],"/TIFs/"))
markers = gsub("\\.tif", "", basename(path))

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

plot_data = tifData

boxplotMarkers = c("PanASyn","pASyn","ParkinPARK2","LRRK2Parkin8", "PARK7DJ1","PINK1PARK6") #PD Disease
#boxplotMarkers = c("PanAmyloidbeta1724", "PHF1Tau","ApoE4","8OHGuano", "pTDP43", "PolyubiK63") #AD Disease

plot_data = tifData[, boxplotMarkers]  

plot_data$SampleType = "Healthy"#___________________________________________________________________________________________change
allTifsData = rbind(allTifsData, plot_data)
#_________________________________________________________________________________________________________________________________________________Run till here for each sample H 

plot_data$SampleType = "AD"#___________________________________________________________________________________________change
allTifsData = rbind(allTifsData, plot_data)
#_________________________________________________________________________________________________________________________________________________Run till here for each sample AD 

plot_data$SampleType = "PD"#________________________________________________________________________________________change
allTifsData = rbind(allTifsData, plot_data)
#_________________________________________________________________________________________________________________________________________________Run till here for each sample  PD

colnames(allTifsData) = gsub("8OHGuano", "X8OHGuano", colnames(allTifsData))

write.csv(allTifsData, paste0(outputFolder, "/allTifsData_PDMarkers.csv") )

#Start running from  here when generated
allTifsData = read.csv("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/MedResScan_TMA/RResults/BoxScatterplotADScaled/allTifsData_ADMarkers.csv")

boxplotMarkers = gsub("8OHGuano", "X8OHGuano", boxplotMarkers)

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

ggplot(scatterPlotdata,aes(x=SampleType, y=value, fill = SampleType)) +
  #geom_boxplot(outlier.alpha = 0) +
  #geom_violin() +
  #geom_point(cex = 2.5, position = position_jitter(height =0 , width = 0.17)) +
  geom_beeswarm (cex= 3.2, size = 3, color = "black", pch =21)+ # to get more functions type ?geom_beeswarm in console
  #geom_quasirandom(cex = 5,size = 5, color = "black", pch =21)+
  scale_fill_manual(values  =c("#606060", "#00AEEF", "#BE1E2D"))+
  #geom_errorbar(stat = "identity", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..)) +
  stat_summary(fun= mean, geom = "crossbar", size=0.3, width = 0.7, color="black")+
  facet_wrap (~variable, scale = "free_y") + 
  # scale_y_continuous(trans=log10_trans(),
  #                    breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = trans_format("log10", math_format(10^.x)))+
  #limits = c(0.00001,0.00006) +
  ylim(c(0,0.00005)) +
  #ylim(c(0.00002,0.00005)) +
  theme_bw(15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename= paste0(outputFolder,"/ScatterBoxPlot_ADMarkers.pdf"), width = 10, height = 8)

#Shapiro-Wilk normality test for the differences and Paired t-Test
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
write.csv(statTestResults, paste0(outputFolder,'/PairedTtest.csv'))

#Shapiro-Wilk normality test for the differences and NonParametric Test
statTestResults1 = data.frame(matrix(NA,length(boxplotMarkers), 4))
colnames(statTestResults1) = c("Normality","Healthy_AD","Healthy_PD","AD_PD")
rownames(statTestResults1) = boxplotMarkers
for (marker in boxplotMarkers){
  res1 = shapiro.test(allTifsData[,marker])
  statTestResults1 [marker,1] = res1$p.value
  # Compute wilcox-test for Healthy_AD
  test_data1 = subset(allTifsData, SampleType %in% c("Healthy", "AD"))
  res1 = wilcox.test(as.formula(paste0(marker,"~SampleType")), test_data1, paired = TRUE)
  statTestResults1 [marker,2] = res1$p.value
  # Compute wilcox-test for Healthy_PD
  test_data1 = subset(allTifsData, SampleType %in% c("Healthy", "PD"))
  res1 = wilcox.test(as.formula(paste0(marker,"~SampleType")), test_data1, paired = TRUE)
  statTestResults1 [marker,3] = res1$p.value
  # Compute wilcox-test for AD_PD
  test_data1 = subset(allTifsData, SampleType %in% c("AD", "PD"))
  res1 = wilcox.test(as.formula(paste0(marker,"~SampleType")), test_data1, paired = TRUE)
  statTestResults1 [marker,4] = res1$p.value
}
write.csv(statTestResults1, paste0(outputFolder,'/NonParametricTest.csv'))
