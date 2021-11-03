###
dataPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/"
outputFolder = paste0(dataPath,"/Pooled")
dir.create(outputFolder, showWarnings= F, recursive = T)

all_data_modMore <- read.csv("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/Pooled/all_data_mod_unique_cell_id.csv")

all_data_mod = data.frame()
listOfPoints = sort(unique(all_data_modMore$Point))
for (filename in listOfPoints) {
  print(filename)
  
  all_data_mod_subset = subset(all_data_modMore, Point == filename) # if doing single run replace SampleNames with Point
  
  for (cell_id in unique(all_data_mod_subset$cellLabelInImage)){
    cell = subset(all_data_mod_subset, cellLabelInImage == cell_id)
    unicell = cell[sample(1:nrow(cell),1),]
    all_data_mod = rbind(all_data_mod, unicell)
  }
  
}
all_data_mod$NumericCode = as.numeric(as.factor(all_data_mod$MantisPopulation))

##to check for duplicates
sum(duplicated(all_data_mod_subset[,c(41,44)])) # for cellLabelInImage (41) and Point (44)
test = subset(all_data_mod, duplicated(all_data_mod[,c(41,44)]))

all_data_mod$X =  NULL
write.csv(all_data_mod, paste0(outputFolder, "/all_data_mod_unique_cell_id.csv"))
