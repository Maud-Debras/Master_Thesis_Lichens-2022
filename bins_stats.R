#Set the working directory

#Import libraries
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(patchwork)
library(magrittr)
library(scales) 
library(colorspace)
library(ggdark)
library(ggthemes)

# Importing the excel workbook for checking the number of sheets it contains
  # getting sheets in this workbook
  sheets = excel_sheets("Binning-Tables.xlsx")
  length(sheets) #number of sampples
  df_complete = data.frame()

  #loop for all sheets
for (s in sheets)
  {
    df = read_excel("Binning-Tables.xlsx", sheet = s, col_names = TRUE, col_types = NULL)
    df_complete = rbind(df_complete, df)
  }
df_complete <- data.frame(df_complete)
View(df_complete)
colnames(df_complete) <- c("Sample","Binning_software","Bin","Total_length","Number_of_contigs","N50","GC_content","Coverage","CheckM_lineage","CheckM_contamination","CheckM_completeness","EukCC_lineage","EukCC_completeness","EukCC_contamination","GTDB_lineage")
df_complete$Sample <- as.character(df_complete$Sample)
df_complete$Binning_software <- as.character(df_complete$Binning_software)
df_complete$Bin <- as.numeric(df_complete$Bin)
df_complete$Total_length <- as.numeric(df_complete$Total_length)
df_complete$Number_of_contigs <- as.numeric(df_complete$Number_of_contigs)
df_complete$N50 <- as.numeric(df_complete$N50)
df_complete$GC_content <- as.numeric(df_complete$GC_content)
df_complete$Coverage <- as.numeric(df_complete$Coverage)
df_complete$CheckM_lineage <- as.character(df_complete$CheckM_lineage)
df_complete$CheckM_contamination <- as.numeric(df_complete$CheckM_contamination)
df_complete$CheckM_completeness <- as.numeric(df_complete$CheckM_completeness)
df_complete$EukCC_lineage <- as.character(df_complete$EukCC_lineage)
df_complete$EukCC_completeness <- as.numeric(df_complete$EukCC_completeness)
df_complete$EukCC_contamination <- as.numeric(df_complete$EukCC_contamination)
df_complete$GTDB_lineage <- as.character(df_complete$GTDB_lineage)
class(df_complete)

# Count the number of bins constructed by each software
sum(df_complete$Binning_software == "CONCOCT") #3288
sum(df_complete$Binning_software == "METABAT") #807

# Mean of bin per sample for each software
vec1_meanbins <- vector()
vec2_meanbins <- vector()
for (s in sheets)
{
  df_s <- subset(df_complete,df_complete$Sample == s)
  vec1_meanbins <- c(vec1_meanbins, sum(df_s$Binning_software == "CONCOCT"))
  vec2_meanbins <- c(vec2_meanbins, sum(df_s$Binning_software == "METABAT"))
}
mean(vec1_meanbins) #59.78182
mean(vec2_meanbins)  #14.65455
max(vec1_meanbins) #92
which.max(vec1_meanbins) #12
sheets[12] #S25
max(vec2_meanbins) #27
which.max(vec2_meanbins) #25
sheets[25]  #S38

# Subset dataframe for each binning tool
df_concoct <- subset(df_complete,df_complete$Binning_software == "CONCOCT")
df_metabat<- subset(df_complete,df_complete$Binning_software == "METABAT")

#Stats bins
nrow(subset(df_concoct,(df_concoct$CheckM_completeness > 0 | df_concoct$EukCC_completeness > 0))) #928
nrow(subset(df_metabat,(df_metabat$CheckM_completeness > 0 | df_metabat$EukCC_completeness > 0))) #576
nrow(subset(df_concoct,(df_concoct$CheckM_completeness > 90 | df_concoct$EukCC_completeness > 90))) #211
nrow(subset(df_metabat,(df_metabat$CheckM_completeness > 90 | df_metabat$EukCC_completeness > 90))) #157
nrow(subset(df_concoct,(df_concoct$CheckM_completeness > 50 | df_concoct$EukCC_completeness > 50))) #502
nrow(subset(df_metabat,(df_metabat$CheckM_completeness > 50 | df_metabat$EukCC_completeness > 50))) #328
nrow(subset(df_concoct,(df_concoct$CheckM_contamination < 5 & df_concoct$CheckM_completeness > 0) | (df_concoct$EukCC_contamination < 5 & df_concoct$EukCC_completeness > 0))) #704
nrow(subset(df_metabat,(df_metabat$CheckM_contamination < 5 & df_metabat$CheckM_completeness > 0) | (df_metabat$EukCC_contamination < 5 & df_metabat$EukCC_completeness > 0))) #506
nrow(subset(df_concoct,(df_concoct$CheckM_contamination < 10 & df_concoct$CheckM_completeness > 0) | (df_concoct$EukCC_contamination < 10 & df_concoct$EukCC_completeness > 0))) #757
nrow(subset(df_metabat,(df_metabat$CheckM_contamination < 10 & df_metabat$CheckM_completeness > 0) | (df_metabat$EukCC_contamination < 10 & df_metabat$EukCC_completeness > 0)))#533
nrow(subset(df_concoct,(df_concoct$CheckM_contamination > 25 | df_concoct$EukCC_contamination > 25))) #128
nrow(subset(df_metabat,(df_metabat$CheckM_contamination > 25 | df_metabat$EukCC_contamination > 25))) #29
nrow(subset(df_concoct,(df_concoct$CheckM_contamination <= 0 & df_concoct$CheckM_completeness <= 0 & df_concoct$EukCC_contamination <= 0 & df_concoct$EukCC_completeness <= 0))) #2360
nrow(subset(df_metabat,(df_metabat$CheckM_contamination <= 0 & df_metabat$CheckM_completeness <= 0 & df_metabat$EukCC_contamination <= 0 & df_metabat$EukCC_completeness <= 0))) #231

#Number of bins of high quality for each software
nrow(subset(df_concoct,(df_concoct$CheckM_contamination < 5 & df_concoct$CheckM_completeness > 90) | (df_concoct$EukCC_contamination < 5 & df_concoct$EukCC_completeness > 90))) #143
nrow(subset(df_metabat,(df_metabat$CheckM_contamination < 5 & df_metabat$CheckM_completeness > 90) | (df_metabat$EukCC_contamination < 5 & df_metabat$EukCC_completeness > 90))) #138
#Number of bins of medium quality for each software
nrow(subset(df_concoct,(df_concoct$CheckM_contamination < 10 & df_concoct$CheckM_completeness > 70) | (df_concoct$EukCC_contamination < 10 & df_concoct$EukCC_completeness > 70))) #263
nrow(subset(df_metabat,(df_metabat$CheckM_contamination < 10 & df_metabat$CheckM_completeness > 70) | (df_metabat$EukCC_contamination < 10 & df_metabat$EukCC_completeness > 70))) #239
#Number of bins of low quality for each software
nrow(subset(df_concoct,(df_concoct$CheckM_contamination < 15 & df_concoct$CheckM_completeness > 50) | (df_concoct$EukCC_contamination < 15 & df_concoct$EukCC_completeness > 50))) #344
nrow(subset(df_metabat,(df_metabat$CheckM_contamination < 15 & df_metabat$CheckM_completeness > 50) | (df_metabat$EukCC_contamination < 15 & df_metabat$EukCC_completeness > 50))) #293
#Number of bins highly inconsistent (contamination > 25%)
nrow(subset(df_concoct,(df_concoct$CheckM_contamination > 25 | df_concoct$EukCC_contamination > 25))) #128
nrow(subset(df_metabat,(df_metabat$CheckM_contamination > 25 | df_metabat$EukCC_contamination > 25))) #29
            
# Number of  bins  with taxonomic assignation
sum(df_complete$GTDB_lineage != "NA")
sum(df_concoct$GTDB_lineage != "NA") # 593
sum(df_metabat$GTDB_lineage != "NA") # 455
sum(df_concoct$CheckM_lineage != "NA") #3288 
sum(df_metabat$CheckM_lineage != "NA") #807
sum(df_concoct$EukCC_lineage != "NA") #141
sum(df_metabat$EukCC_lineage != "NA") #147

#Subset medium quality for each software
medium_concoct <- subset(df_concoct,(df_concoct$CheckM_contamination < 10 & df_concoct$CheckM_completeness > 70) | (df_concoct$EukCC_contamination < 10 & df_concoct$EukCC_completeness > 70))
medium_metabat <- subset(df_metabat,(df_metabat$CheckM_contamination < 10 & df_metabat$CheckM_completeness > 70) | (df_metabat$EukCC_contamination < 10 & df_metabat$EukCC_completeness > 70))

mean(medium_concoct$N50) #43634.17
mean(medium_metabat$N50) #41314.93
mean(medium_concoct$Coverage)  #38.40615
mean(medium_metabat$Coverage) #41.1791
mean(medium_concoct$Number_of_contigs) #836.2357
mean(medium_metabat$Number_of_contigs) #606.4686
mean(medium_concoct$Total_length) #10070533
mean(medium_metabat$Total_length) #8592542

nrow(subset(medium_concoct,(medium_concoct$CheckM_lineage == "p__Cyanobacteria"))) #54
nrow(subset(medium_metabat,(medium_metabat$CheckM_lineage == "p__Cyanobacteria"))) #51
#nrow(subset(medium_concoct,(medium_concoct$EukCC_lineage == "Eukaryota_Fungi" | medium_concoct$EukCC_lineage == "Eukaryota_Fungi_Ascomycota"))) #48
#nrow(subset(medium_metabat,(medium_metabat$EukCC_lineage == "Eukaryota_Fungi" | medium_metabat$EukCC_lineage == "Eukaryota_Fungi_Ascomycota"))) #36

cyano_concoct_medium <- subset(medium_concoct,(medium_concoct$CheckM_lineage == "p__Cyanobacteria"))
cyano_metabat_medium <- subset(medium_metabat,(medium_metabat$CheckM_lineage == "p__Cyanobacteria"))
#fungi_concoct_medium <- subset(medium_concoct,(medium_concoct$EukCC_lineage == "Eukaryota_Fungi" | medium_concoct$EukCC_lineage == "Eukaryota_Fungi_Ascomycota"))
#fungi_metabat_medium <- subset(medium_metabat,(medium_metabat$EukCC_lineage == "Eukaryota_Fungi" | medium_metabat$EukCC_lineage == "Eukaryota_Fungi_Ascomycota"))
#euk <- subset(df_complete,(df_complete$EukCC_lineage != "NA") & (df_complete$EukCC_completeness > 50) & (df_complete$EukCC_contamination < 20))
#unique(euk$EukCC_lineage)

mean(cyano_concoct_medium$Total_length) #8834335
mean(cyano_metabat_medium$Total_length) #7512401
mean(cyano_concoct_medium$Number_of_contigs) #415.2593
mean(cyano_metabat_medium$Number_of_contigs) #256.3725
mean(cyano_concoct_medium$CheckM_completeness) #98.17815
mean(cyano_metabat_medium$CheckM_completeness) #97.36059
mean(cyano_concoct_medium$CheckM_contamination) #1.01463
mean(cyano_metabat_medium$CheckM_contamination) #0.6809804



############################################ Plots statistics #####################################################

normalize <- function (r) {(r-min(r))/(max(r)-min(r))}

#Completion
low_concoct <- subset(df_concoct,(df_concoct$CheckM_completeness > 50) & (df_concoct$CheckM_contamination < 20))
d1 <- sort(low_concoct$CheckM_completeness,decreasing = TRUE)
d1 <- d1[d1 >0]
r1 <- sort(rank(d1))
names1 = rep("CONCOCT",length(d1))
low_metabat <- subset(df_metabat,(df_metabat$CheckM_completeness > 50) & (df_metabat$CheckM_contamination < 20))
d2 <- sort(low_metabat$CheckM_completeness, decreasing = TRUE)
d2 <- d2[d2>0]
r2 <- sort(rank(d2))
names2 = rep("metaBAT2",length(d2))

#Contamination
d3 <- sort(low_concoct$CheckM_contamination)
d3 <- d3[d3 >0]
r3 <- (rank(d3))
names3 = rep("CONCOCT",length(d3))
d4 <- sort(low_metabat$CheckM_contamination)
d4 <- d4[d4 >0]
r4 <- (rank(d4))
names4 = rep("metaBAT2",length(d4))

#F1-scores for fungal bins 
low_euk_concoct <- subset(df_concoct,(df_concoct$EukCC_lineage != "NA") & (df_concoct$EukCC_completeness > 50) & (df_concoct$EukCC_contamination < 20))
low_euk_metabat <- subset(df_metabat,(df_metabat$EukCC_lineage != "NA") & (df_metabat$EukCC_completeness > 50) & (df_metabat$EukCC_contamination < 20))
fungi_concoct_low <- subset(low_euk_concoct,(low_euk_concoct$EukCC_lineage != "Eukaryota" & low_euk_concoct$EukCC_lineage != "Eukaryota_Evosea_Eumycetozoa")) 
fungi_metabat_low <- subset(low_euk_metabat,(low_euk_metabat$EukCC_lineage != "Eukaryota" & low_euk_metabat$EukCC_lineage != "Eukaryota_Evosea_Eumycetozoa"))
  
  # Function to calculate F1-score
  score_func <- function(contam,compl)
  {
    precision = 100 - contam
    score = (2*((precision*compl)/(precision+compl)))/100
    return(score)
  }

conc_scores <- vector()
for (r in 1:nrow(fungi_concoct_low)){
  conc_scores <- c(conc_scores, score_func(fungi_concoct_low[r,"EukCC_contamination"],fungi_concoct_low[r,"EukCC_completeness"]))
}
d5 <- sort(conc_scores,decreasing = TRUE)
d5 <- d5[d5 >0]
r5 <- sort(rank(d5))
names5 = rep("CONCOCT",length(d5))

met_scores <- vector()
for (r in 1:nrow(fungi_metabat_low)){
  met_scores <- c(met_scores, score_func(fungi_metabat_low[r,"EukCC_contamination"],fungi_metabat_low[r,"EukCC_completeness"]))
}
d6 <- sort(met_scores,decreasing = TRUE)
d6 <- d6[d6 >0]
r6 <- sort(rank(d6))
names6 = rep("metaBAT2", length(d6))

# Plot Distributions

# Completion Distribution
data = c(d1,d2) 
rank = c(r1,r2) 
names = c(names1,names2)
df = data.frame(data,rank,names)
ga1 <- ggplot(df, aes(x=rank, y=data, group=names)) +
  geom_line(aes(color=names)) + 
  geom_point(shape=21, size=2.25, fill="white", aes(color=names)) +
  theme_bw() + scale_x_continuous() +
  labs(x="Descending completion rank",y="Estimated bin completeness (%)",color=NULL) +
  theme(panel.grid.minor = element_line(colour="gray95",size=0.01),
        legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95),
        legend.background = element_rect(colour="black"),
        axis.ticks = element_blank(), axis.text.x = element_blank())

# Contamination Distribution
data = c(d3,d4) 
rank = c(r3,r4) 
names = c(names3,names4)
df = data.frame(data,rank,names)
ga1 <- ggplot(df, aes(x=rank, y=data, group=names)) +
  geom_line(aes(color=names)) + 
  geom_point(shape=21, size=2.25, fill="white", aes(color=names)) +
  theme_bw() + scale_x_continuous() +
  labs(x="Ascending contamination rank",y="Estimated bin contamination (%)",color=NULL) +
  theme(panel.grid.minor = element_line(colour="gray95",size=0.01),
        legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95),
        legend.background = element_rect(colour="black"),
        axis.ticks = element_blank(), axis.text.x = element_blank())

# F1-score Distribution
data = c(d5,d6) 
rank = c(r5,r6) 
names = c(names5,names6) 
df = data.frame(data,rank,names)
ga3 <- ggplot(df, aes(x=rank, y=data, group=names)) +
      geom_line(aes(color=names)) + 
      geom_point(shape=21, size=2.25, fill="white", aes(color=names)) +
      theme_bw() + scale_x_continuous() +
      labs(x="Descending score rank",y="F1-score",color=NULL) +
      theme(panel.grid.minor = element_line(colour="gray95",size=0.01),
        legend.justification = c(0.95, 0.95), legend.position = c(0.95, 0.95),
        legend.background = element_rect(colour="black"),
        axis.ticks = element_blank(), axis.text.x = element_blank())


ga1 + scale_color_tableau() + ga2 + scale_color_tableau() + ga3 + scale_color_tableau()
