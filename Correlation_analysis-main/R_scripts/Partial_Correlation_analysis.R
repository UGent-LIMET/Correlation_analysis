# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve, Eewon Tai
# Maintainer: <limet@ugent.be>
# Script: Part II: partial correlation analysis



##########R Pipeline - Part II: partial correlation analysis##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part II: partial correlation analysis - start!")
options(warn=1)
# Part II: partial correlation analysis

#Source of variableMetadata file 1: #NEVER from pipeline, always need to be centrered, normalized...
#INPUT_VARIABLES1 <- VARIABLEMETADATA_EXTERN 

## data_loading
setwd(path_data_in)

##load library
library("reshape")


### GET CONFOUNDERS FROM USER!!!####
# Ask the user for input
raw_input <- readline(prompt = "Enter confounders separated by commas: ")
# Split the string by commas and remove extra whitespace
confounders <- trimws(unlist(strsplit(raw_input, ",")))
# View the result
print(confounders)


#load 2x VMs == SM_merged (format not as strict so, other load)
scalednormalizedvariableMetadata1 <- load_variableMetadata_other(INPUT_VARIABLES1, COLLUMN_NR_START_SAMPLES1) 
scalednormalizedvariableMetadata2 <- load_variableMetadata_other(INPUT_VARIABLES2, COLLUMN_NR_START_SAMPLES2) 
#note: automatically adds X to samples at this point

####
# Verify whether sample names match in both inputs, if not warn the user
#  and remove samples that occur in only one of both inputs
# Then order the inputs according to sample names so that partial correlations make sense
####

if(CHECK_NAMES){
  # check for and remove samples that occur in only 1 input 
  checknames1 <- colnames(scalednormalizedvariableMetadata1[COLLUMN_NR_START_SAMPLES1:ncol(scalednormalizedvariableMetadata1)])
  checknames2 <- colnames(scalednormalizedvariableMetadata2[COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)])
  names1found <- checknames1 %in% checknames2
  names2found <- checknames2 %in% checknames1
  if(sum(names1found == FALSE) > 0){
    warning(paste0("WARNING: IDs ", paste(checknames1[!names1found], collapse = " - "), " listed in INPUT 1 but not found in INPUT 2. Non-matching data is removed."))  
    scalednormalizedvariableMetadata1 <- scalednormalizedvariableMetadata1[,c(1:(COLLUMN_NR_START_SAMPLES1-1), which(names1found)+COLLUMN_NR_START_SAMPLES1-1)]
  }
  if(sum(names2found == FALSE) > 0){
    warning(paste0("WARNING: IDs ", paste(checknames2[!names2found], collapse = " - "), " listed in INPUT 2 but not found in INPUT 1. Non-matching data is removed."))
    scalednormalizedvariableMetadata2 <- scalednormalizedvariableMetadata2[,c(1:(COLLUMN_NR_START_SAMPLES2-1), which(names2found)+COLLUMN_NR_START_SAMPLES2-1)]
  }  
  if(NAME_REORDER){
    warning("Inputs will be matched on sample name. Make sure names correspond in both inputs.")
    scalednormalizedvariableMetadata1 <- scalednormalizedvariableMetadata1[,c(1:(COLLUMN_NR_START_SAMPLES1-1), COLLUMN_NR_START_SAMPLES1 - 1 + order(colnames(scalednormalizedvariableMetadata1[COLLUMN_NR_START_SAMPLES1:ncol(scalednormalizedvariableMetadata1)])))]
    scalednormalizedvariableMetadata2 <- scalednormalizedvariableMetadata2[,c(1:(COLLUMN_NR_START_SAMPLES2-1), COLLUMN_NR_START_SAMPLES2 - 1 + order(colnames(scalednormalizedvariableMetadata2[COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)])))]
  }
} else {
  if(NAME_REORDER){
    # if reordering requested, but check names not performed, throw an error (not usefl)
    stop("NAME_REORDER set to TRUE, but CHECK_NAMES to FALSE, considered inappropraite, please adapt your Configuration.R")
  }
}


#optional: in case of microbiome data as VM2, perform normalisation and centration automatically during run
if(OTU_SCALING_OPTION == OTU_SCALING){
  
  #check only numeric correct (eg ATCG seq at end)
  if(is.na(sum(scalednormalizedvariableMetadata2[, COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)])) ){
    stop("ERROR: Part II: partial correlation analysis stopped because OTU file contains non-numeric data after COLLUMN_NR_START_SAMPLES2")
  }
  
  #normalize according to Pablo script/literat
  #1. compositional (sums to 1 per sample=col)
  scalednormalizedvariableMetadata2[, COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)] <- apply(scalednormalizedvariableMetadata2[, COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)], 2, function(x) x/sum(x))
  
  #2. CLR transformation + autoscale 
  suppressMessages(library(mixOmics))
  scalednormalizedvariableMetadata2[, COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)] <-logratio.transfo(scalednormalizedvariableMetadata2[, COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)], logratio = "CLR", offset = 1*10^-6)
}


#optional: load VM before all calcs to subset using partial correlation threshold 
if(VARIABLES_SELECTION == SUBSET_VM){
  VM <- load_variableMetadata(INPUT_VARIABLES)
}

#2nd optional: decrease VM1 input using subsetted VM. 
#eg. to run 2nd time with only VM1 compids that remain in SUBSET_VM (with p-values for those if big set)
if(VARIABLES_SELECTION2 == SUBSET_VM1_AFTER_VARIABLES_SELECTION){ 
  VM_selection <- load_variableMetadata_other(INPUT_VARIABLES4, COLLUMN_NR_START_SAMPLES4)
  scalednormalizedvariableMetadata1 <- scalednormalizedvariableMetadata1[scalednormalizedvariableMetadata1$CompID %in% VM_selection$CompID,] 
}  




#transpose + re-add variables names
scalednormalizedMetadata_samples1 <- as.data.frame(t(scalednormalizedvariableMetadata1[,COLLUMN_NR_START_SAMPLES1:ncol(scalednormalizedvariableMetadata1)]))
colnames(scalednormalizedMetadata_samples1) <- scalednormalizedvariableMetadata1$CompID
scalednormalizedMetadata_samples2 <- as.data.frame(t(scalednormalizedvariableMetadata2[,COLLUMN_NR_START_SAMPLES2:ncol(scalednormalizedvariableMetadata2)]))
colnames(scalednormalizedMetadata_samples2) <- scalednormalizedvariableMetadata2$CompID


#keep non-unique sampleNames (i.e., sampleNames that occur in both sets)
keep_samplenames <- intersect(rownames(scalednormalizedMetadata_samples1), rownames(scalednormalizedMetadata_samples2)) #keep duplicates only
scalednormalizedMetadata_samples1 <- scalednormalizedMetadata_samples1[rownames(scalednormalizedMetadata_samples1) %in% keep_samplenames,]
scalednormalizedMetadata_samples2 <- scalednormalizedMetadata_samples2[rownames(scalednormalizedMetadata_samples2) %in% keep_samplenames,]


# cols with value O everywhere before log/pareto are now NA. these were previously replaced by 0
# however, NAs can also occur in metadata input and should not be replaced by 0, as such replacement
# can bias partial correlations. All NA columns will result in NA partial correlations, so their removal is not needed
scaledmatrix_samples1 <- scalednormalizedMetadata_samples1 
#scaledmatrix_samples1 <- data.frame(sapply(scaledmatrix_samples1, function(x) ifelse(is.na(x), 0, x)))
#rownames(scaledmatrix_samples1) <- rownames(scalednormalizedMetadata_samples1)

scaledmatrix_samples2 <- scalednormalizedMetadata_samples2 
#scaledmatrix_samples2 <- data.frame(sapply(scaledmatrix_samples2, function(x) ifelse(is.na(x), 0, x)))
#rownames(scaledmatrix_samples2) <- rownames(scalednormalizedMetadata_samples2)


#check2: unique variablenames ok?
colnames_both <- c(colnames(scaledmatrix_samples1), colnames(scaledmatrix_samples2)) #list each unique combi (outer merge)
#if ((ncol(scaledmatrix_samples1) != length(unique(colnames(scaledmatrix_samples1)))) | (ncol(scaledmatrix_samples2) != length(unique(colnames(scaledmatrix_samples2))))){
#  stop("Error - Partial Correlation analysis stopped because non-unique variablenames in INPUT_VARIABLES")
#} #R auto rename "X100.1"
if (length(colnames_both) != length(unique(colnames_both))){
  stop("Error - Partial Correlation analysis stopped because non-unique variablenames between both INPUT_VARIABLES")
}


## set directory to output
setwd(path_data_out_2)


#check ok (normal distribution) per VM
test <- unlist(c(scaledmatrix_samples1)) #just to see over all samples/otu's which distribution + up/down borders
p1 <- hist(test, main ="")
minxaxis <- min(test, na.rm = TRUE)-1
maxaxis <- max(test, na.rm = TRUE)+1
test <- unlist(c(scaledmatrix_samples2)) #just to see over all samples/otu's which distribution + up/down borders
p2 <- hist(test, main ="")

#plot indiv
png(paste(name_project,'_pcorr_histograms_Vars_set1.png', sep=""), width=7, height=5, units="in", res=150)
plot(p1, col=rgb(1,0,0,1/4), main= "Histogram per set1 variables")
dev.off()
png(paste(name_project,'_pcorr_histograms_Vars_set2.png', sep=""), width=7, height=5, units="in", res=150)
plot(p2, col=rgb(0,0,1,1/4), main= "Histogram per set2 variables")
dev.off()

#both on same plot
if(min(test, na.rm = TRUE) < minxaxis | max(test, na.rm = TRUE) > maxaxis){
  #very diff histogram, 1st plot biggest range/freq
  png(paste(name_project,'_pcorr_histograms_Vars_set1+2_overlay.png', sep=""), width=7, height=5, units="in", res=150)
  plot( p2, col=rgb(0,0,1,1/4), main= "Histogram per variables set1+2", xlim= c(minxaxis, maxaxis))  # first histogram, with range x-axis from 2nd hist to be sure both fit
  plot( p1, col=rgb(1,0,0,1/4), add=T)  # second
  dev.off()
}else{
  png(paste(name_project,'_pcorr_histograms_Vars_set1+2_overlay.png', sep=""), width=7, height=5, units="in", res=150)
  plot( p1, col=rgb(1,0,0,1/4), main= "Histogram per variables set1+2", xlim= c(minxaxis, maxaxis))  # first histogram, with range x-axis from 2nd hist to be sure both fit
  plot( p2, col=rgb(0,0,1,1/4), add=T)  # second
  dev.off()
}




## partial correlation between samples
suppressMessages(library(mixOmics)) #not in R3.6.1, R3.5.3
library(RColorBrewer)
colors<-brewer.pal(11,name="RdYlBu")
pal<-colorRampPalette(colors)


#merge both dfs
merged_df <- cbind(scaledmatrix_samples1, scaledmatrix_samples2)
ncol_merged_df <- ncol(merged_df)


if(P_VALUE_CALCULATION == P_VALUES){
  ## for small datasets (after merge), create matrices, heatmaps, matrices p-values 
  #eg clinical vs antropometric
  
  if(ncol_merged_df <= 20){ 
    
    #####HM SvsS per set######
    ## set 1
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(t(scaledmatrix_samples1), confounders)
    write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, "pcorr_coef_SvsS_set1.txt") 
    
    #p values after FDR adjusted pcorr matrix calc
    Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(t(scaledmatrix_samples1), confounders)
    write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, "pvaluesFDR_pcorr_coef_SvsS_set1.txt")
    
    
    #write all info as list
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_SvsS_set1.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    #mean(partial_correlation_scaledmatrix)
    
    name_plot <- paste(name_project, "_Heatmap_pcorr_SvsS_variableset1", sep="")
    heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                        row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    
    #clustering samples order
    ord <- hclust(dist(partial_correlation_scaledmatrix))$order #info about order in object ord$order
    #default clustering: Cluster method: complete; Distance: euclidean see doc: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
    #ord
    dendrogram <-as.dendrogram(hclust(dist(partial_correlation_scaledmatrix)))
    png(paste(name_project, "_dendogram_clust_pcorr_SvsS_set1.png", sep=""), width=10, height=5, units="in", res=300)
    par(cex=.2)
    plot(dendrogram)
    dev.off()
    
    #info from dendogram, write order to pcorr_info
    comp_names2 <- NULL
    comp_names2$SampleName <- rownames(partial_correlation_scaledmatrix) 
    comp_names2 <- data.frame(comp_names2)
    comp_names2$ord <- ord
    write_dataframe_as_txt_file(comp_names2, paste(name_project, "_dendogram_clust_pcorr_SvsS_set1.txt", sep=""))
    
    
    ## set 2
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(t(scaledmatrix_samples2), confounders)
    write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, "pcorr_coef_SvsS_set2.txt") 
    
    #p values after FDR adjusted corr matrix calc
    Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(t(scaledmatrix_samples2), confounders)
    write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, "pvaluesFDR_pcorr_coef_SvsS_set2.txt")
    
    
    #write all info as list
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_SvsS_set2.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    name_plot <- paste(name_project, "_Heatmap_pcorr_SvsS_variableset2", sep="")
    heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                        row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    gc()
    ###########
    
    
    #####HM SvsVar per set######
    #set 1
    partial_correlation_scaledmatrix <- t(scaledmatrix_samples1) #no partial_correlation, just heatmap S vs V (like in statist, univar analysis scripts)
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    if(ncol(partial_correlation_scaledmatrix) <= 200){ #craches with big dataset
      name_plot <- paste(name_project, "_Heatmap_SvsVar_variableset1", sep="")
      heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                          row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    }
    
    #set 2
    partial_correlation_scaledmatrix <- t(scaledmatrix_samples2) #no partial_correlation, just heatmap S vs V
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    if(ncol(partial_correlation_scaledmatrix) <= 200){ #craches with big dataset
      name_plot <- paste(name_project, "_Heatmap_SvsVar_variableset2", sep="")
      heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                          row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    }
    gc()
    ###########
    
    
    ###########HM SvsS merged###########
    
    #write merged df
    write_matrixTT_as_txt_file(merged_df, "merged_partial_correlationMetadata_sets.txt") 
    
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(t(merged_df), confounders)
    write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, "pcorr_coef_SvsS_set1+2.txt") 
    
    #p values after FDR adjusted corr matrix calc
    Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(t(merged_df), confounders)
    write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, "pvaluesFDR_pcorr_coef_SvsS_set1+2.txt")
    
    #write all info as list
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_SvsS_set1+2.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    name_plot <- paste(name_project, "_Heatmap_pcorr_SvsS_variableset1+2", sep="")
    heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                        row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    
    
    #clustering to check
    test <- NULL
    test <- heatmap_comp$row.names
    test <- as.data.frame(test)
    test$order_hm <- c(1:nrow(test))
    gc()
    #########
    
    
    #####HM VarvsVar per set - intra ######
    ## set 1
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(scaledmatrix_samples1, confounders)
    
    #p values after FDR adjusted corr matrix calc
    try(Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(scaledmatrix_samples1, confounders))
    
    #write in parts because too big!
    if(ncol(partial_correlation_scaledmatrix) >= 10000){ #craches with big dataset
      iter = round(ncol(partial_correlation_scaledmatrix)/10000,0)+1
      iter = c(1:iter)
      start = 1
      for(i in iter){
        #start = start + as.numeric(i)    #from 1
        stop = start + 10000 -1            #to 9999
        if (stop > ncol(partial_correlation_scaledmatrix)){
          stop = ncol(partial_correlation_scaledmatrix)
        }
        #print (start)
        #print(stop)
        
        write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_VarvsVar_set1_", as.character(start), "-",as.character(stop), ".txt")) 
        try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix[,start:stop], paste0("pvaluesFDR_pcorr_coef_VarvsVar_set1_", as.character(start), "-",as.character(stop), ".txt")))
        
        start = stop +1
      }
    }else{
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_VarvsVar_set1.txt")) 
      try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, paste0("pvaluesFDR_pcorr_coef_VarvsVar_set1.txt")))
    }
    
    #write all info as list
    #save <- partial_correlation_scaledmatrix
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_VarvsVar_set1.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    if(ncol(partial_correlation_scaledmatrix) <= 200){ #craches with big dataset
      name_plot <- paste(name_project, "_Heatmap_pcorr_VarvsVar_variableset1", sep="")
      heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                          row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    }
    
    #clustering variables order
    ord <- hclust(dist(partial_correlation_scaledmatrix))$order #info about order in object ord$order
    #default clustering: Cluster method: complete; Distance: euclidean see doc: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
    #ord
    dendrogram <-as.dendrogram(hclust(dist(partial_correlation_scaledmatrix)))
    png(paste(name_project, "_dendogram_clust_pcorr_VarvsVar_variableset1.png", sep=""), width=10, height=5, units="in", res=300)
    par(cex=.2)
    plot(dendrogram)
    dev.off()
    
    #info from dendogram, write order to pcorr_info
    comp_names2 <- NULL
    comp_names2$CompID <- rownames(partial_correlation_scaledmatrix) 
    comp_names2 <- data.frame(comp_names2)
    comp_names2$ord <- ord
    write_dataframe_as_txt_file(comp_names2, paste(name_project, "_dendogram_clust_pcorr_VarvsVar_variableset1.txt", sep=""))
    
    
    ## set 2
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(scaledmatrix_samples2, confounders)
    
    
    #p values after FDR adjusted corr matrix calc
    try(Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(scaledmatrix_samples2, confounders))
    
    #write in parts because too big!
    if(ncol(partial_correlation_scaledmatrix) >= 10000){ #craches with big dataset
      iter = round(ncol(partial_correlation_scaledmatrix)/10000,0)+1
      iter = c(1:iter)
      start = 1
      for(i in iter){
        #start = start + as.numeric(i)    #from 1
        stop = start + 10000 -1            #to 9999
        if (stop > ncol(partial_correlation_scaledmatrix)){
          stop = ncol(partial_correlation_scaledmatrix)
        }
        #print (start)
        #print(stop)
        
        write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_VarvsVar_set2_", as.character(start), "-",as.character(stop), ".txt")) 
        try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix[,start:stop], paste0("pvaluesFDR_pcorr_coef_VarvsVar_set2_", as.character(start), "-",as.character(stop), ".txt")))
        
        start = stop +1
      }
    }else{
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_VarvsVar_set2.txt")) 
      try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, paste0("pvaluesFDR_pcorr_coef_VarvsVar_set2.txt")))
    }
    
    #write all info as list
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_VarvsVar_set2.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    if(ncol(partial_correlation_scaledmatrix) <= 200){ #craches with big dataset
      name_plot <- paste(name_project, "_Heatmap_pcorr_VarvsVar_variableset2", sep="")
      heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                          row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    }
    gc()
    ###########
    
    
    ###########HM VarvsVar merged - intra and inter ###########
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(merged_df, confounders)
    
    #p values after FDR adjusted corr matrix calc
    try(Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(merged_df, confounders))
    
    #write in parts because too big!
    if(ncol(partial_correlation_scaledmatrix) >= 10000){ #craches with big dataset
      iter = round(ncol(partial_correlation_scaledmatrix)/10000,0)+1
      iter = c(1:iter)
      start = 1
      for(i in iter){
        #start = start + as.numeric(i)    #from 1
        stop = start + 10000 -1            #to 9999
        if (stop > ncol(partial_correlation_scaledmatrix)){
          stop = ncol(partial_correlation_scaledmatrix)
        }
        #print (start)
        #print(stop)
        
        write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_VarvsVar_set1+2_", as.character(start), "-",as.character(stop), ".txt")) 
        try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix[,start:stop], paste0("pvaluesFDR_pcorr_coef_VarvsVar_set1+2_", as.character(start), "-",as.character(stop), ".txt")))
        
        start = stop +1
      }
    }else{
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_VarvsVar_set1+2.txt")) 
      try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, paste0("pvaluesFDR_pcorr_coef_VarvsVar_set1+2.txt")))
    }
    
    #write all info as list
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_VarvsVar_set1+2.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    if(ncol(partial_correlation_scaledmatrix) <= 200){ #craches with big dataset
      name_plot <- paste(name_project, "_Heatmap_pcorr_VarvsVar_variableset1+2", sep="")
      heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                          row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    }
    gc()
    #########
    
    
    
    #last one for info_corr ranked and boxplot !!
    ###########HM Var set1 vs Var set2 - inter ###########
    partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(scaledmatrix_samples1, scaledmatrix_samples2, confounders)
    
    #p values after FDR adjusted corr matrix calc
    try(Pvalues_partial_correlation_scaledmatrix <- Pvalues_partial_correlation_matrixFDR(scaledmatrix_samples1, scaledmatrix_samples2, confounders))
    
    #write in parts because too big!
    if(ncol(partial_correlation_scaledmatrix) >= 10000){ #craches with big dataset
      iter = round(ncol(partial_correlation_scaledmatrix)/10000,0)+1
      iter = c(1:iter)
      start = 1
      for(i in iter){
        #start = start + as.numeric(i)    #from 1
        stop = start + 10000 -1            #to 9999
        if (stop > ncol(partial_correlation_scaledmatrix)){
          stop = ncol(partial_correlation_scaledmatrix)
        }
        #print (start)
        #print(stop)
        
        write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt")) 
        try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix[,start:stop], paste0("pvaluesFDR_pcorr_coef_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt")))
        
        start = stop +1
      }
    }else{
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_Varset1vsVarset2.txt")) 
      try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, paste0("pvaluesFDR_pcorr_coef_Varset1vsVarset2.txt")))
    }
    
    #write all info as list
    colnames(partial_correlation_scaledmatrix) <- colnames(scaledmatrix_samples2)
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "rho", "adj_p_value")
    try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_Varset1vsVarset2.txt"))) #still with self-corr and 2directions
    
    #1feat all same value => NAN after corr so add rn nan here as well
    partial_correlation_scaledmatrix[is.nan(partial_correlation_scaledmatrix)] = 0
    partial_correlation_scaledmatrix[is.na(partial_correlation_scaledmatrix)] = 0
    
    if(ncol(partial_correlation_scaledmatrix) <= 200){ #craches with big dataset
      name_plot <- paste(name_project, "_Heatmap_pcorr_Varset1vsVarset2", sep="")
      heatmap_comp <- cim(data.matrix(partial_correlation_scaledmatrix), color=pal(100), save='png', name.save = name_plot,
                          row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
    }
    gc()
    #########
    
    
  }else if(ncol_merged_df > 20){
    
    
    
    
    
    
    
    
    ## for large datasets: no visual, only table with corr coef, p-value per var-var comparison
    #eg metab vs antropometric
    #approx 1/2 h (approx 4000)
    
    ###########HM VarvsVar inter ###########
    
    if(CODE_RUN_MODE == CODE_AUTORUN){
      rm(scalednormalizedMetadata_samples1)
      rm(scalednormalizedMetadata_samples2)
      #rm(scaledmatrix_samples1)
      #rm(scaledmatrix_samples2)
      rm(merged_df)
      gc()
    }
    
    #calculates partial_correlation coeffic (Spearman, FDR) 
    # suppress warning for ties
    options(warn = -1)
    #pcorr_info <- partial_correlation_analysis_loop(t(scaledmatrix_samples1), t(scaledmatrix_samples2), confounders) 
    pcorr_info <- partial_correlation_analysis_loop(scaledmatrix_samples1, scaledmatrix_samples2, confounders) 
    
    # reactivate warnings, print as they occur
    options(warn = 1)
    
    # add compound metadata if requested
    if(COMPOUND_TO_NAME2){
      CompID2 <- scalednormalizedvariableMetadata2$CompID
      newCompID2 <- scalednormalizedvariableMetadata2$MetabolitesName
      pcorr_info$CompID2 <- newCompID2[match(pcorr_info$CompID2, CompID2)]
    }
    
    if(COMPOUND_TO_METADATA1){
      compound_metadata_translation <- read.table(paste0(path_data_in, "/", "Compounds_Filter.txt"), header = TRUE, sep = "\t", dec = ",")
      compound_metadata_translation <- compound_metadata_translation[,which(colnames(compound_metadata_translation) %in% c("CompID", "Calc..MW", "m.z", "RT..min.", "Reference.Ion"))]
      CompID1 <- compound_metadata_translation$CompID
      compound_metadata_translation_new <- compound_metadata_translation[match(pcorr_info$CompID1, CompID1),]
      pcorr_info <- cbind(pcorr_info, compound_metadata_translation_new[, -1])
    }
    
    #write in parts because too big!
    if(ncol(pcorr_info) >= 1000000){ #excel max approx 1mio rows, save in parts
      iter = round(ncol(pcorr_info)/1000000,0)+1
      iter = c(1:iter)
      start = 1
      for(i in iter){
        #start = start + as.numeric(i)    #from 1
        stop = start + 1000000 -1            #to 999999
        if (stop > ncol(pcorr_info)){
          stop = ncol(pcorr_info)
        }
        #print (start)
        #print(stop)
        
        #write all info as list
        write_dataframe_as_txt_file(pcorr_info[,start:stop], paste0("partial_correlation_info_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt"))
        
        start = stop +1
      }
    }else{
      #write all info as list
      write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_Varset1vsVarset2.txt"))
    }
    
    #########
    
  }else{  #NEVER matches; works well but on mulicore sometimes core crashes and you loose ALL...
    
    
    
    
    
    
    
    
    ## for ultra large datasets: no visual, only table with corr coef, p-value per var-var comparison
    #eg metab vs microbiome
    #result cannot be opened in excel, so split
    
    ###########HM VarvsVar inter ###########
    
    library(parallel)
    
    if(CODE_RUN_MODE == CODE_AUTORUN){
      rm(scalednormalizedMetadata_samples1)
      rm(scalednormalizedMetadata_samples2)
      #rm(scaledmatrix_samples1)
      #rm(scaledmatrix_samples2)
      rm(merged_df)
      gc()
      nr_of_cores <- detectCores()-1 #5 
      print(nr_of_cores)
      #detectCores() #in Linux, take nr of rbox ef set 5/6-1 cpu avail (ok for shared pcs). take less than you have
      #https://stackoverflow.com/questions/15668893/r-multicore-mcfork-unable-to-fork-cannot-allocate-memory#32171864
      #not use in GUI: If you want to take advantage of parallel processing, it is advised not to use a GUI, as that interrupts the multithreaded processes between your code and the GUI programme
    }
    if(CODE_RUN_MODE == CODE_DEVELOPMENT){
      nr_of_cores <- 1 #Must be exactly 1 on Windows (which uses the master process)
    }
    
    #calculates partial_correlation coeffic (Spearman, FDR) using parallelization
    #https://www.inwt-statistics.com/read-blog/code-performance-in-r-parallelization.html   #mclapply
    #https://stackoverflow.com/questions/3557652/apply-over-two-data-frames                 #how to iter df1 agains df2
    
    #corr matrix calc using lapply (mclapply is parallel in linux)
    bl <- mclapply(scaledmatrix_samples1, mc.cores = nr_of_cores, function(u){
      mclapply(scaledmatrix_samples2, mc.cores = nr_of_cores, function(v){
        partial_correlation_matrix_samples(u,v, confounders) # Function with column from scaledmatrix_samples1 and column from scaledmatrix_samples2 as inputs
      })
    })
    partial_correlation_scaledmatrix <- matrix(unlist(bl), ncol=ncol(scaledmatrix_samples2), byrow=T)
    print("corr done")
    
    #p values after FDR adjusted corr matrix calc using lapply (mclapply is parallel in linux)
    bl <- mclapply(scaledmatrix_samples1, mc.cores = nr_of_cores, function(u){
      mclapply(scaledmatrix_samples2, mc.cores = nr_of_cores, function(v){
        Pvalues_partial_correlation_matrixFDR(u,v, confounders) # Function with column from scaledmatrix_samples1 and column from scaledmatrix_samples2 as inputs
      })
    })
    Pvalues_partial_correlation_scaledmatrix <- matrix(unlist(bl), ncol=ncol(scaledmatrix_samples2), byrow=T)
    print("corr pvalues done")
    
    #write all info as list
    colnames(partial_correlation_scaledmatrix) <- colnames(scaledmatrix_samples2)
    partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
    partial_correlation_df$CompID <- colnames(scaledmatrix_samples1) #depending on which name (system or nic name), as chosen above
    data_melt <- melt(partial_correlation_df, id = c("CompID")) 
    
    Pvalues_partial_correlation_df <- as.data.frame(Pvalues_partial_correlation_scaledmatrix)
    Pvalues_partial_correlation_df$CompID <- rownames(Pvalues_partial_correlation_df) #depending on which name (system or nic name), as chosen above
    data_melt_Pvalues <- melt(Pvalues_partial_correlation_df, id = c("CompID")) 
    
    pcorr_info <- cbind(data_melt, data_melt_Pvalues$value)
    colnames(pcorr_info) <- c("CompID1", "CompID2", "rho", "adj_p_value")
    
    
    
    #write in parts because too big!
    if(ncol(pcorr_info) >= 1000000){ #excel max approx 1mio rows, save in parts
      iter = round(ncol(pcorr_info)/1000000,0)+1
      iter = c(1:iter)
      start = 1
      for(i in iter){
        #start = start + as.numeric(i)    #from 1
        stop = start + 1000000 -1            #to 999999
        if (stop > ncol(pcorr_info)){
          stop = ncol(pcorr_info)
        }
        #print (start)
        #print(stop)
        
        #write all info as list
        write_dataframe_as_txt_file(pcorr_info[,start:stop], paste0("partial_correlation_info_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt"))
        
        start = stop +1
      }
    }else{
      #write all info as list
      write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_Varset1vsVarset2.txt"))
    }
    
    
    #########
    
  }
}



if(P_VALUE_CALCULATION == NO_P_VALUES){
  
  #####HM VarvsVar per set - intra ######
  #set 1
  partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(scaledmatrix_samples1, confounders)
  
  #write in parts because too big!
  if(ncol(partial_correlation_scaledmatrix) >= 10000){ #craches with big dataset
    iter = round(ncol(partial_correlation_scaledmatrix)/10000,0)+1
    iter = c(1:iter)
    start = 1
    for(i in iter){
      #start = start + as.numeric(i)    #from 1
      stop = start + 10000 -1            #to 9999
      if (stop > ncol(partial_correlation_scaledmatrix)){
        stop = ncol(partial_correlation_scaledmatrix)
      }
      #print (start)
      #print(stop)
      
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_VarvsVar_set1_", as.character(start), "-",as.character(stop), ".txt")) 
      
      start = stop +1
    }
  }else{
    write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_VarvsVar_set1.txt")) 
  }
  
  #write all info as list
  #save <- partial_correlation_scaledmatrix
  partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
  partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
  data_melt <- melt(partial_correlation_df, id = c("CompID")) 
  
  pcorr_info <- data_melt
  colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef")
  try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_VarvsVar_set1.txt"))) #still with self-corr and 2directions
  
  
  
  #set 2
  partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(scaledmatrix_samples2, confounders)
  
  #write in parts because too big!
  if(ncol(partial_correlation_scaledmatrix) >= 10000){ #craches with big dataset
    iter = round(ncol(partial_correlation_scaledmatrix)/10000,0)+1
    iter = c(1:iter)
    start = 1
    for(i in iter){
      #start = start + as.numeric(i)    #from 1
      stop = start + 10000 -1            #to 9999
      if (stop > ncol(partial_correlation_scaledmatrix)){
        stop = ncol(partial_correlation_scaledmatrix)
      }
      #print (start)
      #print(stop)
      
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_VarvsVar_set2_", as.character(start), "-",as.character(stop), ".txt")) 
      
      start = stop +1
    }
  }else{
    write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_VarvsVar_set2.txt")) 
  }
  
  #write all info as list
  partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
  partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
  data_melt <- melt(partial_correlation_df, id = c("CompID")) 
  
  pcorr_info <- data_melt
  colnames(pcorr_info) <- c("CompID1", "CompID2", "Corr_coef")
  try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_VarvsVar_set2.txt"))) #still with self-corr and 2directions
  
  gc()
  ###########
  
  
  ###########HM Var set1 vs Var set2 - inter ###########
  #NO P-values, only inter !!
  #approx 10min for 'large' (ncolmerged=5000)
  
  partial_correlation_scaledmatrix <- partial_correlation_matrix_samples(scaledmatrix_samples1, scaledmatrix_samples2, confounders)
  
  #write in parts because too big!
  if(ncol(partial_correlation_scaledmatrix) >= 1000000){ #craches with big dataset
    iter = round(ncol(partial_correlation_scaledmatrix)/1000000,0)+1
    iter = c(1:iter)
    start = 1
    for(i in iter){
      #start = start + as.numeric(i)    #from 1
      stop = start + 1000000 -1            #to 999999
      if (stop > ncol(partial_correlation_scaledmatrix)){
        stop = ncol(partial_correlation_scaledmatrix)
      }
      #print (start)
      #print(stop)
      
      write_matrixTT_as_txt_file(partial_correlation_scaledmatrix[,start:stop], paste0("pcorr_coef_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt")) 
      #try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix[,start:stop], paste0("pvaluesFDR_pcorr_coef_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt")))
      
      start = stop +1
    }
  }else{
    write_matrixTT_as_txt_file(partial_correlation_scaledmatrix, paste0("pcorr_coef_Varset1vsVarset2.txt")) 
    #try(write_matrixTT_as_txt_file(Pvalues_partial_correlation_scaledmatrix, paste0("pvaluesFDR_pcorr_coef_Varset1vsVarset2.txt")))
  }
  
  #write all info as list
  
  #save <- partial_correlation_scaledmatrix
  colnames(partial_correlation_scaledmatrix) <- colnames(scaledmatrix_samples2)
  partial_correlation_df <- as.data.frame(partial_correlation_scaledmatrix)
  partial_correlation_df$CompID <- rownames(partial_correlation_df) #depending on which name (system or nic name), as chosen above
  data_melt <- melt(partial_correlation_df, id = c("CompID")) 
  
  
  pcorr_info <- data_melt #cbind(data_melt, data_melt_Pvalues$value)
  colnames(pcorr_info) <- c("CompID1", "CompID2", "rho")
  #try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_Varset1vsVarset2.txt"))) #still with self-corr and 2directions
  
  
  #write in parts because too big!
  if(ncol(pcorr_info) >= 1000000){ #excel max approx 1mio rows, save in parts
    iter = round(ncol(pcorr_info)/1000000,0)+1
    iter = c(1:iter)
    start = 1
    for(i in iter){
      #start = start + as.numeric(i)    #from 1
      stop = start + 1000000 -1            #to 999999
      if (stop > ncol(pcorr_info)){
        stop = ncol(pcorr_info)
      }
      #print (start)
      #print(stop)
      
      #write all info as list
      write_dataframe_as_txt_file(pcorr_info[,start:stop], paste0("partial_correlation_info_Varset1vsVarset2_", as.character(start), "-",as.character(stop), ".txt"))
      
      start = stop +1
    }
  }else{
    #write all info as list
    write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_Varset1vsVarset2.txt"))
  }
  
}


##rho, p as mumeric!
pcorr_info$rho <- as.numeric(pcorr_info$rho)
if(P_VALUE_CALCULATION == P_VALUES){
  pcorr_info$adj_p_value <- as.numeric(pcorr_info$adj_p_value)
}


##if needed, run part script manual
#pcorr_info <- read.table("partial_correlation_info_Varset1vsVarset2.txt", header=TRUE, sep="\t")    #reserve, if already ran code
#pcorr_info_select <- pcorr_info                                                              #if already selection done in file
#pcorr_info_select <- pcorr_info_top                                                          #if want to keep only TOP_RANK in VM selection


## selection VMs according to selection rho and p-value thresholds
if(VARIABLES_SELECTION == SUBSET_VM){
  
  if(P_VALUE_CALCULATION == P_VALUES){
    pcorr_info_select <- pcorr_info[(abs(pcorr_info$rho) >= PARTIAL_CORRELATION_COEFFICIENT_THRESHOLD) & (pcorr_info$adj_p_value < P_VALUE_THRESHOLD),] 
  }
  if(P_VALUE_CALCULATION == NO_P_VALUES){
    pcorr_info_select <- pcorr_info[abs(pcorr_info$rho) >= PARTIAL_CORRELATION_COEFFICIENT_THRESHOLD,] 
  }
  
  #VM1: only for metabolomics df
  #subset vm only good variables (rows in CompIDs)
  variables_selection <- pcorr_info_select$CompID1 
  if(class(scalednormalizedvariableMetadata1$CompID) == "integer"){
    variables_selection <- as.character(substr(variables_selection, 2, length(variables_selection)))    #rm "X" only if integers, must be to match the VM rules
  }
  variableMetadata_selection <- VM[as.character(VM$CompID) %in% variables_selection,]
  
  #rm "X" in samplenales before write as well
  samplenames <- substring(colnames(variableMetadata_selection[COLLUMN_NR_START_SAMPLES:ncol(variableMetadata_selection)]), 2)
  colnames(variableMetadata_selection) <- c(colnames(variableMetadata_selection[1:(COLLUMN_NR_START_SAMPLES-1)]), samplenames)
  
  #write vm selection
  #write_dataframe_as_txt_file(variables_selection, paste0("pcorr_variables_selection.txt"))
  try({write_dataframe_as_txt_file(variableMetadata_selection, paste0("pcorr_VM_selection.txt"))}) #try if nr rows == 0
}



## ranking of best Rho
#+1 and -1 both BEST, so first make ABS(rho) !!!
#independant if p-values present of not, rank Rho for write ranked file + boxplot all + dotplot top 100

#sort based on order run samples
pcorr_info$ABSrho <- abs(as.numeric(pcorr_info$rho))
pcorr_info <- pcorr_info[order(pcorr_info$ABSrho, decreasing = TRUE),]
try(write_dataframe_as_txt_file(pcorr_info, paste0("partial_correlation_info_Varset1vsVarset2_ABSranked.txt")))
#fyi: ok if too large to open in excel, since other also saved and here ranking is more important


## dotplot top 100 Rho
#top 100, dependent on config retain
TOP_RANK <- round(TOP_RANK, 0)
if(TOP_RANK > nrow(pcorr_info)){
  TOP_RANK <- nrow(pcorr_info)
}
pcorr_info_top <- pcorr_info[1:TOP_RANK,]

dotplot_corr <- plot_dotplot(pcorr_info_top)
png(paste(name_project, "_dotplot_top", as.character(TOP_RANK), "_pcorr_Rho_Varset1vsVarset2.png", sep=""), width=7, height=5, units="in", res=150)
plot(dotplot_corr)
dev.off()



## boxplot all Rho
#boxplot(pcorr_info$ABSrho) #base plot

boxplot_corr <- plot_boxplot(pcorr_info)
png(paste(name_project, "_pcorr_boxplot_Rho_Varset1vsVarset2.png", sep=""), width=7, height=5, units="in", res=150)
try(plot(boxplot_corr))
dev.off()




print("R pipeline - Part II: partial_correlation analysis - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################

