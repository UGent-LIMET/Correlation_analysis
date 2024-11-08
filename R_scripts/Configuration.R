## @knitr INFO
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Configuration


##########Global_settings##########

## @knitr settings
## options
RUN_CODE <- 'run this part of the pipeline'
DONT_RUN_CODE <- 'skip this part of the pipeline, keeps order: pre-processing - targeted analysis - statistical analysis - annotation'

## Adjustments
#Project name:
EXPERIMENT <- 'experiment_name' #structured and short, see READ_ME
POLARITY <- "positive" #{"positive", "negative"} #needed this format for annotation to work
# file_conversion and pre-processing need to be performed seperate
# only choose "both" when no pre-processing needed (eg. merge both ionisationmodes)
USER_COMMENT <- "Tutorial comment" #Add info about experiment, eg. explain (Multiple)Comparisons, to include in reports

COMPOUND_TO_METADATA1 <- TRUE # whether to add (TRUE) or not (FALSE) metadata based on the compound ID in the output
COMPOUND_METADATA1 <- "Compounds_Filter.txt" # should contain column names "Calc. MW", "m/z", "RT (min)" "Reference Ion"

COMPOUND_TO_NAME2 <- TRUE # whether to convert (TRUE) or not (FALSE) the compound ID to its name instead in the output
                           # this is taken from the INPUT_VARIABLES2 object (see below), and
                           # should contain a column name "MetabolitesName"

CHECK_NAMES <- TRUE # whether to remove (TRUE) or keep (FALSE) samples that occur in input 1 but not in input 2, and vice versa
NAME_REORDER <- TRUE # whether to reorder (TRUE) or keep the order (FALSE) of columns/rows in input 1 and 2 according to sample names
                    #   assumes CHECK_NAMES is TRUE, reordering is not useful if discrepant sample names are present

RUN_PART_CORRELATION <- RUN_CODE

#
#####################



##########Correlation statistics##########
if(RUN_PART_CORRELATION == RUN_CODE){
  
  ## optionsÒ
  #in case of microbiome data as VM2, optional to perform normalisation and centration automatically during run
  OTU_SCALING <- 'Performes centered log ratio transformation of OTU table'
  NO_OTU_SCALING <- 'no transformations performed'
  
  #Adjusted p-value (FDR) calculation, takes longer to run
  NO_P_VALUES <- "skip the p-value calculation"
  P_VALUES <- "perform the p-value calculation"
  
  #subset VM (variables selection) for next analysis
  SUBSET_VM <- 'give in adjustments name of VM that needs to be subsetted with only highly correlated variables'
  NO_SUBSET_VM <- 'no VM needed, no subsetting of VM'

  #reduce VM1 input file using subset VM from previous analysis
  SUBSET_VM1_AFTER_VARIABLES_SELECTION <- 'use VM with only highly correlated variables to subset the input metabolomics VM1 file, eg if p-values takes too long to run on whole dataset'
  NO_SUBSET_VM1 <- 'no subsetting of VM1 using required VM file'
  
  
  ## Adjustments
  # variableMetadata file 1: metabolomics
  #NEVER from pipeline, always need to be centrered, normalized...
  INPUT_VARIABLES1 <- 'VM1_part_metabolomics.txt'  #'name.txt' of file.
  COLLUMN_NR_START_SAMPLES1 <- 2
  
  # variableMetadata file 2: other
  #NEVER from pipeline, always need to be centrered, normalized...
  INPUT_VARIABLES2 <- 'VM2_part_source.txt'  #'name.txt' of file. 
  COLLUMN_NR_START_SAMPLES2 <- 4
  
  #optional: in case of microbiome data as VM2, perform normalisation and centration automatically during run
  OTU_SCALING_OPTION <- NO_OTU_SCALING
  
  #Adjusted p-value calculation:
  P_VALUE_CALCULATION <- P_VALUES
  
  #top correlation coeficents (rho) to retain in dotplot
  TOP_RANK <- 100                 #default 100
  
  #selection variables to retain in VM with the same CompIDs as 1st variableMetadata file
  VARIABLES_SELECTION <- NO_SUBSET_VM
  INPUT_VARIABLES <- 'VariableMetadata_metabolomics.txt'  #'name.txt' of file.
  COLLUMN_NR_START_SAMPLES <- 21
  CORRELATION_COEFFICIENT_THRESHOLD <- 0.5      #default 0.5: will retain x < -0.49 and x > 0.49
  P_VALUE_THRESHOLD <- 0.05                     #default 0.05; only used if p-values are calculated
  
  #reduce VM1 using subsetted VM file
  #eg. to run 2nd time with only VM1 compids that remain in SUBSET_VM (if p-values takes too long to run on whole dataset)
  VARIABLES_SELECTION2 <- NO_SUBSET_VM1
  INPUT_VARIABLES4 <- 'VM_selection.txt'  #'name.txt' of file.
  COLLUMN_NR_START_SAMPLES4 <- 20
  
}
#
#####################



