# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Main


##########Configuration##########
#Working directory

## options
PATH_USER <- 'C:/Users/eetai/OneDrive - UGent/Desktop/limet_lab_materials/FFQ_correlations/Correlation_analysis-edit/Correlation_analysis-main/' # example path user
  
CODE_AUTORUN <- 'run code in ternminal automatically'
CODE_DEVELOPMENT <- 'run code manually in Rstudio for development locally'

## Adjustments
PATH <- 'C:/Users/eetai/OneDrive - UGent/Desktop/limet_lab_materials/FFQ_correlations/Correlation_analysis-edit/Correlation_analysis-main/Targeted_FFQ/'

CODE_RUN_MODE <- CODE_AUTORUN

#
#####################


##########Sources##########
#Source functions
path_R_scripts <- file.path(PATH_USER, 'R_scripts')
setwd(path_R_scripts)
source('data_loading.R')
source('data_statistics.R')
source('data_converting.R')
source('data_writing.R')
source('data_plot.R')

if (CODE_RUN_MODE == CODE_AUTORUN){
  #Recognize projectname from Rscript commando
  name_project = "241119_FAME_FFQ_Corr" #commandArgs(trailingOnly=TRUE) #name is given in "Rscript main.r projectname" as argument
  # test if there is one argument: if not, return an error
  if (length(args) == 0) {
    stop("Your projectname must be supplied", call.=FALSE)
  } 
  if (length(args) > 1) {
    stop("Only one projectname must be supplied", call.=FALSE)
  }
}
if (CODE_RUN_MODE == CODE_DEVELOPMENT){
  #give projectname in configuration.r (for code in development @Rstudio)
  name_project <- EXPERIMENT
}

#Source configuration and set input folder
path_data_in <- file.path(PATH, 'Input', name_project) #directory must exist!
setwd(path_data_in)
source('Configuration.R')

#make output directory for project:
dir.create(file.path(PATH, 'Output_20260330', name_project))
path_data_out <- file.path(PATH, 'Output_20260330', name_project)

#
#####################

 
##########R Pipeline - Main##########
print(Sys.time())
start_time_total <- Sys.time()
print("R pipeline - start!")

#Project name:
stopifnot(name_project == EXPERIMENT)
print(name_project)

#Run selected modules
if(RUN_PART_CORRELATION == RUN_CODE){
  setwd(path_R_scripts)
  source('Correlation_analysis.R')
  source('Partial_Correlation_analysis.R')
}
if(RUN_PART_CORRELATION == DONT_RUN_CODE){
  print('no correlation analysis was performed')
}

print(name_project)
print("R pipeline - done!")
print(Sys.time())
end_time_total <- Sys.time()
print(end_time_total - start_time_total)
#
####################
