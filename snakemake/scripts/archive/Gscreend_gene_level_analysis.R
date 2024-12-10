#' ---
#' title: "Gscreend_gene_level"
#' author: "Fer A"
#' date: "2024-04-04"
#' output: html_document
#' ---
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library (tidyverse)
library (readr)
library (dplyr)
library(stringr)
library(gscreend)
library (SummarizedExperiment)
library(tidyr)
library(ggplot2)

#' 
#' 1. Getting the tables with all the data (counts)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pooled_screen_1_counts <- read_csv ("~/Documents/Sanger MPhil/Sym_101Project/Bioinformatic Analysis/gscreend_package/Gene_level/Input data/pooled_screen_1_counts.csv")


#'  
#'  #Check point (number of transcripts that we have) n=734
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
length(unique(pooled_screen_1_counts$Gene))

#' 
#' 
#' 2. Put them in order (Same experiments per TP and replicates) and calculate the Mean of the control counts
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
master_table<- pooled_screen_1_counts %>% 
  rowwise() %>% 
  mutate(Mean_Sym_pool=mean(c(`Sym-pool-1-1`, `Sym-pool-1-2`, `Sym-pool-1-3`)),
         Mean_Sym_pool=round(Mean_Sym_pool,0), 
         .after=`Sym-pool-1-3`)


#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
master_table_1<-master_table %>% 
  pivot_longer(`sym-01-01-TP1`:`sym-03-03-TP6`)

separate_table<-separate_wider_delim(master_table_1, name, "-", 
                     names = c("Project", "Experiment", "Replicate", "TP")) %>% 
  arrange(Experiment, TP)

Master_table_arranged_data<-unite(separate_table, col="Name", Project, Experiment, Replicate, TP, sep="-") %>% 
  pivot_wider(names_from = Name, values_from = value)

#' 
#' 3. For running Gscreend (gene level)
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

list<- list()
k<-1 
exps<- c("01", "02", "03")
for (exp_i in exps){
  for(tp_j in c("TP1", "TP2", "TP3", "TP4", "TP5", "TP6")){
    print(paste0("sym-", exp_i, "-01-", tp_j))
counts_matrix<-cbind(Master_table_arranged_data$Mean_Sym_pool, Master_table_arranged_data[[paste0("sym-", exp_i, "-01-", tp_j)]], Master_table_arranged_data[[paste0("sym-", exp_i, "-02-", tp_j)]], Master_table_arranged_data[[paste0("sym-", exp_i, "-03-", tp_j)]])
  
raw_Data<- data.frame(sgRNA_id=Master_table_arranged_data$unique_guide_name, gene=Master_table_arranged_data$Gene)

col_Data<- data.frame(samplename=c("library", "R1", "R2", "R3"), 
                      timepoint=c("T0", "T1", "T1", "T1"))

se_sym_01_TP1<-SummarizedExperiment(assays=list(counts=counts_matrix), rowData=raw_Data, colData=col_Data)

pse_sym_01_TP1<- createPoolScreenExp(se_sym_01_TP1)

pse_an_sym_01_TP1<- RunGscreend(pse_sym_01_TP1)

df<- data.frame(ResultsTable(pse_an_sym_01_TP1, direction = "negative")) %>% 
   mutate(exp=exp_i, tp=tp_j)
list[[k]]<- df
k<-k+1
  
    }
}

#Master table with all the results
all_data_df<-list_rbind(list) %>% 
  pivot_longer(fdr:lfc, names_to = "variable") %>% 
  mutate(sample_name=paste0(variable, "_", exp, "_", tp))


#' 
#' 6. Formatting the master table 
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for(exp_i in exps){

df_results <- filter(all_data_df, exp == exp_i) |>
  select(-exp, -tp, -variable) |>
  pivot_wider(names_from = sample_name, values_from = value)

#Save the table for each experiment
write_csv(df_results, paste0('Gscreendresults_GL_exp_', exp_i, ".csv"))

}

#' 
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

#' 
#' 
#' 
#' 
#' 
#' 
#' 
