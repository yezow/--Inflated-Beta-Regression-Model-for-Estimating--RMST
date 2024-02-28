current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
###############################
#set your model parameters for setting 1
#pi model parameters
beta_0_pi<- -1 #-1 for setting 1 and 0.5 for setting 2
beta_1_pi<- c(1,2,-1.5) #c(1,2,-1.5) for setting 1 and c(0,0,0) for setting 2
#mu model parameters
alpha_0<- -2
alpha_1<-c(1.2,2)

#number of subjects
nsubj<-100

source('Biometrical_sim.R')
write.csv(results_data_table1[-nrow(results_data_table1),],file = "Table_1_results_100_setting1.csv")
write.csv(results_data_table2,file = "Table_2_results_100_setting1.csv")

#number of subjects
nsubj<-500

source('Biometrical_sim.R')
write.csv(results_data_table1[-nrow(results_data_table1),],file = "Table_1_results_500_setting1.csv")
write.csv(results_data_table2,file = "Table_2_results_500_setting1.csv")

#number of subjects
nsubj<-1000

source('Biometrical_sim.R')
write.csv(results_data_table1[-nrow(results_data_table1),],file = "Table_1_results_1000_setting1.csv")
write.csv(results_data_table2,file = "Table_2_results_1000_setting1.csv")

#number of subjects
nsubj<-1500

source('Biometrical_sim.R')
write.csv(results_data_table1[-nrow(results_data_table1),],file = "Table_1_results_1500_setting1.csv")
write.csv(results_data_table2,file = "Table_2_results_1500_setting1.csv")

###############################
#set your model parameters for setting 2
#pi model parameters
beta_0_pi<- 0.5 
beta_1_pi<- c(0,0,0) 
#mu model parameters
alpha_0<- -2
alpha_1<-c(1.2,2)

#number of subjects
nsubj<-500

source('Biometrical_sim.R')
write.csv(results_data_table2,file = "Table_2_results_500_setting2.csv")

