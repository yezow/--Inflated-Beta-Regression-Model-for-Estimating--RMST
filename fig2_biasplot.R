current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
###for 100 subjects
nsubj<-100
set.seed(123)
source("source_figure2.R")
restricted_mean_CI_RMST_pseudo_empirical
restricted_mean_CI_MI_empirical_indep
restricted_mean_CI_EM_empirical_indep

x <- seq(1:nsubj)

RMST_diff <- Pseudo_RMST_mean_tau_t_est - Mean_tau_T
EM_diff<- Mean_tau_T_est_EM-Mean_tau_T
MI_diff <- Mean_tau_T_est_MI-Mean_tau_T
data_mean_tau_T<-data.frame(Mean_tau_T)
data1_EM_100<-data.frame(x,RMST_diff,EM_diff)
data2_EM_100 <- data.frame(which_out_CI,EM_diff[which_out_CI])
data3_RMST_100 <- data.frame(RMST_which_out_CI,RMST_diff[RMST_which_out_CI])
data1_MI_100<-data.frame(x,RMST_diff,MI_diff)
data2_MI_100 <- data.frame(which_out_CI_MI,MI_diff[which_out_CI_MI])

###for 500 subjects
nsubj<-500
set.seed(2)
source("source_figure2.R")
restricted_mean_CI_RMST_pseudo_empirical
restricted_mean_CI_MI_empirical_indep
restricted_mean_CI_EM_empirical_indep

x <- seq(1:nsubj)

RMST_diff <- Pseudo_RMST_mean_tau_t_est - Mean_tau_T
EM_diff<- Mean_tau_T_est_EM-Mean_tau_T
MI_diff <- Mean_tau_T_est_MI-Mean_tau_T
data_mean_tau_T<-data.frame(Mean_tau_T)
data1_EM_500<-data.frame(x,RMST_diff,EM_diff)
data2_EM_500 <- data.frame(which_out_CI,EM_diff[which_out_CI])
data3_RMST_500 <- data.frame(RMST_which_out_CI,RMST_diff[RMST_which_out_CI])
data1_MI_500<-data.frame(x,RMST_diff,MI_diff)
data2_MI_500 <- data.frame(which_out_CI_MI,MI_diff[which_out_CI_MI])

###for 1000 subjects
nsubj<-1000
set.seed(12)
source("source_figure2.R")
restricted_mean_CI_RMST_pseudo_empirical
restricted_mean_CI_MI_empirical_indep
restricted_mean_CI_EM_empirical_indep

x <- seq(1:nsubj)

RMST_diff <- Pseudo_RMST_mean_tau_t_est - Mean_tau_T
EM_diff<- Mean_tau_T_est_EM-Mean_tau_T
MI_diff <- Mean_tau_T_est_MI-Mean_tau_T
data_mean_tau_T<-data.frame(Mean_tau_T)
data1_EM_1000<-data.frame(x,RMST_diff,EM_diff)
data2_EM_1000 <- data.frame(which_out_CI,EM_diff[which_out_CI])
data3_RMST_1000 <- data.frame(RMST_which_out_CI,RMST_diff[RMST_which_out_CI])
data1_MI_1000<-data.frame(x,RMST_diff,MI_diff)
data2_MI_1000 <- data.frame(which_out_CI_MI,MI_diff[which_out_CI_MI])

###for 1500 subjects
nsubj<-1500
set.seed(12)
source("source_figure2.R")
restricted_mean_CI_RMST_pseudo_empirical
restricted_mean_CI_MI_empirical_indep
restricted_mean_CI_EM_empirical_indep

x <- seq(1:nsubj)

RMST_diff <- Pseudo_RMST_mean_tau_t_est - Mean_tau_T
EM_diff<- Mean_tau_T_est_EM-Mean_tau_T
MI_diff <- Mean_tau_T_est_MI-Mean_tau_T
data_mean_tau_T<-data.frame(Mean_tau_T)
data1_EM_1500<-data.frame(x,RMST_diff,EM_diff)
data2_EM_1500 <- data.frame(which_out_CI,EM_diff[which_out_CI])
data3_RMST_1500 <- data.frame(RMST_which_out_CI,RMST_diff[RMST_which_out_CI])
data1_MI_1500<-data.frame(x,RMST_diff,MI_diff)
data2_MI_1500 <- data.frame(which_out_CI_MI,MI_diff[which_out_CI_MI])
# write.csv(data1_EM,"/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper1\ result/paper1_bias_plots/1500_bias_data1_EM.csv")
# write.csv(data2_EM,"/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper1\ result/paper1_bias_plots/1500_bias_data2_EM.csv")
# write.csv(data3_RMST,"/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper1\ result/paper1_bias_plots/1500_bias_data_RMST.csv")
# write.csv(data1_MI,"/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper1\ result/paper1_bias_plots/1500_bias_data1_MI.csv")
# write.csv(data2_MI,"/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper1\ result/paper1_bias_plots/1500_bias_data2_MI.csv")
# write.csv(data_mean_tau_T,"/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper1\ result/paper1_bias_plots/1500_bias_data_mean_tau_T.csv")

##################################
#plots
##################################
library(ggplot2)
library(gridExtra)
library(cowplot)
data1<-data1_EM_100
data2<-data2_EM_100
data3<-data3_RMST_100

nsubj <- 100
x_value <- c(data1$x[-data2$which_out_CI],data1$x[-data3$RMST_which_out_CI],data2$which_out_CI,data3$RMST_which_out_CI)
y_value <- c(data1$EM_diff[-data2$which_out_CI],data1$RMST_diff[-data3$RMST_which_out_CI],data2$EM_diff.which_out_CI., data3$RMST_diff.RMST_which_out_CI.)
group <- c(rep("a-IBR model",nsubj-nrow(data2)),rep("a-RMST model",nsubj-nrow(data3)),rep("a-IBR outliers",nrow(data2)),rep("a-RMST outliers",nrow(data3)))
data_bias_plot <- data.frame(x_value,y_value,group)
data_bias_plot$group<-factor(data_bias_plot$group,levels=c("a-IBR model","a-IBR outliers","a-RMST model","a-RMST outliers"))

bias_100_1 <- ggplot(data_bias_plot, aes(x=x_value, y=y_value, group=group)) +
  geom_point(aes(shape=group, color=group),size=3.2)+
  scale_shape_manual(values=c(16, 17, 15,8))+
  scale_color_manual(values = c('yellowgreen',"orange2","steelblue3",'violetred2'))+
  theme(legend.title = element_blank(),legend.position="top",text = element_text(size=24))+
  #theme(legend.title = element_blank(),text = element_text(size=15),legend.position="none")+
  labs(x="Subjects",y=expression("Difference between fitted and true "* tau *"-RMST")) 

bias_100 <- ggplot(data_bias_plot, aes(x=x_value, y=y_value, group=group)) +
  geom_point(aes(shape=group, color=group),size=3.2)+
  scale_shape_manual(values=c(16, 17, 15,8))+
  scale_color_manual(values = c('yellowgreen',"orange2","steelblue3",'violetred2'))+
  theme(legend.title = element_blank(),text = element_text(size=19),legend.position="none",axis.text = element_text(size = 19))+
  labs(x="Subjects",y=expression("Difference between fitted and true "* tau *"-RMST"))  

data1<-data1_EM_500
data2<-data2_EM_500
data3<-data3_RMST_500

nsubj <- 500
x_value <- c(data1$x[-data2$which_out_CI],data1$x[-data3$RMST_which_out_CI],data2$which_out_CI,data3$RMST_which_out_CI)
y_value <- c(data1$EM_diff[-data2$which_out_CI],data1$RMST_diff[-data3$RMST_which_out_CI],data2$EM_diff.which_out_CI., data3$RMST_diff.RMST_which_out_CI.)
group <- c(rep("a-IBR model",nsubj-nrow(data2)),rep("a-RMST model",nsubj-nrow(data3)),rep("a-IBR outliers",nrow(data2)),rep("a-RMST outliers",nrow(data3)))
data_bias_plot <- data.frame(x_value,y_value,group)
data_bias_plot$group<-factor(data_bias_plot$group,levels=c("a-IBR model","a-IBR outliers","a-RMST model","a-RMST outliers"))

bias_500 <- ggplot(data_bias_plot, aes(x=x_value, y=y_value, group=group)) +
  geom_point(aes(shape=group, color=group),size=3.2)+
  scale_shape_manual(values=c(16, 17, 15,8))+
  scale_color_manual(values = c('yellowgreen',"orange2","steelblue3",'violetred2'))+
  theme(legend.title = element_blank(),text = element_text(size=19),legend.position="none",axis.text = element_text(size = 19))+
  labs(x="Subjects",y=expression("Difference between fitted and true "* tau *"-RMST"))  

data1<-data1_EM_1000
data2<-data2_EM_1000
data3<-data3_RMST_1000

nsubj <- 1000
x_value <- c(data1$x[-data2$which_out_CI],data1$x[-data3$RMST_which_out_CI],data2$which_out_CI,data3$RMST_which_out_CI)
y_value <- c(data1$EM_diff[-data2$which_out_CI],data1$RMST_diff[-data3$RMST_which_out_CI],data2$EM_diff.which_out_CI., data3$RMST_diff.RMST_which_out_CI.)
group <- c(rep("a-IBR model",nsubj-nrow(data2)),rep("a-RMST model",nsubj-nrow(data3)),rep("a-IBR outliers",nrow(data2)),rep("a-RMST outliers",nrow(data3)))
data_bias_plot <- data.frame(x_value,y_value,group)
data_bias_plot$group<-factor(data_bias_plot$group,levels=c("a-IBR model","a-IBR outliers","a-RMST model","a-RMST outliers"))

bias_1000 <- ggplot(data_bias_plot, aes(x=x_value, y=y_value, group=group)) +
  geom_point(aes(shape=group, color=group),size=3.2)+
  scale_shape_manual(values=c(16, 17, 15,8))+
  scale_color_manual(values = c('yellowgreen',"orange2","steelblue3",'violetred2'))+
  theme(legend.title = element_blank(),text = element_text(size=19),legend.position="none",axis.text = element_text(size = 19))+
  labs(x="Subjects",y=expression("Difference between fitted and true "* tau *"-RMST"))  


data1<-data1_EM_1500
data2<-data2_EM_1500
data3<-data3_RMST_1500
nsubj <- 1500
x_value <- c(data1$x[-data2$which_out_CI],data1$x[-data3$RMST_which_out_CI],data2$which_out_CI,data3$RMST_which_out_CI)
y_value <- c(data1$EM_diff[-data2$which_out_CI],data1$RMST_diff[-data3$RMST_which_out_CI],data2$EM_diff.which_out_CI., data3$RMST_diff.RMST_which_out_CI.)
group <- c(rep("a-IBR model",nsubj-nrow(data2)),rep("a-RMST model",nsubj-nrow(data3)),rep("a-IBR outliers",nrow(data2)),rep("a-RMST outliers",nrow(data3)))
data_bias_plot <- data.frame(x_value,y_value,group)
data_bias_plot$group<-factor(data_bias_plot$group,levels=c("a-IBR model","a-IBR outliers","a-RMST model","a-RMST outliers"))

bias_1500 <- ggplot(data_bias_plot, aes(x=x_value, y=y_value, group=group)) +
  geom_point(aes(shape=group, color=group),size=3.2)+
  scale_shape_manual(values=c(16, 17, 15,8))+
  scale_color_manual(values = c('yellowgreen',"orange2","steelblue3",'violetred2'))+
  theme(legend.title = element_blank(),text = element_text(size=19),legend.position="none",axis.text = element_text(size = 19))+
  labs(x="Subjects",y=expression("Difference between fitted and true "* tau *"-RMST"))  

legend <- get_legend(bias_100_1)
dev.off()
grid.arrange(legend,bias_100, bias_500, bias_1000, bias_1500, ncol=2, nrow=3,layout_matrix = rbind(c(1,1), c(2,3),c(4,5)))

