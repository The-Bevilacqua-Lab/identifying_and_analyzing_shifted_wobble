library(readxl)
library(ggplot2)
library(dplyr)

library(pROC)
library(tidyr)
library(tidyverse)

library(ggpubr)
library(rstatix)
library(ggvenn)

################################# PLOTS FOR PAPER #####################################

ecoli <- read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Lab/Research/Computational/GU_regs_kd/dms_reactivities.xlsx", sheet = "ecoli", na = "nan")
hsapi <- read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Lab/Research/Computational/GU_regs_kd/dms_reactivities.xlsx", sheet = "hsapien", na = "nan")
scere <- read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Lab/Research/Computational/GU_regs_kd/dms_reactivities.xlsx", sheet = "yeast", na = "nan")

#define register 3 dataframes
ecoliUr3 <- ecoli[c(1206, 797, 3966), ]
ecoliGr3 <- ecoli[c(1219, 833, 3958), ]
hsapiUr3 <- hsapi[c(2642), ]
hsapiGr3 <- hsapi[c(2494), ]
scereUr3 <- scere[c(4954, 2520), ]
scereGr3 <- scere[c(5031, 2572), ]




gu_paired_overlap <- function(reactdf,Ur3data,Gr3data,ypos,binnum=50){
  ##subset paired Us and Gs
  Udf <- subset(reactdf, sequence == "U" & pairing == "1")
  Gdf <- subset(reactdf, sequence == "G" & pairing == "1")
  ##calculate bin width
  Ubinw = diff(range(Udf$DMS, na.rm = TRUE))/binnum
  Gbinw = diff(range(Gdf$DMS, na.rm = TRUE))/binnum
  ##define data for lines
  Ur3_lines <- data.frame(
    x = as.numeric(Ur3data$DMS),  # Use the column name directly
    label = paste(Ur3data[[1]], Ur3data[[2]]))
  Gr3_lines <- data.frame(
    x = as.numeric(Gr3data$DMS),  # Use the column name directly
    label = paste(Gr3data[[1]], Gr3data[[2]]))
  ##plot
  ggplot(Udf,aes(x= DMS))+
    geom_histogram(binwidth = Ubinw,aes(y=..density..), fill = "green3", alpha = 0.5)+
    geom_histogram(binwidth = Gbinw,data=Gdf,aes(x=DMS,y=..density..),fill = "gold", alpha = 0.5)+
    #geom_density()+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "Histogram", x = "DMS Reactivity", y = "count")+
    geom_vline(data = Ur3_lines, aes(xintercept = x), color = "green3") +
    geom_vline(data = Gr3_lines, aes(xintercept = x), color = "gold") +
    geom_text(data = Ur3_lines, aes(x = x, y = ypos, label = label),
              vjust = 0, hjust = 0, angle = 90, color = "green4", size = 3.5)+
    geom_text(data = Gr3_lines, aes(x = x, y = 3, label = label),
              vjust = 0, hjust = 0, angle = 90, color = "gold3", size = 3.5)+
    coord_cartesian(xlim=c(-0.05, 1.01),ylim=c(0,10))+
    scale_x_continuous(breaks=c(0.00,0.25, 0.50, 0.75, 1.00))
}

ecoliGU_overlap.plot <- gu_paired_overlap(ecoli,ecoliUr3,ecoliGr3,ypos = 6.5)
ecoliGU_overlap.plot

gu_all_overlap <- function(reactdf,Ur3data,Gr3data,ypos,binnum=50){
  ##subset paired Us and Gs
  Udf <- subset(reactdf, sequence == "U")
  Gdf <- subset(reactdf, sequence == "G")
  ##calculate bin width
  Ubinw = diff(range(Udf$DMS, na.rm = TRUE))/binnum
  Gbinw = diff(range(Gdf$DMS, na.rm = TRUE))/binnum
  ##define data for lines
  Ur3_lines <- data.frame(
    x = as.numeric(Ur3data$DMS),  # Use the column name directly
    label = paste(Ur3data[[1]], Ur3data[[2]]))
  Gr3_lines <- data.frame(
    x = as.numeric(Gr3data$DMS),  # Use the column name directly
    label = paste(Gr3data[[1]], Gr3data[[2]]))
  ##plot
  ggplot(Udf,aes(x= DMS))+
    geom_histogram(binwidth = Ubinw,aes(y=..density..), fill = "green3", alpha = 0.5)+
    geom_histogram(binwidth = Gbinw,data=Gdf,aes(x=DMS,y=..density..),fill = "gold", alpha = 0.5)+
    #geom_density()+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(x = "DMS Reactivity", y = "count")+
    geom_vline(data = Ur3_lines, aes(xintercept = x), color = "green3") +
    geom_vline(data = Gr3_lines, aes(xintercept = x), color = "gold") +
    geom_text(data = Ur3_lines, aes(x = x, y = ypos, label = label),
              vjust = 0, hjust = 0, angle = 90, color = "green4", size = 3.5)+
    geom_text(data = Gr3_lines, aes(x = x, y = 5, label = label),
              vjust = 0, hjust = 0, angle = 90, color = "gold3", size = 3.5)+
    coord_cartesian(xlim=c(-0.1, 1.01),ylim=c(0,8))+
    scale_x_continuous(breaks=c(0.00,0.25, 0.50, 0.75, 1.00))
}

ecoliGU_all_overlap.plot <- gu_all_overlap(ecoli,ecoliUr3,ecoliGr3,ypos = 6)
ecoliGU_all_overlap.plot
hsapiGU_all_overlap.plot <- gu_all_overlap(hsapi,hsapiUr3,hsapiGr3,ypos = 5)
hsapiGU_all_overlap.plot
scereGU_all_overlap.plot <- gu_all_overlap(scere,scereUr3,scereGr3,ypos = 5)
scereGU_all_overlap.plot

#just to check because the yeast plot looks weird for g and u
acgu_all_overlap <- function(reactdf,binnum=50){
  ##subset paired Us and Gs
  Udf <- subset(reactdf, sequence == "U")
  Gdf <- subset(reactdf, sequence == "G")
  Adf <- subset(reactdf, sequence == "A")
  Cdf <- subset(reactdf, sequence == "C")
  #calculate binwidth
  binw = 1.1/binnum
  ##plot
  ggplot(Udf,aes(x= DMS))+
    geom_histogram(binwidth = binw,aes(y=..density..), fill = "green3", alpha = 0.5)+
    geom_histogram(binwidth = binw,data=Gdf,aes(x=DMS,y=..density..),fill = "gold", alpha = 0.5)+
    geom_histogram(binwidth = binw,data=Adf,aes(x=DMS,y=..density..),fill = "red", alpha = 0.5)+
    geom_histogram(binwidth = binw,data=Cdf,aes(x=DMS,y=..density..),fill = "cadetblue1", alpha = 0.5)+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "Histogram", x = "DMS Reactivity", y = "count")+
    coord_cartesian(xlim=c(-0.1, 1.01),ylim=c(0,10))+
    scale_x_continuous(breaks=c(0.00,0.25, 0.50, 0.75, 1.00))
}

scere_check.plot <- acgu_all_overlap(scere)
scere_check.plot
ecoli_check.plot <- acgu_all_overlap(ecoli)
ecoli_check.plot
hsapi_check.plot <- acgu_all_overlap(hsapi)
hsapi_check.plot


pval <- function(fulldf,Ur3df,Gr3df) {
  dfU <- subset(fulldf, sequence == "U")
  # Calculate mean and standard deviation for the DMS column
  Umean_val <- mean(dfU$DMS, na.rm = TRUE)
  Ustdev_val <- sd(dfU$DMS, na.rm = TRUE)
  # Calculate z-scores
  Ur3df$z_score <- (Ur3df$DMS - Umean_val) / Ustdev_val
  # Calculate one-tailed p-values for values greater than the mean
  Ur3df$pval <- 1 - pnorm(Ur3df$z_score)
  #repeat for G
  dfG <- subset(fulldf, sequence == "G")
  Gmean_val <- mean(dfG$DMS, na.rm = TRUE)
  Gstdev_val <- sd(dfG$DMS, na.rm = TRUE)
  Gr3df$z_score <- (Gr3df$DMS - Gmean_val) / Gstdev_val
  Gr3df$pval <- 1 - pnorm(Gr3df$z_score)
  #combine dataframes
  r3df <- rbind(Ur3df,Gr3df)
  print(r3df)
}

ecoli_stat <- pval(ecoli,ecoliUr3,ecoliGr3)
hsapi_stat <- pval(hsapi,hsapiUr3,hsapiGr3)
scere_stat <- pval(scere,scereUr3,scereGr3)

#######################################################################################

############################### OTHER PLOTS I TRIED ###################################

ecolirRNA_react <- read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Lab/Research/Computational/GU_regs_kd/ecolirRNA_react.xlsx", 
                              col_types = c("text", "text","numeric", "text", "numeric", "text", "text", "numeric", "numeric", "numeric", "text"), na = "nan")

ecolirRNA_react <- read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Lab/Research/Computational/GU_regs_kd/ecolirRNA_react.xlsx", na = "nan")
##for some reason ETC is not numeric, check the formatting


#subset Gs and Us
Gdf <- subset(ecolirRNA_react, seq == "G")  #1440 points
Udf <- subset(ecolirRNA_react, seq == "U")  #923 points

Gpaired <- subset(Gdf, pairstate == 0)  #1001
Upaired <- subset(Udf, pairstate == 0)  #564

Gwob <- subset(Gdf, pair_type == "GU")  #191 points
Uwob <- subset(Udf, pair_type == "GU")  #192 points

Gterm <- subset(Gdf, structure == "terminal") #370
Uterm <- subset(Udf, structure == "terminal") #151

Ur3 <- subset(Uwob, color == "3")
Gr3 <- subset(Gwob, color == "3")

#ggplot scatterplot
Ucand <- subset(Uwob, DMS > quantile(Uwob$DMS,0.95) & ETC < quantile(Uwob$ETC,0.10))
Gcand <- subset(Gwob, DMS > quantile(Gwob$DMS,0.95) & ETC < quantile(Gwob$ETC,0.10))

gu_scatter <- function(alldata,r3data,labdata,col){
  ggplot(alldata, aes(x=ETC, y = DMS, fill = color))+
    #geom_jitter(size = 0.6, alpha = 1)+
    geom_jitter(color = "black", shape = 1, alpha = 0.6,size = 2)+
    geom_jitter(data = r3data, aes(x = ETC, y = DMS), color = col, size = 3, shape = 16) +
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    xlab("ETC Reactivity")+
    ylab("DMS Reactivity")+
    geom_text(data = labdata, aes(label = paste(subunit, nt), vjust = 1, hjust = 0, color = "blue", size = 2, fontface = "bold")) # Label points with y > 5
}

G_scatter.plot <- gu_scatter(Gwob,Gr3,Gcand,"gold")
G_scatter.plot
U_scatter.plot <- gu_scatter(Uwob,Ur3,Ucand,"green3")
U_scatter.plot

G_scatter.plot <- gu_scatter(Gdf,Gr3,Gcand,"gold")
G_scatter.plot
U_scatter.plot <- gu_scatter(Udf,Ur3,Ucand,"green3")
U_scatter.plot

G_scatter.plot <- gu_scatter(Gpaired,Gr3,Gcand,"gold")
G_scatter.plot
U_scatter.plot <- gu_scatter(Upaired,Ur3,Ucand,"green3")
U_scatter.plot




#histogram type thing to see where our points of interest fall in terms of other reactivities
gu_dotplot <- function(df,col){
  ggplot(df, aes(x = DMS))+
    #geom_jitter(size = 0.6, alpha = 1)+
    geom_dotplot(aes(fill = factor(color)), alpha = 1,size = 1)+
    scale_fill_manual(values = c("black",col))+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    ylab("DMS Reactivity")
}



#make a dotplot histogram with r3 Gs and Us colored
gu_dotplot <- function (df, col, binnum){
  binw= diff(range(df$DMS))/binnum
  ggplot(df, aes(x=DMS, color = color, fill = color))+
    geom_dotplot(binwidth = binw, stackgroups = TRUE, method = "histodot")+
    scale_fill_manual(values = c("gray",col))+
    scale_color_manual(values = c("gray",col))+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "Dotplot", x = "DMS", y = "count")
}
    
Gdot.plot <- gu_dotplot(Gwob,"gold",50)
Gdot.plot
Udot.plot <- gu_dotplot(Uwob,"green3",60)
Udot.plot



#make a density plot with lines for r3 Gs and Us



gu_hist <- function (df, labdata, col, binnum = 50){
  binw= diff(range(df$DMS, na.rm = TRUE))/binnum
  lab1 = unname(unlist(labdata[1,c(2,3,8)]))
  lab2 = unname(unlist(labdata[2,c(2,3,8)]))
  lab3 = unname(unlist(labdata[3,c(2,3,8)]))
  ggplot(df, aes(x=DMS))+
    geom_histogram(binwidth = binw, aes(y=..density..))+
    #geom_density()+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "Histogram", x = "DMS Reactivity", y = "count")+
    geom_vline(xintercept = c(as.numeric(lab1[3]), as.numeric(lab2[3]), as.numeric(lab3[3])), color = "green3")+
    geom_text(label = paste(lab1[1],lab1[2]), y = 42, x = as.numeric(lab1[3]), vjust = 0, hjust=1, angle = 90, color = "green4", size = 3.5)+
    geom_text(label = paste(lab2[1],lab2[2]), y = 42, x = as.numeric(lab2[3]), vjust = 0, hjust=1, angle = 90, color = "green4", size = 3.5)+
    geom_text(label = paste(lab3[1],lab3[2]), y = 42, x = as.numeric(lab3[3]), vjust = 0, hjust=1, angle = 90, color = "green4", size = 3.5)+
    coord_cartesian(xlim = c(-0.25,1))
}

Udens.plot <- gu_hist(Upaired,Ur3)
Udens.plot

###### THIS IS MY FAVORITE PLOT TO SHOW THIS TTREND ##########
gu_histwindz <- function (df, labdata, col1, col2, ypos, binnum = 50){
  binw= diff(range(df$DMSwindz, na.rm = TRUE))/binnum
  lab1 = unname(unlist(labdata[1,c(2,3,8)]))
  lab2 = unname(unlist(labdata[2,c(2,3,8)]))
  lab3 = unname(unlist(labdata[3,c(2,3,8)]))
  mean = mean(df$DMS, na.rm = TRUE)
  median = median(df$DMS, na.rm=TRUE)
  print (mean)
  ggplot(df, aes(x=DMSwindz))+
    geom_histogram(binwidth = binw,aes(y=..density..))+
    #geom_density()+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "Histogram", x = "DMS Reactivity", y = "count")+
    geom_vline(xintercept = c(as.numeric(lab1[3]), as.numeric(lab2[3]), as.numeric(lab3[3])), color = col1)+
    geom_text(label = paste(lab1[1],lab1[2]), y = ypos, x = as.numeric(lab1[3]), vjust = 0, hjust=0, angle = 90, color = col2, size = 3.5)+
    geom_text(label = paste(lab2[1],lab2[2]), y = ypos, x = as.numeric(lab2[3]), vjust = 0, hjust=0, angle = 90, color = col2, size = 3.5)+
    geom_text(label = paste(lab3[1],lab3[2]), y = ypos, x = as.numeric(lab3[3]), vjust = 0, hjust=0, angle = 90, color = col2, size = 3.5)+
    coord_cartesian(xlim=c(-0.25, 1.01),ylim=c(0,10))+
    scale_x_continuous(breaks=c(-0.25,0.00,0.25, 0.50, 0.75, 1.00))
}
Upaired <- Upaired %>% mutate(DMSwindz= ifelse (DMS > 1.2,1.2,DMS))
Uhist.plot <- gu_histwindz(Upaired,Ur3,"green3","green4", ypos = 6.5)
Uhist.plot

Gpaired <- Gpaired %>% mutate(DMSwindz= ifelse (DMS > 1.2,1.2,DMS))
Ghist.plot <- gu_histwindz(Gpaired,Gr3, "gold", "gold3", ypos = 6.5)
Ghist.plot

gu_windz_overlap <- function(Udf,Gdf,Ulabdata,Glabdata,ypos,binnum=50){
  Ubinw= diff(range(Udf$DMSwindz, na.rm = TRUE))/binnum
  Ulab1 = unname(unlist(Ulabdata[1,c(2,3,8)]))
  Ulab2 = unname(unlist(Ulabdata[2,c(2,3,8)]))
  Ulab3 = unname(unlist(Ulabdata[3,c(2,3,8)]))
  Gbinw= diff(range(Gdf$DMSwindz, na.rm = TRUE))/binnum
  Glab1 = unname(unlist(Glabdata[1,c(2,3,8)]))
  Glab2 = unname(unlist(Glabdata[2,c(2,3,8)]))
  Glab3 = unname(unlist(Glabdata[3,c(2,3,8)]))
  ggplot(Udf,aes(x= DMSwindz))+
    geom_histogram(binwidth = Ubinw,aes(y=..density..), fill = "green3", alpha = 0.5)+
    geom_histogram(binwidth = Gbinw,data=Gdf,aes(x=DMSwindz,y=..density..),fill = "gold", alpha = 0.5)+
    #geom_density()+
    theme_classic(base_size = 22, base_line_size = 0.5)+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
    theme(text = element_text(size = 12, family = "Arial", color = "black"))+
    theme(strip.background = element_rect(fill="white", colour="white", linewidth =1))+
    theme(axis.text = element_text(color = "black"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")+
    labs(title = "Histogram", x = "DMS Reactivity", y = "count")+
    geom_vline(xintercept = c(as.numeric(Ulab1[3]), as.numeric(Ulab2[3]), as.numeric(Ulab3[3])), color = "green3")+
    geom_vline(xintercept = c(as.numeric(Glab1[3]), as.numeric(Glab2[3]), as.numeric(Glab3[3])), color = "gold")+
    geom_text(label = paste(Ulab1[1],Ulab1[2]), y = ypos, x = as.numeric(Ulab1[3]), vjust = 0, hjust=0, angle = 90, color = "green4", size = 3.5)+
    geom_text(label = paste(Ulab2[1],Ulab2[2]), y = ypos, x = as.numeric(Ulab2[3]), vjust = 0, hjust=0, angle = 90, color = "green4", size = 3.5)+
    geom_text(label = paste(Ulab3[1],Ulab3[2]), y = ypos, x = as.numeric(Ulab3[3]), vjust = 0, hjust=0, angle = 90, color = "green4", size = 3.5)+
    geom_text(label = paste(Glab1[1],Glab1[2]), y = 3, x = as.numeric(Glab1[3]), vjust = 0, hjust=0, angle = 90, color = "gold3", size = 3.5)+
    geom_text(label = paste(Glab2[1],Glab2[2]), y = 3, x = as.numeric(Glab2[3]), vjust = 0, hjust=0, angle = 90, color = "gold3", size = 3.5)+
    geom_text(label = paste(Glab3[1],Glab3[2]), y = 3, x = as.numeric(Glab3[3]), vjust = 0, hjust=0, angle = 90, color = "gold3", size = 3.5)+
    coord_cartesian(xlim=c(-0.25, 1.01),ylim=c(0,10))+
    scale_x_continuous(breaks=c(-0.25,0.00,0.25, 0.50, 0.75, 1.00))
}

GU_overlap.plot <- gu_windz_overlap(Upaired,Gpaired,Ur3,Gr3,ypos = 6.5)
GU_overlap.plot






##empirical p-value calculations (how many standard deviations from the mean)

pval <- function(df) {
  # Subset rows where color == "3"
  r3df <- subset(df, df$color == "3")
  newdf <- r3df[, c("rna", "nt", "DMS")]
  # Calculate mean and standard deviation for the DMS column
  mean_val <- mean(df$DMS, na.rm = TRUE)
  stdev_val <- sd(df$DMS, na.rm = TRUE)
  # Calculate z-scores
  newdf$z_score <- (newdf$DMS - mean_val) / stdev_val
  # Calculate one-tailed p-values for values greater than the mean
  newdf$pval <- 1 - pnorm(newdf$z_score)
  print(newdf)
}

stat <- pval(Upaired)
stat2 <- pval(Gpaired)



### get the percentile of each of the reactivities

#add percentile row
Upairedecdf <- ecdf(Upaired$DMS)
Ur3$paired_percentile <- Upairedecdf(Ur3$DMS)
Gpairedecdf <- ecdf(Gpaired$DMS)
Gr3$paired_percentile <- Gpairedecdf(Ur3$DMS)
Uwobecdf <- ecdf(Uwob$DMS)
Gwobecdf <- ecdf(Gwob$DMS)
Ur3$wobble_percentile <- Uwobecdf(Ur3$DMS)
Gr3$wobble_percentile <- Gwobecdf(Gr3$DMS)

Uallecdf <- ecdf(Udf$DMS)
Gallecdf <- ecdf(Gdf$DMS)
Ur3$all_percentile <- Uallecdf(Ur3$DMS)
Gr3$all_percentile <- Gallecdf(Gr3$DMS)


r3comb <- bind_rows(Ur3,Gr3) %>% arrange(rna,nt) %>% select(-color,-ETC,-SHAPE,-subunit,-pairstate,-pair_type)

write.csv(r3comb, file = "/Users/catherinedouds/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Lab/Research/Computational/GU_regs_kd/figures/reactivity_percentile.csv", row.names = FALSE)



