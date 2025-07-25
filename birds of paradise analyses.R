library(picante)
library(geiger)
library(nlme)
library(car)
library(JNplots)
library(phylolm)
library(lmtest)
library(piecewiseSEM)
library(phytools)
library(scales)

setwd("~/Desktop/birds of paradise/")
data <- read.csv("Paradisaeidae.tail_data.csv",header=T,na.strings = c('NA'))
tree <- read.nexus("Ligon.et.al._UltrametricTree")

#### PHYLOGENY ####
species <- data$species
prunned.tree<-drop.tip(tree, setdiff(tree$tip.label, species))
rownames(data) <- data$species
new <- match.phylo.data(prunned.tree, data)
new_phy <- new$phy
new_data <- new$data
name.check(new_phy, new_data)

#### DATA PREPARATION ####
c <- 5
while(c<13){
  new_data[,c] <- as.numeric(as.character(new_data[,c]))
  c <- c+1
}

# averages
new_data$avg_wing <- (new_data$wing_F+new_data$wing_M)/2
new_data$avg_tail_max <- (new_data$tail_max_F+new_data$tail_max_M)/2
new_data$avg_tarsus <- (new_data$tarsus_M+new_data$tarsus_F)/2
new_data$avg_weight <- (new_data$weight_F+new_data$weight_M)/2

# SSD
new_data$SSD_weight <- (new_data$weight_M/new_data$weight_F)-1
new_data$SSD_tarsus <- (new_data$tarsus_M/new_data$tarsus_F)-1

c <- 16
while(c<20){
  new_data[,c] <- as.numeric(as.character((new_data[,c])))
  c <- c+1
}

##### log10 transform #####
c <- 5
while(c<13){
  new_data[,c] <- log10(new_data[,c])
  c <- c+1
}
c <- 16
while(c<20){
  new_data[,c] <- log10(new_data[,c])
  c <- c+1
}

##### only complete cases in dataset #####
new_data_complete <- new_data[complete.cases(new_data$SSD_weight),]
species2 <- rownames(new_data_complete)
prunned.tree<-drop.tip(tree, setdiff(tree$tip.label, species2))
rownames(new_data_complete) <- new_data_complete$species
new <- match.phylo.data(prunned.tree, new_data_complete)
new_phy_complete <- new$phy
new_data_complete <- new$data
name.check(new_phy_complete, new_data_complete)

##### numeric-character #####
c <- 5
while(c<13){
  new_data_complete[,c] <- as.numeric(as.character((new_data_complete[,c])))
  c <- c+1
}
c <- 16
while(c<22){
  new_data_complete[,c] <- as.numeric(as.character((new_data_complete[,c])))
  c <- c+1
}

new_data_complete$AVONET_average_Wt <- as.numeric(as.character(new_data_complete$AVONET_average_Wt))

c <- 1
while(c<length(new_data_complete[,1])+1){
  if(new_data_complete[c,"category"]=="wired"){
    new_data_complete[c,"category"] <- "long"
    c <- c+1
  }else{
    c <- c+1
  }
}

##### estimation of relative tail dimorphism #####
new_data_complete$rel_tail_length_M <- new_data_complete$tail_max_M-0.33*new_data_complete$weight_M
new_data_complete$rel_tail_length_F <- new_data_complete$tail_max_F-0.33*new_data_complete$weight_F
new_data_complete$rel_tail_dim <- new_data_complete$rel_tail_length_M-new_data_complete$rel_tail_length_F

##### datasets for each category #####

new_data_short <- new_data_complete[which(new_data_complete$category=='short'),]
new_data_long <- new_data_complete[-which(new_data_complete$category=='short'),]

new_data_long$rel_tail_dim <- as.numeric(as.character(new_data_long$rel_tail_dim))
new_data_short$rel_tail_dim <- as.numeric(as.character(new_data_short$rel_tail_dim))

##### datasets for each category for sex comparisons #####

tail_max <- c(new_data_complete$tail_max_M,new_data_complete$tail_max_F)
weight <- c(new_data_complete$weight_M,new_data_complete$weight_F)
sex <- (c(rep("M",34),rep("F",34)))
category <- as.character(c(new_data_complete$category,new_data_complete$category))
new_data_sex <- cbind(tail_max,weight,sex,category)
new_data_sex <- as.data.frame(new_data_sex)

c <- 1
while(c<length(new_data_sex[,1])+1){
  if(new_data_sex[c,"category"]=="wired"){
    new_data_sex[c,"category"] <- "long"
    c <- c+1
  }else{
    c <- c+1
  }
}

new_data_sex$sex <- as.factor(new_data_sex$sex)
new_data_sex$category <- as.factor(new_data_sex$category)

new_data_sex_long <- new_data_sex[which(new_data_sex$category=="long"),]
new_data_sex_short <- new_data_sex[which(new_data_sex$category=="short"),]

#### TESTS ####
par(pty='s')

##### RR in body size (Figure 5A) #####
plot(new_data_complete$avg_weight,new_data_complete$SSD_weight,pch=16,
     xlab='(average species weight)',ylab='SSD',xaxt='n')
axis(1, at = c(log10(100),log10(200),log10(300),log10(400)), las=2,
     labels = c(100,200,300,400))
new_data_complete$Clade <- as.factor(new_data_complete$Clade)

# model
mod <- gls(SSD_weight~(avg_weight)*relevel(Clade,ref="Core"),data=new_data_complete, correlation = corPagel(1,new_phy_complete))
summary(mod)

abline(a=-1.2105356,b=0.6707531)
abline(a=-1.2105356+1.7691271,b=0.6707531-0.8492235,lty=2)

##### RR in relative tail length (Figure 5B) #####
plot(new_data_complete$avg_weight,new_data_complete$rel_tail_dim,pch=16,
     xlab='(average species weight)',ylab='relative tail length dimorphism',xaxt='n')
axis(1, at = c(log10(100),log10(200),log10(300),log10(400)), las=2,
     labels = c(100,200,300,400))
abline(0,0,lty=2)

new_data_complete$category <- as.factor(new_data_complete$category)

# model
mod <- gls(rel_tail_dim~avg_weight+relevel(category,ref="short"),data=new_data_complete, correlation = corPagel(1,new_phy_complete))
summary(mod)

abline(a=mean(new_data_short$rel_tail_dim),b=0,col="grey")
abline(a=mean(new_data_long$rel_tail_dim),b=0,col="black")

##### Allometries per group (Figure 6A) #####

plot(new_data_short$weight_F,new_data_short$tail_max_F,col="black",pch=16, ylim=c(1,3),xlim = c(1.5,3),
     ylab = "", xlab ="",bty="n",xaxt="n",yaxt="n")
points(new_data_short$weight_M,new_data_short$tail_max_M,col="black", bg="grey",pch=21)
axis(1, at = c(log10(5),log10(25),log10(50),log10(100),log10(200),log10(400),log10(800)), las=2,
     labels = c(5,25,50,100,200,400,800),cex.axis=0.8)
axis(2, at = c(log10(5),log10(25),log10(50),log10(100),log10(200),log10(400),log10(800)), las=2,
     labels = c(5,25,50,100,200,400,800),cex.axis=0.8)
grid(lty=2,lwd=1,col="lightgrey")

new_data_sex_short$tail_max <- as.numeric(as.character(new_data_sex_short$tail_max))
new_data_sex_short$weight <- as.numeric(as.character(new_data_sex_short$weight))

# model

# interaction not significant
mod <- lm(new_data_sex_short$tail_max ~ new_data_sex_short$weight * new_data_sex_short$sex)
summary(mod)

# only additive terms
mod <- lm(new_data_sex_short$tail_max ~ new_data_sex_short$weight + new_data_sex_short$sex)
summary(mod)

abline(a=1.17236, b=0.39403, col='black', lwd=2)
abline(a=1.17236-0.04077, b=0.39403, col='grey', lwd=2, lty=2)

##### Allometries per group (Figure 6B) #####

plot(new_data_long$weight_F,new_data_long$tail_max_F, col="black",pch=16, ylim=c(1,3),xlim=c(1.5,3),bty='n',xlab="",ylab="",
     xaxt="n",yaxt="n")
points(new_data_long$weight_M,new_data_long$tail_max_M, col="black", bg="grey", pch=21)
axis(1, at = c(log10(5),log10(25),log10(50),log10(100),log10(200),log10(400),log10(800)), las=2,
     labels = c(5,25,50,100,200,400,800),cex.axis=0.8)
axis(2, at = c(log10(5),log10(25),log10(50),log10(100),log10(200),log10(400),log10(800)), las=2,
     labels = c(5,25,50,100,200,400,800),cex.axis=0.8)
grid(lty=2,lwd=1,col="lightgrey")

new_data_sex_long$tail_max <- as.numeric(as.character(new_data_sex_long$tail_max))
new_data_sex_long$weight <- as.numeric(as.character(new_data_sex_long$weight))

# model

# interaction not significant
mod <- lm(new_data_sex_long$tail_max ~ new_data_sex_long$weight * new_data_sex_long$sex)
summary(mod)

# only additive terms
mod <- lm(new_data_sex_long$tail_max ~ new_data_sex_long$weight + new_data_sex_long$sex)
summary(mod)

abline(a=0.15971, b=0.95412, col='black', lwd=2)
abline(a=0.15971+0.33667, b=0.95412, col='grey', lwd=2, lty=2)


# sex-specific models
plot(new_data_complete$weight_M,new_data_complete$tail_max_M, col=as.factor(new_data_complete$category))
mod <- gls(tail_max_M~weight_M*category, data = new_data_complete, correlation = corPagel(1,new_phy_complete))
summary(mod)
abline(a=0.7565754, b=0.8518058)
abline(a=0.7565754+0.7801535, b=0.8518058-0.6303717, col='red')

plot(new_data_complete$weight_F,new_data_complete$tail_max_F, col=as.factor(new_data_complete$category))
mod <- gls(tail_max_F~weight_F*category, data = new_data_complete, correlation = corPagel(1,new_phy_complete))
summary(mod)
abline(a=0.1340421, b=0.9949092)
abline(a=0.1340421+1.2951372, b=0.9949092-0.7157575, col='red')

##### evol.vcv (Table S2) #####

group <- as.factor(new_data_complete$category)
#group <- as.factor(new_data_complete_SSD_weight$mating_system)
names(group) <- rownames(new_data_complete)

c <- 1
w1 <- c()
m1_rate_11 <- c()
m1_rate_12 <- c()
m1_rate_21 <- c()
m1_rate_22 <- c()
m1_r1 <- c()
m1_r2 <- c()
w2 <- c()
m2_rate_11 <- c()
m2_rate_12 <- c()
m2_rate_21 <- c()
m2_rate_22 <- c()
m2_r1 <- c()
m2_r2 <- c()
w2b <- c()
m2b_rate_11 <- c()
m2b_rate_12 <- c()
m2b_rate_21 <- c()
m2b_rate_22 <- c()
m2b_r1 <- c()
m2b_r2 <- c()
w2c <- c()
m2c_rate_11 <- c()
m2c_rate_12 <- c()
m2c_rate_21 <- c()
m2c_rate_22 <- c()
m2c_r1 <- c()
m2c_r2 <- c()
w3 <- c()
m3_rate_11 <- c()
m3_rate_12 <- c()
m3_rate_21 <- c()
m3_rate_22 <- c()
m3_r1 <- c()
m3_r2 <- c()
w3b <- c()
m3b_rate_11 <- c()
m3b_rate_12 <- c()
m3b_rate_21 <- c()
m3b_rate_22 <- c()
m3b_r1 <- c()
m3b_r2 <- c()
w3c <- c()
m3c_rate_11 <- c()
m3c_rate_12 <- c()
m3c_rate_21 <- c()
m3c_rate_22 <- c()
m3c_r1 <- c()
m3c_r2 <- c()
w4 <- c()
m4_rate_11 <- c()
m4_rate_12 <- c()
m4_rate_21 <- c()
m4_rate_22 <- c()
m4_r1 <- c()
m4_r2 <- c()
while(c<101){
  sm <- make.simmap(new_phy_complete,group)
  fit <- evolvcv.lite(sm,new_data_complete[,22:23],models=c("all models"))
  model1_AIC <- as.numeric(as.character(fit$model1[6]))
  m1_rate_11[c] <- fit$model1$R[1,1]
  m1_rate_12[c] <- NA
  m1_rate_21[c] <- fit$model1$R[2,2]
  m1_rate_22[c] <- NA
  m1_r1[c] <- fit$model1$R[1,2]/sqrt(fit$model1$R[1,1]*fit$model1$R[2,2])
  m1_r2[c] <- NA
  model2_AIC <- as.numeric(as.character(fit$model2[6]))
  m2_rate_11[c] <- fit$model2$R$l[1,1]
  m2_rate_12[c] <- fit$model2$R$s[1,1]
  m2_rate_21[c] <- fit$model2$R$l[2,2]
  m2_rate_22[c] <- fit$model2$R$s[2,2]
  m2_r1[c] <- fit$model2$R$l[1,2]/sqrt(fit$model2$R$l[1,1]*fit$model2$R$l[2,2])
  m2_r2[c] <- NA
  model2b_AIC <- as.numeric(as.character(fit$model2b[6]))
  m2b_rate_11[c] <- fit$model2b$R$l[1,1]
  m2b_rate_12[c] <- fit$model2b$R$s[1,1]
  m2b_rate_21[c] <- fit$model2b$R$l[2,2]
  m2b_rate_22[c] <- NA
  m2b_r1[c] <- fit$model2b$R$l[1,2]/sqrt(fit$model2b$R$l[1,1]*fit$model2b$R$l[2,2])
  m2b_r2[c] <- NA
  model2c_AIC <- as.numeric(as.character(fit$model2c[6]))
  m2c_rate_11[c] <- fit$model2c$R$l[1,1]
  m2c_rate_12[c] <- NA
  m2c_rate_21[c] <- fit$model2c$R$l[2,2]
  m2c_rate_22[c] <- fit$model2c$R$s[2,2]
  m2c_r1[c] <- fit$model2c$R$l[1,2]/sqrt(fit$model2c$R$l[1,1]*fit$model2c$R$l[2,2])
  m2c_r2[c] <- NA
  model3_AIC <- as.numeric(as.character(fit$model3[6]))
  m3_rate_11[c] <- fit$model3$R$l[1,1]
  m3_rate_12[c] <- NA
  m3_rate_21[c] <- fit$model3$R$l[2,2]
  m3_rate_22[c] <- NA
  m3_r1[c] <- fit$model3$R$l[1,2]/sqrt(fit$model3$R$l[1,1]*fit$model3$R$l[2,2])
  m3_r2[c] <- fit$model3$R$s[1,2]/sqrt(fit$model3$R$s[1,1]*fit$model3$R$s[2,2])
  model3b_AIC <- as.numeric(as.character(fit$model3b[6]))
  m3b_rate_11[c] <- fit$model3b$R$l[1,1]
  m3b_rate_12[c] <- fit$model3b$R$s[1,1]
  m3b_rate_21[c] <- fit$model3b$R$l[2,2]
  m3b_rate_22[c] <- NA
  m3b_r1[c] <- fit$model3b$R$l[1,2]/sqrt(fit$model3b$R$l[1,1]*fit$model3b$R$l[2,2])
  m3b_r2[c] <- fit$model3b$R$s[1,2]/sqrt(fit$model3b$R$s[1,1]*fit$model3b$R$s[2,2])
  model3c_AIC <- as.numeric(as.character(fit$model3c[6]))
  m3c_rate_11[c] <- fit$model3c$R$l[1,1]
  m3c_rate_12[c] <- NA
  m3c_rate_21[c] <- fit$model3c$R$l[2,2]
  m3c_rate_22[c] <- fit$model3c$R$s[2,2]
  m3c_r1[c] <- fit$model3c$R$l[1,2]/sqrt(fit$model3c$R$l[1,1]*fit$model3c$R$l[2,2])
  m3c_r2[c] <- fit$model3c$R$s[1,2]/sqrt(fit$model3c$R$s[1,1]*fit$model3c$R$s[2,2])
  model4_AIC <- as.numeric(as.character(fit$model4[6]))
  m4_rate_11[c] <- fit$model4$R$l[1,1]
  m4_rate_12[c] <- fit$model4$R$s[1,1]
  m4_rate_21[c] <- fit$model4$R$l[2,2]
  m4_rate_22[c] <- fit$model4$R$s[2,2]
  m4_r1[c] <- fit$model4$R$l[1,2]/sqrt(fit$model4$R$l[1,1]*fit$model4$R$l[2,2])
  m4_r2[c] <- fit$model4$R$s[1,2]/sqrt(fit$model4$R$s[1,1]*fit$model4$R$s[2,2])
  Y <- c(model1_AIC,model2_AIC,model2b_AIC,model2c_AIC,model3_AIC,model3b_AIC,model3c_AIC,model4_AIC)
  #Y <- c(model1_AIC,model2b_AIC,model3_AIC,model3b_AIC)
  W <- aic.w(Y)
  w1[c] <- W[1]
  w2[c] <- W[2]
  w2b[c] <- W[3]
  w2c[c] <- W[4]
  w3[c] <- W[5]
  w3b[c] <- W[6]
  w3c[c] <- W[7]
  w4[c] <- W[8]
  c <- c+1
}

mean(w1)
mean(w2)
mean(w2b)
mean(w2c)
mean(w3)
mean(w3b)
mean(w3c)
mean(w4)

models_output <- cbind(w1,m1_rate_11,m1_rate_12,m1_rate_21,m1_rate_22,m1_r1,m1_r2,
                       w2,m2_rate_11,m2_rate_12,m2_rate_21,m2_rate_22,m2_r1,m2_r2,
                       w2b,m2b_rate_11,m2b_rate_12,m2b_rate_21,m2b_rate_22,m2b_r1,m2b_r2,
                       w2c,m2c_rate_11,m2c_rate_12,m2c_rate_21,m2c_rate_22,m2c_r1,m2c_r2,
                       w3,m3_rate_11,m3_rate_12,m3_rate_21,m3_rate_22,m3_r1,m3_r2,
                       w3b,m3b_rate_11,m3b_rate_12,m3b_rate_21,m3b_rate_22,m3b_r1,m3b_r2,
                       w3c,m3c_rate_11,m3c_rate_12,m3c_rate_21,m3c_rate_22,m3c_r1,m3c_r2,
                       w4,m4_rate_11,m4_rate_12,m4_rate_21,m4_rate_22,m4_r1,m4_r2)

write.csv(models_output,"evol.vcv_iter.csv")


##### evol.vcv (Figure 7) #####

models_output <- read.csv("~/Desktop/birds of paradise/evol.vcv_iterations/evol.vcv_iter.csv")
models_output <- as.data.frame(models_output)

cols <- c("#66c2a5")
hist(models_output$w2,main="",xlim=c(0,1), ylim=c(0,15), breaks = 30, col=alpha(cols,0.6))
cols <- c("dodgerblue")
hist(models_output$w4,main="",breaks = 30, col=alpha(cols,0.6), add=T)
hist(models_output$w1,main="",breaks = 2, col="#b3b3b3", add=T)
hist(models_output$w2b,main="",breaks = 2, col="#b3b3b3", add=T)
hist(models_output$w2c,main="",breaks = 5, col="#b3b3b3", add=T)
hist(models_output$w3,main="",breaks = 3, col="#b3b3b3", add=T)
hist(models_output$w3b,main="",breaks = 2, col="#b3b3b3", add=T)
hist(models_output$w3c,main="",breaks = 2, col="#b3b3b3", add=T)

# 7B
cols <- c("#ffca33")
hist(models_output$m2_rate_11,main="",xlim=c(0,0.05), ylim=c(0,40), breaks = 40, col=alpha(cols,0.6))
cols <- c("#b3b3b3")
hist(models_output$m2_rate_12, col=alpha(cols,0.05), add=T, breaks = 1)
cols <- c("#e78ac3")
hist(models_output$m2_rate_21,main="",breaks = 10, col=alpha(cols,0.6), add=T)
cols <- c("#b3b3b3")
hist(models_output$m2_rate_22,main="", breaks = 1, col=alpha(cols,0.6), add=T)

#7c
cols <- c("#ffca33")
hist(models_output$m4_rate_11,main="",xlim=c(0,0.05), ylim=c(0,40), breaks = 20, col=alpha(cols,0.6))
cols <- c("#b3b3b3")
hist(models_output$m4_rate_12,main="",breaks = 1, col=alpha(cols,0.6), add=T)
cols <- c("#e78ac3")
hist(models_output$m4_rate_21,main="",breaks = 5, col=alpha(cols,0.6), add=T)
cols <- c("#b3b3b3")
hist(models_output$m4_rate_22,main="",breaks = 1, col=alpha(cols,0.6), add=T)

#7D
cols <- c("#b3b3b3")
hist(models_output$m2_r1,main="", xlim=c(0.7,1), ylim=c(0,80), breaks = 10, col=alpha(cols,0.6))

#7E
cols <- c("black")
hist(models_output$m4_r1,main="", xlim=c(0.7,1), ylim=c(0,80), breaks = 20, col=alpha(cols,0.6))
cols <- c("white")
hist(models_output$m4_r2,main="", breaks = 5, col=alpha(cols,0.6),add=T)

