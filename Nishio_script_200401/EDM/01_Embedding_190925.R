
# Create output folder
dir.create("01_Embedding_Out")

# Load library
library(rEDM); packageVersion("rEDM") # v0.6.9

# Load data
temperature <- read.csv("data/2012-2014temp.csv",header=T,sep=",")
H3K27me3 <- read.csv("data/H3K27me3_Omoide-river_all.csv",header=T,sep=",")
H3K4me3 <- read.csv("data/H3K4me3_Omoide-river_all.csv",header=T,sep=",")
RNA <- read.csv("data/RNA_Omoide-river.csv",header=T,sep=",")

# Compile
samplepoint <- c(269,297,311,326,339,353,374,388,402,416,430,444,458,472,486,499,514,528,549,563,577,591,605,619,633,647,661,675,689,703,717,738,752,766,780,794,810,822,836,849,864,878,892,905,920,934,948,962,976,990)
temperature2 <- temperature[samplepoint,]
sp <- smooth.spline(1:50, temperature2$temp, spar=0.5)
x <- 1:50
pred.temp <- predict(sp, x)$y
plot(temperature2$temp,type="o")
lines(pred.temp,col="red")

d <- cbind(pred.temp, 
       H3K27me3$K27.I.STM, H3K27me3$K27.II.STM, H3K27me3$K27.III.STM, H3K27me3$K27.IV.STM,
		   H3K27me3$K27.V.STM, H3K27me3$K27.VI.STM, H3K27me3$K27.VII.STM, H3K27me3$K27.VIII.STM,
		   H3K4me3$K4.I.ACT2, H3K4me3$K4.II.ACT2, H3K4me3$K4.VI.ACT2, H3K4me3$K4.VII.ACT2,
		   H3K4me3$K4.VIII.ACT2, RNA$FLC.ACT2.log10mean)
colnames(d) <- c("Temp","K27_I", "K27_II", "K27_III", "K27_IV", "K27_V", "K27_VI", "K27_VII", "K27_VIII", 
				 "K4_I", "K4_II", "K4_VI", "K4_VII", "K4_VIII", "RNA")
d <- as.data.frame(apply(d, 2, function(x)((x-mean(x))/sd(x))))

# Visualization (brief check)
par(las = 1) # Change y-axis label
plot(d[,"RNA"], type = "o")

# Load function
source("functions/S01_HelperFunctions.R")

# Simplex projection
lib_type_used <- "full"
Etemp <- TestE(d$Temp, d_name="Temp", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)

Ek27I <- TestE(d$K27_I, d_name="K27_I", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27II <- TestE(d$K27_II, d_name="K27_II", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27III <- TestE(d$K27_III, d_name="K27_III", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27IV <- TestE(d$K27_IV, d_name="K27_IV", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27V <- TestE(d$K27_V, d_name="K27_V", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27VI <- TestE(d$K27_VI, d_name="K27_VI", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27VII <- TestE(d$K27_VII, d_name="K27_VII", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek27VIII <- TestE(d$K27_VIII, d_name="K27_VIII", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)

Ek4I <- TestE(d$K4_I, d_name="K4_I", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek4II <- TestE(d$K4_II, d_name="K4_II", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek4VI <- TestE(d$K4_VI, d_name="K4_VI", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek4VII <- TestE(d$K4_VII, d_name="K4_VII", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)
Ek4VIII <- TestE(d$K4_VIII, d_name="K4_VIII", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)

Erna <- TestE(d$RNA, d_name="RNA", E_range = 1:24, lib_type = lib_type_used, E_only = T, visualize=T)




# Test of embedding dimention (K27_II)
pdf("01_Embedding_Out/Etest_K27_II.pdf",width=3.6,height=1.1)
par(mar=c(1.5,1.5,1,5))
par(mgp=c(0.5,0.1,0))
par(ps=7)
par(xpd=T)
test <- simplex(d$K27_II, c(1, length(d$K27_II)), c(1, length(d$K27_II)), num_neighbors = "e+1", E = 1:24, silent = T, stats_only = F)
plot(test$model_output[[1]]$time, test$model_output[[1]]$obs, 
     type="p",cex=0.5,tcl=0.1,xlim=c(0,50),ylim=c(-1,2.2),axes=F,ann=F)
lines(test$model_output[[2]]$time, test$model_output[[2]]$pred, col="black")
lines(test$model_output[[8]]$time, test$model_output[[8]]$pred, col="red")
lines(test$model_output[[17]]$time, test$model_output[[17]]$pred, col="blue")
axis(side=1,at=seq(0,50,10),label=F,tcl=-0.1)
mtext(side=1,at=seq(0,50,10),seq(0,50,10),line=-0.25)
mtext(side=1,"Time points",line=0.25)  
axis(side=2,at=seq(-1,2,1),label=F,tcl=-0.1)
mtext(side=2,at=seq(-1,2,1),seq(-1,2,1),line=0.15,las=1)
mtext(side=2,"K27_II",line=0.5)
box()
legend(par()$usr[2]+0.2, par()$usr[4]+0.7,
       legend=c("Observed",
                expression(paste("Predicted (", italic(E),"=2)")),
                expression(paste("Predicted (", italic(E),"=8)")),
                expression(paste("Predicted (", italic(E),"=17)"))),
       pch=c(1,NA,NA,NA),lty=c(0,1,1,1),pt.cex=0.5,bty="n",x.intersp=0.1,
       y.intersp=0.8,seg.len=0.7,col=c("black","black","red","blue"))
dev.off()

Ek27II <- 2



# Test of embedding dimention (K4_II)
pdf("01_Embedding_Out/Etest_K4_II.pdf",width=3.6,height=1.1)
par(mar=c(1.5,1.5,1,5))
par(mgp=c(0.5,0.1,0))
par(ps=7)
par(xpd=T)
test <- simplex(d$K4_II, c(1, length(d$K4_II)), c(1, length(d$K4_II)), num_neighbors = "e+1", E = 1:24, silent = T, stats_only = F)
plot(test$model_output[[1]]$time, test$model_output[[1]]$obs, 
     type="p",cex=0.5,tcl=0.1,xlim=c(0,50),ylim=c(-2,1.5),axes=F,ann=F)
lines(test$model_output[[4]]$time, test$model_output[[4]]$pred, col="black")
lines(test$model_output[[18]]$time, test$model_output[[18]]$pred, col="red")
axis(side=1,at=seq(0,50,10),label=F,tcl=-0.1)
mtext(side=1,at=seq(0,50,10),seq(0,50,10),line=-0.25)
mtext(side=1,"Time points",line=0.25)  
axis(side=2,at=seq(-2,1,1),label=F,tcl=-0.1)
mtext(side=2,at=seq(-2,1,1),seq(-2,1,1),line=0.15,las=1)
mtext(side=2,"K4_II",line=0.5)
box()
legend(par()$usr[2]+0.2, par()$usr[4]+0.7,
       legend=c("Observed",
                expression(paste("Predicted (", italic(E),"=4)")),
                expression(paste("Predicted (", italic(E),"=18)"))),
       pch=c(1,NA,NA),lty=c(0,1,1),pt.cex=0.5,bty="n",x.intersp=0.1,
       y.intersp=0.8,seg.len=0.7,col=c("black","black","red"))
dev.off()

Ek4II <- 4



# Test of embedding dimention (RNA)
pdf("01_Embedding_Out/Etest_RNA.pdf",width=3.6,height=1.1)
par(mar=c(1.5,1.5,1,5))
par(mgp=c(0.5,0.1,0))
par(ps=7)
par(xpd=T)
test <- simplex(d$RNA, c(1, length(d$RNA)), c(1, length(d$RNA)), num_neighbors = "e+1", E = 1:24, silent = T, stats_only = F)
plot(test$model_output[[1]]$time, test$model_output[[1]]$obs, 
     type="p",cex=0.5,tcl=0.1,xlim=c(0,50),ylim=c(-3,1.5),axes=F,ann=F)
lines(test$model_output[[4]]$time, test$model_output[[4]]$pred, col="black")
lines(test$model_output[[19]]$time, test$model_output[[19]]$pred, col="red")

axis(side=1,at=seq(0,50,10),label=F,tcl=-0.1)
mtext(side=1,at=seq(0,50,10),seq(0,50,10),line=-0.25)
mtext(side=1,"Time points",line=0.25)  
axis(side=2,at=seq(-3,1,1),label=F,tcl=-0.1)
mtext(side=2,at=seq(-3,1,1),seq(-3,1,1),line=0.15,las=1)
mtext(side=2,"RNA",line=0.7)
box()
legend(par()$usr[2]+0.2, par()$usr[4]+0.7,
       legend=c("Observed",
                expression(paste("Predicted (", italic(E),"=4)")),
                expression(paste("Predicted (", italic(E),"=19)"))),
       pch=c(1,NA,NA),lty=c(0,1,1),pt.cex=0.5,bty="n",x.intersp=0.1,
       y.intersp=0.8,seg.len=0.7,col=c("black","black","red"))
dev.off()

Erna <- 4



# Save output
save.image("01_Embedding_Out/01_Embedding_Out_temp190925.RData")

