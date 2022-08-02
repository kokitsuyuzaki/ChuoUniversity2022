# Package Download
list.of.packages <- c("devtools", "gplots")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages,
        repos="https://cloud.r-project.org/", type="source")
}
devtools::install_github("vqv/ggbiplot")

# Package Loading
library("gplots")
library("ggbiplot")

# Data Preprocessing
nests <- rbind(
	c(5, 36, 7.6, 0, 0, 0), # DIF071
	c(5.6, 78, 0, 6.7, 57, 10.2), # DIF077
	c(5.6, 74, 2.2, 9.355, 0, 0), # DIF070
	c(23.2, 51, 0, 4.55, 105, 15.1), # DIF065
	c(4.8, 42, 5.5, 9.5, 120, 28.4), # DIF081
	c(7.62, 90, 0, 0, 17, 18.1), # DIF013
	c(17.78, 14, 43, 8.55, 95, 0), # DIF012
	c(1.27, 90, 7, 0, 55, 9.6), # DIF017
	c(11.43, 58, 0, 6.85, 79, 17),# DIF039
	c(15.24, 22, 0, 5.9, 195, 17), # DIF038
	c(15.24, 89, 0, 7.2, 80, 12), # DIF024
	c(10.16, 77, 0, 0, 29, 10)) # DIF028

rownames(nests) <- c(
	"DIF071_pre", "DIF077_pre", "DIF070_pre", "DIF065_pre", # Pre-monsoon
	"DIF081_mon", "DIF013_mon", "DIF012_mon", "DIF017_mon", # Monsoon
	"DIF039_pos", "DIF038_pos", "DIF024_pos", "DIF028_pos") # Post-monsoon

colnames(nests) <- c("Te", "Angle", "Ts", "ChamberDepth", "ChamberVolume", "Depth")

label <- c(rep("Pre-monsoon", 4), rep("Monsoon", 4), rep("Post-monsoon", 4))

# Baloon plot
balloonplot(t(as.table(nests)), main ="Ant Nests", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)

# PCA w/o scaling
res1 <- prcomp(nests, center=TRUE, scale=FALSE)

ggbiplot(pcobj = res1, choices = 1:2, obs.scale = 1, var.scale = 1,
         groups = label, ellipse = TRUE, circle = TRUE) +
  scale_colour_manual(values = c("#D95F02", "#7570B3", "#1B9E77")) +
  theme(legend.direction = "horizontal", legend.position = "top")

# PCA w/ scaling
res2 <- prcomp(nests, center=TRUE, scale=TRUE)

ggbiplot(pcobj = res2, choices = 1:2, obs.scale = 1, var.scale = 1,
         groups = label, ellipse = TRUE, circle = TRUE) +
  scale_colour_manual(values = c("#D95F02", "#7570B3", "#1B9E77")) +
  theme(legend.direction = "horizontal", legend.position = "top")

# PCA w/o scaling + log transformation
res3 <- prcomp(log10(nests+1), center=TRUE, scale=FALSE)

ggbiplot(pcobj = res1, choices = 1:2, obs.scale = 1, var.scale = 1,
         groups = label, ellipse = TRUE, circle = TRUE) +
  scale_colour_manual(values = c("#D95F02", "#7570B3", "#1B9E77")) +
  theme(legend.direction = "horizontal", legend.position = "top")

# PCA w/ scaling + log transformation
res4 <- prcomp(log10(nests+1), center=TRUE, scale=TRUE)

ggbiplot(pcobj = res2, choices = 1:2, obs.scale = 1, var.scale = 1,
         groups = label, ellipse = TRUE, circle = TRUE) +
  scale_colour_manual(values = c("#D95F02", "#7570B3", "#1B9E77")) +
  theme(legend.direction = "horizontal", legend.position = "top")

# Session Information
sessionInfo()
