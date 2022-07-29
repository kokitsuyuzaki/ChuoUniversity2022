# Package Download
list.of.packages <- c("rTensor", "nnTensor", "iTensor",
    "einsum", "RColorBrewer", "TeachingDemos")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages,
        repos="https://cloud.r-project.org/", type="source")
}

# Package Loading
library("rTensor")
library("nnTensor")
library("iTensor")
library("einsum")
library("RColorBrewer")
library("TeachingDemos")

# Function Definition
plot_coef <- function(coef_matrix, colvec){
    options(repr.plot.width=12, repr.plot.height=6)
    layout(rbind(1:5, 6:10))
    for(i in seq(10)){
        tmp <- as.matrix(coef_matrix[,i])
        dim(tmp) <- c(92, 112)
        image(tmp, col=colvec, main=i)
    }
}

plot_score <- function(score_matrix, facedata, colvec){
    options(repr.plot.width=12, repr.plot.height=12)
    plot(score_matrix[,1:2], pch=16, cex=2,
        axes=FALSE, xlab="", ylab="")
    abline(v=0)
    abline(h=0)
    for(i in seq(10)){
        subplot(image(facedata@data[,,i,1], col=colvec,
            main=i, axes=FALSE),
            x=score_matrix[i,1], y=score_matrix[i,2],
            size=c(0.45, 0.45))
    }
}

# Setting
options(repr.plot.width=12, repr.plot.height=6)
colvec <- rev(brewer.pal(9, "Greys"))

# Data Download
facedata <- load_orl()
str(facedata)

# Plot data (1st angle)
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    image(facedata@data[,,i,1], col=colvec, main=i)
}

# Plot data (2nd angle)
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    image(facedata@data[,,i,2], col=colvec, main=i)
}

# Use 1st Slice (matrix)
matdata <- cs_unfold(facedata[,,1:10,1], m=3)@data

# Average Face
options(repr.plot.width=6, repr.plot.height=6)
average_face <- rowMeans(matdata)
dim(average_face) <- c(92, 112)
image(average_face, col=colvec, main="Average Face")

# PCA
res_pca <- prcomp(t(matdata), center=TRUE, scale=FALSE)

## Plot PCA Coefficient
plot_coef(res_pca$rotation, colvec)

## Plot PCA Score
plot_score(res_pca$x, facedata, colvec)

# SVD
res_svd <- svd(matdata)

## Plot SVD Coefficient
plot_coef(res_svd$u, colvec)

## Plot SVD Score
plot_score(res_svd$v, facedata, colvec)

# NMF
res_nmf <- NMF(matdata, J=10)

## Plot NMF Coefficient
plot_coef(res_nmf$U, colvec)

## Plot NMF Score
plot_score(res_nmf$V, facedata, colvec)

# ICA
res_ica <- ICA(matdata, J=10)

## Plot ICA Coefficient
plot_coef(res_ica$S, colvec)

## Plot ICA Score
plot_score(res_ica$A, facedata, colvec)

# CP (Tensor Decomposition)
norm_facetensor <- facedata[,,1:10,1]
for(i in seq(10)){
    average_facedata <- facedata[,,i,1]@data
    average_facedata[] <- mean(facedata@data[,,i,1])
    norm_facetensor[,,i] <- facedata[,,i,1]@data - average_facedata
}
res_cp <- cp(norm_facetensor, 10)

## Plot CP Coefficient
options(repr.plot.width=12, repr.plot.height=6)
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    coef_cp <- outer(res_cp$U[[1]][, i], res_cp$U[[2]][, i])
    image(coef_cp, col=colvec, main=i)
}

## Plot CP Score
plot_score(res_cp$U[[3]], facedata, colvec)

## Plot Reconstructed Tensor (CP)
rec_cp <- einsum('il,jl,kl->ijk',
    res_cp$U[[1]], res_cp$U[[2]], res_cp$U[[3]])
options(repr.plot.width=12, repr.plot.height=6)
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    image(rec_cp[,,i], col=colvec)
}

# Tucker/HOOI (Tensor Decomposition)
res_tucker <- tucker(norm_facetensor, c(10,10,10))

## Plot Tucker/HOOI Coefficient
layout(rbind(1:5, 6:10))
for(i in seq(3)){
    for(j in seq(3)){
        coef_tucker <- outer(res_tucker$U[[1]][, i], res_tucker$U[[2]][, j])
        image(coef_tucker, col=colvec, main=paste0(i, " x ", j))
    }
}

## Plot Tucker/HOOI Score
plot_score(res_tucker$U[[3]], facedata, colvec)

## Plot Reconstructed Tensor (Tucker)
rec_tucker <- recTensor(res_tucker$Z, res_tucker$U, reverse=TRUE)@data
options(repr.plot.width=12, repr.plot.height=6)
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    image(rec_tucker[,,i], col=colvec)
}

# Session Information
sessionInfo()