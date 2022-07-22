# Package Download
list.of.packages <- c("rTensor", "nnTensor", "iTensor",
    "einsum", "RColorBrewer", "fields")
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
library("fields")

# Data Download
facedata <- load_orl()
str(facedata)

# Plot data (1st angle)
colvec <- brewer.pal(9, "Blues")
layout(rbind(1:3, 4:6))
for(i in seq(6)){
    image(facedata@data[,,i,1], col=colvec, main=i)
}

# Plot data (2nd angle)
layout(rbind(1:3, 4:6))
for(i in seq(6)){
    image(facedata@data[,,i,2], col=colvec, main=i)
}

# Use 1st Slice (matrix)
matdata <- cs_unfold(facedata[,,,1], m=3)@data

# Average Face
average_face <- rowMeans(matdata)
dim(average_face) <- c(92, 112)
image(average_face, col=colvec, main="Average Face")

# PCA
res_pca <- prcomp(t(matdata), center=TRUE, scale=FALSE)

layout(rbind(1:4, 5:8))
for(i in seq(8)){
    score_pca <- as.matrix(res_pca$rotation[,i])
    dim(score_pca) <- c(92, 112)
    image(score_pca, col=colvec, main=i)
}

# SVD
res_svd <- svd(matdata)

layout(rbind(1:4, 5:8))
for(i in seq(8)){
    score_svd <- as.matrix(res_svd$u[,i])
    dim(score_svd) <- c(92, 112)
    image(score_svd, col=colvec, main=i)
}

# NMF
res_nmf <- NMF(matdata,  J=8)

layout(rbind(1:4, 5:8))
for(i in seq(8)){
    score_nmf <- res_nmf$U[,i]
    dim(score_nmf) <- c(92, 112)
    image(score_nmf, col=colvec, main=i)
}

# ICA
res_ica <- ICA(matdata, J=8)

layout(rbind(1:4, 5:8))
for(i in seq(8)){
    score_ica <- res_ica$S[, i]
    dim(score_ica) <- c(92, 112)
    image(score_ica, col=colvec, main=i)
}

# CP (Tensor Decomposition)
norm_facetensor <- facedata[,,,1]
for(i in seq(dim(facedata[,,,1])[3])){
    average_facedata <- facedata[,,i,1]@data
    average_facedata[] <- mean(facedata@data[,,i,1])
    norm_facetensor[,,i] <- facedata[,,i,1]@data - average_facedata
}
res_cp <- cp(norm_facetensor, 20)

layout(rbind(1:5, 6:10, 11:15, 16:20))
for(i in seq(20)){
    score_cp <- outer(res_cp$U[[1]][, i], res_cp$U[[2]][, i])
    image(score_cp, col=colvec, main=i)
}

## Reconstructed Tensor (CP)
rec_cp <- einsum('il,jl,kl->ijk',
    res_cp$U[[1]], res_cp$U[[2]], res_cp$U[[3]])
layout(rbind(1:5, 6:10, 11:15, 16:20))
for(i in seq(20)){
    image(rec_cp[,,i], col=colvec)
}

# Tucker/HOOI (Tensor Decomposition)
res_tucker <- tucker(facedata[,,,1], c(3,3,10))

layout(rbind(1:5, 6:10))
for(i in seq(3)){
    for(j in seq(3)){
        score_tucker <- outer(res_tucker$U[[1]][, i], res_tucker$U[[2]][, j])
        image(score_tucker, col=colvec, main=paste0(i, " x ", j))
    }
}

## Reconstructed Tensor (Tucker)
rec_tucker <- recTensor(res_tucker$Z, res_tucker$U, reverse=TRUE)@data
layout(rbind(1:5, 6:10, 11:15, 16:20))
for(i in seq(20)){
    image(rec_tucker[,,i], col=colvec)
}

# Session Information
sessionInfo()