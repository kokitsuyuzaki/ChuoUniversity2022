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

# Plot data
colvec <- brewer.pal(9, "Blues")
image(facedata@data[,,1,1], col=colvec)

# 1st Slice
matdata <- cs_unfold(facedata[,,,1], m=3)@data

# SVD
res_svd <- svd(matdata)

# Plot Component 1 - 10
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    score_svd <- as.matrix(res_svd$u[,i])
    dim(score_svd) <- c(92, 112)
    image(score_svd, col=colvec, main=i)
}

# NMF
res_nmf <- NMF(matdata,  J=10)

# Plot Component 1 - 10
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    score_nmf <- res_nmf$U[,i]
    dim(score_nmf) <- c(92, 112)
    image(score_nmf, col=colvec, main=i)
}

# ICA
res_ica <- ICA(matdata, J=10)

# Plot Component 1 - 10
layout(rbind(1:5, 6:10))
for(i in seq(10)){
    score_ica <- res_ica$S[, i]
    dim(score_ica) <- c(92, 112)
    image(score_ica, col=colvec, main=i)
}

# CP (Tensor Decomposition)
res_cp <- cp(facedata[,,,1], 20)

# Plot Component 1 - 20
layout(rbind(1:5, 6:10, 11:15, 16:20))
for(i in seq(20)){
    score_cp <- outer(res_cp$U[[1]][, i], res_cp$U[[2]][, i])
    image(score_cp, col=colvec, main=i)
}

# Reconstructed Tensor
rec_cp <- einsum('il,jl,kl->ijk',
    res_cp$U[[1]], res_cp$U[[2]], res_cp$U[[3]])
layout(rbind(1:5, 6:10, 11:15, 16:20))
for(i in seq(20)){
    image(rec_cp[,,i], col=colvec)
}

# Tucker/HOOI (Tensor Decomposition)
res_tucker <- tucker(facedata[,,,1], c(3,3,10))

# Plot Components
layout(rbind(1:5, 6:10))
for(i in seq(3)){
    for(j in seq(3)){
        score_tucker <- outer(res_tucker$U[[1]][, i], res_tucker$U[[2]][, j])
        image(score_tucker, col=colvec, main=paste0(i, " x ", j))
    }
}

# Reconstructed Tensor
rec_tucker <- recTensor(res_tucker$Z, res_tucker$U, reverse=TRUE)@data
layout(rbind(1:5, 6:10, 11:15, 16:20))
for(i in seq(20)){
    image(rec_tucker[,,i], col=colvec)
}
