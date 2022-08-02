# Package Download
list.of.packages <- "nnTensor"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages,
        repos="https://cloud.r-project.org/", type="source")
}

# Package Loading
library("nnTensor")

# Toy data
data <- toyModel("NMF")

# NMF (Frobenius)
res1 <- NMF(data, J=5, algorithm="Frobenius", num.iter=30, viz=TRUE)

# NMF (KL)
res2 <- NMF(data, J=5, algorithm="KL", num.iter=30, viz=TRUE)

# NMF (IS)
res3 <- NMF(data, J=5, algorithm="IS", num.iter=30, viz=TRUE)

# Session Information
sessionInfo()