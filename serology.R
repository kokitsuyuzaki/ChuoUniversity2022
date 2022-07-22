# Package Download
list.of.packages <- c("rTensor", "ggplot2", "reshape2", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages,
        repos="https://cloud.r-project.org/", type="source")
}

# Package Loading
library("rTensor")
library("ggplot2")
library("reshape2")
library("RColorBrewer")

# Function Definition
load_serology <- function(){
	# Download from Zenodo Server
	td <- tempdir()
    tempfile <- paste0(td, "/meyer-lab-systemsSerology-fd9ef61.zip")
    download.file("https://zenodo.org/record/5184449/files/meyer-lab/systemsSerology-v1.0.zip?download=1", tempfile)
    # Preprocessing
    unzip(tempfile, exdir=td)
    csvfile <- paste0(td, "/meyer-lab-systemsSerology-fd9ef61/syserol/data/ZoharCovData.csv")
    data <- read.csv(csvfile, row.names=1, header=TRUE)
    serology <- t(t(as.matrix(data[, 23:140])) - unlist(data[443, 23:140]))
    serology <- serology[1:438, ]
    # Data frame -> Array
    receptor <- c("IgG1", "IgG2", "IgG3", "IgA1", "IgA2", "IgM", "FcRalpha", "FcR2A", "FcR2B", "FcR3A", "FcR3B")
    antigen <- c("S", "RBD", "N", "S1", "S2", "S1.Trimer")
    arr <- array(0, dim=c(11, 6, 438))
    dimnames(arr) <- list(
        receptor = receptor,
        antigen = antigen,
        samples = rownames(serology))
    for(i in receptor){
        for(j in antigen){
            arr[i, j, ] <- serology[, paste0(i, "_", j)]
        }
    }
    # Log Transformation
    arr[which(arr < 0)] <- 10
    arr <- log10(arr)
    # Centering
    for(k in seq_len(dim(arr)[3])){
        arr[,,k] <- arr[,,k] - mean(arr[,,k])
    }
    # Array -> Tensor
    as.tensor(arr)
}

# Tensor Data
covid19 <- load_serology()

# Tucker Decomposition
J = 2
set.seed(123456)
res_tucker <- hosvd(covid19, ranks=rep(J, length=3))

# Receptor Patterns
## Preprocessing
pattern_r <- res_tucker$U[[1]]
rownames(pattern_r) <- dimnames(covid19@data)$receptor
colnames(pattern_r) <- paste("Component", seq(J))
df_r <- melt(pattern_r)
colnames(df_r) <- c("Receptor", "Component", "Value")
df_r$Receptor <- factor(df_r$Receptor, level=rev(unique(df_r$Receptor)))

## Plot Receptor Patterns
g_r <- ggplot(df_r, aes(x=Component, y=Receptor, fill = Value))
g_r <- g_r + geom_tile()
g_r <- g_r +  scale_fill_gradientn("value", colours = brewer.pal(11, "PiYG"))
g_r <- g_r + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g_r <- g_r + theme(text = element_text(size=30))
g_r

# Antigen Patterns
## Preprocessing
pattern_a <- res_tucker$U[[2]]
rownames(pattern_a) <- dimnames(covid19@data)$antigen
colnames(pattern_a) <- paste("Component", seq(J))
df_a <- melt(pattern_a)
colnames(df_a) <- c("Antigen", "Component", "Value")
df_a$Antigen <- factor(df_a$Antigen, level=rev(unique(df_a$Antigen)))

## Plot Antigen Patterns
g_r <- ggplot(df_a, aes(x=Component, y=Antigen, fill = Value))
g_r <- g_r + geom_tile()
g_r <- g_r +  scale_fill_gradientn("value", colours = brewer.pal(11, "PiYG"))
g_r <- g_r + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g_r <- g_r + theme(text = element_text(size=30))
g_r

# Sample Patterns
## Preprocessing
pattern_s <- res_tucker$U[[3]]
rownames(pattern_s) <- dimnames(covid19@data)$samples
colnames(pattern_s) <- paste("Component", seq(J))
df_s <- melt(pattern_s)
colnames(df_s) <- c("Sample", "Component", "Value")
df_s$Sample <- factor(df_s$Sample, level=rev(unique(df_s$Sample)))

## Plot Sample Patterns
g_s <- ggplot(df_s, aes(x=Component, y=Sample, fill = Value))
g_s <- g_s + geom_tile()
g_s <- g_s +  scale_fill_gradientn("value", colours = brewer.pal(11, "PiYG"))
g_s <- g_s + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g_s <- g_s + theme(axis.text.y = element_blank())
g_s <- g_s + theme(text = element_text(size=30))
g_s

# Session Information
sessionInfo()