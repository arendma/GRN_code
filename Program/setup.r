#setup script for phsophoprot
packs = c("BiocManager", "ggplot2", "ggfortify", "amap", "pheatmap", "parallel", "elasticnet", "GeneNet", "igraph")
bcpacks = c("edgeR", "GENIE3", "minet")

for (pack in packs) {
  if(!(requireNamespace(pack, quietly=TRUE))) {
    install.packages(pack)
  }
}

for (bcpack in bcpacks) {
  if(!(requireNamespace(bcpack, quietly=TRUE))) {
    library(BiocManager)
    BiocManager::install(bcpack)
  }
}