#setup script for phsophoprot
packs = c("BiocManager", "readxl", "colorspace", "writexl", "ggplot2")
bcpacks = c("limma", "GO.db", "topGO")

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