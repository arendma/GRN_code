# Code to generate a Gener regulatory netowrk for Chlamydomonas reinhardtii

The code included in the ../Program subfolder recreates the gene regulatory network based on the read count data from 

- *phot* Mutant sreen
- acetate timeline
- [Merchant 2015 diurnal multiomics study](https://academic.oup.com/plcell/article/27/10/2743/6096622?login=true)
- [Merchant 2019 diurnal multiomics study](https://www.pnas.org/content/116/6/2374.abstract)

## Steps

### 1. Setup

Install R according to the [installation guidelines](https://www.r-project.org/) and then source `../Program/setup.r`

Before running the code make sure to have set 

  options(stringsAsFactors = FALSE)

### 2. Aggregating and processing read counts - some info graphics (tested on R.3.6.1)

The script `2-CountProcessing.R` reads the raw read count data from `..\Data` and aggregates it into one large matrix with genes as rows and samples as columns. Next the normalization steps of the edgeR wirkflow using TMM and voom normalization are applied to this matrix to derive at estimates for the log2 count values.
In the second part of the script PCA plots and a heatmap of correlation values are generated to give a first overview of the data.

### 3. network from elastic net regression (tested on R4.2.2)
The script `3-Netreg_elnet_TF_G_me1519am17v2.R` zscore transformes the log values gene wise and aggregates replicates to median values if applicable. 
Then it uses a implementation fo the elastic net algorithm to generate a gene regulatory network **(Runtime on Core-i7 16GB RAM was 5 days)**

### 4. network from other GRN inference methods (tested on R4.2.2)

The script `4-Netreg_bignet_TF_G_me1519am17v2` zscore transformes the log values gene wise and aggregates replicates to median values if applicable. It then generates GRNs inferred using the approaches
- [GENIE3](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012776)
- [Graphical Gaussian Models](https://www.degruyter.com/document/doi/10.2202/1544-6115.1175/html) 
- [CLR](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050008)
- [ARACNE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7)
- [DECONVOLUTION](https://www.nature.com/articles/nbt.2635)
- [SILENCING](https://www.nature.com/articles/nbt.2601)

During later checks it turned out that ARACNE and SILENCING did not produce results of usabel quality. 
The script then aggregates the obtaine GRNs together with the output from `3-Netreg_elnet_TF_G_me1519am17v2.R`and saves them into one R object.
**(Runtime on Core-i7 16GB RAM was 12h)**


### 5. Aggregating network predictions into consensus (tested on R.3.6.1)
The script `5-Gen_consensus.R` uses the previously generated GRNs and aggregates them into a consensus network considering a maximum of 10% of all possible edges per approach. One consensus of all approaches is generated, and a second consensus omitting GRNS of ARACNE and Silencing is generated. The final consensus networks are exported as edge lists in a tab seperated text file: The first column is the node out node, second column the in node after that a column with the rank of the edge in each integrated approach is added and a final column that gives the mean of all ranks. (Ranks of edges not included in a certain approach are imputed as 10% of all possible edges +1). Similar to the process described in [Marbach et al.](https://www.nature.com/articles/nmeth.2016)



