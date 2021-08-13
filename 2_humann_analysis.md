Preperation of HUMAnN data for functional analyses
================

  - [Load in packages and path to
    folder](#load-in-packages-and-path-to-folder)
  - [Prepare and subset the pathway
    data](#prepare-and-subset-the-pathway-data)
  - [Prepare and subset the EC data](#prepare-and-subset-the-ec-data)
  - [Prepare and subset the KO data](#prepare-and-subset-the-ko-data)

# Load in packages and path to folder

``` r
library(tidyverse)
library(phyloseq)
library(vegan)
library(caret)
path = "C:/Users/Tom/Documents/10. semester/V1/data/"
```

# Prepare and subset the pathway data

Adjust the HUMAnN pathway tables:

``` r
#Pathway count
path_count <- read.csv(paste0(path, "humann3_processed_output/Joined_tables/", "Joined_pathabundance.tsv"), sep = "\t", strip.white = T, stringsAsFactors = F, row.names = 1)

#Change J-ID:
colnames(path_count) <- substr(colnames(path_count), start = 1, stop = 6)

#Pull out pathways (not divided by taxonomy)
path_count_unstratified <- path_count[grep("\\|", rownames(path_count), value=T, invert=T),]

#Transpose dataframe -> participants as rows
path_count_unstratified <- t(path_count_unstratified)
dim(path_count_unstratified) #209 472
```

    ## [1] 209 472

``` r
path_count_unstratified <- as.data.frame(path_count_unstratified)
path_count_unstratified$J_ID <- rownames(path_count_unstratified)

#Save
write.table(path_count_unstratified, file=paste0(path, "humann3_processed_output/processed_tables/", "pathway_count_unstratified.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

#Pathway relab
path_relab <- read.csv(paste0(path, "humann3_processed_output/Joined_tables/", "Joined_pathabundance_relab.tsv"), sep = "\t", strip.white = T, stringsAsFactors = F, row.names = 1)

#Change J-ID:
colnames(path_relab) <- substr(colnames(path_relab), start = 1, stop = 6)

#Pull out pathways (not divided by taxonomy)
path_relab_unstratified <- path_relab[grep("\\|", rownames(path_relab), value=T, invert=T),]

#Transpose dataframe -> participants as rows
path_relab_unstratified <- t(path_relab_unstratified)
dim(path_relab_unstratified) #209 472
```

    ## [1] 209 472

``` r
path_relab_unstratified <- as.data.frame(path_relab_unstratified)
path_relab_unstratified$J_ID <- rownames(path_relab_unstratified)

#Save
write.table(path_relab_unstratified, file=paste0(path, "humann3_processed_output/processed_tables/", "pathway_relab_unstratified.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
```

Subset the HUMAnN data based on highest mean and variance for increased
power:

``` r
#Subset based on highest mean and variance
step1.dat <- path_relab_unstratified[,-c(473)][ ,colMeans(path_relab_unstratified[,-c(473)]) >= quantile( colMeans(path_relab_unstratified[,-c(473)]), 0.5) ]
dim(step1.dat) # 209 236
```

    ## [1] 209 236

``` r
# subset to highest mean variance
step2.dat <- apply(step1.dat, 2, var)
step2.dat <- step1.dat[ ,step2.dat >= quantile( step2.dat, 0.50) ]
dim(step2.dat) # 209 118
```

    ## [1] 209 118

``` r
df = cbind(path_relab_unstratified$J_ID, step2.dat)
colnames(df)[1] <- "J_ID"
df[1:4,1:4]
```

    ##          J_ID UNMAPPED UNINTEGRATED
    ## J28338 J28338 0.231806     0.718005
    ## J28339 J28339 0.202328     0.755933
    ## J28340 J28340 0.170398     0.783120
    ## J28341 J28341 0.185062     0.768049
    ##        1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis
    ## J28338                                          0.000528869
    ## J28339                                          0.000380408
    ## J28340                                          0.000484036
    ## J28341                                          0.000402268

``` r
dim(df)#119
```

    ## [1] 209 119

``` r
#Save
write.table(df, file=paste0(path, "humann3_processed_output/processed_tables/", "pathway_relab_unstratified_mean50_var50.txt"), col.names = T, row.names = F, quote = F, sep = '\t') 

#Select the same pathways in the count dataset
count_subset <- path_count_unstratified %>% dplyr::select(colnames(df))
count_subset[1:4,1:4]
```

    ##          J_ID  UNMAPPED UNINTEGRATED
    ## J28338 J28338 1310841.3      4060250
    ## J28339 J28339  603051.2      2253102
    ## J28340 J28340  827579.5      3803422
    ## J28341 J28341 1500899.6      6229059
    ##        1CMET2-PWY: N10-formyl-tetrahydrofolate biosynthesis
    ## J28338                                             2990.707
    ## J28339                                             1133.829
    ## J28340                                             2350.841
    ## J28341                                             3262.488

``` r
dim(count_subset)#119
```

    ## [1] 209 119

``` r
#Save
write.table(count_subset, file=paste0(path, "humann3_processed_output/processed_tables/", "pathway_count_unstratified_mean50_var50.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
```

# Prepare and subset the EC data

Read in HUMAnN EC tables:

``` r
#Pathway abundance
gene_abun <- read.csv(paste0(path, "humann3_processed_output/Joined_tables/", "Joined_genefamilies_relab_EC_rename.tsv"), sep = "\t", strip.white = T, stringsAsFactors = F, row.names = 1)

#Change J-ID:
colnames(gene_abun) <- substr(colnames(gene_abun), start = 1, stop = 6)

#Pull out pathways (not divided by taxonomy)
gene_abun_unstratified <- gene_abun[grep("\\|", rownames(gene_abun), value=T, invert=T),]

#Transpose dataframe -> participants as rows
gene_abun_unstratified <- t(gene_abun_unstratified)
dim(gene_abun_unstratified) #209 2427
```

    ## [1]  209 2427

``` r
gene_abun_unstratified <- as.data.frame(gene_abun_unstratified)
gene_abun_unstratified$J_ID <- rownames(gene_abun_unstratified)

#Save
write.table(gene_abun_unstratified, file=paste0(path, "humann3_processed_output/processed_tables/", "EC_relab_unstratified.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
```

Subset EC dataframe based on highest mean and variance:

``` r
only_path <- gene_abun_unstratified %>% dplyr::select(c(3:2427))
df <- cbind(gene_abun_unstratified$J_ID,gene_abun_unstratified$UNMAPPED,gene_abun_unstratified$UNGROUPED, only_path)

#Subset based on highest mean and variance
# subset to highest mean abundance
step1.dat <- df[,-c(1)][ ,colMeans(df[,-c(1)]) >= quantile( colMeans(df[,-c(1)]), 0.75) ]
dim(step1.dat)# 209 607
```

    ## [1] 209 607

``` r
# subset to highest mean variance
step2.dat <- apply(step1.dat, 2, var)
step2.dat <- step1.dat[ ,step2.dat >= quantile( step2.dat, 0.75) ]
dim(step2.dat) # 209 152
```

    ## [1] 209 152

``` r
df = cbind(df[,1], step2.dat)
colnames(df)[1] <- "J_ID"
colnames(df)[2] <- "UNMAPPED"
colnames(df)[3] <- "UNINTEGRATED"
df[1:4,1:4]
```

    ##          J_ID UNMAPPED UNINTEGRATED 1.1.1.1: Alcohol dehydrogenase
    ## J28338 J28338 0.231806    0.6797106                    3.31807e-05
    ## J28339 J28339 0.202328    0.7164588                    3.98004e-05
    ## J28340 J28340 0.170398    0.7460171                    4.85055e-05
    ## J28341 J28341 0.185062    0.7291379                    6.17249e-05

``` r
#Save
write.table(df, file=paste0(path, "humann3_processed_output/processed_tables/", "EC_relab_unstratified_mean75_var75.txt"), col.names = T, row.names = F, quote = F, sep = '\t')  
```

# Prepare and subset the KO data

Read in HUMAnN KO tables:

``` r
gene_abun <- read.csv(paste0(path, "humann3_processed_output/Joined_tables/", "Joined_genefamilies_relab_KO_rename.tsv"), sep = "\t", strip.white = T, stringsAsFactors = F, row.names = 1)

#Change J-ID:
colnames(gene_abun) <- substr(colnames(gene_abun), start = 1, stop = 6)

#Pull out pathways (not divided by taxonomy)
gene_abun_unstratified <- gene_abun[grep("\\|", rownames(gene_abun), value=T, invert=T),]

#Transpose dataframe -> participants as rows
gene_abun_unstratified <- t(gene_abun_unstratified)
dim(gene_abun_unstratified) #209 6177
```

    ## [1]  209 6177

``` r
gene_abun_unstratified <- as.data.frame(gene_abun_unstratified)
gene_abun_unstratified$J_ID <- rownames(gene_abun_unstratified)

#Save
write.table(gene_abun_unstratified, file=paste0(path, "humann3_processed_output/processed_tables/", "KO_relab_unstratified.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
```

Subset KO dataframe based on highest mean and variance:

``` r
only_path <- gene_abun_unstratified %>% dplyr::select(c(3:6177))
df <- cbind(gene_abun_unstratified$J_ID,gene_abun_unstratified$UNMAPPED,gene_abun_unstratified$UNGROUPED, only_path)

# subset to highest mean abundance
step1.dat <- df[,-c(1)][ ,colMeans(df[,-c(1)]) >= quantile( colMeans(df[,-c(1)]), 0.75) ]
dim(step1.dat) #209 1545
```

    ## [1]  209 1545

``` r
# subset to highest mean variance
step2.dat <- apply(step1.dat, 2, var)
step2.dat <- step1.dat[ ,step2.dat >= quantile( step2.dat, 0.75) ]
dim(step2.dat) # 209 387
```

    ## [1] 209 387

``` r
df = cbind(df[,1], step2.dat)
colnames(df)[1] <- "J_ID"
colnames(df)[2] <- "UNMAPPED"
colnames(df)[3] <- "UNINTEGRATED"
df[1:4,1:4]
```

    ##          J_ID UNMAPPED UNINTEGRATED
    ## J28338 J28338 0.231806    0.7333464
    ## J28339 J28339 0.202328    0.7417751
    ## J28340 J28340 0.170398    0.7657038
    ## J28341 J28341 0.185062    0.7615598
    ##        K00014: shikimate dehydrogenase [EC:1.1.1.25]
    ## J28338                                   5.61485e-05
    ## J28339                                   6.91865e-05
    ## J28340                                   9.16378e-05
    ## J28341                                   6.77431e-05

``` r
#Save
write.table(df, file=paste0(path, "humann3_processed_output/processed_tables/", "KO_relab_unstratified_mean75_var75.txt"), col.names = T, row.names = F, quote = F, sep = '\t') 
```
