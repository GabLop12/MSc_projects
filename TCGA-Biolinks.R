# R script to download data from TCGA using TCGAbiolinks R package

setwd("A:/My-workshops/cancer-databases/TCGAbiolinks/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

##version.string R version 4.2.2 (2021-11-01)
packageVersion("TCGAbiolinks")
##[1] '2.26.0'
packageVersion('SummarizedExperiment')
#[1] ‘1.28.0’

##error 'TCGAbiolinks' in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]): namespace 'rlang' 1.0.6 is already loaded, but >= 1.1.0 is required
##packageVersion("rlang")
##remove.packages("rlang")
###install.packages("rlang")
##packageVersion("rlang")

#BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

# get a list of projects
gdcprojects <- getGDCprojects()
TCGAbiolinks:::getProjectSummary('TCGA-BRCA')

###DOWNLOADING OF GENE EXPRESSION DATA

# building a broad query
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)


# build a query to download TCGA-BRCA gene expression data for specific samples------------
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c('TCGA-B6-A0RI-01A-11R-A057-13', 'TCGA-BH-A1FU-11A-23R-A14D-07','TCGA-AO-A03U-01B-21R-A10J-07'))

getResults(query_TCGA)

# build a query to download TCGA-BRCA gene expression data for all BRCA samples------------
query_TCGA_all_samples <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')

# download data - GDCdownload
GDCdownload(query_TCGA)

# prepare data
tcga_brca_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'fpkm_unstrand')

###Downloading of clinical data
library(TCGAbiolinks)
query_clinical <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query_clinical)
clinical.BCRtab.all <- GDCprepare(query_clinical)
names(clinical.BCRtab.all)



###DOWNLOADING OF COPY NUMBER DATA
query_copy_number <- GDCquery(
  project = "TCGA-ACC",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",              
  access = "open"
)
GDCdownload(query_copy_number)
data <- GDCprepare(query_copy_number)
copynumber_matrix <- assay(data, 'copy_number')


###DOWNLOADING OF METHYLATION DATA

# build a query to retrieve DNA methylation data --------------
query_methly <- GDCquery(project = 'TCGA-GBM',
                         data.category = 'DNA Methylation',
                         platform = 'Illumina Human Methylation 27',
                         access = 'open',
                         data.type = 'Methylation Beta Value',
                         barcode = c('TCGA-19-0962-01B-01D-0521-05', 'TCGA-06-0137-01A-01D-0218-05'))

output_query_methyl <- getResults(query_methly)

GDCdownload(query_methly)
#dna.meth <- GDCprepare(query_methly, summarizedExperiment = TRUE)
#assay(dna.meth)  


###DOWNLOADING OF MUTATION DATA

# download mutation data from TCGA ----------------------
query_mutation <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Simple Nucleotide Variation',
                           access = 'open',
                           barcode = c('TCGA-LL-A73Y-01A-11D-A33E-09','TCGA-LL-A73Y-10B-01D-A33H-09',
                                       'TCGA-E9-A1NH-01A-11D-A14G-09','TCGA-E9-A1NH-11A-33D-A14G-09'))

output_query_mutation <- getResults(query_mutation)

GDCdownload(query_mutation)

maf <- GDCprepare(query_mutation, summarizedExperiment = TRUE)


