
## ---- ORG_DownloadData
#download organoid data

## Velasco2019 study

## 3 month organoids
#PGP1 Line, 3 month, Batch 1 (of 2)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_PGP1.3mon.txt.gz",
              "data/ORG/Velasco2019_PGP1-3mo-Batch1.txt.gz")
gunzip("data/ORG/Velasco2019_PGP1-3mo-Batch1.txt.gz")
#PGP1 line, 3 month, Batch 2 (of 2)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_PGP1.3mon.batch2.txt.gz",
              "data/ORG/Velasco2019_PGP1-3mo-Batch2.txt.gz")
gunzip("data/ORG/Velasco2019_PGP1-3mo-Batch2.txt.gz")
#HUES66 line, 3 month, Batch 1 (of 1)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_HUES66.3mon.txt.gz",
              "data/ORG/Velasco2019_HUES66-3mo-Batch1.txt.gz")
gunzip("data/ORG/Velasco2019_HUES66-3mo-Batch1.txt.gz")

## 6 month organoids
#PGP1 line, 6 month, Batch 1 (of 3)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_PGP1.6mon.txt.gz",
              "data/ORG/Velasco2019_PGP1-6mo-Batch1.txt.gz")
gunzip("data/ORG/Velasco2019_PGP1-6mo-Batch1.txt.gz")
#PGP1 line, 6 month, Batch 2 (of 3)
##this looks missing
#PGP1 line, 6 month, Batch 3 (of 3)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_PGP1.6mon.batch3.txt.gz",
              "data/ORG/Velasco2019_PGP1-6mo-Batch3.txt.gz")
gunzip("data/ORG/Velasco2019_PGP1-6mo-Batch3.txt.gz")
#GM line, 6 month, Batch 1 (of 1)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_GM.6mon.txt.gz",
         "data/ORG/Velasco2019_GM-6mo-Batch1.txt.gz")
gunzip("data/ORG/Velasco2019_GM-6mo-Batch1.txt.gz")
#11a, 6 month, Batch 1 (of 1)
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129519/suppl/GSE129519_expression_11a.6mon.txt.gz",
         "data/ORG/Velasco2019_11a-6mo-Batch1.txt.gz")
gunzip("data/ORG/Velasco2019_11a-6mo-Batch1.txt.gz")


## Giandomenico 2019
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124174/suppl/GSE124174_org75_seuratdata.txt.gz",
              "data/ORG/Giandomenico2019_75d.txt.gz")
gunzip("data/ORG/Giandomenico2019_H9-75d-Batch1.txt.gz")
