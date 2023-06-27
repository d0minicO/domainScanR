## script with example usages of domainScanR

# load libraries
library(tidyverse)
library(magrittr)
library(httr)
library(viridis)



# load the domainScanR function from github
devtools::source_url("https://github.com/d0minicO/domainScanR/blob/main/domainScanR.R?raw=TRUE")


#############

## selection of gene names
names =
  c(
    "CUL3",
    "CUL4B",
    "CUL7",
    "DDRGK1",
    "DDX21",
    "DDX23",
    "DDX54",
    "DEPDC1B",
    "DGCR8",
    "DHX57",
    "ECT2"
    )

domainScanR(input=names)

ggsave(filename = paste0("domainScanR_example_gene_name.pdf"),
       width=8,
       height=4)

## uniprot codes
ids =
  c(
    "P19474",
    "O15164",
    "P14373",
    "Q15654",
    "P62256",
    "Q9H9Y6",
    "P30876",
    "P13727",
    "Q9Y2Y8",
    "O43900",
    "Q6VN20",
    "Q96S59",
    "Q9H871"
    )

dat = domainScanR(input=ids,
            data_type="uniprot")


# save the plot
ggsave(plot=dat[[1]],
       filename = "domainScanR_example_uniprot_IDs.pdf",
       width=8,
       height=4)

# save the table
dat[[2]] %>%
  write.table(file="domainScanR_example_uniprot_IDs_table.tsv",
              quote=F,
              sep="\t",
              row.names=F)