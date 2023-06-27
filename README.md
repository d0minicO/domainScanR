# domainScanR
## _R function for quantifying enrichment of domains in a set of proteins_


domainScanR is a custom R function that takes gene names or Uniprot IDs and outputs an enrichment plot and a table quantifying enrichment of protein domains

---
# Install
### Part 1 -- install package dependencies if you have not already got them
```sh

install.packages("tidyverse")
install.packages("magrittr")
install.packages("httr")
install.packages("viridis")

```
Check you can load all these libraries successfully. If you see no errors then you are good to go!
```sh
library(tidyverse)
library(magrittr)
library(httr)
library(viridis)
```

### Part 2 -- source domainScanR from github
run this single line to load the funtion directly from github
```sh
devtools::source_url("https://github.com/d0minicO/domainScanR/blob/main/domainScanR.R?raw=TRUE")
```

# Arguments
 - input: the list of genes or proteins to search for enrichment of domains in. A character vector of gene names OR uniprot IDs to search
 - data_type (optional): whether gene names or uniprot IDs have been provided. One of c("gene", "uniprot)
 --- Default is to use gene names (data_type="gene")
  
  - background (optional): which background to use, default is whole proteome, or provide your own list. either "proteome" or your own character vecor of Genes / IDs.
  --- Default is to use whole proteome (background="proteome")
  
  - stat_test (optional): # statistical test to use. One of c("Fisher", "Chi", "Hypergeometric")
  --- Default is to use Fisher test (most stringent) (stat_test="Fisher")
  
  - p_adj (optional): p-adjustment method to use. One of p.adjust.methods: c("holm","hochberg","hommel","bonferroni,"BH","BY","fdr","none")
  --- Default is to use "BH" (Benjamini & Hochberg method) (p_adj="BH)
  
  - thresh (optional): threshold for p-adjusted filtering. Can be any numeric value
  --- Default is 0.05  (thresh=0.05)
  
  - to_plot (optional): The top N domains to plot, ranked by p value. Can be any numeric value or use NA to turn off plotting
  --- Default is to plot top 10 (to_plot=10)
  
  - plot_name (optional): character vector name for the title of the plot

# Example usage

Basic usage

```sh
## with gene names
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


## with uniprot codes
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

domainScanR(input=ids,
            data_type="uniprot")

# save the plot if you want to
ggsave(filename = paste0("domainScanR_example_uniprot_IDs.pdf"),
       width=8,
       height=4)
			
```


# Output
By Default, the function returns a list
- First item is the plot
- Second item is a tibble with the full table of result of the statistical test
--- NOTE: if to_plot=NA is used, just the table is returned

[![Example image](https://ibb.co/z7y4Mwy)](https://github.com/d0minicO/domainScanR/blob/main/example_output/domainScanR_example_uniprot_IDs.png)

#### example table

| domain_name | interpro_id | count_input | total_input | freq_input | count_bkg | freq_bkg | total_bkg | domainRatio | p_adjusted | stat_test | Genes | IDs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CRA domain | IPR013144 | 3 | 13 | 23.08 | 6 | 0.029 | 20407 | 0.5 | 3.38E-07 | Fisher | RANBP10, RANBP9, RMND5A | Q6VN20, Q96S59, Q9H871 |
| CTLH/CRA C-terminal to LisH motif domain | IPR024964 | 3 | 13 | 23.08 | 6 | 0.029 | 20407 | 0.5 | 3.38E-07 | Fisher | RANBP10, RANBP9, RMND5A | Q6VN20, Q96S59, Q9H871 |
| CTLH, C-terminal LisH motif | IPR006595 | 3 | 13 | 23.08 | 10 | 0.049 | 20407 | 0.3 | 7.66E-07 | Fisher | RANBP10, RANBP9, RMND5A | Q6VN20, Q96S59, Q9H871 |
| B30.2/SPRY domain | IPR001870 | 4 | 13 | 30.77 | 96 | 0.470 | 20407 | 0.0417 | 2.99E-06 | Fisher | TRIM21, TRIM27, RANBP10, RANBP9 | P19474, P14373, Q6VN20, Q96S59 |
| SPRY domain | IPR003877 | 4 | 13 | 30.77 | 92 | 0.450 | 20407 | 0.0435 | 2.99E-06 | Fisher | TRIM21, TRIM27, RANBP10, RANBP9 | P19474, P14373, Q6VN20, Q96S59 |
| LIS1 homology motif | IPR006594 | 3 | 13 | 23.08 | 28 | 0.137 | 20407 | 0.107 | 5.98E-06 | Fisher | RANBP10, RANBP9, RMND5A | Q6VN20, Q96S59, Q9H871 |
| Eosinophil major basic protein, C-type lectin-like domain | IPR033816 | 2 | 13 | 15.38 | 2 | 0.0098 | 20407 | 1 | 1.12E-05 | Fisher | PRG2, PRG3 | P13727, Q9Y2Y8 |
| Ran binding protein 9/10, SPRY domain | IPR035782 | 2 | 13 | 15.38 | 2 | 0.0098 | 20407 | 1 | 1.12E-05 | Fisher | RANBP10, RANBP9 | Q6VN20, Q96S59 |


  
### Column Descriptions:
- domain_name  : The name of the protein domain.
- interpro_id  : The InterPro identifier associated with the protein domain.
- count_input  : The number of proteins with the domain in the input list of protein IDs / Gene names.
- total_input  : The total number of proteins in the input list of protein IDs / Gene names.
- freq_input   : The frequency (percentage) of proteins with the domain in the input list of protein IDs / Gene names.
- count_bkg    : The number of proteins with the domain in the background list of protein IDs / Gene names.
- freq_bkg     : The frequency (percentage) of proteins with the domain in the background list of protein IDs / Gene names.
- total_bkg    : The total number of proteins in the background list of protein IDs / Gene names.
- domainRatio  : The ratio of the proteins identified with the domain compared to the total number of proteins with that domain.
- p_adjusted   : The p-value from the statistical test, adjusted for multiple comparisons.
- stat_test    : The statistical test used to compare the frequency of the domain in the input list and the background list.
- Genes        : The gene names associated with the domain in the input list of protein IDs / Gene names.
- IDs          : The Uniprot protein IDs associated with the domain in the input list of protein IDs / Gene names.
  


---

# How does it work?

domainScanR uses the extensive interpro database (https://www.ebi.ac.uk/interpro/) to retreive domain annotations for the reviewed human proteome available on uniprot (2023-04-19, UP000005640)
Statistical tests are done in a similar manner to GO terms enrichment. Either Fisher test, Chi-square test, or Hypergeometric tests are supported. Custom backgrounds are supported. P values are by default adjusted using Benjamini & Hochberg method, though other methods can be applied.

---

# Troubleshooting

dominic.owens ..at?? utoronto.ca