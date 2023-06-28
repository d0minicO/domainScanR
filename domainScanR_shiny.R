
### function to quantify enrichment of protein domains in a list of proteins
## MODIFIED TO RUN ON SHINY APP
## PASSES THE INTERPRO DATABASE AS AN ARGUMENT RATHER THAN DYNAMICALLY LOAD FROM GITHUB

domainScanR_shiny <- function(input,
                        interpro,
                        data_type="gene",
                        background = "proteome",
                        stat_test = "Fisher",
                        p_adj = "BH",
                        thresh = 0.05,
                        to_plot = 10,
                        plot_name=NULL
                        ){
  
  ### function to quantify enrichment of protein domains in a list of proteins
  
  # required packages
  require(tidyverse)
  require(magrittr)
  require(viridis)
  require(httr)
  
  
  ###################
  #### Arguments ####
  ###################
  
  # input: the list of genes or proteins to search for enrichment of domains in
  ## a character vector of gene names OR uniprot IDs to search
  
  # interpro: the tibble containing data from interpro.Rds (needed to load from the shiny/data directory)
  
  # data_type (optional): whether gene names or uniprot IDs have been provided
  ## one of c("gene", "uniprot)
  ## default is to use gene names (data_type="gene")
  
  # background (optional): which background to use, default is whole proteome, or provide your own list
  ## either "proteome" or your own character vecor of Genes / IDs
  ## default is to use whole proteome (background="proteome")
  
  # stat_test (optional): # statistical test to use
  # one of c("Fisher", "Chi", "Hypergeometric")
  ## default is to use Fisher test (most stringent) (stat_test="Fisher")
  
  # p_adj (optional): p-adjustment method to use
  # one of p.adjust.methods: c("holm","hochberg","hommel","bonferroni,"BH","BY","fdr","none")
  ## default is to use "BH" (Benjamini & Hochberg method) (p_adj="BH)
  
  # thresh (optional): threshold for p-adjusted filtering
  # can be any numeric value
  ## default is 0.05  (thresh=0.05)
  
  # to_plot (optional):
  # can be any numeric value or use NA to turn off plotting
  ## default is to plot top 10 (to_plot=10)
  
  # plot_name (optional): character vector name for the title of the plot
  
  
  ################
  #### Return ####
  ################
  
  # default is to return a list
  # first item is the plot
  # second item is a tibble with the full table of result of the statistical test
  # NOTE: if to_plot=NA is used, just the table is returned
  
  # Column Descriptions:
  # domain_name  : The name of the protein domain.
  # interpro_id  : The InterPro identifier associated with the protein domain.
  # count_input  : The number of proteins with the domain in the input list of protein IDs / Gene names.
  # total_input  : The total number of proteins in the input list of protein IDs / Gene names.
  # freq_input   : The frequency (percentage) of proteins with the domain in the input list of protein IDs / Gene names.
  # count_bkg    : The number of proteins with the domain in the background list of protein IDs / Gene names.
  # freq_bkg     : The frequency (percentage) of proteins with the domain in the background list of protein IDs / Gene names.
  # total_bkg    : The total number of proteins in the background list of protein IDs / Gene names.
  # domainRatio  : The ratio of the proteins identified with the domain compared to the total number of proteins with that domain.
  # p_adjusted   : The p-value from the statistical test, adjusted for multiple comparisons.
  # stat_test    : The statistical test used to compare the frequency of the domain in the input list and the background list.
  # Genes        : The gene names associated with the domain in the input list of protein IDs / Gene names.
  # IDs          : The Uniprot protein IDs associated with the domain in the input list of protein IDs / Gene names.
  
  
  ###################
  #### LOAD DATA ####
  ###################
  
  ### load interpro annotated human proteome file from github
  #cat("Loading interpro database from github... \n")
  
  cat("Loading interpro database... \n")
  
  # Define the URL of the .Rds file on GitHub
  #url = "https://github.com/d0minicO/interactR/blob/main/data/Interpro.Rds?raw=TRUE"
  
  # Download and read in the .Rds file
  # GET the data
  #response = GET(url)
  
  # Write the content to a temporary file
  #temp = tempfile(fileext = ".Rds")
  #writeBin(content(response), temp)
  
  # Read the .Rds file
  #interpro = readRDS(temp)
  
  if(nrow(interpro)==20407){
    ### load interpro annotated human proteome file from github
    cat("Database looks good! \n")
  } else {
    warning("Interpro databse seems off... sorry! Contact dominic.owens at utoronto.ca\n")
    return()
  }
  
  
  # Check if the file exists before trying to remove it
  #if (file.exists(temp)) {
  #  # Remove the temporary file
  #  file.remove(temp)
  #  if (file.exists(temp)) {
  #    message("Database temp file removal failed")
  #  } else {
  #    cat("Database temp file successfully removed\n")
  #  }
  #} else {
  #  message("Database file does not exist...")
  #}
  
  
  
  ########################
  ## MAIN FUNCTION WORK ##
  ########################
  
  
  #####
  # 1 #
  #####
  
  ## drop any NAs as this causes too many matches at step 3
  input = input[!is.na(input)]
  
  ## construct a df with colname matching the data type
  ## check to make sure data_type is either "gene" or "uniprot"
  
  if(data_type=="gene"){
    cat("Using gene name \n")
    input = tibble(Gene=input)
  } else if (data_type=="uniprot"){
    input = tibble(ID=input)
  } else {
    warning("data_type must be set to either\ngene (HGNC approved gene symbol) OR \nuniprot (uniprot accession ID)")
    return()
  }
  
  
  #####
  # 2 #
  #####
  
  ## chose the background set and if custom then intersect with full proteome and report matches
  if(background =="proteome"){
    bkg=interpro
    cat("Uniprot proteome with", nrow(bkg), "entries used for background\n")
  } else if (background !="proteome"& data_type=="uniprot"){
    cat("user-provided custom background, expecting Uniprot IDs as data_type should be same for background and input\n")
    bkg=
      interpro %>%
      filter(ID %in% background)
    # number of background ids matched
    cat(nrow(bkg), "backround IDs found in interpro database out of", length(background),"\n")
  } else if (background !="proteome"& data_type=="gene"){
    cat("user-provided custom background, expecting gene names as data_type should be same for background and input\n")
    bkg=
      interpro %>%
      filter(Gene %in% background)
    # number of background ids matched
    cat(nrow(bkg), "backround Genes found in interpro database out of", length(background),"\n")
  }
  
  
  #####
  # 3 #
  #####
  
  ## intersect input with annotated proteome]
  input_data = left_join(input,interpro)
  
  # report match lengths
  
  if(data_type=="gene" & nrow(input_data)>1){
    cat(nrow(input_data), "input genes found in interpro database out of", nrow(input),"\n")
  } else if (data_type=="uniprot" & nrow(input_data)>0){
    cat(nrow(input_data), "input IDs found in interpro database out of", nrow(input),"\n")
  } else {
    warning("Not enough matches to interpro database, check input data_type","\n")
    return()
  }
  
  
  #input_data %>%
  #  group_by(Gene) %>%
  #  dplyr::count() %>%
  #  arrange(desc(n))
  
  #####
  # 4 #
  #####
  
  ## calulate domain frequencies in input/ background
  cat("Calculating domain frequencies \n")
  
  domFreqCount <- function(data){
    counts =
      data %>%
      mutate(split=strsplit(domains,";")) %>%
      unnest(cols = c(split)) %>%
      drop_na() %>%
      group_by(split) %>%
      dplyr::count() %>%
      ungroup() %>%
      arrange(desc(n)) %>%
      separate(split,into=c("domain_name","interpro_id",NA),sep=":::") %>%
      mutate(domain_name=trimws(domain_name)) %>%
      dplyr::rename(protein_count_with_dom=n)
    return(counts)
  }
  
  input_counts = domFreqCount(input_data)
  bkg_counts = domFreqCount(bkg)
  
  
  #####
  # 5 #
  #####
  
  ## now calculate frequencies of total proteins in input and background
  input_counts %<>%
    mutate(total_input=nrow(input_data)) %>%
    mutate(freq_input=protein_count_with_dom*100/total_input)
  
  bkg_counts %<>%
    mutate(total_bkg=nrow(bkg)) %>%
    mutate(freq_bkg=protein_count_with_dom*100/total_bkg)
  
  
  
  #####
  # 6 #
  #####
  
  ## join to both to construct a table for stats
  ## seperate on over and under represented
  ## just keep over represented domains (frequency above background)
  comb = 
    full_join(input_counts,bkg_counts,by="interpro_id") %>%
    drop_na() %>%
    dplyr::rename(domain_name = domain_name.x,
                  count_input = protein_count_with_dom.x,
                  count_bkg = protein_count_with_dom.y) %>%
    dplyr::select(domain_name,
                  interpro_id,
                  count_input,
                  total_input,
                  freq_input,
                  count_bkg,
                  freq_bkg,
                  total_bkg) %>%
    filter(freq_input>freq_bkg)
  
  
  #####
  # 8 #
  #####
  
  ## now loop through each overrepresented domain found and perform fischer, chi-square, and hypergeometric tests
  cat("Now performing statistical tests \n")
  
  
  domains = unique(comb$interpro_id)
  d = domains[1]
  out.df = tibble()
  for(d in domains){
    
    #cat(d,"\n")
    
    temp =
      comb %>% 
      filter(interpro_id==!!d)
    
    ## get numbers for a contingency table
    # for the interactors
    input_with = temp$count_input
    input_without = temp$total_input-input_with
    
    # for the background
    bkg_with = temp$count_bkg
    bkg_without = temp$total_bkg-bkg_with
    
    # construct a contingency table
    cont_table =
      data.frame(input = c(input_with,input_without),
                 background = c(bkg_with,bkg_without))
    
    # chisq test p val
    # suppressing  Chi-squared approximation may be incorrect warnings
    res=suppressWarnings({
      chisq.test(cont_table)
    })
    chi_p = res$p.value
    
    # fisher test
    res=fisher.test(cont_table)
    fish_p = res$p.value
    
    # hypergeometric test
    #M = Total number of genes
    M = temp$total_bkg
    #n = Total number of DE genes
    n = temp$total_input
    #N = Total number of genes with a specific GO term
    N = temp$count_bkg
    #x = Number of DE genes with the GO term
    x = temp$count_input
    
    # perform hypergeometric test
    hyper_pvalue = 1 - phyper(x - 1, N, M - N, n)
    
    
    ## construct an output df
    ## calculate gene ratio
    temp %<>%
      mutate(Fisher_p=fish_p,
             Chi_p=chi_p,
             Hypergeometric_p=hyper_pvalue,
             domainRatio = count_input/count_bkg)
    
    out.df %<>% rbind.data.frame(temp)
    
  }
  
  
  
  
  #####
  # 9 #
  #####
  
  ## adjust p-values and filtering
  cat("Adjust p-values and filtering based on",stat_test,"test \n")
  
  ## adjust all the p values
  out.df$Fisher_p.adj = p.adjust(out.df$Fisher_p,method=p_adj)
  out.df$Chi_p.adj = p.adjust(out.df$Chi_p,method=p_adj)
  out.df$Hypergeometric_p.adj = p.adjust(out.df$Hypergeometric_p,method=p_adj)
  
  
  # filter based on this test p value
  
  if(stat_test=="Fisher"){
    out.df %<>%
      filter(Fisher_p.adj<thresh) %>%
      mutate(p_adjusted=Fisher_p.adj,
             stat_test=stat_test)
  } else if (stat_test=="Chi"){
    out.df %<>%
      filter(Chi_p.adj<thresh) %>%
      mutate(p_adjusted=Chi_p.adj,
             stat_test=stat_test)
  } else if (stat_test=="Hypergeometric"){
    out.df %<>%
      filter(Hypergeometric_p.adj<thresh) %>%
      mutate(p_adjusted=Hypergeometric_p.adj,
             stat_test=stat_test)
  } else {
    message("Incorrect statistical test name used, must be either Fisher, Chi, or Hypergeometric \n")
    return()
  }
  
  # clean up the table
  out.df %<>%
    dplyr::select(1:8,domainRatio,p_adjusted,stat_test)
  
  
  #####
  # 9 #
  #####
  
  ## adjust p-values and filtering
  cat("Reporting proteins with enriched domains \n")
  
  ## return the proteins found in each category as a comma separated column
  out.df =
    input_data %>%
    mutate(split=strsplit(domains,";")) %>%
    unnest(cols = c(split)) %>%
    drop_na() %>%
    tidyr::separate(split,into=c("domain_name","interpro_id",NA),sep=":::") %>%
    mutate(domain_name=trimws(domain_name)) %>%
    dplyr::select(Gene,ID,domain_name,interpro_id) %>%
    filter(interpro_id %in% out.df$interpro_id) %>%
    right_join(out.df,by = c("domain_name", "interpro_id")) %>%
    group_by(domain_name,across(-c(Gene,ID))) %>%
    summarise(Genes = paste(Gene, collapse = ", "),
              IDs = paste(ID, collapse = ", "),.groups = "drop") %>%
    arrange(p_adjusted)
  
  
  
  ######
  # 10 #
  ######
  
  ## plot the top X enriched domains
  if(!is.na(to_plot)){
    cat("Plotting circle plots of top", to_plot, "enriched domains \n")
    
    temp =
      out.df %>%
      dplyr::slice(1:to_plot) %>%
      mutate(name=paste0(domain_name," (",interpro_id,")"))
    
    temp$name=factor(temp$name,levels=rev(temp$name))
    
    p =
      ggplot(temp,aes(x=-log10(p_adjusted),y=name,fill=domainRatio*100,size=freq_input)) +
      geom_segment(aes(x = 0, y = name, xend = -log10(p_adjusted), yend = name), color = "grey",size=0.1,linetype="dashed")+
      geom_point(shape=21)+
      labs(y="Interpro domain",x=paste0("-log10(",stat_test,".p.adjusted)"),size="% of proteins\nwith domain",fill="% of domain-\ncontaining proteins\nfound")+
      scale_fill_viridis()+
      theme_bw()+
      theme(panel.grid=element_blank(),
            line=element_line(size=0.1),
            legend.key.size = unit(2,"mm"))
    
    ## add a title if requested
    if (!is.null(plot_name)){
      p = 
        p+
        ggtitle(plot_name)
    }
    
  }
  
  
  
  ######
  # 11 #
  ######
  
  ## construct a list to return
  ## either a list of table and plot
  ## or just the table
  
  if(!is.na(to_plot)){
    cat("Returning plot and table as a list \n")
    
    out_data = list(p,out.df)
    
  } else if (is.na(to_plot)){
    cat("Returning table only \n")
    out_data = out.df
  }
  
  return(out_data)

}
