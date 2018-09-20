#' @title wcd
#' @description Calculate weighted co-publication distance given a list of PubMed IDs (PMIDs) and an annotation file containing GeneID and PMID links.
#' @param pmids A one-column tibble containing the total list of PMID accessions in the topic under consideration. A test file is provided in the package (test_pmids) for format requirement.
#' @param annot A tibble of GeneID x PMID relationships, cite counts, and pub year. A  test file containing 100,000 relationships is provided in the package (test_annotation) for format requirement.
#' @param year integer() An integer specifying the cutoff year. DEFAULT is now
#' @param weighted logical() Whether to calculate weighted co-publication factor (T) or non-weighted co-publication factor (F). DEFAULT is T
#' @keywords semantic similarity
#' @export
#' @examples
#' wcd(pmids = test_pmids, annot = test_annotation)


wcd<- function(pmids,
               annot,
               year = as.integer(format(Sys.Date(), "%Y")),
               weighted = TRUE
               ){

  require(dplyr)

  # Logistic transform
  lt <- function(a, b, c, x){
    a/((1 + (b)*exp(-c*x)))
  }

  # Weibull distribution
  wt <- function(lambda, k, x){
    (k/lambda) * ((x/lambda)^(k-1)) * exp(-(x/lambda)^k)
  }

  if(!is.numeric(year) | year < 1950 | year > as.integer(format(Sys.Date(), "%Y"))){stop("Please enter a year between 1950 and present.")}

  # If weighted, calculate the immediacy and impact of the publications, otherwise each publication weights 1.00
  if(weighted){
    annot_small_cit <- annot %>% dplyr::distinct(Citations)
    annot_small_cit$cit <- lapply(annot_small_cit$Citations, function(x) lt(1,6,2,log10(x+1))) %>% as.numeric()
    annot <- annot %>% dplyr::left_join(annot_small_cit)

    annot_small_year <- annot %>% dplyr::distinct(Year)
    annot_small_year$imm = lapply(annot_small_year$Year, function(x) wt(1,1.25,(year-x+1)/10)) %>% as.numeric()
    annot <- annot %>% dplyr::left_join(annot_small_year)

    annot$w = (1 * annot$cit) + (1 * annot$imm) + 1
  } else {
    annot$w = 1
  }

  ## Subsetting Gene2Pubmed file and retaining only the PMIDs that appear in the Keyword searches (Set of Publications with term)
  annot. <- dplyr::filter(annot, PubMed_ID %in% pmids$pmid, Year <= year)

  ## If no qualifying left, exit
  if(nrow(annot.)==0){stop("No qualifying publication found.")}



  ## Count the frequency of each gene in the subsetted Gene2Pubmed file
  termCount <- dplyr::group_by(annot., GeneID) %>% dplyr::summarize(Term_count = sum(w)) %>% dplyr::arrange(-Term_count) %>% dplyr::ungroup()

  # GETTING THE CONSTANTS ASSOCIATED WITH THE QUERY
  # Tallying the total number of gene-linked publications in the corpus
  nP <- annot %>%  dplyr::filter(Year <= year) %>% dplyr::distinct(PubMed_ID, .keep_all=T) %>% dplyr::select(w) %>% sum()

  # Tallying the total number of gene-linked publications in the search term
  nPt <- annot. %>% dplyr::distinct(PubMed_ID, .keep_all=T) %>% dplyr::select(w) %>% sum()

  # CALCULATING SEMANTIC SIMILARITY
  # This is the linked publication counts for a particular gene (any search term)
  totalCount <- annot %>% dplyr::filter(Year <= year) %>%
                          dplyr::group_by(GeneID) %>%
                          dplyr::summarize(Total_count = sum(w)) %>%
                          dplyr::arrange(-Total_count) %>%
                          dplyr::ungroup()
  termCount <- termCount %>% dplyr::left_join(totalCount)

  # Get WCD numerator
  termCount <- termCount %>% dplyr::mutate(WCD_n = ifelse(nPt > Total_count, log10(nPt), log10(Total_count)) - log10(Term_count) )
  # Get WCD denominator
  termCount <- termCount %>% dplyr::mutate(WCD_d = log10(nP) - ifelse(nPt < Total_count, log10(nPt), log10(Total_count)))
  # Get WCD
  termCount <- termCount %>% dplyr::mutate(WCD = (WCD_n/WCD_d))
  termCount <- termCount %>% dplyr::mutate(Z = (WCD-mean(termCount$WCD))/sd(termCount$WCD)) %>%
                             dplyr::mutate(P = pnorm(Z)) %>%
                             dplyr::arrange(P)

  return(termCount)

}


#' @title pmid
#' @description Retrieves a list of PMIDs from the NCBI EUtils API based on a search term
#' @param q character() A character string specifying search topic
#' @keywords PubMed query
#' @export
#' @examples
#' pmid("fibrosis NOT cystic")


pmid <- function(q){

  require(dplyr)
  require(httr)
  require(foreach)
  require(doParallel)
  require(iterators)
  require(snow)
  require(xml2)

  cl <- makeCluster(4, type='SOCK', outfile="")
  registerDoParallel(cl)


  # 2018-04-24 Fix hyphen
  q <- gsub("-", " ", q)

  # 2017-08-25 sanitize all other pucntuations
  q <- gsub("[[:punct:]]", "", q)

  # 2017-09-06 change space to + sign after sanitizing other punctuations, not before, otherwise queries with more than one words won't work
  q <- gsub(" ","+", q)

  pub_count_address <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=", q, "&rettype=count")

  # NewL: 2018-04-20: use httr and xml2 to get pub count
  pc <- httr::RETRY("GET", pub_count_address)
  total_pub_count <- xml2::xml_contents(xml2::xml_children(httr::content(pc))) %>% as.character() %>% as.integer()

  # Quit if no PMIDs
  #if(total_pub_count == 0){return(data.frame(pmid = character()))}
  if(total_pub_count == 0){stop("The query resulted in no retrievable publication.")}

  ret_interval <- min(99999, total_pub_count)
  ret_times <- total_pub_count %/% ret_interval


  # Alternative 2018-04-20: now using httr and xml2 to retrieve and read the results, with parallelization and automatic retry.

  row_vec = seq(from=0, to=ret_times, length=ret_times+1)
  all_pmids <- foreach (this_row=row_vec, .combine='c', .packages=c('httr')) %dopar%{

    pmid_address <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",q,"&retmax=",
                           ret_interval, "&retstart=", this_row*ret_interval)
    # Retrieve the address content, with retry
    ret = httr::RETRY("GET", pmid_address)
    # Read the retrieved content as XML, then read the contents of the XML children nodes
    ret_content = xml2::xml_contents(xml2::xml_children(httr::content(ret)))
    # Read the contents of the Id fields again to get the PMID names
    as.character(xml2::xml_contents(ret_content))
  }

  # Convert to data.frame for backward compatibility
  all_pmids <- tibble(pmid=all_pmids)

  stopCluster(cl)

  #all_pmids <- mutate(all_pmids, pmid = as.integer(pmid))
  return(all_pmids)
}

