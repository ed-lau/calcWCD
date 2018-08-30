#' @title wcd
#' @description Calculate weighted co-publication distance given a list of PubMed IDs (PMIDs) and an annotation file containing GeneID and PMID links.
#' @param pmids A one-column dataframe containing the total list of PMID accessions in the topic under consideration. A test file is provided in test_pmids.Rds for format requirement.
#' @param annotation A tibble of GeneID x PMID relationships, cite counts, and pub year. A  test file containing 100,000 relationships is provided in test_annotation.Rds for format requirement.
#' @keywords semantic similarity
#' @export
#' @examples
#' wcd(pmids = test_pmids, annotation = test_annotation)



wcd<- function(pmids, annotation){

  require(dplyr)

  logistic_transform <- function(a, b, c, x){a/((1 + (b)*exp(-c*x)))}
  weibull_transform <- function(lambda, k, x){(k/lambda) * ((x/lambda)^(k-1)) * exp(-(x/lambda)^k)}

  annotation$cit = lapply(annotation$Citations, function(x) logistic_transform(1,6,2,log10(x+1))) %>% as.numeric()
  annotation$imm = lapply(annotation$Year, function(x) weibull_transform(1,1.25,(as.integer(format(Sys.Date(), "%Y"))-x+1)/10))  %>% as.numeric()
  annotation$w = (1 * annotation$cit) + (1 * annotation$imm) + 1

  ## Subsetting Gene2Pubmed file and retaining only the PMIDs that appear in the Keyword searches (Set of Publications with term)
  annotation. <- dplyr::filter(annotation, PubMed_ID %in% pmids$pmid) #search_term[,1])

  ## If no gene, return null
  if(nrow(annotation.)==0){return(NULL)}

  ## Count the frequency of each gene in the subsetted Gene2Pubmed file
  term_count_table <- dplyr::group_by(annotation., GeneID) %>% dplyr::summarize(Term_count = sum(w)) %>% dplyr::arrange(-Term_count) %>% ungroup()

  ################### GETTING THE CONSTANTS ASSOCIATED WITH THE QUERY #################
  # Tallying the total number of gene-linked publications in the corpus
  total_pubs_tax <- annotation %>% dplyr::distinct(PubMed_ID, .keep_all=T) %>% dplyr::select(w) %>% sum()

  # Tallying the total number of gene-linked publications in the search term
  total_pubs_tax_term <- annotation. %>% dplyr::distinct(PubMed_ID, .keep_all=T) %>% dplyr::select(w) %>% sum()

  ################### CALCULATING SEMANTIC SIMILARITY ############################
  # This is the linked publication counts for a particular gene (any search term)
  total_count_table <- dplyr::group_by(annotation, GeneID) %>% dplyr::summarize(Total_count = sum(w)) %>% dplyr::arrange(-Total_count) %>% dplyr::ungroup()
  term_count_table <- dplyr::left_join(term_count_table, total_count_table)

  # Get WCD numerator
  term_count_table <- dplyr::mutate(term_count_table, WCD_n = ifelse(total_pubs_tax_term > Total_count, log10(total_pubs_tax_term), log10(Total_count)) - log10(Term_count) )
  # Get WCD denominator
  term_count_table <- dplyr::mutate(term_count_table, WCD_d = log10(total_pubs_tax) - ifelse(total_pubs_tax_term < Total_count, log10(total_pubs_tax_term), log10(Total_count)) )
  # Get WCD
  term_count_table <- dplyr::mutate(term_count_table, WCD = (WCD_n/WCD_d))
  term_count_table <- dplyr::mutate(term_count_table, Z = (WCD-mean(term_count_table$WCD))/sd(term_count_table$WCD))
  term_count_table <- dplyr::mutate(term_count_table, P = pnorm(Z))
  term_count_table <- dplyr::mutate(term_count_table, WCD = round(WCD,3), Z = round(Z,2), P = signif(P,2)) %>% dplyr::arrange(P)

  return(term_count_table)

}
