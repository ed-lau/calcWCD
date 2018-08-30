#' @title wcd
#' @description Calculate weighted co-publication distance given a list of PubMed IDs (PMIDs) and an annotation file containing GeneID and PMID links.
#' @param pmids A one-column dataframe containing the total list of PMID accessions in the topic under consideration. A test file is provided in test_pmids.Rds for format requirement.
#' @param annot A tibble of GeneID x PMID relationships, cite counts, and pub year. A  test file containing 100,000 relationships is provided in test_annotation.Rds for format requirement.
#' @param year integer() An integer specifying the cutoff year. DEFAULT is now
#' @keywords semantic similarity
#' @export
#' @examples
#' wcd(pmids = test_pmids, annot = test_annotation)


wcd<- function(pmids,
               annot,
               year = as.integer(format(Sys.Date(), "%Y"))
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

  annot$cit = lapply(annot$Citations, function(x) lt(1,6,2,log10(x+1))) %>% as.numeric()
  annot$imm = lapply(annot$Year, function(x) wt(1,1.25,(year-x+1)/10)) %>% as.numeric()
  annot$w = (1 * annot$cit) + (1 * annot$imm) + 1

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
