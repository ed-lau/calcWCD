#' @title wcd
#' @description Calculate weighted co-publication distance given a list of PubMed IDs (PMIDs) and an annotation file containing GeneID and PMID links.
#' @param pmids A one-column tibble containing the total list of PMID accessions in the topic under consideration. A test file is provided in the package (test_pmids) for format requirement.
#' @param annot A tibble of GeneID x PMID relationships, cite counts, and pub year. A  test file containing 100,000 relationships is provided in the package (test_annotation) for format requirement.
#' @param year integer() An integer specifying the cutoff year. DEFAULT is now
#' @param weighted logical() Whether to calculate weighted co-publication factor (T) or non-weighted co-publication factor (F). DEFAULT is T
#' @keywords semantic similarity
#' @returns weighted correlation table
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
    # annot_small_cit <- annot %>% dplyr::distinct(Citations)
    # annot_small_cit$cit <- lapply(annot_small_cit$Citations, function(x) lt(1,6,2,log10(x+1))) %>% as.numeric()
    # annot <- annot %>% dplyr::left_join(annot_small_cit)
    #
    # annot_small_year <- annot %>% dplyr::distinct(Year)
    # annot_small_year$imm = lapply(annot_small_year$Year, function(x) wt(1,1.25,(year-x+1)/10)) %>% as.numeric()
    # annot <- annot %>% dplyr::left_join(annot_small_year)

    annot <- annot %>% dplyr::mutate(cit = lt(1, 6, 2, log10(Citations + 1)),
                                     imm = wt(1, 1.25, (year - Year + 1)/ 10))

    annot$w = (1 * annot$cit) + (1 * annot$imm) + 1
  } else {
    annot$w = 1
  }

  ## Subsetting Gene2Pubmed file and retaining only the PMIDs that appear in the Keyword searches (Set of Publications with term)
  annot. <- dplyr::filter(annot, PubMed_ID %in% pmids$pmid, Year <= year)

  ## If no qualifying left, exit
  if(nrow(annot.)==0){
    print("No qualifying publication found.")
    null_df <- data.frame(GeneID=NA,
                          Term_count=NA,
                          Total_count=NA,
                          WCD_n=NA,
                          WCD_d=NA,
                          WCD=NA,
                          Z=NA,
                          P=NA)

    return (null_df) }



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
#' @description Retrieves a list of PMIDs from the EuropePMC REST API based on a search term.
#' Note that as of 2023, NCBI EUtils E Search only returns the first 10,000 entries.
#' Hence we now use the EuropePMC API which returns 1000 results per page but allows multiple page retrieval
#' Because of the slower speed we now add a max_retrieval option
#' @param q character() A character string specifying search topic
#' @param max_retrieval integer() An integer specifying the maximum number of retrieved PMIDs
#' @keywords PubMed query
#' @returns A tibble containing one column named pmid, with the list of retrieved PMIDs.
#' @export
#' @examples
#' pmid("fibrosis NOT cystic", max_retrieval=1000)


pmid <- function(q,
                 max_retrieval = Inf){

  require(dplyr)
  require(httr)
  # require(foreach)
  # require(doParallel)
  # require(iterators)
  # require(snow)
  # require(xml2)

  # cl <- makeCluster(4, type='SOCK', outfile="")
  # registerDoParallel(cl)


  # 2018-04-24 Fix hyphen
  q <- gsub("-", " ", q)

  # 2017-08-25 sanitize all other pucntuations
  q <- gsub("[[:punct:]]", "", q)

  # 2023-07-13 for the European PMC Rest API, replace space with %20 in the URL
  q <- gsub(" ","%20", q)

  # 2017-09-06 change space to + sign after sanitizing other punctuations, not before, otherwise queries with more than one words won't work
  # q <- gsub(" ","+", q)

  # pub_count_address <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=", q, "&rettype=count")
  #
  # # NewL: 2018-04-20: use httr and xml2 to get pub count
  # pc <- httr::RETRY("GET", pub_count_address)
  # total_pub_count <- xml2::xml_contents(xml2::xml_children(httr::content(pc))) %>% as.character() %>% as.integer()
  #
  # print(total_pub_count)
  #
  # # Quit if no PMIDs
  # #if(total_pub_count == 0){return(data.frame(pmid = character()))}
  # if(total_pub_count == 0){stop("The query resulted in no retrievable publication.")}
  #
  # ret_interval <- min(9999, total_pub_count)
  # ret_times <- total_pub_count %/% ret_interval


  # Alternative 2018-04-20: now using httr and xml2 to retrieve and read the results, with parallelization and automatic retry.
  # row_vec = seq(from=0, to=ret_times, length=ret_times+1)
  # all_pmids <- foreach (this_row=row_vec, .combine='c', .packages=c('httr')) %dopar%{
  #
  #   pmid_address <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",q,"&retmax=",
  #                          ret_interval, "&retstart=", this_row*ret_interval, '&apikey=', apikey)
  #   print(pmid_address)
  #   # Retrieve the address content, with retry
  #   ret = httr::RETRY("GET", pmid_address)
  #   # Read the retrieved content as XML, then read the contents of the XML children nodes
  #   ret_content = xml2::xml_contents(xml2::xml_children(httr::content(ret)))
  #   # Read the contents of the Id fields again to get the PMID names
  #   as.character(xml2::xml_contents(ret_content))
  # }
  # stopCluster(cl)

  #



  read_result <- function(ret, lis){
    # """
    # Get PMIDs from the returned page
    # ret: retrieved content
    # lis: total list
    # """
    if(length(ret_content$resultList$result) > 0){
      for(i in 1:length(ret$resultList$result)){
          lis <- c(lis, ret$resultList$result[[i]]$pmid)

      }
    }
    return (lis)
  }


  total_list <- c()
  # Initial cursor
  cur <- "*"
  # Initial URL
  next_url <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=", q, "&resultType=idlist&synonym=FALSE&cursorMark=", cur, "&pageSize=1000&format=json&sort_cited:y")


  # Loop through the pages while there is a "next_url" and we haven't downloaded enough PMIDs to reach the max_retrieval yet.
  while(!is.null(next_url) & length(total_list) < max_retrieval){
    ret = httr::RETRY("GET", gsub(" ", "%20", next_url))
    ret_content = httr::content(ret, type  = "application/json")

    pc <- ret_content$hitCount
    if(pc == 0){ stop("The query resulted in no retrievable publication.")}

    total_list <- read_result(ret_content, total_list)

    print(paste0(length(unique(total_list)), " of ", pc))

    next_url <- ret_content$nextPageUrl

  }



  # Convert to data.frame for backward compatibility
  all_pmids <- tibble(pmid=total_list)
  print(all_pmids)


  #all_pmids <- mutate(all_pmids, pmid = as.integer(pmid))
  return(all_pmids)
}

# Time comparison
# system.time(pmid("heart failure", max_retrieval=10000))
# system.time(europepmc::epmc_search("heart failure", limit = 10000))

# Test
# gene2pubmed_test <- readRDS("data/g2p_pub_comb_cit_precalc_9606_2018-04-13_small.Rds")
# gene2pubmed_test <- readRDS("data/g2p_pub_comb_cit_temp_9606_2023-07-23_small.Rds")
# pmid_test <- pmid("Heart Failure", max_retrieval = 10000)
# wcd(pmid_test, gene2pubmed_test)
