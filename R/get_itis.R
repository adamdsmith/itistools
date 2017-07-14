#' Retrieve specific taxon information from ITIS based on scientific name query
#'
#' @param scientific_names string of scientific names with which to query ITIS
#'
#' @return tibble containing the queried scientific name, the valid ITIS scientific
#'  name if the taxon is of 'Species' rank, the ITIS common name (if present), and
#'  the ITIS taxon rank of the queried scientific name

#' @examples
#' \dontrun{
#' sn <- c("Thryothorus ludovicianus", "Cardinalis cardinalis", "Bidens alba",
#'         "Pachydiplax longipennis", "Cladonia evansii", "Ischnura hastata",
#'         "Smilax bona-nox", "Plantae", "Malvaviscus", "Anaxyrus", "Gulo gulo")
#' get_itis(sn)
#' }

get_itis <- function(scientific_names) {

  scientific_names <- unique(scientific_names)

  # Have to split lengthy requests so API can handle it
  if (length(scientific_names) > 100) {
    group_sn <- cut(seq_along(scientific_names), ceiling(length(scientific_names)/100), labels = FALSE)
  } else group_sn <- rep(1, length(scientific_names))

  itis <- lapply(unique(group_sn), function(i) {
    tmp_sn <- scientific_names[which(group_sn == i)]
    sci_query <- paste0('nameWOInd:(', paste(shQuote(tmp_sn), collapse = " "), ')')
    invisible(solrium::solr_connect("http://services.itis.gov/", verbose = FALSE))
    itis <- solrium::solr_search(q = sci_query,
                                 fl = c('tsn', 'nameWOInd', 'usage', 'rank', 'acceptedTSN',
                                        'vernacular', 'hierarchySoFarWRanks'),
                                 rows = length(tmp_sn) + 20) %>% # allow room for multiple returned matches
      group_by(nameWOInd) %>%
      # If multiple matches, preferentially keep valid/accepted, if available, or first record
      slice(max(which(usage %in% c("accepted", "valid")), 1)) %>%
      ungroup()
  })
  itis <- bind_rows(itis)

  # Save unmatched scientific names to add in later
  unmatched <- scientific_names[which(!(scientific_names %in% itis$nameWOInd))]

  # Add *missing* vernacular if not present..
  if (!("vernacular" %in% colnames(itis))) itis$vernacular <- NA_character_

  # Simplify ITIS data.frame
  itis <- mutate(itis,
                 # Accepted scientific name
                 valid_sci_name = retrieve_sci_name(hierarchySoFarWRanks),
                 # Return most common common name in ITIS...and capitalize it
                 itis_com_name = Cap(get_vernac(vernacular)),
                 # Get rank if no longer species after correcting TSN
                 itis_taxon_rank = retrieve_rank(hierarchySoFarWRanks)) %>%
    select(sci_name = nameWOInd, valid_sci_name, itis_com_name, itis_taxon_rank = rank)

  # Retrieve common names of changed scientific names, if missing...
  needs_com_name <- itis %>%
    filter(!identical(sci_name, valid_sci_name),
           is.na(itis_com_name)) %>%
    pull(valid_sci_name)

  # Again splitting, if necessary, to keep API happy
  if (length(needs_com_name) > 100) {
    group_cn <- cut(seq_along(needs_com_name), ceiling(length(needs_com_name)/100), labels = FALSE)
  } else group_cn <- rep(1, length(needs_com_name))

  fix_cn <- lapply(unique(group_cn), function(i) {
    tmp_cn <- needs_com_name[which(group_cn == i)]
    sci_query <- paste0('nameWOInd:(', paste(shQuote(tmp_cn), collapse = " "), ')')
    invisible(solrium::solr_connect("http://services.itis.gov/", verbose = FALSE))
    itis <- solrium::solr_search(q = sci_query,
                                 fl = c('nameWOInd', 'vernacular'),
                                 rows = length(tmp_cn) + 20) %>% # allow room for multiple returned matches
      group_by(nameWOInd) %>%
      slice(1) %>% ungroup() %>%
      mutate(itis_com_name = Cap(get_vernac(vernacular))) %>%
      select(valid_sci_name = nameWOInd, itis_com_name)
  })
  fix_cn <- bind_rows(fix_cn)
  itis$itis_com_name[match(fix_cn$valid_sci_name, itis$valid_sci_name)] <- fix_cn$itis_com_name

  # If necessary, add in unmatched records
  if (length(unmatched) > 0)
    itis <- bind_rows(itis, data.frame(sci_name = unmatched,
                                       stringsAsFactors = FALSE))

  itis

}
