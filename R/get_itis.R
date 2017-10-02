#' Retrieve specific taxon information from ITIS based on scientific name query
#'
#' @param scientific_names string of scientific names with which to query ITIS
#'
#' @return tibble containing the queried scientific name, the valid ITIS scientific
#'  name if the taxon is of 'Species' rank, the ITIS common name (if present), and
#'  the ITIS taxon rank of the queried scientific name
#'
#' @export
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
    tmp_sn <- solrium::solr_search(q = sci_query,
                                 fl = c('tsn', 'nameWOInd', 'usage', 'rank', 'acceptedTSN',
                                        'vernacular', 'hierarchySoFarWRanks'),
                                 rows = length(tmp_sn) + 20) # allow room for multiple returned matches
    if (nrow(tmp_sn) > 0) {
      tmp_sn <- tmp_sn %>%
        group_by(.data$nameWOInd) %>%
        # If multiple matches, preferentially keep valid/accepted, if available, or first record
        slice(max(which(.data$usage %in% c("accepted", "valid")), 1)) %>%
        ungroup()
    } else tibble()
  })

  itis <- bind_rows(itis)

  if (nrow(itis) > 0) {

    # Save unmatched scientific names to add in later
    unmatched <- scientific_names[which(!(scientific_names %in% itis$nameWOInd))]

    # Add *missing* vernacular if not present..
    if (!("vernacular" %in% colnames(itis))) itis$vernacular <- NA_character_

    # Simplify ITIS data.frame
    itis <- mutate(itis,
                   # Taxon class
                   class = retrieve_class(.data$hierarchySoFarWRanks),
                   # Accepted scientific name
                   valid_sci_name = retrieve_sci_name(.data$hierarchySoFarWRanks),
                   # Return most common common name in ITIS...and capitalize it
                   itis_com_name = Cap(get_vernac(.data$vernacular)),
                   # Get rank if no longer species after correcting TSN
                   itis_taxon_rank = retrieve_rank(.data$hierarchySoFarWRanks)) %>%
      select(sci_name = .data$nameWOInd, .data$valid_sci_name, .data$itis_com_name,
             .data$class, .data$itis_taxon_rank)

    # Retrieve common names and class of changed scientific names, if missing...
    needs_com_name <- itis %>%
      filter(!identical(.data$sci_name, .data$valid_sci_name),
             !is.na(.data$valid_sci_name),
             is.na(.data$itis_com_name)) %>%
      pull(.data$valid_sci_name) %>% unique()

    # Again splitting, if necessary, to keep API happy
    if (length(needs_com_name) > 100) {
      group_cn <- cut(seq_along(needs_com_name), ceiling(length(needs_com_name)/100), labels = FALSE)
    } else group_cn <- rep(1, length(needs_com_name))

    if (length(group_cn) > 0) {
      fix_cn <- lapply(unique(group_cn), function(i) {
        tmp_cn <- needs_com_name[which(group_cn == i)]
        sci_query <- paste0('nameWOInd:(', paste(shQuote(tmp_cn), collapse = " "), ')')
        invisible(solrium::solr_connect("http://services.itis.gov/", verbose = FALSE))
        tmp_cn <- solrium::solr_search(q = sci_query,
                                       fl = c('nameWOInd', 'vernacular'),
                                       # allow room for multiple returned matches
                                       rows = length(tmp_cn) + 20,
                                       callopts = httr::timeout(timeout)) %>%
          # Ensure vernacular column is present..
          bind_rows(data.frame(vernacular = NA_character_,
                               stringsAsFactors = FALSE)) %>%
          group_by(.data$nameWOInd) %>%
          slice(1) %>% ungroup() %>%
          mutate(itis_com_name = Cap(get_vernac(.data$vernacular)),
                 # Taxon class
                 class = retrieve_class(.data$hierarchySoFarWRanks)) %>%
          select(valid_sci_name = .data$nameWOInd, .data$itis_com_name,
                 .data$class) %>%
          filter(!is.na(.data$itis_com_name))
      })

      fix_cn <- bind_rows(fix_cn)
      keep <- anti_join(itis, fix_cn, by = "valid_sci_name")
      update <- semi_join(itis, fix_cn, by = "valid_sci_name") %>%
        left_join(fix_cn, by = "valid_sci_name") %>%
        mutate(itis_com_name = .data$itis_com_name.y) %>%
        select(.data$sci_name, .data$valid_sci_name,
               .data$itis_com_name, .data$class, .data$itis_taxon_rank)
      itis <- bind_rows(keep, update)
    }

    # If necessary, add in unmatched records
    if (length(unmatched) > 0)
      itis <- bind_rows(itis, data.frame(sci_name = unmatched,
                                         stringsAsFactors = FALSE))

    itis

  } else tibble()

}
