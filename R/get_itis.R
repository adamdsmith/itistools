#' Retrieve specific taxon information from ITIS based on scientific name query
#'
#' @param scientific_names string of scientific names with which to query ITIS
#' @param timeout integer indicating time, in seconds, to allow for HTTP requests to
#'  process. Default is 30 seconds, which should be more than adequate for most
#'  requests.
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

get_itis <- function(scientific_names, timeout = 20L) {

  scientific_names <- Cap(unique(scientific_names), "first")
  save_sn <- scientific_names
  # Remove blanks
  scientific_names <- scientific_names[!is_missing(scientific_names)]

  # Have to split lengthy requests so API can handle it
  if (length(scientific_names) > 50)
    group_sn <- cut(seq_along(scientific_names),
                    ceiling(length(scientific_names)/50), labels = FALSE)
  else
    group_sn <- rep(1, length(scientific_names))

  for (i in 1:3) { # Try up to 3 times to set up SOLR connection
    con <- try(solrium::SolrClient$new(host = "services.itis.gov",
                                       scheme = "https",
                                       port = NULL), silent = TRUE)
    if (!inherits(con, "error") || i == 3) break
    Sys.sleep(stats::runif(1, 5, 10))
  }

  if (inherits(con, "error")) {
    stop("ITIS connection failed.")
  }

  itis <- lapply(unique(group_sn), function(i) {
    tmp_sn <- scientific_names[which(group_sn == i)]
    # Encode for ITIS query
    tmp_sn <- gsub(' ', '\\\\%20', tmp_sn)
    sci_query <- paste0('nameWOInd:',
                        paste(tmp_sn,
                              collapse = utils::URLencode(" OR nameWOInd:")))
    tmp_sn <- con$search(params = list(q = sci_query,
                                       fl = c('tsn', 'nameWOInd', 'usage',
                                              'rank', 'acceptedTSN',
                                              'vernacular', 'hierarchySoFarWRanks')),
                         # allow room for multiple returned matches
                         rows = length(tmp_sn) + 20,
                         callopts = list(timeout = timeout))
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
    unmatched <- save_sn[which(!(save_sn %in% itis$nameWOInd))]

    # Ensure vernacular column is present..
    itis <- bind_rows(itis, data.frame(vernacular = character(0),
                                       stringsAsFactors = FALSE))

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
    if (length(needs_com_name) > 50)
      group_cn <- cut(seq_along(needs_com_name), ceiling(length(needs_com_name)/50), labels = FALSE)
    else
      group_cn <- rep(1, length(needs_com_name))

    if (length(group_cn) > 0) {
      fix_cn <- lapply(unique(group_cn), function(i) {

        tmp_cn <- needs_com_name[which(group_cn == i)]
        # Encode for ITIS query
        tmp_cn <- gsub(' ', '\\\\%20', tmp_cn)
        sci_query <- paste0('nameWOInd:',
                            paste(tmp_cn,
                                  collapse = utils::URLencode(" OR nameWOInd:")))
        tmp_cn <- con$search(params = list(q = sci_query,
                                           fl = c('nameWOInd', 'vernacular',
                                                  'hierarchySoFarWRanks')),
                             # allow room for multiple returned matches
                             rows = length(tmp_cn) + 20,
                             callopts = list(timeout = timeout))
        if (nrow(tmp_cn) > 0) {
          tmp_cn <- tmp_cn %>%
          # Ensure vernacular column is present..
          bind_rows(data.frame(vernacular = character(0),
                               stringsAsFactors = FALSE)) %>%
          group_by(.data$nameWOInd) %>%
          slice(1) %>% ungroup() %>%
          mutate(itis_com_name = Cap(get_vernac(.data$vernacular)),
                 # Taxon class
                 class = retrieve_class(.data$hierarchySoFarWRanks)) %>%
          select(valid_sci_name = .data$nameWOInd, .data$itis_com_name,
                 .data$class) %>%
          filter(!is.na(.data$itis_com_name))
        } else {
          tibble(valid_sci_name = character(0),
                 itis_com_name = character(0),
                 class = character(0))
        }
      })

      fix_cn <- bind_rows(fix_cn)

      # Handle valid scientific names that for various reasons failed to retrieve from ITIS
      # Typically due to hybrid name formatting (e.g., Gaillardia X grandiflora)
      # TODO: Fix this?
      failed_valid_sn <- needs_com_name[!needs_com_name %in% fix_cn$valid_sci_name]
      if (length(failed_valid_sn) > 0) {
        failed_valid_sn <- itis %>%
          filter(valid_sci_name %in% failed_valid_sn) %>%
          select(valid_sci_name, itis_com_name, class)
        fix_cn <- bind_rows(fix_cn, failed_valid_sn)
      }

      keep <- anti_join(itis, fix_cn, by = "valid_sci_name")
      update <- semi_join(itis, fix_cn, by = "valid_sci_name") %>%
        left_join(fix_cn, by = "valid_sci_name") %>%
        select(.data$sci_name, .data$valid_sci_name,
               itis_com_name = .data$itis_com_name.y,
               class = .data$class.y,
               .data$itis_taxon_rank)
      itis <- bind_rows(keep, update)
    }

    # If necessary, add in unmatched records
    if (length(unmatched) > 0)
      itis <- bind_rows(itis, data.frame(sci_name = unmatched,
                                         stringsAsFactors = FALSE))

    itis

  } else tibble()

}
