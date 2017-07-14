#' Options for word capitalization
#'
#' @param string string of words to capitalize
#' @param words whether to capitalize \code{all} words (title case; default) or
#'  only the \code{first}. \code{first} DOES NOT preserve capitalization.
#'
#' @return string of equal length to \code{string} with capitalization applied
#' @export
#' @examples
#' # Title case
#' Cap("put me in title case please")
#'
#' # First word only
#' Cap("only the first word this time.", "first")
#'
#' # Existing capitalization is NOT preserved
#' Cap("Since Oklahoma is not the first word, it will not be capitalized.", "first")
#'
#' # Vectorization is accommodated
#' Cap(c("title case please", "and me too"))
#'
#' # Handles missing values
#' Cap(c("title case please", NA_character_, "and me too"))

Cap <- function(string, words = c("all", "first")) {
  words <- match.arg(words)
  isna <- is.na(string)
  string <- tolower(string)
  if (identical(words, "all")) {
    s <- strsplit(string, " ")
    s <- sapply(s, function(i) {
      paste(toupper(substring(i, 1,1)), substring(i, 2), sep="", collapse=" ")
    })
  } else {
    s <- paste0(toupper(substr(string, 1, 1)),
                substr(string, 2, nchar(string)))
  }
  s[isna] <- NA_character_
  s
}

#' Check for various forms of missing strings beyond NA and ""
#'
#' @param string string to evaluate for missing content. Missing content is defined
#'  as NA_character_, empty strings, and strings comprised entirely of white space
#'
#' @return logical of equal length to \code{string} indicating whether content is
#'  missing

is_missing <- function(string) {
  is.na(string) | nchar(string) == 0 | grepl("^ +$", string)
}

#' Extract vernacular
#'
#' @param string string returned in 'vernacular' field of
#'  \code{\link[solrium]{solr_search}} to \url{http://services.itis.gov/}. It
#'  is a very specific format...
#'
#' @return string of equal length to \code{string} returning the first vernacular
#'  English name, in all lower case

get_vernac <- function(string) {
  string <- sapply(string, function(i) {
    eng <- grepl("\\$English\\$", i)
    dol_signs <- gregexpr("\\$", i)[[1]]
    if (eng) {
      end <- regexpr("\\$English\\$", i)[[1]] - 1
      strt <- dol_signs[which(dol_signs == regexpr("\\$English\\$", i)[[1]]) - 1] + 1
      string <- substr(i, strt, end)
    } else {
      string <- substr(i, 2, dol_signs[2] - 1)
    }
  })
  tolower(unname(string))
}

#' Extract taxon rank
#'
#' @param string string returned in 'hierarchySoFarWRanks' field of
#'  \code{\link[solrium]{solr_search}} to \url{http://services.itis.gov/}. It
#'  is a very specific format...
#'
#' @return string of equal length to \code{string} returning the lowest taxonomic
#'  rank associated with that taxon

retrieve_rank <- function(string) {
  out <- sapply(string, function(i) {
    tmp <- strsplit(i, "$", fixed = TRUE)[[1]]
    tmp <- tmp[length(tmp)]
    sub(":.*$", "", tmp)
  })
  unname(out)
}

#' Extract valid scientific name, if present rank
#'
#' @param string string returned in 'hierarchySoFarWRanks' field of
#'  \code{\link[solrium]{solr_search}} to \url{http://services.itis.gov/}. It
#'  is a very specific format...
#'
#' @return string of equal length to \code{string} returning the valid ITIS
#'  scientific name, if present

retrieve_sci_name <- function(string) {
  ifelse(grepl("\\$Species:", string),
       sub(".*Species: *(.*?) *\\$.*$", "\\1", string),
       NA_character_)
}
