# Functions for importing data

#========================================================================>
# Misc functions

# Check specified filename glob
check_glob <- function(filename) {

  files <- Sys.glob(filename)
  if (length(files) == 0) {
    msg = sprintf('No files matching "%s"', filename)
    stop(msg)
  }

  return(files)
}

# Match instrument type
match_instrument <- function(instrument) {

  if (length(instrument) > 1) {
    msg = 'Instrument name must be provided as a single string'
    stop(msg)
  }

  # Trimming terminal whitespace
  instrument = str_trim(instrument)

  patterns = c('z2' = '[zZ]2')

  matches = str_detect(instrument, patterns)

  valid_names <- c('Z2')
  valid_string <- paste(valid_names, collapse = ', ')

  if (sum(matches) > 1) {
    msg = sprintf('Ambiguous instrument name. Supported instruments:\n\t%s',
                   valid_string)
    stop(msg)
  }
  else if (sum(matches) == 0) {
    msg = sprintf('Invalid instrument name. Supported instruments:\n\t%s',
                   valid_string)
    stop(msg)
  }
  else {
    return(names(patterns)[matches])
  }
}

#========================================================================>
# Importing

#------------------------------------------------------------------------
#' Import sample description metadata from text output from 
#' the Coulter Z2 multisizer
#'
#' Imports sample description metadata from text output from a Coulter Z2
#' multisizer and converts it into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#
#' @return A data.frame with columns of acquisition parameters.
#'
#' @examples
#' # To do
#' @export

import_z2_meta <- function(filename) {

  files <- check_glob(filename) 
  n.files <- length(files)

  # Initializing output
  options(stringsAsFactors = FALSE)

  # Setting up vector of regex
  patterns <- c('sample.id' = '(?<=\nSampleID=).*',
                'start.time' = '(?<=\nStartTime=).*',
                'n.bins' = '(?<=nBins=).*',
                'aperture' = '(?<=\nAper=)\\d+',
                'kd' = '(?<=\nKd=).*',
                'volume' = '(?<=\nVol=).*',
                'pa.gain' = '(?<=\nPAGn=).*',
                'mn.gain' = '(?<=\nMnGn=).*',
                'current' = '(?<=\nCur=).*',
                'counts' = '(?<=ThCnts0=).*')
  
  n.patterns <- length(patterns)
  col.names <- c(names(patterns[-n.patterns]), 'count.all', 'count.upper')

  out <- as.data.frame(matrix('', nrow = n.files, ncol = (n.patterns + 1)))
  colnames(out) <- col.names

  # Looping through each file
  for (i in 1:n.files) {

    # Reading file as string
    d.string <- readChar(files[i], file.info(files[i])$size)

    # Checking file type
    file.type <- str_extract(d.string, '(?<=ftype=).*')
    file.type <- tolower(str_replace_all(file.type, '\\s+', ''))
    if (is.na(file.type) || file.type != '300multisizer') {
      msg <- 'Incorrect file type.'
      error(msg)
    }
    
    # Finding split between numeric and metadata 
    meta <- str_split(d.string, '\r\n(?=\\[#Extra\\])')[[1]][1]

    # Matching required patterns
    matches <- str_match(meta, patterns)

    # Getting rid of spaces
    matches <- str_trim(matches)

    # Trimming time
    matches[2] <- paste(str_split(matches[2], '\\s+')[[1]][-1], collapse = ' ')

    # Parsing counts
    counts <- str_split(matches[n.patterns], ',')[[1]]
    matches <- c(matches[-n.patterns], counts[2], counts[3])

    # Attaching matches
    out[i,] <- matches
  }
  
  out$filename <- basename(files)
  out <- out %>%
           select(filename, everything()) %>%
           mutate(sample.id = as.numeric(sample.id),
                  start.time = ymd_hms(start.time, tz='EST'),
                  aperture = as.numeric(aperture),
                  kd = as.numeric(kd),
                  volume = as.numeric(volume),
                  pa.gain = as.numeric(pa.gain),
                  mn.gain = as.numeric(mn.gain),
                  current = as.numeric(current),
                  count.all = as.numeric(count.all),
                  count.upper = as.numeric(count.upper),
                  count.between = count.all - count.upper)

  return(out)
}

#------------------------------------------------------------------------
#' Import numeric data from text output from Coulter Z2 multisizer
#'
#' Imports numeric text output from a Coulter Z2 multisizer and 
#' converts it into long format.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#
#' @return A data.frame with the following columns:  
#'
#' @examples
#' # To do
#' @export

import_z2_data <- function(filename) {

  files <- check_glob(filename) 
  n.files <- length(files)

  # Initializing output
  options(stringsAsFactors = FALSE)

  out <- list()

  # Looping through each file
  for (i in 1:n.files) {

    # Reading file as string
    d.string <- readChar(files[i], file.info(files[i])$size)

    # Checking file type
    file.type <- str_extract(d.string, '(?<=ftype=).*')
    file.type <- tolower(str_replace_all(file.type, '\\s+', ''))
    if (is.na(file.type) || file.type != '300multisizer') {
      msg <- 'Incorrect file type.'
      error(msg)
    }
    
    # Finding split between numeric and metadata 
    d <- str_split(d.string, '\r\n(?=\\[#Extra\\])')[[1]][2]

    # Matching required patterns
    pattern <- '(?<=\\[#Bindiam\\]\r\n).*(?=\r\n\\[Binunits\\])'
    diam <- str_match(d, regex(pattern, dotall = TRUE))
    pattern <- '(?<=\\[#Binheight\\]\r\n).*(?=\r\n\\[end\\])'
    height <- str_match(d, regex(pattern, dotall = TRUE))

    diam <- as.numeric(str_trim(str_split(diam, '\r\n')[[1]]))
    height <- as.numeric(str_trim(str_split(height, '\r\n')[[1]]))

    n <- length(diam)

    out[[basename(files[i])]] <- data.frame(lower = diam[1:(n-1)],
                                            upper = diam[2:n],
                                            count = height[1:(n-1)])
  }

  out <- melt(out, id = c('lower','upper', 'count')) %>%
           select(filename = L1, everything())

  return(out)
}

#========================================================================>
# Generic functions

#------------------------------------------------------------------------
#' Import sample data from Coulter output.
#'
#' Imports sample data from exported output and converts it 
#' into long format. Import procedure depends on specified instrument. 
#' [This function is a placeholder as only Z2 import is currently
#' supported]. Calls import_Z2_data() and/or import_Z2_meta() as
#' necessary.
#'
#' @param filename Filename or pattern for matching multiple filenames.
#' @param instrumet String specifying instrument. Currently limited to Z2.
#' @param data.type One of either "data", "meta", or "full" to specify data 
#'                  to import.
#
#' @return A data.frame.
#'
#' @examples
#' # To do
#' @export

import_coulter <- function(filename, instrument = 'z2', data.type = 'data') {

  instrument = match_instrument(instrument)

  meta.functions <- c('z2' = 'import_z2_meta')
  data.functions <- c('z2' = 'import_z2_data')

  if (!instrument %in% names(meta.functions)) {
    msg <- 'Input for this instrument is not supported'
    stop(msg)
  }

  if (!data.type %in% c('data', 'meta', 'full')) {
    msg <- 'data.type must be one of "data", "meta", or "full"'
    error(msg)
  }

  if (data.type == 'meta') {
    out <- do.call(meta.functions[instrument], c(list(filename = filename)))
  } else if (data.type == 'data') {
    out <- do.call(data.functions[instrument], c(list(filename = filename)))
  } else if (data.type == 'full') {
    m <- do.call(meta.functions[instrument], c(list(filename = filename)))
    d <- do.call(data.functions[instrument], c(list(filename = filename)))
    out <- left_join(d, m, by = c('filename'))
  }

  return(out)
} 
