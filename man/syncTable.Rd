\name{syncTable}
\alias{syncTable}
\title{Sync whole genome annotation tables from UCSC's server}
\usage{
  syncTable(genome, table, local, url = url.default)
}
\arguments{
  \item{genome}{The UCSC string specific to the genome to
  be downloaded (e.g. "hg19", "hg19", "mm10", etc)}

  \item{table}{The UCSC string specific for the table to
  sync (e.g. "knownGene", "kgXref", etc)}

  \item{local}{The root of a local directory where all
  files downloaded from UCSC are stored}

  \item{url}{The root of the remote rsync URL to download
  UCSC data from (set by default to
  \code{rsync://hgdownload.cse.ucsc.edu/goldenPath/})}
}
\value{
  Prints status of the sync to the screen. No output
  object.
}
\description{
  This function will create the necessary directory tree
  under \code{local}, download, and extract the text file
  data and schema for given table from the given genome. If
  the data already exists, the function will only
  re-download and re-extract if the modified dates are
  different between the local and remote copies.
}

