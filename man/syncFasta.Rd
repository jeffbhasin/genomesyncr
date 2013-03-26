\name{syncFasta}
\alias{syncFasta}
\title{Sync whole genome FASTA from UCSC's server}
\usage{
  syncFasta(genome = "", local, url = url.default)
}
\arguments{
  \item{genome}{The UCSC string specific to the genome to
  be downloaded (e.g. "hg19", "hg19", "mm10", etc)}

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
  under \code{local}, download, and extract the whole
  genome FASTA for the genome given as \code{genome}. If
  the FASTA archive already exists, the function will only
  re-download and re-extract if the modified dates are
  different between the local and remote copies.
}

