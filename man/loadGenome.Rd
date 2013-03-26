\name{loadGenome}
\alias{loadGenome}
\title{Load downloaded whole genome FASTA into a GenomicRanges object}
\usage{
  loadGenome(genome, local)
}
\arguments{
  \item{genome}{The UCSC string specific to the genome to
  be downloaded (e.g. "hg19", "hg19", "mm10", etc)}

  \item{local}{The root of a local directory where all
  files downloaded from UCSC are stored}
}
\value{
  A \code{DNAStringSet} object with one sequence per
  chromosome in UCSC's FASTA release for the given genome.
}
\description{
  After downloading the FASTA for a genome using
  \code{\link{syncFasta}}, this function loads all
  chromosomes into a \code{DNAStringSet} object.
}

