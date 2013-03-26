\name{loadChr}
\alias{loadChr}
\title{Load downloaded chromosome into a GenomicRanges object}
\usage{
  loadChr(genome, chr, local)
}
\arguments{
  \item{genome}{The UCSC string specific to the genome to
  be downloaded (e.g. "hg19", "hg19", "mm10", etc)}

  \item{chr}{The string specific to the chromosome to load
  (e.g. "chr1", "chr22", "chr_M", "chr1_random", etc)}

  \item{local}{The root of a local directory where all
  files downloaded from UCSC are stored}
}
\value{
  A \code{DNAStringSet} object with one sequence for the
  chromosome from UCSC's FASTA release for the given
  genome.
}
\description{
  After downloading the FASTA for a genome using
  \code{\link{syncFasta}}, this function loads a given
  chromosome only into a \code{DNAStringSet} object.
}

