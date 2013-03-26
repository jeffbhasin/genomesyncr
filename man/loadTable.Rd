\name{loadTable}
\alias{loadTable}
\title{Load downloaded annotation table into a data.frame object}
\usage{
  loadTable(genome, table, local)
}
\arguments{
  \item{genome}{The UCSC string specific to the genome to
  be downloaded (e.g. "hg19", "hg19", "mm10", etc)}

  \item{table}{The UCSC string specific for the table to
  sync (e.g. "knownGene", "kgXref", etc)}

  \item{local}{The root of a local directory where all
  files downloaded from UCSC are stored}
}
\value{
  A data.frame object with UCSC's column headers for the
  specified table.
}
\description{
  After downloading the table for a genome using
  \code{\link{syncTable}}, this function loads all data
  into a data.frame object. The headers for the table are
  extracted from the accompanying SQL schema file.
}

