# #############################################################################
# Download, update, and parse UCSC genome browser tables
# Author: Jeffrey Bhasin <jeffb@case.edu>
# Created: 2013-03-26
# #############################################################################

# =============================================================================
# Packages and Globals

#library(GenomicRanges)
#library(Biostrings)
#library(stringr)

url.default <- "rsync://hgdownload.cse.ucsc.edu/goldenPath/"
# =============================================================================

# =============================================================================
# Sequence Functions

# -----------------------------------------------------------------------------
#' Sync whole genome FASTA from UCSC's server
#'
#' This function will create the necessary directory tree under \code{local}, download, and extract the whole genome FASTA for the genome given as \code{genome}. If the FASTA archive already exists, the function will only re-download and re-extract if the modified dates are different between the local and remote copies.
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param local The root of a local directory where all files downloaded from UCSC are stored
#' @param url The root of the remote rsync URL to download UCSC data from (set by default to \code{rsync://hgdownload.cse.ucsc.edu/goldenPath/})
#' @return Prints status of the sync to the screen. No output object.
#' @export
syncFasta <- function(genome="",local, url=url.default)
{
	# Only works for genomes with "chromFa.tar.gz" files available

	# Account for the fact that on hg18 it is a .zip file
	ext <- ".tar.gz"
	if(genome=="hg18") {ext <- ".zip"}

	# Full url to rsync from
	url.dl <- paste(url, genome , "/bigZips/chromFa",ext, sep="")

	# Create genome directory if it does not exist
	local.dir1 <- paste(local, genome, sep="")
	dir.create(local.dir1, showWarnings = FALSE)

	# Create bigZips directory if it does not exist
	local.dir2 <- paste(local, genome, "/bigZips/", sep="")
	dir.create(local.dir2, showWarnings = FALSE)

	# Full path of file to save locally
	local.file <- paste(local, genome, "/bigZips/chromFa",ext, sep="")

	# Save file modified time - this is how we'll check if rsync updated the file or not
	file.mtime1 <- file.info(local.file)$mtime

	print(paste("Checking local file ",local.file, " against remote file ", url.dl,sep=""))
	system(paste("rsync -avzhP ", url.dl," ",local.file, sep=""))

	# Recheck file mod time
	file.mtime2 <- file.info(local.file)$mtime

	if((file.mtime1 != file.mtime2) | (is.na(file.mtime1)))
	{
		print(paste("New version downloaded: ", file.mtime1, " -> ", file.mtime2,sep=""))

		# Extract the new contents
		extract.dir <- paste(local.dir2, "chromFa/", sep="")
		print(paste("Extracting new genome to ", extract.dir, sep=""))
		unlink(extract.dir)
		dir.create(extract.dir, showWarnings = FALSE)
		
		if(genome=="hg18")
		{
			system(paste("unzip ", local.file, " -d ", extract.dir,sep=""))
		} else
		{
			system(paste("tar -zxvf ", local.file, " -C ", extract.dir,sep=""))
		}
	} else
	{
		print(paste("Version retained: ", file.mtime1, " == ", file.mtime2,sep=""))
	}
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Load downloaded whole genome FASTA into a GenomicRanges object
#'
#' After downloading the FASTA for a genome using \code{\link{syncFasta}}, this function loads all chromosomes into a \code{DNAStringSet} object.
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param local The root of a local directory where all files downloaded from UCSC are stored
#' @return A \code{DNAStringSet} object with one sequence per chromosome in UCSC's FASTA release for the given genome.
#' @export
loadGenome <- function(genome, local)
{
	local.dir <- paste(local, genome, "/bigZips/chromFa", sep="")

	local.fasta <- list.files(local.dir,full.names=TRUE)

	readDNAStringSet(local.fasta, format="fasta")
}
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
#' Load downloaded chromosome into a GenomicRanges object
#'
#' After downloading the FASTA for a genome using \code{\link{syncFasta}}, this function loads a given chromosome only into a \code{DNAStringSet} object.
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param chr The string specific to the chromosome to load (e.g. "chr1", "chr22", "chr_M", "chr1_random", etc)
#' @param local The root of a local directory where all files downloaded from UCSC are stored
#' @return A \code{DNAStringSet} object with one sequence for the chromosome from UCSC's FASTA release for the given genome.
#' @export
loadChr <- function(genome, chr, local)
{
	local.fasta <- paste(local, genome, "/bigZips/chromFa/", chr, ".fa",sep="")

	readDNAStringSet(local.fasta, format="fasta")
}
# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Annotation Functions

# -----------------------------------------------------------------------------
#' Sync whole genome annotation tables from UCSC's server
#'
#' This function will create the necessary directory tree under \code{local}, download, and extract the text file data and schema for given table from the given genome. If the data already exists, the function will only re-download and re-extract if the modified dates are different between the local and remote copies.
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param table The UCSC string specific for the table to sync (e.g. "knownGene", "kgXref", etc)
#' @param local The root of a local directory where all files downloaded from UCSC are stored
#' @param url The root of the remote rsync URL to download UCSC data from (set by default to \code{rsync://hgdownload.cse.ucsc.edu/goldenPath/})
#' @return Prints status of the sync to the screen. No output object.
#' @export
syncTable <- function(genome,table,local,url=url.default)
{
	# Full url to rsync from - the data flatfile and the SQL schema header
	url.dl.txt <- paste(url, genome , "/database/",table,".txt.gz",sep="")
	url.dl.sql <- paste(url, genome , "/database/",table,".sql",sep="")

	# Create genome directory if it does not exist
	local.dir1 <- paste(local, genome, sep="")
	dir.create(local.dir1, showWarnings = FALSE)

	# Create database directory if it does not exist
	local.dir2 <- paste(local, genome, "/database/", sep="")
	dir.create(local.dir2, showWarnings = FALSE)

	# Full path of files to save locally
	local.file.txt <- paste(local, genome, "/database/", table, ".txt.gz", sep="")
	local.file.sql <- paste(local, genome, "/database/", table, ".sql", sep="")

	# Save file modified time - this is how we'll check if rsync updated the file or not
	file.mtime1.txt <- file.info(local.file.txt)$mtime

	print(paste("Checking local file ",local.file.txt, " against remote file ", url.dl.txt,sep=""))
	system(paste("rsync -avzhP ", url.dl.txt," ",local.file.txt, sep=""))

	# Recheck file mod time
	file.mtime2.txt <- file.info(local.file.txt)$mtime

	# Save file modified time - this is how we'll check if rsync updated the file or not
	file.mtime1.sql <- file.info(local.file.sql)$mtime

	print(paste("Checking local file ",local.file.sql, " against remote file ", url.dl.sql,sep=""))
	system(paste("rsync -avzhP ", url.dl.sql," ",local.file.sql, sep=""))

	# Recheck file mod time
	file.mtime2.sql <- file.info(local.file.sql)$mtime

	if( ((file.mtime1.txt != file.mtime2.txt) | (is.na(file.mtime1.txt))) | ((file.mtime1.sql != file.mtime2.sql) | (is.na(file.mtime1.sql))) )
	{
		print(paste("New TXT.GZ (data) version downloaded: ", file.mtime1.txt, " -> ", file.mtime2.txt,sep=""))
		print(paste("New SQL (schema) version downloaded: ", file.mtime1.sql, " -> ", file.mtime2.sql,sep=""))

		# Extract the new contents
		local.file.txt.un <- paste(local, genome, "/database/", table, ".txt", sep="")
		print(paste("Extracting new TXT file to ", local.file.txt.un, sep=""))
		system(paste("gunzip -c ", local.file.txt, " > ", local.file.txt.un ,sep=""))
	} else
	{
		print(paste("TXT.GZ (data) version retained: ", file.mtime1.txt, " == ", file.mtime2.txt,sep=""))
		print(paste("SQL (schema) version retained: ", file.mtime1.sql, " == ", file.mtime2.sql,sep=""))
	}

}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Load downloaded annotation table into a data.frame object
#'
#' After downloading the table for a genome using \code{\link{syncTable}}, this function loads all data into a data.frame object. The headers for the table are extracted from the accompanying SQL schema file.
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param table The UCSC string specific for the table to sync (e.g. "knownGene", "kgXref", etc)
#' @param local The root of a local directory where all files downloaded from UCSC are stored
#' @return A data.frame object with UCSC's column headers for the specified table.
#' @export
loadTable <- function(genome, table, local)
{
	local.file.txt.un <- paste(local, genome, "/database/", table, ".txt", sep="")
	local.file.sql <- paste(local, genome, "/database/", table, ".sql", sep="")

	# First load in the table data
	txt <- read.table(file=local.file.txt.un, comment.char="", header=FALSE, stringsAsFactors=FALSE, sep="\t", quote="")
	txt.nCols <- ncol(txt)

	# Parse SQL schema to get the row names
	# Parsing by the backticks worked, until I noticed some of the older tables didn't have them
	# Now using a line by line reader
	#sql <- readChar(local.file.sql, file.info(local.file.sql)$size)

	# Want all strings encased in backquotes
	#ex <- unlist(str_extract_all(sql,"`.*?`"))
	#rp <- str_replace_all(ex,"`","")

	# Assuming first 3 are the table name, then taking the rest equal to the number of columns
	#mycols <- rp[(1:txt.nCols)+3]

	mycols <- getHeadersFromSQL(local.file.sql)

	# Set table with these row names
	names(txt) <- mycols

	txt
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Extract column names from UCSC's SQL schema files
#'
#' Parse UCSC's SQL files to extract column headers for the data files.
#' @param sql.file Path to the SQL file
#' @return A vector of column names.
getHeadersFromSQL <- function(sql.file)
{
	cols <- c()

	con  <- file(sql.file, open = "r")

	extract <- FALSE
	while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
	{
		linevec <- unlist(str_split(oneLine," "))
		#print(linevec)
		if((extract==TRUE)&((linevec[3]=="KEY")|(linevec[3]=="PRIMARY")|(linevec[3]=="UNIQUE")))
		{
			extract <- FALSE
		}
		if(extract==TRUE)
		{
			#print(linevec[3])
			colname <- linevec[3]
			colname <- str_replace_all(linevec[3],"`","")
			cols <- c(cols, colname)
		}
		if(linevec[1]=="CREATE")
		{
			extract <- TRUE
		}
	} 

	close(con)
}
# -----------------------------------------------------------------------------


# =============================================================================
