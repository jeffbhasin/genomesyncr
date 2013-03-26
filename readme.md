# genomesyncr: UCSC Genome Browser Data Synchronizer for R

## Goals
* Provide reproducible interface to access most commonly used sequence and annotation data from UCSC
* Download and sync table flatfiles from UCSC's FTP
* Keep UCSC's folder structure: A partial local mirror of the tables you need
* Provide some faculty for parsing the tables based on their schemas and creating GenomicRanges objects from them
* Provide some common joins on KnownGene to add gene symbol, cytoBand, etc (may be duplicating GenomicFeatures)

## Usage
In development.

## Proposed Workflow
* Set directory of goldenPath mirror
* One function to load sequence data into a GRanges object for any species/build
	* Download if not already available
* One function to load annotation data into a GRanges object
	* Check if it exists in our local mirror
	* If it does, check with rsync that it hasn't been updated
	* If updated/never downloaded, pull down new version with rsync
	* Load and parse the table into a GRanges