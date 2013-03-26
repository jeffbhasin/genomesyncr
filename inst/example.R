library(genomesyncr)

path.local <- "~/Documents/myucsc/"

# Fetch some genomes
syncFasta(genome="mm9", path.local)
syncFasta(genome="mm10", path.local)
syncFasta(genome="hg18", path.local)
syncFasta(genome="hg19", path.local)

# Load genomic data
hg18 <- loadGenome(genome="hg18", path.local)
hg18.22 <- loadChr("hg18","chr22", path.local)
mm10.1 <- loadChr("mm10", "chr1", path.local)

# Fetch some tables
syncTable(genome="hg19", table="knownGene", path.local)
kg <- loadTable(genome="hg19", table="knownGene", path.local)