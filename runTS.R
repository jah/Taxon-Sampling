#-----------------------------------------------------------------------------#
# Parameters for the analysis, uncomment and edit them if you prefer to leave
# it automated for future uses.

idsFile <- "inputs/test1"
multifasta <- "inputs/multifasta"
taxondir <- "../KOMODO2/taxdump"

m <- 100
method <- "diversity"
randomize <- "no"
replacement <- "no"
ignoreIDs <- NULL
requireIDs <- NULL
ignoreNonLeafID <- NULL
outFile <- "output"

#-----------------------------------------------------------------------------#

source("TaxonSampling.R")

#library("CHNOSZ")
#library("ape")

nodes <- suppressMessages(getnodes(taxondir))
countIDs <- TS_TaxonomyData(idsFile, nodes)

nodes <- Simplify_Nodes(nodes, countIDs)
outputIDs <- TS_Algorithm(1, m, nodes, countIDs, method, randomize,
                          replacement, ignoreIDs, requireIDs,
                          ignoreNonLeafID)

WriteFasta(idsFile, multifasta, outputIDs, outFile)
