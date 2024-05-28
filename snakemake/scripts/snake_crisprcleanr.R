suppressMessages(library('CRISPRcleanR'))

# Load the data
cell_line <- snakemake@wildcards['line']
libraryAnnotation <- read.table(snakemake@input[['guides']], header = TRUE, sep = ",", row.names = 1)
counts_file <- snakemake@input[['counts']]

outdir <- snakemake@params[['outdir']]
min.ngenes <- snakemake@params[['min_ngenes']]  # default=3
min.width <- snakemake@params[['min_width']]  # default=2

## Loading, median-normalizing and computing fold-changes
normANDfcs <- ccr.NormfoldChanges(
    filename = counts_file,
    display = TRUE,
    min_reads = 30,
    ncontrols = 3,  # We have 3 plasmid pools for controls
    EXPname = cell_line,
    saveToFig = TRUE,
    outdir = outdir,
    libraryAnnotation = libraryAnnotation
)

## Genome-sorting the fold changes
gwSortedFCs <- ccr.logFCs2chromPos(
    foldchanges = normANDfcs$logFCs,
    libraryAnnotation = libraryAnnotation
)

## Identifying and correcting biased sgRNAs' fold changes
correctedFCs_and_segments <- suppressWarnings(ccr.GWclean(
    gwSortedFCs = gwSortedFCs,
    display = TRUE,
    label = cell_line,
    saveTO = outdir,
    min.ngenes=min.ngenes,
    min.width=min.width,
))

## Write the corrected fold changes to file
write.csv(correctedFCs_and_segments$corrected_logFCs, file = snakemake@output[['corrlogFC']])
write.csv(normANDfcs$logFCs, file = snakemake@output[['logFC']])
write.csv(correctedFCs_and_segments$segments, file = snakemake@output[['segments']])
write.csv(correctedFCs_and_segments$SORTED_sgRNAs, file = snakemake@output[['sortedsgRNAs']])

print("Done")
