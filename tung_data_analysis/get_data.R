##### This file is written based on the Github folder (https://github.com/jdblischak/singleCellSeq) that the original authors for the data provided. You need to first download their folder to get the data
## need to put your own directory for the data folder 
DATADIR = "singleCellSeq/data/"

anno <- read.table(paste(DATADIR, "annotation.txt", sep = ""), 
                   header = TRUE,
                   stringsAsFactors = FALSE)
molecules <- read.table(psate(DATADIR, "molecules.txt", sep = ""), 
                        header = TRUE,
                        stringsAsFactors = FALSE)
quality_single_cells <- scan(paste(DATADIR, "quality-single-cells.txt", sep = ""),
                             what = "character")

## Prepare single cell molecule data following the same code from the original authors
molecules_single <- molecules[, colnames(molecules) %in% quality_single_cells]
anno_single <- anno[anno$sample_id %in% quality_single_cells, ]
## remove one replicate with bad quality
molecules_single <- molecules_single[, !(anno_single$individual == 19098 & anno_single$batch == 2)]
anno_single <- anno_single[!(anno_single$individual == 19098 & anno_single$batch == 2), ]

expressed_single <- rowSums(molecules_single) > 0
molecules_single <- molecules_single[expressed_single, ]
overexpressed_genes <- rownames(molecules_single)[apply(molecules_single, 1,
                                                        function(x) any(x >= 1024))]
molecules_single <- molecules_single[!(rownames(molecules_single) %in% overexpressed_genes), ]

ercc_rows_molecules_single <- grep("ERCC", rownames(molecules_single))
bio.g <- molecules_single[-ercc_rows_molecules_single, ]

## the data is included in the folder 
saveRDS(bio.g, file = "data.rds")


