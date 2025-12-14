.libPaths("/root/RnD/Aneufinder/R_libs/")
library(AneuFinder)

args <- commandArgs(trailingOnly=TRUE)
bam_file <- args[1]
output_dir <- args[2]

Aneufinder(inputfolder = dirname(bam_file),
           outputfolder = output_dir,
           assembly = "hg19",
           binsizes = c(1e6),
           stepsizes = c(1e6),
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           method = "HMM",
           correction.method = "GC",
           GC.BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
           num.trials = 10)
