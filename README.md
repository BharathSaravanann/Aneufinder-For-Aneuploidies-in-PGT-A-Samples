This is a tools for Finding aneuploidies in PGT-A samples, Using HMM model it will find the reads in per chr area to find is there a gain or loss of chromosomes

This tools will conclude whether is that a Monosomy (1-somy), Trisomy (3-somy), Nullisomy (0-somy), Normal (2-somy) based on its state using HMM model

To Run Aneufiner

Rscript Aneufinder.R Bam Output_location

To Generate stats from Model

Rscript Report.R Model/RData Output_location
