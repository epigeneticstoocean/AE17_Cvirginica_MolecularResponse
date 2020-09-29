### Script checking number of OA-induced DML in mantle tissue (this study) that overlap with DML in gonadal tissue

# Directory path
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/")

### DML files from Venkataraman et al. (2020)
# Github repo : https://github.com/epigeneticstoocean/paper-gonad-meth
# Bed file with 598 DML 
ym <- read.delim("results/DNAm/2019-04-05-DML-Destrand-5x-Locations.tab",
                 sep="\t",header = FALSE)
# Annotation file
ym_annot <- read.csv("results/DNAm/2019-06-20-Master-DML-Annotation.csv")

## DML from mantle tissue
dml <- read.csv("results/manuscript/Supp/Supplemental_TableS3.3__DMLlist.csv")
# Extract only OA-induced DML
dml_OA <- dml[dml$Comparison == "Treatment" | dml$Comparison == "Treatment - Day 80" | dml$Comparison == "Treatment - Day 9",]

# Look at overlapping DML among tissues
out <- paste0(dml_OA$CpG.Position,"_",dml_OA$X) %in% paste0(ym$V1,"_",ym$V2)
dml_overlap <- dml_OA[out,]
# Retrieve Annotations from Venkataraman et al. (2020) data to confirm they are functionally identical
out <-  paste0(ym_annot$chr,"_",ym_annot$start) %in% paste0(dml_overlap$CpG.Position,"_",as.numeric(as.character(dml_overlap$X))+1)
ym_annot_overlap <- ym_annot[out,]
