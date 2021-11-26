#-------------------------------------------------------------------------------
#                           NSCLC data from UCSC Xena
#-------------------------------------------------------------------------------

# Clinical data
LUNG_clinicalMatrix <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/LUNG_clinicalMatrix", TRUE, "\t"))
rownames(LUNG_clinicalMatrix) <- LUNG_clinicalMatrix[, 1]
LUNG_clinicalMatrix <- LUNG_clinicalMatrix[, -1]
clinical <- LUNG_clinicalMatrix[, c("age_at_initial_pathologic_diagnosis", 
                                    "pathologic_T", "number_pack_years_smoked", 
                                    "pre_bronchodilator_fev1_percent", "gender",
                                    "X_primary_disease")]
colnames(clinical) <- c("Age", "StageT", "PYS", "FEV1", "Gender", "LUAD")
clinical[, "StageT"] <- as.numeric(sapply(1:nrow(clinical), function(x) 
    unlist(strsplit(clinical[, "StageT"][x], ""))[2]))
clinical[which(clinical[, "StageT"] == 1), "StageT"] <- 0
clinical[which(clinical[, "StageT"] == 2 | clinical[, "StageT"] == 3 | 
                   clinical[, "StageT"] == 4), "StageT"] <- 1
clinical[which(clinical[, "LUAD"] == "lung adenocarcinoma"), "LUAD"] <- 1
clinical[which(clinical[, "LUAD"] == "lung squamous cell carcinoma"), "LUAD"] <- 0
clinical[which(clinical[, "Gender"] == "FEMALE"), "Gender"] <- 0
clinical[which(clinical[, "Gender"] == "MALE"), "Gender"] <- 1
clinical <- apply(clinical, 2, function(x) as.numeric(x))
rownames(clinical) <- rownames(LUNG_clinicalMatrix)
clinical <- clinical[- which(apply(clinical, 1, function(x) any(is.na(x)))), ]

# Gene expression data
HiSeqV2_LUNG <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/HiSeqV2_LUNG", TRUE, "\t"))
rownames(HiSeqV2_LUNG) <- HiSeqV2_LUNG[, 1]
HiSeqV2_LUNG <- HiSeqV2_LUNG[, -1]
colnames(HiSeqV2_LUNG) <- sapply(1:ncol(HiSeqV2_LUNG), function(x) 
    paste(unlist(strsplit(colnames(HiSeqV2_LUNG)[x], "\\.")), collapse = "-"))

# Pack sample with both gene and clinical information available
patient.id <- intersect(colnames(HiSeqV2_LUNG), rownames(clinical))
patient.id <- sort(patient.id)
patient.id <- patient.id[-which(duplicated(sapply(1:length(patient.id), function(x)
    paste(unlist(strsplit(patient.id[x], "-"))[1:3], collapse = "-"))))]
HiSeqV2_LUNG <- HiSeqV2_LUNG[, patient.id]
clinical <- clinical[patient.id, ]
gene <- apply(HiSeqV2_LUNG, 1, as.numeric)
rownames(gene) <- patient.id
gene <- gene[, -which(apply(gene, 2, function(x) (sum(x == 0) / length(x)) > 0.3))]

# Marginal screening
pv <- sapply(1:ncol(gene), function(x) 
    summary(lm(clinical[, "FEV1"] ~
                   gene[, x] + clinical[, c("Age", "PYS", "StageT", "LUAD",
                                            "Gender")]))$coefficients["gene[, x]", "Pr(>|t|)"])
names(pv) <- colnames(gene)
gene.reduced <- gene[, order(pv, decreasing = FALSE)[1:300]]

write.csv(clinical, "./RealDataAnalysis/RealDataAnalysisData/UCSC_NSCLC_FEV1_clinical.csv")
write.csv(gene.reduced, "./RealDataAnalysis/RealDataAnalysisData/UCSC_NSCLC_FEV1_gene_screen.csv")


#-------------------------------------------------------------------------------
#                            LGG data from UCSC Xena
#------------------------------------------------------------------------------- 

library(survival)

# Clinical data
LGG_clinicalMatrix <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/LGG_clinicalMatrix", TRUE, "\t"))
rownames(LGG_clinicalMatrix) <- LGG_clinicalMatrix[, 1]
LGG_clinicalMatrix <- LGG_clinicalMatrix[, -1]
clinical <- LGG_clinicalMatrix[, c("age_at_initial_pathologic_diagnosis", 
                                   "gender", "neoplasm_histologic_grade")]
colnames(clinical) <- c("Age", "Gender", "Grade")
clinical[which(clinical[, "Gender"] == "FEMALE"), "Gender"] <- 0
clinical[which(clinical[, "Gender"] == "MALE"), "Gender"] <- 1
clinical[which(clinical[, "Grade"] == "G2"), "Grade"] <- 0
clinical[which(clinical[, "Grade"] == "G3"), "Grade"] <- 1
clinical <- apply(clinical, 2, function(x) as.numeric(x))
rownames(clinical) <- rownames(LGG_clinicalMatrix)
clinical <- clinical[- which(apply(clinical, 1, function(x) any(is.na(x)))), ]

# Survival data
LGG_survival <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/LGG_survival.txt", TRUE, "\t"))
rownames(LGG_survival) <- LGG_survival[, 1]
LGG_survival <- LGG_survival[, - (match(c("sample", "X_PATIENT"), 
                                        colnames(LGG_survival)))]
survival <- Surv(as.numeric(LGG_survival[, "OS.time"]), 
                 as.numeric(LGG_survival[, "OS"]))
names(survival) <- rownames(LGG_survival)
survival <- survival[-which(apply(survival, 1, function(x) 
    any(is.na(x) | x[1] == 0))), ]

# Protein expression data
RPPA_LGG <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/RPPA_LGG", TRUE, "\t"))
rownames(RPPA_LGG) <- RPPA_LGG[, 1]
RPPA_LGG <- RPPA_LGG[, -1]

# Gene expression data
HiSeqV2_LGG <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/HiSeqV2_LGG", TRUE, "\t"))
rownames(HiSeqV2_LGG) <- HiSeqV2_LGG[, 1]
HiSeqV2_LGG <- HiSeqV2_LGG[, -1]

colnames(RPPA_LGG) <- sapply(1:ncol(RPPA_LGG), function(x) 
    paste(unlist(strsplit(colnames(RPPA_LGG)[x], "\\.")), collapse = "-"))
colnames(HiSeqV2_LGG) <- sapply(1:ncol(HiSeqV2_LGG), function(x) 
    paste(unlist(strsplit(colnames(HiSeqV2_LGG)[x], "\\.")), collapse = "-"))

patient.id <- intersect(intersect(names(survival), rownames(clinical)), 
                        intersect(colnames(RPPA_LGG), colnames(HiSeqV2_LGG)))

survival <- survival[patient.id, ]
clinical <- clinical[patient.id, ]
protein <- RPPA_LGG[, patient.id]
gene <- HiSeqV2_LGG[, patient.id]

survival <- survival[-which(sapply(1:length(patient.id), function(x)
    unlist(strsplit(patient.id[x], "-"))[4]) == "02"), ]
clinical <- clinical[-which(sapply(1:length(patient.id), function(x)
    unlist(strsplit(patient.id[x], "-"))[4]) == "02"), ]
protein <- protein[, -which(sapply(1:length(patient.id), function(x)
    unlist(strsplit(patient.id[x], "-"))[4]) == "02")]
gene <- gene[, -which(sapply(1:length(patient.id), function(x)
    unlist(strsplit(patient.id[x], "-"))[4]) == "02")]
patient.id <- names(survival)

gene <- apply(gene, 1, as.numeric)
rownames(gene) <- patient.id
gene <- gene[, - which(apply(gene, 2, function(x) (sum(x == 0)/length(x)) > 0.3))]
protein <- apply(protein, 1, as.numeric)
rownames(protein) <- patient.id
protein <- protein[, - which(apply(protein, 2, function(x) any(is.na(x))))]

write.csv(clinical, "./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_clinical.csv")
write.csv(survival, "./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_survival.csv")
write.csv(protein, "./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_protein.csv")
write.csv(gene, "./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_gene.csv")
