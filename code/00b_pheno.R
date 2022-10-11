#00b_pheno.R

####----- read pheno file

pheno <- as.data.frame(read_excel("../suppl/pheno/pSS_Pheno.xlsx"))

pheno <- pheno[, c("Single cell_Sample_id_publ",
                   "Age symptom onset",
                   "Age diagnosis",
                   "First symptoms",
                   "ANA",
                   "SSA Ro52",
                   "SSA Ro60",
                   "SSA",
                   "SSB",
                   "SSA and SSB",
                   "RNP",
                   "Sm",
                   "RF",
                   "Focus score Greenspan",
                   "Anemia Hb<120g/L",
                   "Leucopenia <4.0x109/L",
                   "Lymphopenia <1.0x109/L",
                   "P-IgG >15g/L",
                   "Raynaud",
                   "Arthritis",
                   "Purpura",
                   "Major salivary gland swelling",
                   "Hypothyreoidism",
                   "Treatment at B cell sampling, HCQ",
                   "Prednisolone")]

colnames(pheno)[1] <- "orig.ident"

pheno <- pheno %>% add_row(orig.ident = c("C001", "C002", "C003", "C004")) %>% arrange(orig.ident)
colnames(pheno) <- gsub(" Hb<120g/L", "", colnames(pheno))
colnames(pheno) <- gsub(" <4.0x109/L", "", colnames(pheno))
colnames(pheno) <- gsub(" <1.0x109/L", "", colnames(pheno))
colnames(pheno) <- gsub(" >15g/L", "", colnames(pheno))
colnames(pheno) <- gsub(", HCQ", " HCQ", colnames(pheno))
colnames(pheno) <- gsub(" ", "_", colnames(pheno))

#split first symptoms column
pheno$arthralgia_firstSymptom <- ifelse(grepl("arthralgia", pheno$First_symptoms), "YES", "NO")
pheno$arthritis_firstSymptom <- ifelse(grepl("arthritis", pheno$First_symptoms), "YES", "NO")
pheno$hemiparesis_firstSymptom <- ifelse(grepl("hemiparesis", pheno$First_symptoms), "YES", "NO")
pheno$myalgia_firstSymptom <- ifelse(grepl("myalgia", pheno$First_symptoms), "YES", "NO")
pheno$sicca_firstSymptom <- ifelse(grepl("sicc", pheno$First_symptoms), "YES", "NO")
pheno$erythema_nodosum_firstSymptom <- ifelse(grepl("erythema", pheno$First_symptoms), "YES", "NO")
pheno$major_salivary_gland_swelling_firstSymptom <- ifelse(grepl("salivary", pheno$First_symptoms), "YES", "NO")
pheno$Raynaud_firstSymptom <- ifelse(grepl("Raynaud", pheno$First_symptoms), "YES", "NO")
pheno$pleuritis_entesitis_tendinitis_firstSymptom <- ifelse(grepl("pleuritis", pheno$First_symptoms), "YES", "NO")
#pheno <- pheno[, c(-4)] #remove First_symptoms column

pheno[, 4:ncol(pheno)][is.na(pheno[, 4:ncol(pheno)])] <- "NO"

#add patient group
pheno$patient_group1 <- paste0(ifelse(pheno$SSA == "YES", "SSA", ""),
                              "_",
                              ifelse(pheno$SSB == "YES", "SSB", ""))
pat_annot <- c(SSAB = "SSA_SSB", SSA = "SSA_", DNEG = "NA_NA", DNEG = "_")
pheno$patient_group1 <- names(pat_annot)[match(pheno$patient_group1, pat_annot)]
pheno$patient_group1[grep("^C00", pheno$orig.ident) ] <- "CTRL"
pheno$patient_group1 <- factor(pheno$patient_group1, levels = c("CTRL", "DNEG", "SSA", "SSAB"))

pat_annot <- c(CTRL = "CTRL", `SSA-` = "DNEG", `SSA+` = "SSA", `SSAB` = "SSAB")
pheno$patient_group <- names(pat_annot)[match(pheno$patient_group1, pat_annot)]
pheno$patient_group <- factor(pheno$patient_group, levels = names(pat_annot))
rownames(pheno) <- pheno$orig.ident

saveRDS(pheno, paste0("../results/pheno_", Sys.Date(), ".rds"))


