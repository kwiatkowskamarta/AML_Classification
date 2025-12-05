# SKRYPT 01: Data Download (TCGA-LAML)
# ------------------------------------------------------------------------------
# Pobranie danych klinicznych i transkryptomicznych (RNA-seq)
# Źródło: GDC (Genomic Data Commons) via TCGAbiolinks
# ------------------------------------------------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)

message(">>> ROZPOCZYNANIE POBIERANIA DANYCH DLA PROJEKTU TCGA-LAML <<<")

# Utworzenie folderu na pliki tymczasowe
if(!dir.exists("GDCdata")) dir.create("GDCdata")

# --- 1. DANE KLINICZNE (Etykiety / Y) ---
message("1. Pobieranie danych klinicznych.")

query_clin <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

# Pobieranie do folderu GDCdata
GDCdownload(query_clin, directory = "GDCdata")
clinical_data <- GDCprepare(query_clin, directory = "GDCdata")

# Główną tabela pacjentów
clinical_patient <- clinical_data$clinical_patient_laml

message(paste("   Pobrano dane kliniczne dla", nrow(clinical_patient), "pacjentów."))

# --- 2. DANE RNA-SEQ (Cechy / X) ---
message("2. Pobieranie danych RNA-seq (Gene Expression Quantification).")

query_rna <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
)

# Trybu 'api' i mniejsze paczek (chunks), aby uniknąć błędów sieciowych
GDCdownload(query_rna, 
            method = "api", 
            files.per.chunk = 10,
            directory = "GDCdata")

data_rna <- GDCprepare(query_rna, directory = "GDCdata")

message(paste("   Pobrano macierz ekspresji:", nrow(data_rna), "genów x", ncol(data_rna), "próbek."))

# --- 3. ZAPIS SUROWYCH DANYCH ---
message("3. Zapisywanie danych do folderu data/raw/.")

if(!dir.exists("data/raw")) dir.create("data/raw", recursive = TRUE)

save(data_rna, clinical_patient, file = "data/raw/TCGA_LAML_Raw.RData")

message(">>> Dane zostały zapisane w pliku 'data/raw/TCGA_LAML_Raw.RData'. <<<")