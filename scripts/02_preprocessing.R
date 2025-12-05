# SKRYPT 02: Preprocessing (Czyszczenie i Normalizacja)
# ------------------------------------------------------------------------------
# Cel: Przygotowanie tabeli "ML-ready" z surowych danych.
# Input: data/raw/TCGA_LAML_Raw.RData
# Output: data/processed/TCGA_LAML_Cleaned.RData
# ------------------------------------------------------------------------------

library(SummarizedExperiment)
library(dplyr)
library(edgeR) 

message(">>> PRZETWARZANIE DANYCH <<<")

# 1. WCZYTANIE SUROWYCH DANYCH
load("data/raw/TCGA_LAML_Raw.RData")

# Wyciągamy macierz (geny x próbki)
counts_matrix <- assay(data_rna)
message("Wczytano macierz: ", nrow(counts_matrix), " genów x ", ncol(counts_matrix), " próbek.")

# 2. PRZYGOTOWANIE DANYCH KLINICZNYCH (TARGET)
target_col <- "fab_category"

# Sprawdzenie czy kolumna istnieje (zabezpieczenie)
if(!target_col %in% colnames(clinical_patient)) {
  stop(paste("BŁĄD: Nie znaleziono kolumny", target_col))
}

message("Czyszczenie danych klinicznych.")

# Wybieramy ID pacjenta i Target, a następnie czyścimy śmieci
clin_clean <- clinical_patient %>%
  select(bcr_patient_barcode, target = all_of(target_col)) %>%
  # Krok 1: Usunięcie metadanych (wiersze nagłówkowe z pliku Biotab)
  filter(!grepl("CDE_ID", target)) %>%            # Usuwa wiersz z ID
  filter(!grepl("morphology_code", target)) %>%   # Usuwa ten wiersz z długim napisem
  # Krok 2: Usunięcie braków diagnozy
  filter(target != "Not Classified") %>%
  filter(!is.na(target)) %>%
  # Krok 3: Ujednolicenie nazw (M0 Undifferentiated -> M0)
  mutate(target = ifelse(target == "M0 Undifferentiated", "M0", target))

message("Liczba pacjentów po wyczyszczeniu diagnozy: ", nrow(clin_clean))
print(table(clin_clean$target)) # Podgląd dla pewności

# 3. PAROWANIE PRÓBEK (GENY <-> PACJENCI)
# ID w RNA (np. TCGA-AB-2805-03A...) są dłuższe niż w klinice (TCGA-AB-2805).
rna_barcodes <- colnames(counts_matrix)
rna_patient_ids <- substr(rna_barcodes, 1, 12) 

# Wspólne ID
common_patients <- intersect(clin_clean$bcr_patient_barcode, rna_patient_ids)
message("Liczba pacjentów do sparowania (Match): ", length(common_patients))

if(length(common_patients) < 50) {
  stop("UWAGA: Zbyt mało wspólnych pacjentów! Sprawdź formaty ID.")
}

# Filtracja obu tabele
# a) Kliniczne - tylko ci, co mają RNA
clin_final <- clin_clean %>% 
  filter(bcr_patient_barcode %in% common_patients) %>%
  arrange(bcr_patient_barcode) # Sortujemy alfabetycznie

# b) Genetyczne - tylko ci, co mają diagnozę
match_indices <- match(common_patients, substr(colnames(counts_matrix), 1, 12))
counts_final <- counts_matrix[, match_indices]

# Ostateczne sprawdzenie kolejności
if(!all(substr(colnames(counts_final), 1, 12) == clin_final$bcr_patient_barcode)) {
  stop("BŁĄD KRYTYCZNY: Kolejność pacjentów się nie zgadza!")
}

# 4. FILTRACJA GENÓW (SZUM)
# Usuwamy geny, które mają bardzo niską ekspresję
# Zasada: zostawiamy geny, które mają > 10 zliczeń u przynajmniej 20% pacjentów
keep_genes <- rowSums(counts_final > 10) >= (0.2 * ncol(counts_final))
counts_filtered <- counts_final[keep_genes, ]

message("Filtracja genów: zredukowano z ", nrow(counts_final), " do ", nrow(counts_filtered))

# 5. NORMALIZACJA (LOG2 CPM)
message("Normalizacja danych (Log2 CPM).")
dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge) # TMM normalization
cpm_log <- cpm(dge, log = TRUE, prior.count = 1)

# 6. TRANSPOZYCJA (DLA ML)
# Zamieniamy miejscami: Wiersze=Pacjenci, Kolumny=Geny
ml_matrix <- t(cpm_log) 
ml_df <- as.data.frame(ml_matrix)

# Dodajemy kolumnę Target (Y)
ml_df$Target <- factor(clin_final$target)

# Czyścimy nazwy kolumn (geny)
colnames(ml_df) <- make.names(colnames(ml_df))

message("Gotowa tabela do ML: ", nrow(ml_df), " wierszy x ", ncol(ml_df), " kolumn.")

# 7. ZAPIS
if(!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(ml_df, file = "data/processed/TCGA_LAML_Cleaned.RData")

message(">>> Dane gotowe do uczenia maszynowego. <<<")