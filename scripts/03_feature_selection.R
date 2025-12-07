# SKRYPT 03: Feature Selection (Boruta) - FULL VERSION
# ------------------------------------------------------------------------------
# Cel: Identyfikacja biomarkerów genetycznych bez wstępnej redukcji.
# ------------------------------------------------------------------------------

library(Boruta)
library(dplyr)
library(caret)

set.seed(42)

message(">>> SELEKCJA CECH (PEŁNA ANALIZA BORUTA) <<<")

# 1. WCZYTANIE DANYCH
load("data/processed/TCGA_LAML_Cleaned.RData")

message("Wczytano dane.")
message("   Liczba pacjentów: ", nrow(ml_df))
message("   Liczba genów do sprawdzenia: ", ncol(ml_df) - 1)

# 2. PRZYGOTOWANIE DANYCH (X i Y)
target <- ml_df$Target
features <- ml_df %>% select(-Target)

# 3. KONFIGURACJA I URUCHOMIENIE BORUTY
message("\nRozpoczecie obliczen algorytmem Boruta.")

boruta_output <- Boruta(
  x = features,
  y = target,
  doTrace = 2,
  maxRuns = 200 
)

message("\n--- KONIEC OBLICZEŃ ---")

# 4. DIAGNOSTYKA WYNIKÓW
print(boruta_output)

if(any(boruta_output$finalDecision == "Tentative")) {
  message("Znaleziono atrybuty niepewne (Tentative).")
  boruta_output <- TentativeRoughFix(boruta_output) # jeśli po 200 rundach metoda Boruta nadal nie jest pewna co do niektórych genów, zmuszamy ją do podjęcia decyzji (TentativeRoughFix).
  print(boruta_output)
}

# 5. EKSTRAKCJA KLUCZOWYCH GENÓW
# tylko te, które zostały potwierdzone (Confirmed)
final_genes <- getSelectedAttributes(boruta_output, withTentative = FALSE)
message("\nLiczba zidentyfikowanych biomarkerów: ", length(final_genes))

# Wyświetlenie Top 10 najważniejszych (wg importance)
imps <- attStats(boruta_output)
imps_confirmed <- imps[imps$decision == "Confirmed", ]
top_genes <- row.names(imps_confirmed)[order(imps_confirmed$meanImp, decreasing = TRUE)]

message("Top 10 najważniejszych genów:")
print(head(top_genes, 10))

# 6. ZAPIS WYNIKÓW
ml_data_final <- ml_df %>% 
  select(all_of(final_genes), Target)

if(!dir.exists("results/models")) dir.create("results/models", recursive = TRUE)

save(boruta_output, ml_data_final, final_genes, file = "results/models/Boruta_Results.RData")

message(">>> Wyniki zapisane w results/models/Boruta_Results.RData <<<")
