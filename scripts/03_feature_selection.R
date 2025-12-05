# SKRYPT 03: Feature Selection (Boruta) - FULL VERSION
# ------------------------------------------------------------------------------
# Cel: Identyfikacja biomarkerów genetycznych bez wstępnej redukcji.
# Metoda: Algorytm Boruta (All-relevant feature selection).
#         Analiza wszystkich ~27,000 genów.
# ------------------------------------------------------------------------------

library(Boruta)
library(dplyr)
library(caret)

# Ustawienie ziarna losowości dla powtarzalności wyników (Reproducibility)
set.seed(42)

message(">>> SELEKCJA CECH (PEŁNA ANALIZA BORUTA) <<<")

# 1. WCZYTANIE DANYCH
load("data/processed/TCGA_LAML_Cleaned.RData")

message("Wczytano dane.")
message("   Liczba pacjentów: ", nrow(ml_df))
message("   Liczba genów do sprawdzenia: ", ncol(ml_df) - 1)

# 2. PRZYGOTOWANIE DANYCH (X i Y)
# Rozdzielamy macierz na cechy i cel. Używamy składni x/y zamiast formuły,
# ponieważ przy 27k kolumnach formuła R (Target ~ .) jest bardzo wolna.
target <- ml_df$Target
features <- ml_df %>% select(-Target)

# 3. KONFIGURACJA I URUCHOMIENIE BORUTY
message("\nRozpoczynam obliczenia algorytmem Boruta.")
message("To jest pełna analiza. Może potrwać od 30 min do kilku godzin.")
message("Postęp będzie wyświetlany poniżej.")

# maxRuns = 200 -> Zwiększamy liczbę iteracji, bo przy tylu genach algorytm
# potrzebuje więcej czasu, żeby zdecydować o "niepewnych" (Tentative) atrybutach.
# doTrace = 2 -> Szczegółowy podgląd postępu
boruta_output <- Boruta(
  x = features,
  y = target,
  doTrace = 2,
  maxRuns = 200 
)

message("\n--- KONIEC OBLICZEŃ ---")

# 4. DIAGNOSTYKA WYNIKÓW
print(boruta_output)

# Naprawa atrybutów "Tentative" (Niepewnych)
# Jeśli po 200 rundach Boruta nadal nie jest pewna co do niektórych genów,
# zmuszamy ją do podjęcia decyzji (TentativeRoughFix).
if(any(boruta_output$finalDecision == "Tentative")) {
  message("Znaleziono atrybuty niepewne (Tentative). Dokonuję ostatecznej decyzji...")
  boruta_output <- TentativeRoughFix(boruta_output)
  print(boruta_output)
}

# 5. EKSTRAKCJA KLUCZOWYCH GENÓW
# Wyciągamy tylko te, które zostały potwierdzone (Confirmed)
final_genes <- getSelectedAttributes(boruta_output, withTentative = FALSE)
message("\nLiczba zidentyfikowanych biomarkerów: ", length(final_genes))

# Wyświetlenie Top 10 najważniejszych (wg Importance)
# (To przyda się do tekstu pracy)
imps <- attStats(boruta_output)
imps_confirmed <- imps[imps$decision == "Confirmed", ]
top_genes <- row.names(imps_confirmed)[order(imps_confirmed$meanImp, decreasing = TRUE)]

message("Top 10 najważniejszych genów:")
print(head(top_genes, 10))

# 6. ZAPIS WYNIKÓW
# Zapisujemy pełny obiekt Boruta (do wykresów) oraz odchudzony dataset do modeli
ml_data_final <- ml_df %>% 
  select(all_of(final_genes), Target)

if(!dir.exists("results/models")) dir.create("results/models", recursive = TRUE)

save(boruta_output, ml_data_final, final_genes, file = "results/models/Boruta_Results.RData")

message(">>> Wyniki zapisane w results/models/Boruta_Results.RData <<<")