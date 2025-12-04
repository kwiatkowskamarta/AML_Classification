# SKRYPT 03: Feature Selection (Boruta)
# ------------------------------------------
# OPIS: Wybiera najważniejsze geny (biomarkery) spośród 20,000 dostępnych.
# INPUT: data/processed/TCGA_LAML_Cleaned.RData
# OUTPUT: results/selected_features.rds

# 1. Wczytanie przetworzonych danych

# 2. Wstępna redukcja (opcjonalnie: usunięcie skorelowanych zmiennych)

# 3. Uruchomienie algorytmu Boruta (to może potrwać!)
#    - Trenowanie na całym zbiorze lub tylko treningowym

# 4. Wyświetlenie wyników (wykres ważności cech)

# 5. Zapisanie listy wybranych genów (np. top 50-100)
