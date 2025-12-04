# SKRYPT 02: Data Preprocessing & Cleaning
# ------------------------------------------
# OPIS: Czyści dane, filtruje brakujące próbki i normalizuje ekspresję genów.
# INPUT: data/raw/TCGA_LAML_Raw.RData
# OUTPUT: data/processed/TCGA_LAML_Cleaned.RData

# 1. Wczytanie surowych danych

# 2. Przygotowanie macierzy X (geny) i wektora Y (podtypy)
#    - Usunięcie pacjentów bez zdefiniowanego podtypu (klasyfikacja FAB)
#    - Transpozycja macierzy (Pacjenci jako wiersze)

# 3. Filtracja genów
#    - Usunięcie genów o zerowej wariancji (stałych)
#    - Usunięcie genów o bardzo niskiej ekspresji

# 4. Normalizacja danych
#    - Zastosowanie transformacji logarytmicznej (log2) lub VST

# 5. Zapisanie gotowego datasetu do folderu data/processed/
