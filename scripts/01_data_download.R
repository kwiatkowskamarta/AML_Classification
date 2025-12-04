# SKRYPT 01: Data Download
# ------------------------------------------
# OPIS: Pobiera surowe dane RNA-seq i kliniczne z bazy TCGA dla projektu LAML.
# INPUT: Brak (pobieranie z API)
# OUTPUT: Plik data/raw/TCGA_LAML_Raw.RData

# 1. Ładowanie bibliotek (TCGAbiolinks)

# 2. Definicja zapytania (GDCquery) - dane kliniczne

# 3. Pobranie i przygotowanie danych klinicznych

# 4. Definicja zapytania (GDCquery) - dane transkryptomiczne (STAR - Counts)

# 5. Pobranie i przygotowanie danych RNA-seq

# 6. Zapisanie surowych obiektów na dysk (data/raw/)
