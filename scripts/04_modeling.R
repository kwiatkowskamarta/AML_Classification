# SKRYPT 04: Machine Learning Models
# ------------------------------------------
# OPIS: Trenuje i porównuje modele klasyfikacyjne (RF, SVM, kNN).
# INPUT: data/processed/TCGA_LAML_Cleaned.RData + results/selected_features.rds
# OUTPUT: results/model_performance.csv

# 1. Wczytanie danych i wyfiltrowanie tylko wybranych genów

# 2. Podział na zbiór treningowy i testowy (np. 75/25)
#    - Ustawienie ziarna losowości (set.seed) dla powtarzalności

# 3. Konfiguracja walidacji krzyżowej (np. 5-fold CV)

# 4. Trening modeli:
#    - Model 1: k-NN (Baseline)
#    - Model 2: SVM (Linear/Radial)
#    - Model 3: Random Forest

# 5. Predykcja na zbiorze testowym

# 6. Ewaluacja i wizualizacja
#    - Macierze pomyłek (Confusion Matrix)
#    - Porównanie metryk (Accuracy, Kappa)

# 7. Zapisanie wyników
