# Porównawcza analiza skuteczności algorytmów ML w klasyfikacji białaczki AML (TCGA)

## O projekcie
Repozytorium zawiera kod źródłowy oraz dokumentację do pracy magisterskiej pt.:
> **"Porównawcza analiza skuteczności algorytmów uczenia maszynowego i metod statystycznych w klasyfikacji podtypów molekularnych ostrej białaczki szpikowej (AML) na podstawie danych RNA-seq."**

Głównym celem projektu jest sprawdzenie, czy algorytmy uczenia maszynowego (Random Forest, SVM, XGBoost) pozwalają na skuteczniejszą klasyfikację podtypów białaczki niż klasyczne metody statystyczne, przy wykorzystaniu danych transkryptomicznych z bazy **The Cancer Genome Atlas (TCGA)**.
---

## Uruchamianie projektu:
### Wymagania wstępne
* R (wersja 4.0 lub nowsza)
* RStudio
* Połączenie z internetem (do pobrania danych TCGA)

### Instrukcja krok po kroku

1.  **Sklonuj repozytorium:**
    ```bash
    git clone [https://github.com/mkwiatk9/AML_Classification_Thesis.git](https://github.com/mkwiatk9/AML_Classification_Thesis.git)
    ```
2.  **Otwórz projekt:**
    Kliknij dwukrotnie plik `AML_Classification_Thesis.Rproj` w folderze projektu.

3.  **Zainstaluj zależności:**
    Uruchom przygotowany skrypt instalacyjny, który automatycznie skonfiguruje środowisko `renv` i pobierze wymagane biblioteki:
    ```r
    source("install_dependencies.R")
    ```

4.  **Uruchom Pipeline:**
    Skrypty należy uruchamiać w kolejności numerycznej:

    * `scripts/01_data_download.R` - Pobiera surowe dane z TCGA (może zająć kilka minut).
    * `scripts/02_preprocessing.R` - Czyści dane, normalizuje i filtruje próbki.
    * `scripts/03_feature_selection.R` - Wybiera najważniejsze geny (Boruta).
    * `scripts/04_modeling.R` - Trenuje modele i generuje wyniki.

---

## Struktura Projektu

```text
├── data/
│   ├── raw/               # Surowe dane z TCGA (ignorowane przez git)
│   └── processed/         # Wyczyszczone dane gotowe do modeli
├── scripts/
│   ├── 01_data_download.R    # Pobieranie danych (TCGAbiolinks)
│   ├── 02_preprocessing.R    # Normalizacja i czyszczenie
│   ├── 03_feature_selection.R # Selekcja cech (Boruta)
│   └── 04_modeling.R         # Trening i walidacja modeli (caret)
├── results/               # Wyniki analiz, wykresy, zapisane modele
├── renv/                  # Pliki konfiguracyjne środowiska R
├── .gitignore             # Pliki ignorowane przez Git (np. duże dane)
├── install_dependencies.R # Skrypt instalacyjny dla nowych użytkowników
├── renv.lock              # Dokładna lista wersji pakietów (dla reprodukowalności)
└── README.md              # Dokumentacja projektu