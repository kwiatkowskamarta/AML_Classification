# Klasyfikacja podtypów białaczki AML z wykorzystaniem uczenia maszyowego

Projekt realizowany w ramach pracy magisterskiej. Celem jest opracowanie narzędzia bioinformatycznego do automatycznej klasyfikacji podtypów ostrej białaczki szpikowej (AML) na podstawie profilu ekspresji genów (RNA-seq).

## Główne funkcjonalności

Projekt realizuje pełne przetwarzanie danych:
1.  **Pobieranie Danych:** Automatyczna integracja z bazą NCI Genomic Data Commons (GDC) poprzez API, pobierająca dane projektu TCGA-LAML.
2.  **Preprocessing:** Zaawansowane czyszczenie danych klinicznych, parowanie próbek i normalizacja ekspresji genów (Log2 CPM).
3.  **Selekcja Biomarkerów:** Wykorzystanie algorytmu Boruta do redukcji wymiarowości z 26,000 genów do 150 kluczowych biomarkerów.
4.  **Machine Learning:** Trening i walidacja modeli klasyfikacyjnych (Random Forest vs k-NN vs SVM) z wykorzystaniem powtarzanej walidacji krzyżowej (Repeated Cross-Validation).
5.  **Interaktywne GUI:** Aplikacja webowa (R Shiny) umożliwiająca lekarzom i badaczom wizualizację wyników oraz symulację diagnostyczną.

## Instrukcja Uruchomienia

Projekt jest w pełni reprodukowalny dzięki wykorzystaniu pakietu `renv`.

### Wymagania
* R (wersja 4.0+)
* RStudio
* Git

### Krok po kroku

1.  **Sklonuj repozytorium:**
    ```bash
    git clone [https://github.com/twoj-nick/AML_Classification.git](https://github.com/twoj-nick/AML_Classification.git)
    ```
2.  **Otwórz projekt:** Kliknij plik `AML_Classification.Rproj`.
3.  **Zainstaluj biblioteki:**
    Uruchom przygotowany skrypt, który skonfiguruje środowisko:
    ```r
    source("install_dependencies.R")
    ```
4.  **Uruchom analizę (Opcjonalnie):**
    ```r
    source("scripts/00_run_pipeline.R")
    ```
5.  **Uruchom Aplikację GUI:**
    Aby zobaczyć gotowy dashboard z wynikami:
    ```r
    shiny::runApp("scripts/05_GUI.R")
    ```

## Struktura Projektu

* `data/` - Dane.
* `scripts/`
    * `01_data_download.R` - Pobieranie danych z TCGA.
    * `02_preprocessing.R` - Czyszczenie i normalizacja.
    * `03_feature_selection.R` - Algorytm Boruta.
    * `04_modeling.R` - Trening Random Forest/k-NN.
    * `05_GUI.R` - Kod aplikacji Shiny Dashboard.
* `results/` - Zapisane modele (.RData) i wykresy.
* `renv.lock` - Plik blokady wersji pakietów (dla reprodukowalności).

## Autorka
**Marta Kwiatkowska**
Praca Magisterska
Politechnika Warszawska
