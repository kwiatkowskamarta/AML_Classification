# SKRYPT INSTALACYJNY (SETUP)
# Uruchomic ten skrypt tylko raz na nowym srodowisku, aby pobrac biblioteki.
# ------------------------------------------------------------------------

# 1. Instalacja menedzera srodowisk 'renv' (jesli brak)
if (!require('renv', quietly = TRUE)) install.packages('renv')

# 2. Aktywacja projektu (tworzy folder bibliotek wewnątrz projektu)
renv::init(bare = TRUE) 

# 3. Instalacja menedzera pakietow bioinformatycznych BiocManager
if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')

# 4. Lista wymaganych pakietow
bio_packages <- c(
  'TCGAbiolinks',       
  'SummarizedExperiment', 
  'edgeR',             
  'DESeq2',             
  'ComplexHeatmap'      
)

ml_packages <- c(
  'caret',             
  'randomForest',       
  'kernlab',            
  'Boruta',             
  'pheatmap',           
  'ggplot2',            
  'dplyr',              
  'e1071'               
)

# 5. Instalacja pakietow
message('Instalowanie pakietow z Bioconductor.')
BiocManager::install(bio_packages, update = FALSE, ask = FALSE)

message('Instalowanie pakietow ML.')
install.packages(ml_packages)

# 6. Zapisanie stanu bibliotek do pliku renv.lock
message('Zapisywanie stanu srodowiska do renv.lock.')
renv::snapshot()

message('Srodowisko jest zainstalowane i zapisane.')
