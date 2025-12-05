# SKRYPT 04: Modeling (Wersja Fail-Safe)
# ------------------------------------------------------------------------------
# Cel: Udowodnienie wyższości Random Forest.
# ------------------------------------------------------------------------------

library(caret)
library(randomForest)
library(kernlab)
library(pheatmap)
library(dplyr)

set.seed(777)

message(">>> START: ETAP MODELOWANIA <<<")

# 1. WCZYTANIE DANYCH
load("results/models/Boruta_Results.RData")

# 2. FILTRACJA I PODZIAŁ
# Usuwamy klasy < 10 próbek
class_counts <- table(ml_data_final$Target)
classes_to_keep <- names(class_counts)[class_counts >= 10]

ml_data_filtered <- ml_data_final %>%
  filter(Target %in% classes_to_keep) %>%
  mutate(Target = droplevels(Target))

message("Klasy w analizie: ", paste(levels(ml_data_filtered$Target), collapse=", "))

train_index <- createDataPartition(ml_data_filtered$Target, p = 0.75, list = FALSE)
train_set <- ml_data_filtered[train_index, ]
test_set  <- ml_data_filtered[-train_index, ]

# 3. KONFIGURACJA TRENINGU
fit_control <- trainControl(
  method = "repeatedcv",
  number = 3,           # 3-fold 
  repeats = 5,          
  classProbs = TRUE,
  savePredictions = "final"
)

message("\n>>> ROZPOCZYNAM TRENING MODELI <<<")

# --- MODEL 1: Random Forest  ---
message("1. Trenowanie Random Forest.")
set.seed(123)
model_rf <- train(
  Target ~ ., data = train_set, 
  method = "rf", 
  trControl = fit_control, 
  ntree = 500,
  tuneLength = 5
)
message("   -> Random Forest: done!")

# --- MODEL 2: k-NN ---
message("2. Trenowanie k-NN...")
set.seed(123)
model_knn <- train(
  Target ~ ., data = train_set, 
  method = "knn", 
  trControl = fit_control, 
  preProcess = c("center", "scale"),
  tuneLength = 10
)
message("   -> k-NN: done!")

# --- MODEL 3: SVM Linear ---
message("3. Trenowanie SVM (Linear)...")
set.seed(123)
model_svm_linear <- tryCatch({
  train(
    Target ~ ., data = train_set, 
    method = "svmLinear", 
    trControl = fit_control, 
    preProcess = c("center", "scale"),
    tuneGrid = expand.grid(C = 1) # Sztywny parametr (ułatwia zbieżność)
  )
}, error = function(e) { 
  message("   ⚠️ OSTRZEŻENIE: SVM Linear nie zbiegł się (Błąd matematyczny). Pomijam go.")
  return(NULL) 
})

if(!is.null(model_svm_linear)) message("   -> SVM: done!")

message("--- TRENING ZAKOŃCZONY ---")

# 4. ZESTAWIENIE WYNIKÓW
# Tworzymy listę tylko z tych modeli, które się udały
models_list <- list(Random_Forest = model_rf, kNN = model_knn)
if(!is.null(model_svm_linear)) models_list$SVM <- model_svm_linear

results <- resamples(models_list)

message("\n--- WYNIKI TRENINGOWE (Cross-Validation) ---")
print(summary(results))

# 5. EWALUACJA TESTOWA
message("\n--- WERYFIKACJA NA ZBIORZE TESTOWYM ---")

evaluate_model <- function(model, name) {
  if(is.null(model)) return(NULL)
  pred <- predict(model, test_set)
  cm <- confusionMatrix(pred, test_set$Target)
  acc <- round(cm$overall['Accuracy'] * 100, 2)
  kappa <- round(cm$overall['Kappa'], 2)
  message(paste0(name, ": Accuracy = ", acc, "%, Kappa = ", kappa))
  return(cm)
}

# Najważniejszy wynik!
cm_rf <- evaluate_model(model_rf, "Random Forest")
cm_knn <- evaluate_model(model_knn, "k-NN")
if(!is.null(model_svm_linear)) cm_svm <- evaluate_model(model_svm_linear, "SVM")

# 6. ZAPIS
if(!dir.exists("results/models")) dir.create("results/models", recursive = TRUE)
save(models_list, results, cm_rf, file = "results/models/Final_Models.RData")

message("\n>>> Analiza zakończona. <<<")