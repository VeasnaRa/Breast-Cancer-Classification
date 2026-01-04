library(ggplot2)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggcorrplot)
library(naniar)
library(visdat)
library(rstatix)
library(DescTools)
library(car)
library(limma)
library(survival)
library(pheatmap)
library(diptest) 
library(forcats)
library(glmnet)
library(caTools)
library(pROC)
library(gridExtra)
library(smotefamily)

## Function helper for models fitting
fit_all_models <- function(X_train_input, X_test_input, Y_train_input, Y_test_input
                           , model_name = "Model") {
  
  cat("\n=== FITTING MODELS:", model_name, "===\n")
  cat("Total features:", ncol(X_train_input), "\n")
  cat("Train samples:", nrow(X_train_input), "\n")
  cat("Test samples:", nrow(X_test_input), "\n\n")
  
  results <- list()
  
  # --- 1. Logistic Regression ---
  cat("1. Logistic Regression...\n")
  train_df <- data.frame(X_train_input, Y = Y_train_input)
  test_df <- data.frame(X_test_input, Y = Y_test_input)
  
  model_logistic <- glm(Y ~ ., data = train_df, family = binomial())
  
  pred_train_log <- predict(model_logistic, train_df, type = "response")
  pred_test_log <- predict(model_logistic, test_df, type = "response")
  
  results$logistic <- list(
    model = model_logistic
    , n_features = ncol(X_train_input)
    , train_auc = as.numeric(auc(Y_train_input, pred_train_log))
    , test_auc = as.numeric(auc(Y_test_input, pred_test_log))
    , test_acc = mean((pred_test_log > 0.5) == Y_test_input)
    , pred_test = pred_test_log
  )
  
  # --- 2. Ridge ---
  cat("2. Ridge Regression...\n")
  set.seed(42)
  cv_ridge <- cv.glmnet(X_train_input, Y_train_input, family = "binomial"
                        , alpha = 0, nfolds = 10, standardize = FALSE)
  
  pred_train_ridge <- predict(cv_ridge, X_train_input, s = cv_ridge$lambda.1se, type = "response")
  pred_test_ridge <- predict(cv_ridge, X_test_input, s = cv_ridge$lambda.1se, type = "response")
  
  results$ridge <- list(
    model = cv_ridge
    , lambda = cv_ridge$lambda.1se
    , n_features = sum(coef(cv_ridge, s = cv_ridge$lambda.1se)[-1] != 0)
    , train_auc = as.numeric(auc(Y_train_input, as.vector(pred_train_ridge)))
    , test_auc = as.numeric(auc(Y_test_input, as.vector(pred_test_ridge)))
    , test_acc = mean((pred_test_ridge > 0.5) == Y_test_input)
    , pred_test = as.vector(pred_test_ridge)
  )
  
  # --- 3. ElasticNet ---
  cat("3. ElasticNet...\n")
  set.seed(42)
  cv_elasticnet <- cv.glmnet(X_train_input, Y_train_input, family = "binomial"
                             , alpha = 0.5, nfolds = 10, standardize = FALSE)
  
  pred_train_elastic <- predict(cv_elasticnet, X_train_input, s = cv_elasticnet$lambda.1se, type = "response")
  pred_test_elastic <- predict(cv_elasticnet, X_test_input, s = cv_elasticnet$lambda.1se, type = "response")
  
  results$elasticnet <- list(
    model = cv_elasticnet
    , lambda = cv_elasticnet$lambda.1se
    , n_features = sum(coef(cv_elasticnet, s = cv_elasticnet$lambda.1se)[-1] != 0)
    , train_auc = as.numeric(auc(Y_train_input, as.vector(pred_train_elastic)))
    , test_auc = as.numeric(auc(Y_test_input, as.vector(pred_test_elastic)))
    , test_acc = mean((pred_test_elastic > 0.5) == Y_test_input)
    , pred_test = as.vector(pred_test_elastic)
  )
  
  # --- 4. Lasso ---
  cat("4. Standard Lasso...\n")
  set.seed(42)
  cv_lasso <- cv.glmnet(X_train_input, Y_train_input, family = "binomial"
                        , alpha = 1, nfolds = 10, standardize = FALSE)
  
  pred_train_lasso <- predict(cv_lasso, X_train_input, s = cv_lasso$lambda.1se, type = "response")
  pred_test_lasso <- predict(cv_lasso, X_test_input, s = cv_lasso$lambda.1se, type = "response")
  
  results$lasso <- list(
    model = cv_lasso
    , lambda = cv_lasso$lambda.1se
    , n_features = sum(coef(cv_lasso, s = cv_lasso$lambda.1se)[-1] != 0)
    , train_auc = as.numeric(auc(Y_train_input, as.vector(pred_train_lasso)))
    , test_auc = as.numeric(auc(Y_test_input, as.vector(pred_test_lasso)))
    , test_acc = mean((pred_test_lasso > 0.5) == Y_test_input)
    , pred_test = as.vector(pred_test_lasso)
  )
  
  # --- 5. Adaptive Lasso ---
  cat("5. Adaptive Lasso...\n")
  ridge_coefs <- coef(cv_ridge, s = cv_ridge$lambda.1se)[-1]
  adaptive_weights <- 1 / (abs(ridge_coefs) + 0.01)
  
  set.seed(42)
  cv_adaptive <- cv.glmnet(X_train_input, Y_train_input, family = "binomial"
                           , alpha = 1, penalty.factor = adaptive_weights
                           , nfolds = 10, standardize = FALSE)
  
  pred_train_adaptive <- predict(cv_adaptive, X_train_input, s = cv_adaptive$lambda.1se, type = "response")
  pred_test_adaptive <- predict(cv_adaptive, X_test_input, s = cv_adaptive$lambda.1se, type = "response")
  
  results$adaptive <- list(
    model = cv_adaptive
    , lambda = cv_adaptive$lambda.1se
    , n_features = sum(coef(cv_adaptive, s = cv_adaptive$lambda.1se)[-1] != 0)
    , train_auc = as.numeric(auc(Y_train_input, as.vector(pred_train_adaptive)))
    , test_auc = as.numeric(auc(Y_test_input, as.vector(pred_test_adaptive)))
    , test_acc = mean((pred_test_adaptive > 0.5) == Y_test_input)
    , pred_test = as.vector(pred_test_adaptive)
  )
  
  # --- 6. uniLasso ---
  cat("6. uniLasso...\n")
  n_features <- ncol(X_train_input)
  uni_coefs <- numeric(n_features)
  
  for (i in 1:n_features) {
    X_uni <- X_train_input[, i, drop = FALSE]
    uni_df <- data.frame(X = X_uni[, 1], Y = Y_train_input)
    
    tryCatch({
      uni_model <- glm(Y ~ X, data = uni_df, family = binomial())
      uni_coefs[i] <- coef(uni_model)[2]
    }, error = function(e) {
      uni_coefs[i] <<- 0
    })
  }
  
  uni_coefs[is.na(uni_coefs)] <- 0
  uni_weights <- 1 / (abs(uni_coefs) + 0.01)
  uni_weights[is.infinite(uni_weights)] <- 1
  
  set.seed(42)
  cv_unilasso <- cv.glmnet(X_train_input, Y_train_input, family = "binomial"
                           , alpha = 1, penalty.factor = uni_weights
                           , nfolds = 10, standardize = FALSE)
  
  pred_train_unilasso <- predict(cv_unilasso, X_train_input, s = cv_unilasso$lambda.1se, type = "response")
  pred_test_unilasso <- predict(cv_unilasso, X_test_input, s = cv_unilasso$lambda.1se, type = "response")
  
  results$unilasso <- list(
    model = cv_unilasso
    , lambda = cv_unilasso$lambda.1se
    , n_features = sum(coef(cv_unilasso, s = cv_unilasso$lambda.1se)[-1] != 0)
    , train_auc = as.numeric(auc(Y_train_input, as.vector(pred_train_unilasso)))
    , test_auc = as.numeric(auc(Y_test_input, as.vector(pred_test_unilasso)))
    , test_acc = mean((pred_test_unilasso > 0.5) == Y_test_input)
    , pred_test = as.vector(pred_test_unilasso)
  )
  
  # Summary table
  comparison_df <- data.frame(
    Feature_Set = model_name
    , Model = c("Logistic", "Ridge", "ElasticNet", "Lasso", "Adaptive Lasso", "uniLasso")
    , Features = c(results$logistic$n_features
                   , results$ridge$n_features
                   , results$elasticnet$n_features
                   , results$lasso$n_features
                   , results$adaptive$n_features
                   , results$unilasso$n_features)
    , Train_AUC = c(results$logistic$train_auc
                    , results$ridge$train_auc
                    , results$elasticnet$train_auc
                    , results$lasso$train_auc
                    , results$adaptive$train_auc
                    , results$unilasso$train_auc)
    , Test_AUC = c(results$logistic$test_auc
                   , results$ridge$test_auc
                   , results$elasticnet$test_auc
                   , results$lasso$test_auc
                   , results$adaptive$test_auc
                   , results$unilasso$test_auc)
    , Test_Accuracy = c(results$logistic$test_acc
                        , results$ridge$test_acc
                        , results$elasticnet$test_acc
                        , results$lasso$test_acc
                        , results$adaptive$test_acc
                        , results$unilasso$test_acc)
  )
  
  results$summary <- comparison_df
  
  # Export metrics to CSV
  csv_filename <- paste0(gsub(" ", "_", model_name), "_all_models_metrics.csv")
  if (!dir.exists("model_metrics")) {
    dir.create("model_metrics", recursive = TRUE)
  }
  write.csv(comparison_df, file.path("model_metrics", csv_filename), row.names = FALSE)
  cat("Exported metrics to:", file.path("model_metrics", csv_filename), "\n\n")
  
  return(results)
}

## Functions for graph plotting
plot_model_comparison <- function(combined_df) {
  
  # Plot 1: Test AUC by Feature Set and Model
  p1 <- ggplot(combined_df, aes(x = Feature_Set, y = Test_AUC, fill = Model, group = Model)) +
    geom_bar(stat = "identity", position = "dodge", color = "#023047", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", Test_AUC))
              , position = position_dodge(width = 0.9)
              , vjust = -0.3
              , size = 2.5
              , fontface = "bold") +
    labs(title = "Model Performance Across Feature Sets"
         , subtitle = "Test AUC comparison for different feature selections"
         , x = "Feature Set"
         , y = "Test AUC"
         , fill = "Model") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.title = element_text(face = "bold", color = "#023047")
          , axis.text.x = element_text(angle = 45, hjust = 1)
          , legend.position = "bottom") +
    scale_fill_manual(values = c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500", "#e63946"))
  
  print(p1)
  
  # Plot 2: Feature Selection Count
  p2 <- ggplot(combined_df, aes(x = Feature_Set, y = Features, fill = Model, group = Model)) +
    geom_bar(stat = "identity", position = "dodge", color = "#023047", linewidth = 0.3) +
    geom_text(aes(label = Features)
              , position = position_dodge(width = 0.9)
              , vjust = -0.3
              , size = 2.5
              , fontface = "bold") +
    labs(title = "Feature Selection Across Feature Sets"
         , subtitle = "Number of selected features by model type"
         , x = "Feature Set"
         , y = "Number of Features"
         , fill = "Model") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.title = element_text(face = "bold", color = "#023047")
          , axis.text.x = element_text(angle = 45, hjust = 1)
          , legend.position = "bottom") +
    scale_fill_manual(values = c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500", "#e63946"))
  
  print(p2)
  
  # Plot 3: Test AUC vs Sparsity
  p3 <- ggplot(combined_df, aes(x = Features, y = Test_AUC, color = Model, shape = Feature_Set)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_line(aes(group = Model), linewidth = 1, alpha = 0.5) +
    labs(title = "Model Performance vs Feature Count"
         , subtitle = "Trade-off between sparsity and predictive accuracy"
         , x = "Number of Features"
         , y = "Test AUC"
         , color = "Model"
         , shape = "Feature Set") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.title = element_text(face = "bold", color = "#023047")
          , legend.position = "right") +
    scale_color_manual(values = c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500", "#e63946"))
  
  print(p3)
  
  # Plot 4: Train vs Test AUC (Overfitting Check)
  p4 <- ggplot(combined_df, aes(x = Train_AUC, y = Test_AUC, color = Model, shape = Feature_Set)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#e63946", linewidth = 1) +
    labs(title = "Train vs Test AUC"
         , subtitle = "Points below diagonal indicate overfitting"
         , x = "Train AUC"
         , y = "Test AUC"
         , color = "Model"
         , shape = "Feature Set") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.title = element_text(face = "bold", color = "#023047")
          , legend.position = "right") +
    scale_color_manual(values = c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500", "#e63946")) +
    coord_cartesian(xlim = c(0.5, 1), ylim = c(0.5, 1))
  
  print(p4)
  
  # Plot 5: Heatmap of Test AUC
  p5 <- ggplot(combined_df, aes(x = Feature_Set, y = Model, fill = Test_AUC)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", Test_AUC))
              , color = "white"
              , fontface = "bold"
              , size = 4) +
    scale_fill_gradient2(low = "#023047", mid = "#219ebc", high = "#ffb703"
                         , midpoint = median(combined_df$Test_AUC)
                         , limits = c(min(combined_df$Test_AUC), max(combined_df$Test_AUC))) +
    labs(title = "Model Performance Heatmap"
         , subtitle = "Test AUC across models and feature sets"
         , x = "Feature Set"
         , y = "Model"
         , fill = "Test AUC") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.title = element_text(face = "bold", color = "#023047")
          , axis.text.x = element_text(angle = 45, hjust = 1)
          , legend.position = "right")
  
  print(p5)
  
}

## functions for feature importance plotting
plot_feature_importance <- function(model_results, top_n = 20, feature_set_name = "Model") {
  
  model_names <- c("Logistic", "Ridge", "ElasticNet", "Lasso", "Adaptive", "uniLasso")
  model_keys <- c("logistic", "ridge", "elasticnet", "lasso", "adaptive", "unilasso")
  
  for (i in 1:6) {
    model_name <- model_names[i]
    model_key <- model_keys[i]
    
    if (model_key == "logistic") {
      coefs <- coef(model_results[[model_key]]$model)
    } else {
      coefs <- coef(model_results[[model_key]]$model, s = model_results[[model_key]]$lambda)
    }
    
    coef_df <- as.data.frame(as.matrix(coefs))
    colnames(coef_df) <- "Coefficient"
    coef_df$Feature <- rownames(coef_df)
    coef_df <- coef_df[coef_df$Feature != "(Intercept)", ]
    coef_df <- coef_df[coef_df$Coefficient != 0, ]
    
    if (nrow(coef_df) > 0) {
      coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]
      
      if (nrow(coef_df) > top_n) {
        coef_df <- coef_df[1:top_n, ]
      }
      
      coef_df <- coef_df[order(coef_df$Coefficient), ]
      coef_df$Feature <- factor(coef_df$Feature, levels = coef_df$Feature)
      coef_df$Direction <- ifelse(coef_df$Coefficient > 0, "Positive", "Negative")
      
      p <- ggplot(coef_df, aes(x = Feature, y = Coefficient, fill = Direction)) +
        geom_col(color = "#023047", linewidth = 0.3) +
        coord_flip() +
        scale_fill_manual(values = c("Positive" = "#219ebc", "Negative" = "#e63946")) +
        labs(title = paste(model_name, "-", feature_set_name)
             , subtitle = paste("Top", nrow(coef_df), "non-zero features")
             , x = NULL
             , y = "Coefficient Value"
             , fill = NULL) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
              , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
              , axis.title = element_text(face = "bold", color = "#023047")
              , legend.position = "bottom")
      
      print(p)
      cat(sprintf("%s: %d non-zero features\n", model_name, nrow(coef_df)))
      
    } else {
      p <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5
                 , label = "No features selected"
                 , size = 5, color = "#555555") +
        labs(title = paste(model_name, "-", feature_set_name)
             , subtitle = "All coefficients = 0") +
        theme_void() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
              , plot.subtitle = element_text(hjust = 0.5, color = "#555555"))
      
      print(p)
      cat(sprintf("%s: 0 non-zero features\n", model_name))
    }
  }
  
}

## Function for exporting metrics to CSV
export_metrics_to_csv <- function(results_object, filename, output_dir = ".") {
  
  # Check if results object has a summary component
  if (is.null(results_object$summary)) {
    cat("Warning: No summary found in results object for", filename, "\n")
    return(NULL)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created directory:", output_dir, "\n")
  }
  
  # Generate full file path
  filepath <- file.path(output_dir, filename)
  
  # Write to CSV
  write.csv(results_object$summary, filepath, row.names = FALSE)
  cat("Exported metrics to:", filepath, "\n")
  
  return(filepath)
}

## Function for exporting all model results
export_all_metrics <- function(model_results_list, output_dir = "model_metrics") {
  
  cat("\n=== EXPORTING ALL METRICS TO CSV ===\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created directory:", output_dir, "\n")
  }
  
  exported_files <- list()
  
  # Loop through all model results and export
  for (model_name in names(model_results_list)) {
    filename <- paste0(model_name, "_metrics.csv")
    filepath <- export_metrics_to_csv(model_results_list[[model_name]]
                                      , filename
                                      , output_dir)
    
    if (!is.null(filepath)) {
      exported_files[[model_name]] <- filepath
    }
  }
  
  cat("\nTotal files exported:", length(exported_files), "\n")
  cat("All metrics saved to directory:", output_dir, "\n")
  
  return(exported_files)
}

# Function for fit single model
fit_single_model_across_features <- function(model_type = "lasso"
                                             , X_train_all, X_test_all
                                             , Y_train, Y_test
                                             , n_clinical
                                             , top_genes_ranked
                                             , gene_sets = c(5000, 1000, 500, 100, 50, 20)) {
  
  cat("\n=== FITTING", toupper(model_type), "ACROSS FEATURE SETS ===\n\n")
  
  results_list <- list()
  feature_set_names <- c("Clinical_Only", paste0("Clinical_TOP", gene_sets))
  
  # 1. Clinical Only
  cat("Fitting Clinical_Only...\n")
  X_train_clinical <- X_train_all[, 1:n_clinical]
  X_test_clinical <- X_test_all[, 1:n_clinical]
  
  result_clinical <- fit_one_model(model_type, X_train_clinical, X_test_clinical
                                   , Y_train, Y_test, "Clinical_Only")
  results_list[[1]] <- result_clinical
  
  # 2. Clinical + Gene sets
  for (i in seq_along(gene_sets)) {
    n_genes <- gene_sets[i]
    feature_set_name <- paste0("Clinical_TOP", n_genes)
    
    cat(sprintf("Fitting %s...\n", feature_set_name))
    
    top_n_genes <- rownames(top_genes_ranked)[1:n_genes]
    gene_indices <- which(colnames(X_train_all) %in% top_n_genes)
    
    X_train_subset <- X_train_all[, c(1:n_clinical, gene_indices)]
    X_test_subset <- X_test_all[, c(1:n_clinical, gene_indices)]
    
    result_subset <- fit_one_model(model_type, X_train_subset, X_test_subset
                                   , Y_train, Y_test, feature_set_name)
    results_list[[i + 1]] <- result_subset
  }
  
  # Combine results into comparison table
  comparison_df <- do.call(rbind, lapply(results_list, function(x) x$summary))
  
  cat("\n=== SUMMARY TABLE ===\n")
  print(comparison_df)
  
  # Export metrics to CSV
  csv_filename <- paste0(model_type, "_across_features_metrics.csv")
  if (!dir.exists("model_metrics")) {
    dir.create("model_metrics", recursive = TRUE)
  }
  write.csv(comparison_df, file.path("model_metrics", csv_filename), row.names = FALSE)
  cat("Exported metrics to:", file.path("model_metrics", csv_filename), "\n")
  
  # Plot comparison
  plot_single_model_comparison(comparison_df, model_type)
  
  # Plot feature importance for each feature set
  for (i in seq_along(results_list)) {
    plot_single_feature_importance(results_list[[i]]$model_obj
                                   , model_type
                                   , feature_set_names[i]
                                   , top_n = 20)
  }
  
  return(list(results = results_list, summary = comparison_df, Y_test = Y_test))
}


# function for fit one model with custom data
fit_one_model <- function(model_type, X_train, X_test, Y_train, Y_test, feature_set_name) {
  
  train_df <- data.frame(X_train, Y = Y_train)
  test_df  <- data.frame(X_test, Y = Y_test)
  
  if (model_type == "logistic") {
    model <- glm(Y ~ ., data = train_df, family = binomial())
    pred_train <- predict(model, train_df, type = "response")
    pred_test <- predict(model, test_df, type = "response")
    n_features <- ncol(X_train)
    
  } else if (model_type == "ridge") {
    set.seed(42)
    model <- cv.glmnet(X_train, Y_train, family = "binomial"
                       , alpha = 0, nfolds = 10, standardize = FALSE)
    pred_train <- as.vector(predict(model, X_train, s = model$lambda.1se, type = "response"))
    pred_test <- as.vector(predict(model, X_test, s = model$lambda.1se, type = "response"))
    n_features <- sum(coef(model, s = model$lambda.1se)[-1] != 0)
    
  } else if (model_type == "elasticnet") {
    set.seed(42)
    model <- cv.glmnet(X_train, Y_train, family = "binomial"
                       , alpha = 0.5, nfolds = 10, standardize = FALSE)
    pred_train <- as.vector(predict(model, X_train, s = model$lambda.1se, type = "response"))
    pred_test <- as.vector(predict(model, X_test, s = model$lambda.1se, type = "response"))
    n_features <- sum(coef(model, s = model$lambda.1se)[-1] != 0)
    
  } else if (model_type == "lasso") {
    set.seed(42)
    model <- cv.glmnet(X_train, Y_train, family = "binomial"
                       , alpha = 1, nfolds = 10, standardize = FALSE)
    pred_train <- as.vector(predict(model, X_train, s = model$lambda.1se, type = "response"))
    pred_test <- as.vector(predict(model, X_test, s = model$lambda.1se, type = "response"))
    n_features <- sum(coef(model, s = model$lambda.1se)[-1] != 0)
    
  } else if (model_type == "adaptive") {
    set.seed(42)
    cv_ridge <- cv.glmnet(X_train, Y_train, family = "binomial"
                          , alpha = 0, nfolds = 10, standardize = FALSE)
    ridge_coefs <- coef(cv_ridge, s = cv_ridge$lambda.1se)[-1]
    adaptive_weights <- 1 / (abs(ridge_coefs) + 0.01)
    
    set.seed(42)
    model <- cv.glmnet(X_train, Y_train, family = "binomial"
                       , alpha = 1, penalty.factor = adaptive_weights
                       , nfolds = 10, standardize = FALSE)
    pred_train <- as.vector(predict(model, X_train, s = model$lambda.1se, type = "response"))
    pred_test <- as.vector(predict(model, X_test, s = model$lambda.1se, type = "response"))
    n_features <- sum(coef(model, s = model$lambda.1se)[-1] != 0)
    
  } else if (model_type == "unilasso") {
    n_feat <- ncol(X_train)
    uni_coefs <- numeric(n_feat)
    
    for (i in 1:n_feat) {
      X_uni <- X_train[, i, drop = FALSE]
      uni_df <- data.frame(X = X_uni[, 1], Y = Y_train)
      tryCatch({
        uni_model <- glm(Y ~ X, data = uni_df, family = binomial())
        uni_coefs[i] <- coef(uni_model)[2]
      }, error = function(e) {
        uni_coefs[i] <<- 0
      })
    }
    
    uni_coefs[is.na(uni_coefs)] <- 0
    uni_weights <- 1 / (abs(uni_coefs) + 0.01)
    uni_weights[is.infinite(uni_weights)] <- 1
    
    set.seed(42)
    model <- cv.glmnet(X_train, Y_train, family = "binomial"
                       , alpha = 1, penalty.factor = uni_weights
                       , nfolds = 10, standardize = FALSE)
    pred_train <- as.vector(predict(model, X_train, s = model$lambda.1se, type = "response"))
    pred_test <- as.vector(predict(model, X_test, s = model$lambda.1se, type = "response"))
    n_features <- sum(coef(model, s = model$lambda.1se)[-1] != 0)
  }
  
  train_auc <- as.numeric(auc(Y_train, pred_train))
  test_auc  <- as.numeric(auc(Y_test, pred_test))
  test_acc  <- mean((pred_test > 0.5) == Y_test)
  
  summary_df <- data.frame(
    Feature_Set = feature_set_name
    , Model = toupper(model_type)
    , Features = n_features
    , Train_AUC = train_auc
    , Test_AUC = test_auc
    , Test_Accuracy = test_acc
  )
  
  return(list(model_obj = model, summary = summary_df, pred_test = pred_test))
}

# Function for single model comparision on different set
plot_single_model_comparison <- function(df, model_type) {
  
  # Plot 1: Test AUC trend
  p1 <- ggplot(df, aes(x = Feature_Set, y = Test_AUC)) +
    geom_bar(stat = "identity", fill = "#219ebc", color = "#023047", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", Test_AUC)), vjust = -0.3
              , fontface = "bold", size = 3) +
    labs(title = paste(toupper(model_type), "- Performance Across Feature Sets")
         , x = "Feature Set"
         , y = "Test AUC") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , axis.title = element_text(face = "bold", color = "#023047")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  
  # Plot 2: Feature count
  p2 <- ggplot(df, aes(x = Feature_Set, y = Features)) +
    geom_bar(stat = "identity", fill = "#ffb703", color = "#023047", linewidth = 0.3) +
    geom_text(aes(label = Features), vjust = -0.3, fontface = "bold", size = 3) +
    labs(title = paste(toupper(model_type), "- Selected Features")
         , x = "Feature Set"
         , y = "Number of Features") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , axis.title = element_text(face = "bold", color = "#023047")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
  
  # Plot 3: Train vs Test AUC
  p3 <- ggplot(df, aes(x = Train_AUC, y = Test_AUC, label = Feature_Set)) +
    geom_point(size = 4, color = "#219ebc") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#e63946") +
    geom_text(vjust = -0.5, size = 3) +
    labs(title = paste(toupper(model_type), "- Train vs Test AUC")
         , x = "Train AUC"
         , y = "Test AUC") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , axis.title = element_text(face = "bold", color = "#023047"))
  
  print(p3)
}

# Function for plot a single model feature importance
plot_single_feature_importance <- function(model_obj, model_type, feature_set_name, top_n = 20) {
  
  if (model_type == "logistic") {
    coefs <- coef(model_obj)
  } else {
    coefs <- coef(model_obj, s = model_obj$lambda.1se)
  }
  
  coef_df <- as.data.frame(as.matrix(coefs))
  colnames(coef_df) <- "Coefficient"
  coef_df$Feature <- rownames(coef_df)
  coef_df <- coef_df[coef_df$Feature != "(Intercept)", ]
  coef_df <- coef_df[coef_df$Coefficient != 0, ]
  
  if (nrow(coef_df) > 0) {
    coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]
    
    if (nrow(coef_df) > top_n) {
      coef_df <- coef_df[1:top_n, ]
    }
    
    coef_df <- coef_df[order(coef_df$Coefficient), ]
    coef_df$Feature <- factor(coef_df$Feature, levels = coef_df$Feature)
    coef_df$Direction <- ifelse(coef_df$Coefficient > 0, "Positive", "Negative")
    
    p <- ggplot(coef_df, aes(x = Feature, y = Coefficient, fill = Direction)) +
      geom_col(color = "#023047", linewidth = 0.3) +
      coord_flip() +
      scale_fill_manual(values = c("Positive" = "#219ebc", "Negative" = "#e63946")) +
      labs(title = paste(toupper(model_type), "-", feature_set_name)
           , subtitle = paste("Top", nrow(coef_df), "non-zero features")
           , x = NULL
           , y = "Coefficient Value"
           , fill = NULL) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
            , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
            , axis.title = element_text(face = "bold", color = "#023047")
            , legend.position = "bottom")
    
    print(p)
    
  } else {
    cat(sprintf("  %s - %s: No features selected\n", toupper(model_type), feature_set_name))
  }
}


## Function for single Model Classificiont Metrics
plot_classification_metrics_single <- function(comparison_results, threshold = 0.5, csv_filename = NULL) {
  
  cat("\n=== CLASSIFICATION METRICS ===\n\n")
  
  results_list <- comparison_results$results
  feature_sets <- sapply(results_list, function(x) x$summary$Feature_Set)
  model_type <- results_list[[1]]$summary$Model[1]
  
  metrics_list <- list()
  
  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    feature_set <- feature_sets[i]
    
    # Get predictions
    Y_test <- comparison_results$Y_test
    pred_prob <- result$pred_test
    pred_class <- ifelse(pred_prob > threshold, 1, 0)
    
    # Confusion matrix
    TP <- sum(pred_class == 1 & Y_test == 1)
    TN <- sum(pred_class == 0 & Y_test == 0)
    FP <- sum(pred_class == 1 & Y_test == 0)
    FN <- sum(pred_class == 0 & Y_test == 1)
    
    # Metrics
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
    recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
    specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
    f1_score <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
    test_auc <- as.numeric(auc(Y_test, pred_prob))
    
    metrics_list[[i]] <- data.frame(
      Feature_Set = feature_set
      , TP = TP
      , TN = TN
      , FP = FP
      , FN = FN
      , Accuracy = accuracy
      , Precision = precision
      , Recall = recall
      , Specificity = specificity
      , F1_Score = f1_score
      , AUC = test_auc
    )
    
    cat(sprintf("%s:\n", feature_set))
    cat(sprintf("  TP=%d TN=%d FP=%d FN=%d\n", TP, TN, FP, FN))
    cat(sprintf("  Accuracy=%.3f Precision=%.3f Recall=%.3f F1=%.3f AUC=%.3f\n\n",
                accuracy, precision, recall, f1_score, test_auc))
  }
  
  metrics_df <- do.call(rbind, metrics_list)
  
  # Plot 1: Confusion Matrix Heatmap
  cm_long <- reshape2::melt(metrics_df[, c("Feature_Set", "TP", "TN", "FP", "FN")], id.vars = "Feature_Set")
  colnames(cm_long) <- c("Feature_Set", "Metric", "Count")
  
  p1 <- ggplot(cm_long, aes(x = Feature_Set, y = Metric, fill = Count)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Count), color = "white", fontface = "bold", size = 4) +
    scale_fill_gradient(low = "#023047", high = "#ffb703") +
    labs(title = paste(model_type, "- Confusion Matrix Across Feature Sets")
         , subtitle = "TP, TN, FP, FN counts"
         , x = NULL
         , y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  
  # Plot 2: Performance Metrics
  perf_long <- reshape2::melt(metrics_df[, c("Feature_Set", "Accuracy", "Precision", "Recall", "F1_Score", "AUC")], id.vars = "Feature_Set")
  colnames(perf_long) <- c("Feature_Set", "Metric", "Value")
  
  p2 <- ggplot(perf_long, aes(x = Feature_Set, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "#023047", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", Value))
              , position = position_dodge(width = 0.9)
              , vjust = -0.3
              , size = 2.2) +
    labs(title = paste(model_type, "- Classification Metrics")
         , subtitle = "Accuracy, Precision, Recall, F1-Score, AUC"
         , x = NULL
         , y = "Score"
         , fill = "Metric") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.text.x = element_text(angle = 45, hjust = 1)
          , legend.position = "bottom") +
    scale_fill_manual(values = c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500")) +
    ylim(0, 1.1)
  
  print(p2)
  
  # Plot 3: F1-Score trend
  p3 <- ggplot(metrics_df, aes(x = Feature_Set, y = F1_Score, group = 1)) +
    geom_line(color = "#219ebc", linewidth = 1.2) +
    geom_point(size = 4, color = "#023047") +
    geom_text(aes(label = sprintf("%.3f", F1_Score)), vjust = -0.8, fontface = "bold") +
    labs(title = paste(model_type, "- F1-Score Across Feature Sets")
         , subtitle = "Trend of model performance"
         , x = NULL
         , y = "F1-Score") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1)
  
  print(p3)
  
  # Plot 4: AUC Trend
  p4 <- ggplot(metrics_df, aes(x = Feature_Set, y = AUC, group = 1)) +
    geom_line(color = "#ffb703", linewidth = 1.2) +
    geom_point(size = 4, color = "#023047") +
    geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.8, fontface = "bold") +
    labs(title = paste(model_type, "- AUC Across Feature Sets")
         , subtitle = "Area Under the ROC Curve"
         , x = NULL
         , y = "AUC") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1)
  
  print(p4)
  
  # Plot 5: Sensitivity vs Specificity
  p5 <- ggplot(metrics_df, aes(x = Specificity, y = Recall, color = Feature_Set, label = Feature_Set)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(vjust = -0.8, size = 3) +
    labs(title = paste(model_type, "- Sensitivity vs Specificity")
         , subtitle = "Trade-off across feature sets"
         , x = "Specificity (True Negative Rate)"
         , y = "Sensitivity (Recall / True Positive Rate)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555")
          , legend.position = "bottom") +
    xlim(0, 1) + ylim(0, 1)
  
  print(p5)
  
  cat("\n=== SUMMARY TABLE ===\n")
  print(metrics_df)
  
  # Export to CSV if filename is provided
  if (!is.null(csv_filename)) {
    if (!dir.exists("model_metrics")) {
      dir.create("model_metrics", recursive = TRUE)
    }
    filepath <- file.path("model_metrics", csv_filename)
    write.csv(metrics_df, filepath, row.names = FALSE)
    cat("\nExported classification metrics to:", filepath, "\n")
  }
  
  return(metrics_df)
}

## Function for Class imbalanced
apply_smote <- function(X_train, Y_train, k = 5) {
  
  cat("\n=== APPLYING SMOTE ===\n")
  cat("Before SMOTE:\n")
  cat("  Alive (0):", sum(Y_train == 0), "\n")
  cat("  Dead (1):", sum(Y_train == 1), "\n")
  cat("  Ratio:", round(sum(Y_train == 0) / sum(Y_train == 1), 2), ":1\n\n")
  
  # Combine X and Y for SMOTE
  train_data <- data.frame(X_train, Class = Y_train)
  
  # Apply SMOTE
  smote_result <- SMOTE(train_data[, -ncol(train_data)], train_data$Class, K = k)
  
  # Extract balanced data
  X_train_balanced <- smote_result$data[, -ncol(smote_result$data)]
  Y_train_balanced <- smote_result$data$class
  
  cat("After SMOTE:\n")
  cat("  Alive (0):", sum(Y_train_balanced == 0), "\n")
  cat("  Dead (1):", sum(Y_train_balanced == 1), "\n")
  cat("  Ratio:", round(sum(Y_train_balanced == 0) / sum(Y_train_balanced == 1), 2), ":1\n")
  cat("  Total samples:", nrow(X_train_balanced), "\n\n")
  
  return(list(
    X_train = as.matrix(X_train_balanced),
    Y_train = as.numeric(Y_train_balanced)
  ))
}

feature_importance <- function(model_obj, model_name = "Model", top_n = 20, gene_names) {
  
  cat("\n", rep("=", 80), "\n", sep="")
  cat("FEATURE IMPORTANCE ANALYSIS:", model_name, "\n")
  cat(rep("=", 80), "\n\n")
  
  # --- 1. Extract Coefficients ---
  coefs <- coef(model_obj, s = model_obj$lambda.1se)
  coef_df <- data.frame(
    Feature = rownames(coefs)[-1]
    , Coefficient = as.vector(coefs)[-1]
  )
  
  # Remove zeros
  coef_df <- coef_df[coef_df$Coefficient != 0, ]
  coef_df$Abs_Coef <- abs(coef_df$Coefficient)
  coef_df <- coef_df[order(-coef_df$Abs_Coef), ]
  
  # --- 2. Classify Features (CORRECTED) ---
  coef_df$Type <- ifelse(coef_df$Feature %in% gene_names, "Genomic", "Clinical")
  
  # --- 3. Add Odds Ratios ---
  coef_df$Odds_Ratio <- exp(coef_df$Coefficient)
  
  # --- 4. Add Direction ---
  coef_df$Direction <- ifelse(coef_df$Coefficient > 0
                              , "Increases Death Risk"
                              , "Decreases Death Risk")
  
  # --- 5. Add Rank ---
  coef_df$Rank <- 1:nrow(coef_df)
  
  # --- 6. Summary Statistics ---
  cat("=== SUMMARY ===\n")
  cat("Total features selected:", nrow(coef_df), "\n")
  cat("  Clinical:", sum(coef_df$Type == "Clinical"), "\n")
  cat("  Genomic:", sum(coef_df$Type == "Genomic"), "\n\n")
  
  cat("Direction:\n")
  cat("  Increases death risk:", sum(coef_df$Direction == "Increases Death Risk"), "\n")
  cat("  Decreases death risk:", sum(coef_df$Direction == "Decreases Death Risk"), "\n\n")
  
  cat("Coefficient range:\n")
  cat("  Min:", round(min(coef_df$Coefficient), 4), "\n")
  cat("  Max:", round(max(coef_df$Coefficient), 4), "\n")
  cat("  Mean (absolute):", round(mean(coef_df$Abs_Coef), 4), "\n\n")
  
  cat("Odds Ratio range:\n")
  cat("  Min:", round(min(coef_df$Odds_Ratio), 4), "\n")
  cat("  Max:", round(max(coef_df$Odds_Ratio), 4), "\n\n")
  
  # --- 7. Top Features Table ---
  cat("=== TOP", top_n, "FEATURES ===\n\n")
  top_features <- head(coef_df, top_n)
  print_df <- top_features[, c("Rank", "Feature", "Type", "Coefficient", "Odds_Ratio", "Direction")]
  print_df$Coefficient <- round(print_df$Coefficient, 4)
  print_df$Odds_Ratio <- round(print_df$Odds_Ratio, 4)
  print(print_df, row.names = FALSE)
  cat("\n")
  
  # --- 8. Top Clinical Features ---
  clinical_df <- coef_df[coef_df$Type == "Clinical", ]
  if(nrow(clinical_df) > 0) {
    cat("=== TOP CLINICAL FEATURES ===\n\n")
    top_clinical <- head(clinical_df, 10)
    print_clinical <- top_clinical[, c("Feature", "Coefficient", "Odds_Ratio", "Direction")]
    print_clinical$Coefficient <- round(print_clinical$Coefficient, 4)
    print_clinical$Odds_Ratio <- round(print_clinical$Odds_Ratio, 4)
    print(print_clinical, row.names = FALSE)
    cat("\n")
  }
  
  # --- 9. Top Genomic Features ---
  genomic_df <- coef_df[coef_df$Type == "Genomic", ]
  if(nrow(genomic_df) > 0) {
    cat("=== TOP GENOMIC FEATURES ===\n\n")
    top_genomic <- head(genomic_df, 10)
    print_genomic <- top_genomic[, c("Feature", "Coefficient", "Odds_Ratio", "Direction")]
    print_genomic$Coefficient <- round(print_genomic$Coefficient, 4)
    print_genomic$Odds_Ratio <- round(print_genomic$Odds_Ratio, 4)
    print(print_genomic, row.names = FALSE)
    cat("\n")
  }
  
  # --- 10. Strongest Protective Features ---
  protective <- coef_df[coef_df$Direction == "Decreases Death Risk", ]
  protective <- protective[order(protective$Coefficient), ]
  top_protective <- head(protective, 5)
  print_protective <- top_protective[, c("Feature", "Type", "Coefficient", "Odds_Ratio")]
  print_protective$Coefficient <- round(print_protective$Coefficient, 4)
  print_protective$Odds_Ratio <- round(print_protective$Odds_Ratio, 4)
  print(print_protective, row.names = FALSE)
  cat("\n")
  
  # --- 11. Strongest Risk Features ---
  risk <- coef_df[coef_df$Direction == "Increases Death Risk", ]
  risk <- risk[order(-risk$Coefficient), ]
  top_risk <- head(risk, 5)
  print_risk <- top_risk[, c("Feature", "Type", "Coefficient", "Odds_Ratio")]
  print_risk$Coefficient <- round(print_risk$Coefficient, 4)
  print_risk$Odds_Ratio <- round(print_risk$Odds_Ratio, 4)
  print(print_risk, row.names = FALSE)
  cat("\n")
  
  # --- 12. Improved Visualizations with ggplot2 ---
  
  # Plot 1: Top 15 Features (Clinical vs Genomic)
  top_plot_df <- head(coef_df, 15)
  top_plot_df$Feature <- factor(top_plot_df$Feature
                                , levels = top_plot_df$Feature[order(top_plot_df$Coefficient)])
  
  coef_range <- range(top_plot_df$Coefficient)
  y_expansion <- max(abs(coef_range)) * 0.3
  y_limits <- c(min(coef_range) - y_expansion, max(coef_range) + y_expansion)
  
  p1 <- ggplot(top_plot_df, aes(x = Feature, y = Coefficient, fill = Type)) +
    geom_col(color = "#023047", linewidth = 0.5, width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Coefficient))
              , hjust = ifelse(top_plot_df$Coefficient > 0, -0.15, 1.15)
              , size = 3.5
              , fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#666666", linewidth = 1) +
    coord_flip(ylim = y_limits, clip = "off") +
    scale_fill_manual(values = c("Clinical" = "#219ebc", "Genomic" = "#fb8500")) +
    labs(title = paste(model_name, "- Top 15 Features")
         , subtitle = "Clinical vs Genomic Features"
         , x = NULL
         , y = "Coefficient Value"
         , fill = "Feature Type") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555", size = 10)
          , axis.title.x = element_text(face = "bold", color = "#023047", size = 10, margin = margin(t = 10))
          , axis.text.y = element_text(face = "bold", size = 9)
          , axis.text.x = element_text(size = 9)
          , legend.position = "bottom"
          , legend.title = element_text(face = "bold", size = 10)
          , legend.text = element_text(size = 9)
          , panel.grid.major.y = element_blank()
          , plot.margin = margin(t = 10, r = 50, b = 10, l = 10))
  
  print(p1)
  
  # Plot 2: Top Clinical Features
  if(nrow(clinical_df) > 0) {
    top_clin_plot <- head(clinical_df, 10)
    top_clin_plot$Feature <- factor(top_clin_plot$Feature
                                    , levels = top_clin_plot$Feature[order(top_clin_plot$Coefficient)])
    top_clin_plot$Color <- ifelse(top_clin_plot$Coefficient > 0, "Risk", "Protective")
    
    p2 <- ggplot(top_clin_plot, aes(x = Feature, y = Coefficient, fill = Color)) +
      geom_col(color = "#023047", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.3f", Coefficient))
                , hjust = ifelse(top_clin_plot$Coefficient > 0, -0.1, 1.1)
                , size = 3
                , fontface = "bold") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#666666", linewidth = 1) +
      coord_flip(ylim = c(min(top_clin_plot$Coefficient) * 1.2, max(top_clin_plot$Coefficient) * 1.2)) +
      scale_fill_manual(values = c("Risk" = "#e63946", "Protective" = "#219ebc")
                        , labels = c("Protective" = "Decreases Death Risk", "Risk" = "Increases Death Risk")) +
      labs(title = "Top 10 Clinical Features"
           , subtitle = "Effect on mortality risk"
           , x = NULL
           , y = "Coefficient Value"
           , fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13, color = "#023047")
            , plot.subtitle = element_text(hjust = 0.5, color = "#555555", size = 10)
            , axis.title.x = element_text(face = "bold", color = "#023047", size = 11)
            , axis.text.y = element_text(face = "bold", size = 10)
            , axis.text.x = element_text(size = 10)
            , legend.position = "bottom"
            , legend.text = element_text(size = 10)
            , panel.grid.major.y = element_blank()
            , plot.margin = margin(10, 40, 10, 10))
    
    print(p2)
  }
  
  # Plot 3: Top Genomic Features
  if(nrow(genomic_df) > 0) {
    top_gen_plot <- head(genomic_df, 10)
    top_gen_plot$Feature <- factor(top_gen_plot$Feature
                                   , levels = top_gen_plot$Feature[order(top_gen_plot$Coefficient)])
    top_gen_plot$Color <- ifelse(top_gen_plot$Coefficient > 0, "Risk", "Protective")
    
    p3 <- ggplot(top_gen_plot, aes(x = Feature, y = Coefficient, fill = Color)) +
      geom_col(color = "#023047", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.3f", Coefficient))
                , hjust = ifelse(top_gen_plot$Coefficient > 0, -0.1, 1.1)
                , size = 3
                , fontface = "bold") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#666666", linewidth = 1) +
      coord_flip(ylim = c(min(top_gen_plot$Coefficient) * 1.2, max(top_gen_plot$Coefficient) * 1.2)) +
      scale_fill_manual(values = c("Risk" = "#e63946", "Protective" = "#219ebc")
                        , labels = c("Protective" = "Decreases Death Risk", "Risk" = "Increases Death Risk")) +
      labs(title = "Top 10 Genomic Features"
           , subtitle = "Gene expression effects on mortality"
           , x = NULL
           , y = "Coefficient Value"
           , fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13, color = "#023047")
            , plot.subtitle = element_text(hjust = 0.5, color = "#555555", size = 10)
            , axis.title.x = element_text(face = "bold", color = "#023047", size = 11)
            , axis.text.y = element_text(face = "bold", size = 10)
            , axis.text.x = element_text(size = 10)
            , legend.position = "bottom"
            , legend.text = element_text(size = 10)
            , panel.grid.major.y = element_blank()
            , plot.margin = margin(10, 40, 10, 10))
    
    print(p3)
  }
  
  # Plot 4: Odds Ratio Distribution
  or_range <- range(coef_df$Odds_Ratio)
  
  p4 <- ggplot(coef_df, aes(x = Odds_Ratio)) +
    geom_histogram(bins = 30, fill = "#8ecae6", color = "white", alpha = 0.9) +
    geom_vline(xintercept = 1, color = "#e63946", linewidth = 1.5, linetype = "dashed") +
    annotate("text", x = 1, y = Inf, label = "OR = 1\n(No Effect)"
             , color = "#e63946", vjust = 2, fontface = "bold", size = 4) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(title = "Odds Ratio Distribution"
         , subtitle = paste0("Range: ", sprintf("%.2f", or_range[1]), " to ", sprintf("%.2f", or_range[2]))
         , x = "Odds Ratio (OR)"
         , y = "Frequency (Number of Features)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13, color = "#023047")
          , plot.subtitle = element_text(hjust = 0.5, color = "#555555", size = 10)
          , axis.title = element_text(face = "bold", color = "#023047", size = 11)
          , axis.text = element_text(size = 10)
          , panel.grid.minor = element_blank())
  
  print(p4)
  
  # --- 13. Export to CSV ---
  csv_filename <- paste0(gsub(" ", "_", gsub("[^[:alnum:] ]", "", model_name)), "_feature_importance.csv")
  if (!dir.exists("model_metrics")) {
    dir.create("model_metrics", recursive = TRUE)
  }
  write.csv(coef_df, file.path("model_metrics", csv_filename), row.names = FALSE)
  cat("\nExported feature importance to:", file.path("model_metrics", csv_filename), "\n")
  
  return(list(
    all_features = coef_df
    , clinical = clinical_df
    , genomic = genomic_df
    , protective = protective
    , risk = risk
  ))
  
}