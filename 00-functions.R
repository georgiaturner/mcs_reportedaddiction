
 
### colour palettes




################################################################################
# Perform MDS functions
################################################################################

### MDS to get the dimensions on the training dataset
perform_mds_train <- function(data, k = 2, remove_cols = NULL) {
  # Optionally remove ID column
  if (!is.null(remove_cols)) {
    data <- dplyr::select(data, -all_of(remove_cols))
  }
  # Compute distance matrix
  dist_matrix <- dist(data)
  # Perform classical (metric) MDS
  mds_result  <- cmdscale(dist_matrix, eig = TRUE, k = k)
}

### Identify the coords for the same dimensions on the test dataset
project_mds_test <- function(train_data, test_data, mds_result, k = 2, remove_cols = NULL) {
  if (!is.null(remove_cols)) {
    train_data <- dplyr::select(train_data, -all_of(remove_cols))
    test_data  <- dplyr::select(test_data,  -all_of(remove_cols))
  }
  
  # Mean center the data based on train (for consistency)
  train_centered <- scale(train_data, center = TRUE, scale = FALSE)
  test_centered  <- scale(test_data,  center = attr(train_centered, "scaled:center"), scale = FALSE)
  
  # Compute inner product matrix (B) from cmdscale eigenvectors/values
  eigvals <- mds_result$eig[1:k]
  X_train <- mds_result$points[, 1:k]
  
  # Approximate using scalar product with centered data
  test_coords <- as.matrix(test_centered) %*% (t(as.matrix(train_centered)) %*% X_train) %*% diag(1 / eigvals)
  
  test_coords
}

flip_mds_dimensions <- function(mds_result, flip_dims = c(1, 2)) {
  mds_result$points[, flip_dims] <- -mds_result$points[, flip_dims]
  return(mds_result)
}


plot_mds_variable_correlations <- function(
    vars_train, vars_test, mds_train, mds_test,
    title_train = "Train", title_test = "Test",
    top_n = NULL,                 # Optional: top N variables by importance
    rename_vars = FALSE,          # Optional: whether to rename variables
    var_name_map = NULL           # Named vector or list: old_name = new_name
) {
  # Load required packages
  require(pheatmap)
  require(gridExtra)
  
  # Step 1: Create MDS data frames with proper column names
  mds_df_train <- as.data.frame(mds_train)
  mds_df_test  <- as.data.frame(mds_test)
  colnames(mds_df_train) <- paste0("Dim", seq_len(ncol(mds_df_train)))
  colnames(mds_df_test)  <- paste0("Dim", seq_len(ncol(mds_df_test)))
  
  # Step 2: Remove zero-variance columns
  vars_train_clean <- vars_train[, apply(vars_train, 2, function(x) sd(x, na.rm = TRUE) != 0)]
  vars_test_clean  <- vars_test[,  apply(vars_test,  2, function(x) sd(x, na.rm = TRUE) != 0)]
  
  mds_train_clean <- mds_df_train[, apply(mds_df_train, 2, function(x) sd(x, na.rm = TRUE) != 0)]
  mds_test_clean  <- mds_df_test[,  apply(mds_df_test,  2, function(x) sd(x, na.rm = TRUE) != 0)]
  
  # Step 3: Compute correlation matrices
  cor_train <- cor(vars_train_clean, mds_train_clean, use = "pairwise.complete.obs")
  cor_test  <- cor(vars_test_clean,  mds_test_clean,  use = "pairwise.complete.obs")
  
  # Step 4: Keep only common rows (variables) and align order
  common_vars <- intersect(rownames(cor_train), rownames(cor_test))
  cor_train_filtered <- cor_train[common_vars, , drop = FALSE]
  cor_test_reordered <- cor_test[common_vars, , drop = FALSE]
  
  # ✅ Step 5: Compute importance and optionally select top N variables
  importance_train <- rowSums(abs(cor_train_filtered))
  importance_test  <- rowSums(abs(cor_test_reordered))
  importance_mean  <- (importance_train + importance_test) / 2
  
  # Sort by importance
  ordered_vars <- names(sort(importance_mean, decreasing = TRUE))
  
  if (!is.null(top_n)) {
    top_vars <- ordered_vars[1:min(top_n, length(ordered_vars))]
    cor_train_filtered <- cor_train_filtered[top_vars, , drop = FALSE]
    cor_test_reordered <- cor_test_reordered[top_vars, , drop = FALSE]
    importance_mean <- importance_mean[top_vars]
  } else {
    # Reorder all by importance even if not filtering
    cor_train_filtered <- cor_train_filtered[ordered_vars, , drop = FALSE]
    cor_test_reordered <- cor_test_reordered[ordered_vars, , drop = FALSE]
    importance_mean <- importance_mean[ordered_vars]
  }
  
  # ✅ Step 6: Optionally rename variables using mapping
  if (rename_vars && !is.null(var_name_map)) {
    # Expect var_name_map as a *named vector or list* where names = old, values = new
    old_names <- rownames(cor_train_filtered)
    new_names <- ifelse(old_names %in% names(var_name_map),
                        var_name_map[old_names],
                        old_names)
    
    rownames(cor_train_filtered) <- new_names
    rownames(cor_test_reordered) <- new_names
  }
  
  # Step 7: Plot heatmaps side by side
  p1 <- pheatmap(cor_train_filtered, silent = TRUE,
                 main = title_train,
                 cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
  
  p2 <- pheatmap(cor_test_reordered, silent = TRUE,
                 main =title_test,
                 cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
  
  # Step 8: Display side-by-side
  grid.newpage()
  grid.draw(gridExtra::arrangeGrob(p1[[4]], p2[[4]], ncol = 2))
  
  # ✅ Step 9: Return numerical results
  return(list(
    cor_train = cor_train_filtered,
    cor_test = cor_test_reordered,
    importance = importance_mean,
    p_train = p1,
    p_test  = p2
  ))
}

fit_and_relabel_gmm_by_wellbeing <- function(mds_points, warwick_wellbeing, G = 2, seed = 3) {
  set.seed(seed)
  
  # Fit GMM model
  gmm_model <- Mclust(mds_points, G = G, verbose = FALSE)
  
  # Create data frame from MDS points and add cluster labels and wellbeing scores
  df <- as.data.frame(mds_points)
  colnames(df) <- paste0("Dim", seq_len(ncol(df)))
  df$OriginalCluster <- gmm_model$classification
  df$warwick_wellbeing <- warwick_wellbeing
  
  # Compute mean wellbeing per cluster
  cluster_means <- aggregate(warwick_wellbeing ~ OriginalCluster, data = df, mean)
  
  # Determine new order: cluster with highest mean wellbeing gets label 1, etc.
  new_order <- cluster_means$OriginalCluster[order(-cluster_means$warwick_wellbeing)]
  label_map <- setNames(seq_along(new_order), new_order)
  
  # Remap classification labels
  gmm_model$classification <- label_map[as.character(gmm_model$classification)]
  gmm_model$z <- gmm_model$z[, new_order]
  
  return(gmm_model)
}


plot_mds_clusters <- function(
    mds_result_points,
    gmm_model,
    title = "MDS with GMM Clustering",
    cluster_colours = NULL  # 👈 new argument
) {
  # Build data frame
  mds_df <- as.data.frame(mds_result_points)
  colnames(mds_df) <- paste0("Dim", 1:ncol(mds_df))  # Automatically handle >2D
  
  # Add original GMM classification
  mds_df$Cluster <- as.factor(gmm_model$classification)
  mds_df$Uncertainty <- gmm_model$uncertainty
  
  # Base plot
  p <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Cluster)) +
    geom_point(aes(alpha = 1 - Uncertainty), size = 3) +
    stat_ellipse(aes(group = Cluster), color = "black", linetype = "dashed", linewidth = 1) +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    theme_minimal() +
    labs(
      title = title,
      x = "Dimension 1",
      y = "Dimension 2",
      color = "Cluster"
    )
  
  # Apply custom cluster colours if provided
  if (!is.null(cluster_colours)) {
    p <- p + scale_color_manual(values = cluster_colours)
  } else {
    p <- p + scale_color_brewer(palette = "Set1")
  }
  
  return(p)
}

################################################################################
# GMM functions
################################################################################

fit_and_relabel_gmm_by_dim1 <- function(mds_points, G = 2, seed = 3) {
  set.seed(seed)
  
  # Fit GMM model
  gmm_model <- Mclust(mds_points, G = G, verbose = FALSE)
  
  # Create data frame from MDS points and add cluster labels
  df <- as.data.frame(mds_points)
  colnames(df) <- paste0("Dim", seq_len(ncol(df)))
  df$OriginalCluster <- gmm_model$classification
  
  # Compute mean Dim1 per cluster
  cluster_means <- aggregate(Dim1 ~ OriginalCluster, data = df, mean)
  
  # Determine new order: cluster with lowest mean Dim1 gets label 1, etc.
  new_order <- cluster_means$OriginalCluster[order(cluster_means$Dim1)]
  label_map <- setNames(seq_along(new_order), new_order)
  
  # Remap classification labels
  gmm_model$classification <- label_map[as.character(gmm_model$classification)]
  gmm_model$z <- gmm_model$z[,new_order]
  
  return(gmm_model)
}

fit_and_relabel_gmm_by_wellbeing <- function(mds_points, warwick_wellbeing, G = 2, seed = 3) {
  set.seed(seed)
  
  # Fit GMM model
  gmm_model <- Mclust(mds_points, G = G, verbose = FALSE)
  
  # Create data frame from MDS points and add cluster labels and wellbeing scores
  df <- as.data.frame(mds_points)
  colnames(df) <- paste0("Dim", seq_len(ncol(df)))
  df$OriginalCluster <- gmm_model$classification
  df$warwick_wellbeing <- warwick_wellbeing
  
  # Compute mean wellbeing per cluster
  cluster_means <- aggregate(warwick_wellbeing ~ OriginalCluster, data = df, mean)
  
  # Determine new order: cluster with highest mean wellbeing gets label 1, etc.
  new_order <- cluster_means$OriginalCluster[order(-cluster_means$warwick_wellbeing)]
  label_map <- setNames(seq_along(new_order), new_order)
  
  # Remap classification labels
  gmm_model$classification <- label_map[as.character(gmm_model$classification)]
  gmm_model$z <- gmm_model$z[, new_order]
  
  return(gmm_model)
}

################################################################################
# elastic net regularised regression functions
################################################################################

################################################################################
# choose alpha and lambda parameters

get_best_alpha_lambda <- function(params_path, whichReg_string) {
  # Step 1: load mse results
  mse_file <- list.files(
    params_path, 
    pattern = paste0("boostrap_results_", whichReg_string, "_df_mse"), 
    full.names = TRUE
  )
  mse_df <- read.csv(mse_file)
  
  # Step 2: find alpha with lowest mean_mse
  best_alpha <- mse_df$alpha[which.min(mse_df$mean_mse)]
  
  # Step 3: load lambda results
  lambda_file <- list.files(
    params_path, 
    pattern = paste0("boostrap_results_", whichReg_string, "_df_lambda_1se"), 
    full.names = TRUE
  )
  lambda_df <- read.csv(lambda_file)
  
  # Step 4: select mean_lambda for best alpha
  best_lambda <- lambda_df$mean_lambda[which.min(abs(lambda_df$alpha - best_alpha))]
  
  # Return results
  return(data.frame(alpha = best_alpha, mean_lambda = best_lambda))
}



fit_elastic_net_CV <- function(x, y, family = family) {
  
  x <- x %>% as.matrix()
  
  # Initialize output
  list.of.fits        <- list()
  results             <- data.frame()
  coef_matrices       <- list()
  selected_predictors <- list()
  
  for (i in 0:10) {
    alpha_val <- i / 10
    fit.name  <- paste0("alpha", alpha_val)
    
    print(fit.name)
    
    # Fit model
    fit <- cv.glmnet(
      x, y,
      type.measure = "mse",
      alpha        = alpha_val,
      family       = family
    )
    
    list.of.fits[[fit.name]] <- fit
    
    
    # Extract the cvm corresponding to lambda.1se
    lambda_1se  <- fit$lambda.1se
    lambda_idx  <- which(fit$lambda == lambda_1se)
    mse         <- fit$cvm[lambda_idx]
    
    # Store MSE
    temp     <- data.frame(alpha = alpha_val, 
                           mse = mse, 
                           fit.name = fit.name,
                           lambda_1se = lambda_1se)
    results  <- rbind(results, temp)
    
    
    # Extract non-zero coefficients (excluding intercept)
    
    coef_matrix    <- as.matrix(coef(fit, s = "lambda.1se"))
    non_zero_coefs <- rownames(coef_matrix)[coef_matrix[, 1] != 0]
    non_zero_coefs <- setdiff(non_zero_coefs, "(Intercept)")
    
    coef_matrices[[fit.name]]       <- coef_matrix
    selected_predictors[[fit.name]] <- non_zero_coefs
  }
  
  return(list(
    fits       = list.of.fits,
    results    = results,
    coef_matrices = coef_matrices,
    predictors    = selected_predictors
  ))
  
}

bootstrap_alpha_selection <- function(x, y, B = 50, family = "gaussian") {
  
  alpha_choices     <- seq(0, 1, by = 0.1)
  alpha_selection   <- numeric(B)
  mse_matrix        <- matrix(NA, nrow = B, ncol = length(alpha_choices))
  lambda_1se_matrix <- matrix(NA, nrow = B, ncol = length(alpha_choices))
  
  alpha_labels <- paste0("alpha", alpha_choices)
  
  colnames(mse_matrix)        <- alpha_labels
  colnames(lambda_1se_matrix) <- alpha_labels
  
  for (b in 1:B) {
    set.seed(1000 + b)  # Make reproducible but vary seeds
    
    #### fit elastic net models
    
    fit <- fit_elastic_net_CV(x, y, family = family)
    
    best_row           <- fit$results[which.min(fit$results$mse), ]
    best_alpha         <- best_row$alpha
    alpha_selection[b] <- best_alpha
    
    # Match each alpha's metrics
   # fit_order <- match(alpha_labels, fit$results$fit.name)
    
    mse_matrix[b, ] <- fit$results$mse
    lambda_1se_matrix[b, ] <- fit$results$lambda_1se
    
    cat("Iteration", b, "best alpha =", best_alpha, "\n")
  }
  
  # Frequency and summary
  alpha_freq <- table(alpha_selection)
  avg_mse    <- colMeans(mse_matrix, na.rm = TRUE)
  sd_mse     <- apply(mse_matrix, 2, sd, na.rm = TRUE)
  
  avg_lambda_1se <- colMeans(lambda_1se_matrix, na.rm = TRUE)
  sd_lambda_1se  <- apply(lambda_1se_matrix, 2, sd, na.rm = TRUE)
  
  return(list(
    alpha_freq = alpha_freq,
    avg_mse = avg_mse,
    sd_mse = sd_mse,
    alpha_selection = alpha_selection,
    mse_matrix = mse_matrix,
    lambda_1se_matrix = lambda_1se_matrix,
    avg_lambda_1se = avg_lambda_1se,
    sd_lambda_1se = sd_lambda_1se
  ))
}





# ===============================
# Functions
# ===============================

# Return cluster report as string
report_clusters <- function(train_gmm, test_gmm, train_data, test_data, gmm_name) {
  k <- length(unique(train_gmm$classification))  # detect number of clusters
  out <- paste0("\n----- ", gmm_name, " -----\n")
  
  # Train
  for (i in 1:k) {
    n_i <- sum(train_gmm$classification == i)
    out <- paste0(out,
                  "Train ", gmm_name, ", cluster ", i, " N: ", n_i,
                  " out of ", nrow(train_data),
                  " (", round(n_i * 100 / nrow(train_data), 2), "%)\n"
    )
  }
  
  out <- paste0(out, "\n")
  
  # Test
  for (i in 1:k) {
    n_i <- sum(test_gmm$classification == i)
    out <- paste0(out,
                  "Test ", gmm_name, ", cluster ", i, " N: ", n_i,
                  " out of ", nrow(test_data),
                  " (", round(n_i * 100 / nrow(test_data), 2), "%)\n"
    )
  }
  paste0(out, "\n")
}

# Return silhouette report as string
report_silhouette <- function(train_gmm, test_gmm, train_mds, test_mds, gmm_name) {
  sil_train_crisp <- mean(silhouette(train_gmm$classification, dist(train_mds$points))[, 3])
  sil_train_fuzzy <- SIL.F(Xca = train_mds$points, U = train_gmm$z, alpha = 1)
  
  sil_test_crisp <- mean(silhouette(test_gmm$classification, dist(test_mds))[, 3])
  sil_test_fuzzy <- SIL.F(Xca = test_mds, U = test_gmm$z, alpha = 1)
  
  paste0(
    "\n===== Silhouette results for ", gmm_name, " =====\n",
    "Train silhouette (crisp): ", round(sil_train_crisp, 4), "\n",
    "Train silhouette (fuzzy): ", round(sil_train_fuzzy, 4), "\n\n",
    "Test silhouette (crisp): ", round(sil_test_crisp, 4), "\n",
    "Test silhouette (fuzzy): ", round(sil_test_fuzzy, 4), "\n\n\n"
  )
}


plot_bootstrap_alpha_results <- function(bootstrap_results, save_path = "./../data_processed/explore_elastic_net_parameters/") {
  
  # Capture object name automatically
  results_name <- deparse(substitute(bootstrap_results))
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  
  alpha_values <- seq(0, 1, by = 0.1)
  
  # -------------------------------
  # Plot 1: Mean MSE vs Alpha
  df_mse <- data.frame(
    alpha = alpha_values,
    mean_mse = bootstrap_results$avg_mse,
    sd_mse = bootstrap_results$sd_mse
  )
  
  p1 <- ggplot(df_mse, aes(x = alpha, y = mean_mse)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_mse - sd_mse, ymax = mean_mse + sd_mse), width = 0.02) +
    labs(
      title = "Mean MSE vs Alpha (with SD error bars)",
      x = "Alpha",
      y = "Mean MSE"
    ) +
    theme_minimal()
  
  print(p1)
  ggsave(file.path(save_path, paste0(timestamp, "_", results_name, "_plot_mse.png")),
         plot = p1, width = 7, height = 5, dpi = 300, bg = "white")
  
  # -------------------------------
  # Plot 2: Frequency of Best Alpha Values
  df_freq <- as.data.frame(bootstrap_results$alpha_freq)
  colnames(df_freq) <- c("alpha", "frequency")
  df_freq$alpha <- as.numeric(as.character(df_freq$alpha))
  
  p2 <- ggplot(df_freq, aes(x = factor(alpha), y = frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = "Frequency of Best Alpha Values Across Bootstraps",
      x = "Alpha",
      y = "Frequency"
    ) +
    theme_minimal()
  
  print(p2)
  ggsave(file.path(save_path, paste0(timestamp, "_", results_name, "_plot_freq.png")),
         plot = p2, width = 7, height = 5, dpi = 300, bg = "white")
  
  # -------------------------------
  # Plot 3: Mean Lambda 1SE vs Alpha
  df_lambda_1se <- data.frame(
    alpha = alpha_values,
    mean_lambda = bootstrap_results$avg_lambda_1se,
    sd_lambda = bootstrap_results$sd_lambda_1se
  )
  
  p3 <- ggplot(df_lambda_1se, aes(x = alpha, y = mean_lambda)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_lambda - sd_lambda, ymax = mean_lambda + sd_lambda), width = 0.02) +
    labs(
      title = "Mean Lambda.1se vs Alpha (with SD error bars)",
      x = "Alpha",
      y = "Mean Lambda.1se"
    ) +
    theme_minimal()
  
  print(p3)
  ggsave(file.path(save_path, paste0(timestamp, "_", results_name, "_plot_lambda.png")),
         plot = p3, width = 7, height = 5, dpi = 300, bg = "white")
  
  # -------------------------------
  # Save DataFrames as CSVs
  write.csv(df_mse,
            file = file.path(save_path, paste0(timestamp, "_", results_name, "_df_mse.csv")),
            fileEncoding = "UTF-8", row.names = FALSE)
  
  write.csv(df_freq,
            file = file.path(save_path, paste0(timestamp, "_", results_name, "_df_freq.csv")),
            fileEncoding = "UTF-8", row.names = FALSE)
  
  write.csv(df_lambda_1se,
            file = file.path(save_path, paste0(timestamp, "_", results_name, "_df_lambda_1se.csv")),
            fileEncoding = "UTF-8", row.names = FALSE)
  
  message("Saved plots and results for ", results_name, " at ", save_path)
  
  return(
    list(
      df_mse = df_mse,
      df_freq = df_freq,
      df_lambda_1se = df_lambda_1se
    )
  )
}

fit_elastic_net_train <- function(x, y, alpha_val = 0.5, family = "gaussian") {
  
  x <- x %>% as.matrix()
  
  # Fit model
  fit <- cv.glmnet(
    x, y,
    type.measure = "mse",
    alpha        = alpha_val,
    family       = family
  )
  
  # Extract non-zero coefficients (excluding intercept)
  
  coef_matrix    <- as.matrix(coef(fit, s = "lambda.1se"))
  non_zero_coefs <- rownames(coef_matrix)[coef_matrix[, 1] != 0]
  non_zero_coefs <- setdiff(non_zero_coefs, "(Intercept)")
  selected_predictors <- non_zero_coefs
  
  return(list(
    fit         = fit,
    coef_matrix = coef_matrix,
    predictors  = selected_predictors
  ))
}

fit_elastic_net <- function(x, y, alpha_val, lambda_val, family = "gaussian") {
  
  x <- x %>% as.matrix()
  
  # Fit elastic net model with fixed alpha and lambda
  fit <- glmnet(
    x, y,
    alpha  = alpha_val,
    lambda = lambda_val,
    family = family
  )
  
  # Extract non-zero coefficients (excluding intercept)
  coef_matrix    <- as.matrix(coef(fit))
  non_zero_coefs <- rownames(coef_matrix)[coef_matrix[, 1] != 0]
  non_zero_coefs <- setdiff(non_zero_coefs, "(Intercept)")
  selected_predictors <- non_zero_coefs
  
  return(list(
    fit         = fit,
    coef_matrix = coef_matrix,
    predictors  = selected_predictors
  ))
}




bootstrap_feature_importance <- function(x, y, alpha_val, lambda_val, family = "gaussian",
                                         n_bootstrap = 1000, seed = 123,
                                         feature_selection_proportion = 0.95,
                                         colours = c("blue", "red"),
                                         which_dataset = "dataset") {
  set.seed(seed)
  x <- as.matrix(x)
  n <- nrow(x)
  
  coef_list <- vector("list", n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    if (i %% 100 == 0 || i == 1 || i == n_bootstrap) {
      message(paste("Bootstrap iteration:", i, "/", n_bootstrap))
    }
    
    # Bootstrap sample
    sample_idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[sample_idx, , drop = FALSE]
    y_boot <- y[sample_idx]
    
    model <- fit_elastic_net(x_boot, y_boot, alpha_val, lambda_val, family)
    coefs <- as.vector(model$coef_matrix)
    names(coefs) <- rownames(model$coef_matrix)
    coef_list[[i]] <- coefs
  }
  
  # Combine coefficients
  coef_df <- bind_rows(lapply(coef_list, function(coef_vec) {
    tibble(feature = names(coef_vec), coefficient = coef_vec)
  }), .id = "bootstrap") %>%
    filter(feature != "(Intercept)")
  
  # Summarize
  summary_df <- coef_df %>%
    group_by(feature) %>%
    summarise(
      mean = mean(coefficient),
      lower = quantile(coefficient, 0.025),
      upper = quantile(coefficient, 0.975),
      prop_nonzero = mean(coefficient != 0),
      .groups = "drop"
    ) %>%
    filter(prop_nonzero >= feature_selection_proportion) %>%
    mutate(direction = ifelse(mean >= 0, "positive", "negative"))
  
  # Plot
  p <- ggplot(summary_df, aes(x = reorder(feature, mean), y = mean, color = direction)) +
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("positive" = colours[1], "negative" = colours[2])) +
    coord_flip() +
    labs(
      title = paste("Bootstrapped Feature Importance\n(alpha = ", alpha_val,
                    ", lambda = ", round(lambda_val, 3), ")\n\n",
                    "Only features selected in\n≥", feature_selection_proportion, 
                    "% of bootstraps\n\n", which_dataset, sep = ""),
      x = "Features",
      y = "Mean Coefficient with 95% CI",
      color = "Direction"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))
  
  return(list(summary = summary_df, plot = p))
}


bootstrap_feature_importance_wrapper <- function(
    x_train, y_train, x_test, y_test,
    alpha_val, lambda_val, family = "gaussian",
    n_bootstrap = 1000, seed = 123,
    feature_selection_proportion = 0.95,
    colours = c("purple", "darkgreen")
) {
  
  # --- Train ---
  res_train <- bootstrap_feature_importance(x_train, y_train,
                                            alpha_val, lambda_val, family,
                                            n_bootstrap, seed,
                                            feature_selection_proportion,
                                            colours, "Train")
  
  # --- Test (all features) ---
  res_test <- bootstrap_feature_importance(x_test, y_test,
                                           alpha_val, lambda_val, family,
                                           n_bootstrap, seed,
                                           feature_selection_proportion,
                                           colours, "Test (all features)")
  
  # --- Test (restricted to train-selected features) ---
  train_features <- res_train$summary$feature
  keep_idx <- colnames(x_test) %in% train_features
  
  if (any(keep_idx)) {
    res_test_restricted <- bootstrap_feature_importance(x_test[, keep_idx, drop = FALSE], y_test,
                                                        alpha_val, lambda_val, family,
                                                        n_bootstrap, seed,
                                                        feature_selection_proportion,
                                                        colours, "Test (restricted to Train features)")
  } else {
    res_test_restricted <- NULL
    message("No predictor variables from training survived into the restricted test set. Returning only Train and Test (all features).")
  }
  
  # Return all three
  return(list(
    train = res_train,
    test_all = res_test,
    test_restricted = res_test_restricted
  ))
}

