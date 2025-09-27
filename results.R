library(ggplot2)
library(caret)
library(pampam)
library(pROC)
library(MLmetrics)
library(glassoFast)

load("C:/Users/tomas/Desktop/Github/tese/glioma-RNASeq-2021-classification.RData")

astro_RNA$y <- 1
gbm_RNA$y <- 2
oligo_RNA$y <- 3
# Merge datasets
data <- rbind(astro_RNA, gbm_RNA, oligo_RNA)
# Response and predictors
y <- as.factor(data$y)
x <- as.matrix(subset(data, select = -c(names, sample.type, y)))

acc <- function(pred, truth) { 
  mean(pred == truth)
}

fit_to_coefmatrix <- function(fit, lambda) {
  coef <- coef(fit, s = lambda)
  ncols <- length(coef)
  nrows <- length(coef[[1]])
  matrix <- matrix(0, nrow = nrows, ncol = ncols)
  for (i in 1:ncols) {
    matrix[, i] <- coef[[i]][,1]
    i <- i + 1
  }
  rownames(matrix) <- rownames(coef[[1]])
  return(matrix)
}

build_W_from_cor <- function(X, cap = 0.95) {
  R <- cor(X)
  R[!is.finite(R)] <- 0
  R <- pmin(pmax(R, -cap), cap)
  
  denom <- 1 - R^2
  
  # Off-diagonal entries: -2 * rho / (1 - rho^2)
  W <- -2 * R / denom
  diag(W) <- 0  # clear diagonals for now
  
  # Diagonal entries: 2 * sum_{s!=i} 1 / (1 - rho_{is}^2)
  diag(W) <- 2 * (rowSums(1 / denom) - 1 / diag(denom))
  
  return(W)  # ensure W is returned as a matrix
}

inv_sqrt_W <- function(W, tau = 1e-3) {
  Wreg <- W + diag(tau, ncol(W))
  e <- eigen(Wreg, symmetric = TRUE) #eigen takes too long for all the data
  d <- pmax(e$values, 1e-10)
  Wi2 <- e$vectors %*% diag(1 / sqrt(d)) %*% t(e$vectors)
  dimnames(Wi2) <- dimnames(W)   # <-- preserve feature names
  
  return(Wi2)
}

weighted_f1 <- function(y, pred){
  # F1 per class
  classes <- levels(y)
  f1_per_class <- sapply(classes, function(cl) {
    F1_Score(y_true = y, y_pred = pred, positive = cl)
  })
  
  # Weighted F1 (weighted by class frequency in y_test)
  class_weights <- table(y) / length(y)
  f1_weighted3 <- sum(f1_per_class * class_weights)
}


# --- Run experiments ---
n_reps <- 100
outdir1 <- "C:/Users/tomas/Desktop/NÃO APAGAR/results"

results_te <- matrix(NA, nrow = n_reps, ncol = 10, dimnames = list(NULL, c("Lasso", "Ridge", "AdaptiveLasso", "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))
results_tr <- matrix(NA, nrow = n_reps, ncol = 10, dimnames = list(NULL, c("Lasso", "Ridge",  "AdaptiveLasso", "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))
results_n_zero <- matrix(NA, nrow = n_reps, ncol = 10, dimnames = list(NULL, c("Lasso", "Ridge", "AdaptiveLasso",  "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))
results_auc <- matrix(NA, nrow = n_reps, ncol = 10, dimnames = list(NULL, c("Lasso", "Ridge", "AdaptiveLasso",  "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))
results_f1 <- matrix(NA, nrow = n_reps, ncol = 10, dimnames = list(NULL, c("Lasso", "Ridge", "AdaptiveLasso",  "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))
results_duration <- matrix(NA, nrow = n_reps, ncol = 10, dimnames = list(NULL, c("Lasso", "Ridge", "AdaptiveLasso",  "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))
results_duration_total <- matrix(NA, nrow = n_reps, ncol = 8, dimnames = list(NULL, c("AdaptiveLasso", "ElasticNet","ElasticNet+Corr","ElasticNet+Corr+Adaptive","Glasso&CorW", "Glasso&CorW&Adaptive", "Glasso&Corwi", "Glasso&Corwi&Adaptive")))

results_cm_l<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_r<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_al<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_enet<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_c<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_c_al<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_gl_W<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_gl_W_al<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_gl_wi<-matrix(NA, nrow = n_reps, ncol = 9)
results_cm_gl_wi_al<-matrix(NA, nrow = n_reps, ncol = 9)

for (rep in 1:n_reps) {
  # Load data
  load("C:/Users/tomas/Desktop/Github/tese/glioma-RNASeq-2021-classification.RData")
  
  astro_RNA$y <- 1
  gbm_RNA$y <- 2
  oligo_RNA$y <- 3
  # Merge datasets
  data <- rbind(astro_RNA, gbm_RNA, oligo_RNA)
  # Response and predictors
  y <- as.factor(data$y)
  x <- as.matrix(subset(data, select = -c(names, sample.type, y)))
  
  # Train/Test Split
  id_tr  <- createDataPartition(y, p = 0.8, list = FALSE)
  Xtr_raw <-x[id_tr, , drop = FALSE]
  Xte_raw <- x[-id_tr, , drop = FALSE]
  y_tr    <- y[id_tr]
  y_te    <- y[-id_tr]
  
  set.seed(rep)  # ensures reproducibility
  foldid <- sample(rep(1:10, length.out = length(y_tr)))  # 10-fold CV
  
  nzv <- nearZeroVar(Xtr_raw)                   # detect
  Xtr_raw_filtered <- Xtr_raw[, -nzv]           # remove
  Xte_raw_filtered <- Xte_raw[, -nzv]
  
  # Standardize
  mu  <- colMeans(Xtr_raw_filtered)
  sdv <- apply(Xtr_raw_filtered, 2, sd); sdv[sdv == 0] <- 1
  Xtr <- scale(Xtr_raw_filtered, center = mu, scale = sdv)
  Xte <- scale(Xte_raw_filtered, center = mu, scale = sdv)
  
  #Lasso
  duration_lasso<-as.numeric(system.time(fit_lasso <- pampam::cv.glmnet(
    Xtr, y_tr, family = "multinomial",
    type.measure = "class", standardize = FALSE,
    alpha = 1, foldid = foldid))[3])
  
  pred_l_te <- predict(fit_lasso, Xte, s="lambda.min", type="class")
  pred_l_tr <- predict(fit_lasso, Xtr, s="lambda.min", type="class")
  acc_l_te <- acc(pred_l_te, y_te)
  acc_l_tr <- acc(pred_l_tr, y_tr)  
  auc_l <- as.numeric(multiclass.roc(y_te, as.vector(pred_l_te))$auc)
  f1_l<-weighted_f1(y_te, pred_l_te)
  cm_l <- as.vector(confusionMatrix(as.factor(pred_l_te), y_te)$table)
  
  coeffs_l <- fit_to_coefmatrix(fit_lasso, fit_lasso$lambda.min)[-1, , drop = FALSE]
  nz_l <- as.matrix(coeffs_l[rowSums(abs(coeffs_l)) != 0, , drop = FALSE])
  nz_l_n <- nrow(nz_l)
  
  #Ridge
  duration_ridge<-as.numeric(system.time(fit_ridge <- pampam::cv.glmnet(
    Xtr, y_tr, family = "multinomial",
    type.measure = "class", standardize = FALSE,
    alpha = 0, foldid = foldid))[3])
  
  pred_r_te <- predict(fit_ridge, Xte, s="lambda.min", type="class")
  pred_r_tr <- predict(fit_ridge, Xtr, s="lambda.min", type="class")
  acc_r_te <- acc(pred_r_te, y_te)
  acc_r_tr <- acc(pred_r_tr, y_tr)  
  auc_r <- as.numeric(multiclass.roc(y_te, as.vector(pred_r_te))$auc)
  f1_r<-weighted_f1(y_te, pred_r_te)
  cm_r <- as.vector(confusionMatrix(as.factor(pred_r_te), y_te)$table)
  
  coeffs_r <- fit_to_coefmatrix(fit_ridge, fit_ridge$lambda.min)[-1, , drop = FALSE]
  nz_r <- as.matrix(coeffs_r[rowSums(abs(coeffs_r)) != 0, , drop = FALSE])
  nz_r_n <- nrow(nz_r)
  
  #Adaptive LASSO
  weights <- 1 / abs(coeffs_r)
  weights <- pmin(weights, quantile(weights, 0.95))
  
  duration_al<-as.numeric(system.time(fit_alasso <- pampam::cv.glmnet(
    Xtr, y_tr, family="multinomial", 
    type.measure = "class", standardize = FALSE,
    alpha=1, penalty.factor = weights,
    foldid = foldid))[3])
  
  pred_al_te <- predict(fit_alasso, Xte, s="lambda.min", type="class")
  pred_al_tr <- predict(fit_alasso, Xtr, s="lambda.min", type="class")
  acc_al_te <- acc(pred_al_te, y_te)
  acc_al_tr <- acc(pred_al_tr, y_tr)  
  auc_al <- as.numeric(multiclass.roc(y_te, as.vector(pred_al_te))$auc)
  f1_al<-weighted_f1(y_te, pred_al_te)
  cm_al <- as.vector(confusionMatrix(as.factor(pred_al_te), y_te)$table)
  
  coeffs_al <- fit_to_coefmatrix(fit_alasso, fit_alasso$lambda.min)[-1, , drop = FALSE]
  nz_al <- as.matrix(coeffs_al[rowSums(abs(coeffs_al)) != 0, , drop = FALSE])
  nz_al_n <- nrow(nz_al)
  
  # Elastic Net
  duration_enet<-as.numeric(system.time(fit_elastic <- pampam::cv.glmnet(
    Xtr, y_tr, family = "multinomial",
    type.measure = "class", standardize = FALSE,
    alpha = 0.13, foldid = foldid))[3])

  pred_enet_te <- predict(fit_elastic, Xte, s="lambda.min", type="class")
  pred_enet_tr <- predict(fit_elastic, Xtr, s="lambda.min", type="class")
  acc_enet_te <- acc(pred_enet_te, y_te)
  acc_enet_tr <- acc(pred_enet_tr, y_tr)  
  auc_enet <- as.numeric(multiclass.roc(y_te, as.vector(pred_enet_te))$auc)
  f1_enet<-weighted_f1(y_te, pred_enet_te)
  cm_enet <- as.vector(confusionMatrix(as.factor(pred_enet_te), y_te)$table)
  
  coeffs_enet <- fit_to_coefmatrix(fit_elastic, fit_elastic$lambda.min)[-1, , drop = FALSE]
  nz_enet <- as.matrix(coeffs_enet[rowSums(abs(coeffs_enet)) != 0, , drop = FALSE])
  nz_enet_names <- rownames(nz_enet)
  nz_enet_n <- nrow(nz_enet)
  
  Xtr_reduced <- Xtr[, nz_enet_names, drop = FALSE]
  Xte_reduced <- Xte[, nz_enet_names, drop = FALSE]

  # Elastic Net + Correlation-Based
  W    <- build_W_from_cor(Xtr_reduced)
  Wi2  <- inv_sqrt_W(W)
  XtrS <- Xtr_reduced %*% Wi2
  XteS <- Xte_reduced %*% Wi2

  duration_c<-as.numeric(system.time(fit_corr <- pampam::cv.glmnet(
    XtrS, y_tr, family="multinomial", 
    type.measure = "class", standardize = FALSE,
    alpha=0, foldid = foldid))[3])
  
  pred_c_te <- predict(fit_corr, XteS, s="lambda.min", type="class")
  pred_c_tr <- predict(fit_corr, XtrS, s="lambda.min", type="class")
  acc_c_te <- acc(pred_c_te, y_te)
  acc_c_tr <- acc(pred_c_tr, y_tr) 
  auc_c <- as.numeric(multiclass.roc(y_te, as.vector(pred_c_te))$auc)
  f1_c<-weighted_f1(y_te, pred_c_te)
  cm_c <- as.vector(confusionMatrix(as.factor(pred_c_te), y_te)$table)
  
  coeffs_corr <- Wi2%*%fit_to_coefmatrix(fit_corr, fit_corr$lambda.min)[-1, , drop = FALSE]
  nz_c_n <- nrow(coeffs_corr)
  
  # Elastic Net + Correlation-Based + Adaptive LASSO
  weights <- 1 / abs(coeffs_corr)
  weights <- pmin(weights, quantile(weights, 0.95))
  
  duration_c_al<-as.numeric(system.time(fit_corr_alasso <- pampam::cv.glmnet(
    Xtr_reduced, y_tr, family="multinomial", 
    type.measure = "class", standardize = FALSE,
    alpha=1, penalty.factor = weights,
    foldid = foldid))[3])
  
  pred_c_al_te <- predict(fit_corr_alasso, Xte_reduced, s="lambda.min", type="class")
  pred_c_al_tr <- predict(fit_corr_alasso, Xtr_reduced, s="lambda.min", type="class")
  acc_c_al_te <- acc(pred_c_al_te, y_te)
  acc_c_al_tr <- acc(pred_c_al_tr, y_tr)  
  auc_c_al <- as.numeric(multiclass.roc(y_te, as.vector(pred_c_al_te))$auc)
  f1_c_al<-weighted_f1(y_te, pred_c_al_te)
  cm_c_al <- as.vector(confusionMatrix(as.factor(pred_c_al_te), y_te)$table)
  
  coeffs_corr_alasso <- fit_to_coefmatrix(fit_corr_alasso, fit_corr_alasso$lambda.min)[-1, , drop = FALSE]
  nz_corr_alasso <- as.matrix(coeffs_corr_alasso[rowSums(abs(coeffs_corr_alasso)) != 0, , drop = FALSE])
  nz_cor_alasso_n <- nrow(nz_corr_alasso)
  
  rm(astro_RNA, data, gbm_RNA, oligo_RNA, unclass_RNA, x, nzv, y)
  gc()
  
  cov_duration<-system.time(cov_matrix<-cov(Xtr))[3]
  graph_duration<-system.time(graph<-glassoFast(cov_matrix, rho = 0.95))[3]
  wi_duration<-system.time(wi<-graph$wi)[3]
  
  rm(graph, cov_matrix)
  gc()
  
  degree <- rowSums(wi != 0) - 1
  selected_features <- which(degree > 0)
  
  wi_sub<-wi[selected_features, selected_features]
  rm(wi)
  gc()
  
  Xtr_reduced <- Xtr[, selected_features]
  Xte_reduced <- Xte[, selected_features]
  
  # Standardize again
  
  mu  <- colMeans(Xtr_reduced)
  sdv <- apply(Xtr_reduced, 2, sd); sdv[sdv == 0] <- 1
  Xtr_reduced <- scale(Xtr_reduced, center = mu, scale = sdv)
  Xte_reduced <- scale(Xte_reduced, center = mu, scale = sdv)
  
  #Graphical Lasso + Correlation_wi
  Wi2  <- inv_sqrt_W(wi_sub)
  XtrS_raw <- Xtr_reduced %*% Wi2
  XteS_raw <- Xte_reduced %*% Wi2
  
  mu  <- colMeans(XtrS_raw)
  sdv <- apply(XtrS_raw, 2, sd); sdv[sdv == 0] <- 1
  XtrS <- scale(Xtr_reduced, center = mu, scale = sdv)
  XteS <- scale(Xte_reduced, center = mu, scale = sdv)
  
  duration_gl_wi<-as.numeric(system.time(fit_gl_wi <- pampam::cv.glmnet(XtrS, y_tr, family="multinomial", 
                                                                        type.measure = "class", standardize = FALSE,
                                                                        alpha=0, foldid = foldid))[3])
  
  pred_gl_wi_te <- predict(fit_gl_wi, XteS, s="lambda.min", type="class")
  pred_gl_wi_tr <- predict(fit_gl_wi, XtrS, s="lambda.min", type="class")
  acc_gl_wi_te <- acc(pred_gl_wi_te, y_te)
  acc_gl_wi_tr <- acc(pred_gl_wi_tr, y_tr) 
  auc_gl_wi <- as.numeric(multiclass.roc(y_te, as.vector(pred_gl_wi_te))$auc)
  f1_gl_wi <-weighted_f1(y_te, pred_gl_wi_te)
  cm_gl_wi <- as.vector(confusionMatrix(as.factor(pred_gl_wi_te), y_te)$table)
  
  coeffs_gl_wi <- Wi2%*%fit_to_coefmatrix(fit_gl_wi, fit_gl_wi$lambda.min)[-1, , drop = FALSE]
  nz_gl_wi_n <- nrow(coeffs_gl_wi)
  
  #Graphical Lasso + Correlation_wi + Adaptive
  weights <- 1 / (abs(coeffs_gl_wi) + 1e-8)
  weights <- pmin(weights, quantile(weights, 0.95))
  
  duration_gl_wi_al<-as.numeric(system.time(fit_gl_wi_al <- pampam::cv.glmnet(Xtr_reduced, y_tr, family="multinomial", 
                                                                              type.measure = "class", standardize = FALSE, 
                                                                              alpha=1, penalty.factor = weights, foldid = foldid))[3])
  
  pred_gl_wi_al_te <- predict(fit_gl_wi_al, Xte_reduced, s="lambda.min", type="class")
  pred_gl_wi_al_tr <- predict(fit_gl_wi_al, Xtr_reduced, s="lambda.min", type="class")
  acc_gl_wi_al_te <- acc(pred_gl_wi_al_te, y_te)
  acc_gl_wi_al_tr <- acc(pred_gl_wi_al_tr, y_tr) 
  auc_gl_wi_al <- as.numeric(multiclass.roc(y_te, as.vector(pred_gl_wi_al_te))$auc)
  f1_gl_wi_al <-weighted_f1(y_te, pred_gl_wi_al_te)
  cm_gl_wi_al <- as.vector(confusionMatrix(as.factor(pred_gl_wi_al_te), y_te)$table)
  
  coeffs_gl_wi_al <- fit_to_coefmatrix(fit_gl_wi_al, fit_gl_wi_al$lambda.min)[-1, , drop = FALSE]
  nz_gl_wi_al <- as.matrix(coeffs_gl_wi_al[rowSums(abs(coeffs_gl_wi_al)) != 0, , drop = FALSE])
  nz_gl_wi_al_n <- nrow(nz_gl_wi_al)
  
  
  #Graphical Lasso + Correlation_W
  W <- build_W_from_cor(Xtr_reduced)
  Wi2  <- inv_sqrt_W(W)
  XtrS_raw <- Xtr_reduced %*% Wi2
  XteS_raw <- Xte_reduced %*% Wi2
  
  mu  <- colMeans(XtrS_raw)
  sdv <- apply(XtrS_raw, 2, sd); sdv[sdv == 0] <- 1
  XtrS <- scale(Xtr_reduced, center = mu, scale = sdv)
  XteS <- scale(Xte_reduced, center = mu, scale = sdv)
  
  duration_gl_W<-as.numeric(system.time(fit_gl_W <- pampam::cv.glmnet(XtrS, y_tr, family="multinomial", 
                                                                      type.measure = "class", standardize = FALSE,
                                                                      alpha=0, foldid = foldid))[3])
  
  pred_gl_W_te <- predict(fit_gl_W, XteS, s="lambda.min", type="class")
  pred_gl_W_tr <- predict(fit_gl_W, XtrS, s="lambda.min", type="class")
  acc_gl_W_te <- acc(pred_gl_W_te, y_te)
  acc_gl_W_tr <- acc(pred_gl_W_tr, y_tr) 
  auc_gl_W <- as.numeric(multiclass.roc(y_te, as.vector(pred_gl_W_te))$auc)
  f1_gl_W <-weighted_f1(y_te, pred_gl_W_te)
  cm_gl_W <- as.vector(confusionMatrix(as.factor(pred_gl_W_te), y_te)$table)
  
  coeffs_gl_W <- Wi2%*%fit_to_coefmatrix(fit_gl_W, fit_gl_W$lambda.min)[-1, , drop = FALSE]
  nz_gl_W_n <- nrow(coeffs_gl_W)
  
  #Graphical Lasso + Correlation_W + Adaptive
  weights <- 1 / (abs(coeffs_gl_W) + 1e-8)
  weights <- pmin(weights, quantile(weights, 0.95))
  
  duration_gl_W_al<-as.numeric(system.time(fit_gl_W_al <- pampam::cv.glmnet(Xtr_reduced, y_tr, family="multinomial", 
                                                                            type.measure = "class", standardize = FALSE,
                                                                            alpha=1, penalty.factor = weights, foldid = foldid))[3])
  
  pred_gl_W_al_te <- predict(fit_gl_W_al, Xte_reduced, s="lambda.min", type="class")
  pred_gl_W_al_tr <- predict(fit_gl_W_al, Xtr_reduced, s="lambda.min", type="class")
  acc_gl_W_al_te <- acc(pred_gl_W_al_te, y_te)
  acc_gl_W_al_tr <- acc(pred_gl_W_al_tr, y_tr) 
  auc_gl_W_al <- as.numeric(multiclass.roc(y_te, as.vector(pred_gl_W_al_te))$auc)
  f1_gl_W_al <-weighted_f1(y_te, pred_gl_W_al_te)
  cm_gl_W_al <- as.vector(confusionMatrix(as.factor(pred_gl_W_al_te), y_te)$table)
  
  coeffs_gl_W_al <- fit_to_coefmatrix(fit_gl_W_al, fit_gl_W_al$lambda.min)[-1, , drop = FALSE]
  nz_gl_W_al <- as.matrix(coeffs_gl_W_al[rowSums(abs(coeffs_gl_W_al)) != 0, , drop = FALSE])
  nz_gl_W_al_n <- nrow(nz_gl_W_al)
  
  
  # Store
  results_te[rep, ] <- c(acc_l_te, acc_r_te, acc_al_te, acc_enet_te, acc_c_te, acc_c_al_te, acc_gl_W_te, acc_gl_W_al_te, acc_gl_wi_te, acc_gl_wi_al_te) #Test Accuracy
  results_tr[rep, ] <- c(acc_l_tr, acc_r_tr, acc_al_tr, acc_enet_tr, acc_c_tr, acc_c_al_tr, acc_gl_W_tr, acc_gl_W_al_tr, acc_gl_wi_tr, acc_gl_wi_al_tr) #Train Accuracy
  results_n_zero[rep, ] <- c(nz_l_n, nz_r_n, nz_al_n, nz_enet_n, nz_c_n, nz_cor_alasso_n, nz_gl_W_n, nz_gl_W_al_n, nz_gl_wi_n, nz_gl_wi_al_n) #Nonzero Coefficients Count
  results_auc[rep, ] <- c(auc_l, auc_r, auc_al, auc_enet, auc_c,auc_c_al, auc_gl_W, auc_gl_W_al, auc_gl_wi, auc_gl_wi_al) #AUC
  results_f1[rep, ] <- c(f1_l, f1_r, f1_al, f1_enet, f1_c, f1_c_al, f1_gl_W, f1_gl_W_al, f1_gl_wi, f1_gl_wi_al) #F1
  results_duration[rep, ] <- c(duration_lasso, duration_ridge, duration_al, duration_enet, duration_c, duration_c_al, duration_gl_W, duration_gl_W_al, duration_gl_wi, duration_gl_wi_al) #Duration (s)
  results_duration_total[rep, ] <- c(duration_ridge+duration_al, duration_enet, duration_enet+duration_c, duration_enet+duration_c+duration_c_al, cov_duration+graph_duration+wi_duration+duration_gl_W,
                                     cov_duration+graph_duration+wi_duration+duration_gl_W+duration_gl_W_al,
                                     cov_duration+graph_duration+wi_duration+duration_gl_wi,
                                     cov_duration+graph_duration+wi_duration+duration_gl_wi+duration_gl_wi_al) #Duration (s)

  results_cm_l[rep, ] <- cm_l
  results_cm_r[rep, ] <- cm_r
  results_cm_al[rep, ] <- cm_al
  results_cm_enet[rep, ] <- cm_enet
  results_cm_c[rep, ] <- cm_c
  results_cm_c_al[rep, ] <- cm_c_al
  results_cm_gl_W[rep, ] <- cm_gl_W
  results_cm_gl_W_al[rep, ] <- cm_gl_W_al
  results_cm_gl_wi[rep, ] <- cm_gl_wi
  results_cm_gl_wi_al[rep, ] <- cm_gl_wi_al
  
  write.csv(as.data.frame(nz_l), file = file.path(outdir1, "coeffs2/lasso", paste0("nz_l_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(nz_r), file = file.path(outdir1, "coeffs2/ridge", paste0("nz_r_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(nz_al), file = file.path(outdir1, "coeffs2/alasso", paste0("nz_al_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(nz_enet), file = file.path(outdir1, "coeffs2/enet", paste0("nz_enet_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(coeffs_corr), file = file.path(outdir1, "coeffs2/enet_c", paste0("nz_enet_c_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(nz_corr_alasso), file = file.path(outdir1, "coeffs2/enet_c_al", paste0("nz_enet_c_al_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(foldid), file = file.path("C:/Users/tomas/Desktop/NÃO APAGAR/foldids2", paste0("foldid_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(coeffs_gl_wi), file = file.path(outdir1, "coeffs2/glasso_wi", paste0("nz_gl_wi_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(nz_gl_wi_al), file = file.path(outdir1, "coeffs2/glasso_wi_al", paste0("nz_gl_wi_al_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(coeffs_gl_W), file = file.path(outdir1, "coeffs2/glasso_W", paste0("nz_gl_W_", rep, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(nz_gl_W_al), file = file.path(outdir1, "coeffs2/glasso_W_al", paste0("nz_gl_W_al_", rep, ".csv")), row.names = TRUE)
  
  if (rep %% 1 == 0) cat(rep, "done\n")
}

# --- Summarize ---
res_summary_te <- apply(results_te, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_tr <- apply(results_tr, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_n_zero <- apply(results_n_zero, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_auc <- apply(results_auc, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_f1 <- apply(results_f1, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_duration <- apply(results_duration, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_duration_total <- apply(results_duration_total, 2, function(v) c(mean=mean(v), sd=sd(v)))

res_summary_cm_l <- apply(results_cm_l, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_r <- apply(results_cm_r, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_al <- apply(results_cm_al, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_enet <- apply(results_cm_enet, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_c <- apply(results_cm_c, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_c_al <- apply(results_cm_c_al, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_gl_W <- apply(results_cm_gl_W, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_gl_W_al <- apply(results_cm_gl_W_al, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_gl_wi <- apply(results_cm_gl_wi, 2, function(v) c(mean=mean(v), sd=sd(v)))
res_summary_cm_gl_wi_al <- apply(results_cm_gl_wi_al, 2, function(v) c(mean=mean(v), sd=sd(v)))

setwd("C:/Users/tomas/Desktop/NÃO APAGAR/results/new correlation")

saveRDS(res_summary_te, "res_summary_te.rds")
saveRDS(res_summary_tr, "res_summary_tr.rds")
saveRDS(res_summary_n_zero, "res_summary_n_zero.rds")
saveRDS(res_summary_auc, "res_summary_auc.rds")
saveRDS(res_summary_f1, "res_summary_f1.rds")
saveRDS(res_summary_duration, "res_summary_duration.rds")
saveRDS(res_summary_duration_total, "res_summary_duration_total.rds")

saveRDS(results_te, "results_te.rds")
saveRDS(results_tr, "results_tr.rds")
saveRDS(results_n_zero, "results_n_zero.rds")
saveRDS(results_auc, "results_auc.rds")
saveRDS(results_f1, "results_f1.rds")
saveRDS(results_duration, "results_duration.rds")
saveRDS(results_duration_total, "results_duration_total.rds")

saveRDS(res_summary_cm_l, "res_summary_cm_l.rds")
saveRDS(res_summary_cm_r, "res_summary_cm_r.rds")
saveRDS(res_summary_cm_al, "res_summary_cm_al.rds")
saveRDS(res_summary_cm_enet, "res_summary_cm_enet.rds")
saveRDS(res_summary_cm_c, "res_summary_cm_c.rds")
saveRDS(res_summary_cm_c_al, "res_summary_cm_c_al.rds")
saveRDS(res_summary_cm_gl_W, "res_summary_cm_gl_W.rds")
saveRDS(res_summary_cm_gl_W_al, "res_summary_cm_gl_W_al.rds")
saveRDS(res_summary_cm_gl_wi, "res_summary_cm_gl_wi.rds")
saveRDS(res_summary_cm_gl_wi_al, "res_summary_cm_gl_wi_al.rds")

saveRDS(results_cm_l, "results_cm_l.rds")
saveRDS(results_cm_r, "results_cm_r.rds")
saveRDS(results_cm_al, "results_cm_al.rds")
saveRDS(results_cm_enet, "results_cm_enet.rds")
saveRDS(results_cm_c, "results_cm_c.rds")
saveRDS(results_cm_c_al, "results_cm_c_al.rds")
saveRDS(results_cm_gl_W, "results_cm_gl_W.rds")
saveRDS(results_cm_gl_W_al, "results_cm_gl_W_al.rds")
saveRDS(results_cm_gl_wi, "results_cm_gl_wi.rds")
saveRDS(results_cm_gl_wi_al, "results_cm_gl_wi_al.rds")
  
  
print(round(res_summary_te,4))
print(round(res_summary_tr,4))
print(round(res_summary_auc,4))
print(round(res_summary_f1,4))
print(round(res_summary_duration,4))
print(round(res_summary_n_zero,4))
print(round(res_summary_duration_total,4))

# Put all objects into a list
summaries <- list(
  TestError   = res_summary_te,
  TrainError  = res_summary_tr,
  AUC         = res_summary_auc,
  F1          = res_summary_f1,
  Duration    = res_summary_duration_total,
  N_zero      = res_summary_n_zero
)

# Extract the "Lasso" column from each and combine
lasso_metrics <- sapply(summaries, function(x) x[ , "Glasso&Corwi&Adaptive"])

# If you also want rownames (mean, sd) preserved:
lasso_metrics <- do.call(cbind, lapply(summaries, function(x) x[, "Glasso&Corwi&Adaptive", drop=FALSE]))
colnames(lasso_metrics) <- names(summaries)

print(round(lasso_metrics, 4))

    
