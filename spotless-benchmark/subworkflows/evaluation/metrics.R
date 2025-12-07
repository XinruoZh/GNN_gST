#!/usr/bin/env Rscript
library(precrec)
library(magrittr)
library(philentropy)
library(stringr)

getConfusionMatrix <- function(known_props, test_props, threshold = 0.05){
  test_props[test_props < threshold] <- 0
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  missing_rows <- which(rowSums(is.na(known_props)) > 0)
  
  if (length(missing_rows) > 0){
    test_props <- test_props[-missing_rows,]
    known_props <- known_props[-missing_rows,]
  }
  for (i in 1:nrow(known_props)){
    for (j in 1:ncol(known_props)){
      if (known_props[i, j] > 0 & test_props[i, j] > 0){
        tp <- tp + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] == 0){
        tn <- tn + 1
      } else if (known_props[i, j] > 0 & test_props[i, j] == 0){
        fn <- fn + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] > 0){
        fp <- fp + 1
      }
    }
  }
  return(list(tp=tp, tn=tn, fn=fn, fp=fp))
}

par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

# --- 1. Load Reference Data ---
ground_truth_data <- readRDS(par$sp_input)

if (class(ground_truth_data) == "list"){ 
  celltype_cols <- !grepl("^name$|^region$|^spot_no$",
                          colnames(ground_truth_data$relative_spot_composition))
}

ncells <- sum(celltype_cols)
known_props <- ground_truth_data$relative_spot_composition[,celltype_cols]

# Clean known_props names
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
known_props <- known_props[,sort(colnames(known_props), method="shell")]

# --- 2. Load Deconvolution Results ---
# Check if file exists
if (!file.exists(par$props_file)) {
    stop(paste("Proportions file not found:", par$props_file))
}

deconv_matrix <- read.table(par$props_file, sep="\t", header=TRUE, check.names=FALSE)

# --- FIX: Handle Row Names (Spot IDs) ---
# If the first column is character/factor, assume it is the Spot ID index
if (!is.numeric(deconv_matrix[,1])) {
    # Check if these look like the spot IDs in known_props
    # (Optional verify, but generally safe to assume first col is index if char)
    rownames(deconv_matrix) <- deconv_matrix[,1]
    deconv_matrix <- deconv_matrix[,-1]
}

# --- FIX: Clean and Sort Columns ---
# Ensure column names are cleaned exactly like known_props to allow matching
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")

# Sort columns alphabetically
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

# Ensure numeric matrix
deconv_matrix <- as.matrix(deconv_matrix)

# --- 3. Align Columns ---
# 1. Identify missing columns in Deconv (add as 0)
missing_in_deconv <- setdiff(colnames(known_props), colnames(deconv_matrix))
if (length(missing_in_deconv) > 0) {
    zero_mat <- matrix(0, nrow=nrow(deconv_matrix), ncol=length(missing_in_deconv))
    colnames(zero_mat) <- missing_in_deconv
    deconv_matrix <- cbind(deconv_matrix, zero_mat)
}

# 2. Identify extra columns in Deconv (add to Known as 0 to avoid dimension mismatch)
# (Ideally shouldn't happen if using same reference, but good for safety)
missing_in_known <- setdiff(colnames(deconv_matrix), colnames(known_props))
if (length(missing_in_known) > 0) {
    zero_mat <- matrix(0, nrow=nrow(known_props), ncol=length(missing_in_known))
    colnames(zero_mat) <- missing_in_known
    known_props <- cbind(known_props, zero_mat)
}

# 3. Final Sort to ensure identical order
common_cols <- sort(colnames(known_props), method="shell")
known_props <- known_props[, common_cols]
deconv_matrix <- deconv_matrix[, common_cols]

# Optional: Remapping logic
if (!is.null(par$remap)){
  conversion <- read.table(par$remap, sep="\t") %>% setNames(c("old_annot", "new_annot"))
  # Check if all conversion columns exist
  present_cols <- colnames(deconv_matrix)
  
  # Only run mapping if columns are present
  deconv_matrix <- sapply(unique(conversion$new_annot), function(new_celltype){
                   old_cols <- conversion$old_annot[conversion$new_annot == new_celltype]
                   old_cols <- old_cols[old_cols %in% present_cols]
                   if(length(old_cols) == 0) return(rep(0, nrow(deconv_matrix)))
                   if(length(old_cols) == 1) return(deconv_matrix[,old_cols])
                   rowSums(deconv_matrix[,old_cols, drop=FALSE])
                   })
                   
  known_props <- sapply(unique(conversion$new_annot), function(new_celltype){
                   old_cols <- conversion$old_annot[conversion$new_annot == new_celltype]
                   old_cols <- old_cols[old_cols %in% present_cols]
                   if(length(old_cols) == 0) return(rep(0, nrow(known_props)))
                   if(length(old_cols) == 1) return(known_props[,old_cols])
                   rowSums(known_props[,old_cols, drop=FALSE])
                   })
}

ncells <- ncol(deconv_matrix)

# --- 4. Calculate Metrics ---

# Correlation and RMSE
# Ensure both are numeric matrices
known_props <- as.matrix(known_props)
deconv_matrix <- as.matrix(deconv_matrix)

corr_spots <- mean(diag(cor(t(known_props), t(deconv_matrix))), na.rm=TRUE)
RMSE <- mean(sqrt(rowSums((known_props-deconv_matrix)**2)/ncells), na.rm=TRUE)

# Classification metrics
conf_mat <- getConfusionMatrix(known_props, deconv_matrix)
accuracy <- round((conf_mat$tp + conf_mat$tn) / (conf_mat$tp + conf_mat$tn + conf_mat$fp + conf_mat$fn), 3)
sensitivity <- round(conf_mat$tp / (conf_mat$tp + conf_mat$fn), 3)
specificity <- round(conf_mat$tn / (conf_mat$tn + conf_mat$fp), 3)
balanced_accuracy <- round((sensitivity + specificity) / 2, 3)
precision <- round(conf_mat$tp / (conf_mat$tp + conf_mat$fp), 3)
F1 <- round((2* precision * sensitivity) / (precision + sensitivity), 3)
F2 <- round((5 * precision * sensitivity) / (4*precision + sensitivity), 3)

# Precrec package
known_props_binary <- ifelse(known_props > 0, "present", "absent") %>%
                      reshape2::melt() %>% dplyr::select(value)

# Area under precision-recall curve
# Wrap in tryCatch as evalmod can fail on edge cases (e.g. all zeros)
tryCatch({
    eval_prc <- evalmod(scores = c(as.matrix(deconv_matrix)), labels=known_props_binary)
    prc <- subset(auc(eval_prc), curvetypes == "PRC")$aucs
    roc <- subset(auc(eval_prc), curvetypes == "ROC")$aucs
}, error = function(e) {
    prc <<- NA
    roc <<- NA
})

# Jensen-shannon divergence
jsd <- suppressMessages(
        sapply(1:nrow(known_props), function(i) {
          JSD(as.matrix(rbind(known_props[i,], deconv_matrix[i,])))
  })) %>% mean(na.rm=TRUE)

metrics <- data.frame("corr"=corr_spots, "RMSE"=RMSE,
                      "accuracy"=accuracy, "balanced_accuracy"=balanced_accuracy,
                      "sensitivity"=sensitivity, "specificity"=specificity,
                      "precision"=precision, "F1"=F1, "F2"=F2, "prc"=prc, "roc"=roc,
                      "jsd"=jsd)

write.table(metrics, file=par$output, row.names=FALSE)
