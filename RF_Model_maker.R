# Load the necessary libraries
library(terra)
library(sf)
library(data.table)
library(randomForest)
library(viridis)
library(ggplot2)
library(vip)
library(pROC)
library(yardstick)
library(pdp)
library(caret)
library(spThin)
library(rpart.plot)
library(reprtree)
library(plotrix)
library(tree)
library(partykit)
library(dplyr)
library(grid)
library(igraph)
library(ggraph)

# Define a function to impute NA values with the mean for a raster layer
impute_na_with_mean <- function(r) {
  vals <- values(r)
  na_mean <- mean(vals, na.rm = TRUE)
  vals[is.na(vals)] <- na_mean
  values(r) <- vals
  return(r)
}

# Function to manually thin the points based on a minimum distance
thin_points_sf <- function(points, min_dist) {
  points_sf <- st_as_sf(points, coords = c("longitude", "latitude"), crs = 4326)
  thinned_sf <- points_sf[1, ]  # Start with the first point
  
  for (i in 2:nrow(points_sf)) {
    point <- points_sf[i, ]
    distances <- st_distance(thinned_sf, point)
    if (all(as.numeric(distances) >= min_dist * 1000)) {  # Convert km to meters
      thinned_sf <- rbind(thinned_sf, point)
    }
  }
  
  return(thinned_sf)
}

# Load presence data
presence_data <- fread("F:\\Modeling\\Stinknet_Updated_RF\\presance\\ONPI Pres.csv")

# Adjust column names based on actual data structure
presence_data_clean <- presence_data[, .(latitude = get("latitude"), longitude = get("longitude"))]
presence_data_clean <- na.omit(presence_data_clean)

# Check if presence_data_clean is empty
if (nrow(presence_data_clean) == 0) {
  stop("Presence data is empty after cleaning.")
}

# Define minimum distance in kilometers
min_dist_km <- 15  # Example: 10 km

# Apply the custom thinning function
presence_data_thinned_sf <- thin_points_sf(presence_data_clean, min_dist_km)

# Convert the thinned sf object back to a data table
presence_data_thinned <- as.data.table(st_coordinates(presence_data_thinned_sf))
setnames(presence_data_thinned, c("X", "Y"), c("longitude", "latitude"))

# Visualize presence data points after thinning
ggplot(presence_data_thinned, aes(x = longitude, y = latitude)) +
  geom_point(color = "red") +
  ggtitle("Presence Data Points After Thinning") +
  xlab("Longitude") +
  ylab("Latitude")

# Convert to sf object
presence_data_sf <- st_as_sf(presence_data_thinned, coords = c("longitude", "latitude"), crs = 4326)


# Load the template raster (Min_Temp)
template_raster <- rast("F:/Modeling/Stinknet_Updated_RF/AZ_Rasters/Min_Temp.tif")
cat("Template raster loaded with resolution:", res(template_raster), "\n")

# Set target CRS
target_crs <- "EPSG:4326"  # WGS 1984

# Load other rasters
raster_directory <- "F:/Modeling/Stinknet_Updated_RF/AZ_Rasters"
raster_files <- list.files(path = raster_directory, pattern = "\\.tif$", full.names = TRUE)
aligned_rasters <- list()

# Process each raster file
for (file in raster_files) {
  cat("Processing", file, "\n")
  r <- rast(file)
  
  # Check if the raster is read correctly
  if (is.null(r)) {
    cat("Failed to load", file, "\n")
    next
  }
  
  # Check for non-NA values in the original raster
  original_values <- values(r)
  if (all(is.na(original_values))) {
    cat("All values are NA in the original raster", file, "\n")
    next
  }
  
  # Project to WGS 1984 if necessary
  if (crs(r) != target_crs) {
    cat("Projecting", file, "to WGS 1984\n")
    r <- project(r, target_crs)
  }
  
  # Resample to match the template raster
  cat("Resampling", file, "to match template raster\n")
  r <- resample(r, template_raster, method = "bilinear")
  
  # Check new resolution
  new_res <- res(r)
  cat("New resolution of", file, ":", new_res, "\n")
  
  # Impute NA values with mean for continuous rasters
  r <- impute_na_with_mean(r)
  
  # Check for non-NA values in the resampled raster
  resampled_values <- values(r)
  if (all(is.na(resampled_values))) {
    cat("All values are NA in the resampled raster", file, "\n")
    next
  }
  
  raster_name <- tools::file_path_sans_ext(basename(file))
  aligned_rasters[[raster_name]] <- r
}

# Combine aligned rasters into a stack
raster_stack <- rast(aligned_rasters)

# Check if the stack was created successfully
summary(raster_stack)

# Extract raster values for presence data
presence_values <- terra::extract(raster_stack, presence_data_thinned_sf)
presence_values <- presence_values[, -1]
presence_df <- as.data.table(presence_data_thinned_sf)
presence_df <- cbind(presence_df, presence_values)
presence_df$response <- 1

# Generate background data
set.seed(17172)
background_points <- spatSample(raster_stack, size = nrow(presence_data_thinned), xy = TRUE, as.df = TRUE)
background_sf <- st_as_sf(background_points, coords = c("x", "y"), crs = 4326)
background_values <- terra::extract(raster_stack, background_sf)
background_values <- background_values[, -1]
background_df <- as.data.table(background_sf)
background_df <- cbind(background_df, background_values)
background_df$response <- 0

# Print summary of background data
cat("Summary of background data:\n")
print(summary(background_df))
cat("\n")

# Ensure presence_df and background_df have the same columns
common_cols <- intersect(names(presence_df), names(background_df))
presence_df <- presence_df[, ..common_cols]
background_df <- background_df[, ..common_cols]

# Combine presence and background data
combined_data <- rbind(presence_df, background_df, fill = TRUE)

# Remove geometry column if it exists
combined_data <- combined_data[, !"geometry", with = FALSE]

# Remove rows with any remaining NA values
combined_data_clean <- na.omit(combined_data)
combined_data_clean$response <- as.factor(combined_data_clean$response)

# Rename levels of the response variable to valid R variable names
levels(combined_data_clean$response) <- make.names(levels(combined_data_clean$response))

# Print summary of combined data
cat("Summary of combined data:\n")
print(summary(combined_data_clean))
cat("\n")

# Define response and predictors
response_var <- "response"
predictors <- setdiff(names(combined_data_clean), response_var)

# Set up cross-validation control
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Define the grid of hyperparameters to search
tune_grid <- expand.grid(mtry = seq(2, min(6, length(predictors)), by = 1))

# Initialize results storage
results <- data.frame(ntree = integer(), mtry = integer(), ROC = numeric())

# Loop over ntree values and train the model for each value
ntree_values <- seq(200, 1200, by = 100)
for (ntree in ntree_values) {
  rf_model_cv <- train(
    response ~ .,
    data = combined_data_clean,
    method = "rf",
    trControl = control,
    tuneGrid = tune_grid,
    metric = "ROC",
    ntree = ntree,
    importance = TRUE
  )
  
  # Extract cross-validated AUC
  cv_results <- rf_model_cv$results
  best_mtry <- cv_results[which.max(cv_results$ROC), "mtry"]
  best_roc <- max(cv_results$ROC)
  
  # Store results
  results <- rbind(results, data.frame(ntree = ntree, mtry = best_mtry, ROC = best_roc))
}

# Print the results
print(results)

# Select the best model based on ROC
best_model_index <- which.max(results$ROC)
best_ntree <- results$ntree[best_model_index]
best_mtry <- results$mtry[best_model_index]

# Train the final model with the best hyperparameters
final_rf_model <- randomForest(
  response ~ .,
  data = combined_data_clean,
  ntree = best_ntree,
  mtry = best_mtry,
  nodesize = 4,  # Adjust nodesize to control tree depth
  maxnodes = 20  # Adjust maxnodes to control tree size
)

# Evaluate the model using OOB error
print(final_rf_model)

# Check variable importance to identify redundant features
varImpPlot(final_rf_model)
# Create a unique filename for the prediction output
output_filename <- paste0("F:\\Modeling\\Stinknet_Updated_RF\\RF_Models\\Model_5", Sys.Date(), ".tif")

# Predict raster using the trained model
predictions <- tryCatch({
  terra::predict(raster_stack, final_rf_model, type = "prob", index = 2, filename=output_filename, na.rm=TRUE, progress="text", overwrite= FALSE)
}, error = function(e) {
  cat("Error in prediction step:", e$message, "\n")
  NULL
})

if (!is.null(predictions)) {
  predicted_raster <- rast(output_filename)
  plot(predicted_raster, main="Predicted Probability of Presence", col=viridis::viridis(100))
} else {
  cat("Prediction step failed. Please check the compatibility of the model and the raster stack.\n")
}

# Define a custom prediction function for permutation importance
predict_function <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, 2]
}

# Calculate permutation importance with AUC as the metric
perm_importance <- vip::vi(
  object = rf_model_cv$finalModel,
  method = "permute",
  target = "response",
  train = combined_data_clean,
  metric = yardstick::roc_auc_vec,
  pred_wrapper = predict_function,
  nsim = 50,  # Number of permutations
  sample_frac = 1,  # Use the whole dataset for each permutation
  smaller_is_better = FALSE
)

# Plot permutation importance using ggplot2
vip::vip(perm_importance, num_features = length(predictors), geom = "col") +
  theme_minimal() +
  ggtitle("Permutation Variable Importance Plot") +
  xlab("Variables") +
  ylab("Importance")

# Create PDPs for all predictor variables
for (var in predictors) {
  pd <- partial(rf_model_cv$finalModel, pred.var = var, train = combined_data_clean, prob = TRUE)
  p <- autoplot(pd) +
    ggtitle(paste("Partial Dependence Plot for", var)) +
    xlab(var) +
    ylab("Predicted Probability of Presence") +
    theme_minimal()
  print(p)
  cat("Press Enter to continue to the next plot...")
  readline()
}


# Plot ROC curve for each fold
for (fold in folds) {
  fold_predictions <- predictions[predictions$Resample == fold, ]
  roc_obj <- roc(fold_predictions$obs, fold_predictions$yes)
  plot(roc_obj, col = "blue", add = TRUE)
}

# Extract the predictions and plot ROC curves for each fold
predictions <- rf_model_cv$pred
folds <- unique(predictions$Resample)

# Initialize plot
plot(0, 0, type = "n", xlab = "1 - Specificity", ylab = "Sensitivity", xlim = c(0, 1), ylim = c(0, 1), main = "ROC Curves for Each Fold")

# Add a legend with the fold names
legend("bottomright", legend = paste("Fold", folds), col = "blue", lty = 1)

# Extract and plot cross-validated AUC
cv_results <- rf_model_cv$results
cv_auc <- cv_results[which.max(cv_results$ROC), "ROC"]

# Generate the overall ROC curve
train_predictions <- predict(rf_model_cv, combined_data_clean, type = "prob")[, 2]
roc_obj <- roc(combined_data_clean$response, train_predictions)
plot(roc_obj, main = "Cross-Validated ROC Curve for Random Forest Model", col = "blue")
mtext(paste("AUC:", round(cv_auc, 3)), side = 1, line = 2.5, at = 0.5, col = "blue")









# Load necessary libraries
library(randomForest)
library(ggplot2)
library(igraph)
library(ggraph)
library(tibble)
library(tidyr)
library(dplyr)
library(grid)

# Ensure igraph and tidyr are loaded
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
library(igraph)

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(tidyr)

# Function to plot a tree
plot_tree <- function(final_model, tree_num) {
  tree <- randomForest::getTree(final_model, k = tree_num, labelVar = TRUE) %>%
    tibble::rownames_to_column("rowname") %>%
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  
  graph_frame <- tree %>%
    select(rowname, `left daughter`, `right daughter`) %>%
    tidyr::pivot_longer(cols = c(`left daughter`, `right daughter`),
                        names_to = "direction", values_to = "to") %>%
    filter(to != 0) %>%
    select(from = rowname, to)
  
  graph <- graph_from_data_frame(graph_frame)
  
  V(graph)$name <- as.character(V(graph))
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`[match(V(graph)$name, tree$rowname)]))
  V(graph)$leaf_label <- as.character(tree$prediction[match(V(graph)$name, tree$rowname)])
  V(graph)$split <- as.character(round(tree$`split point`[match(V(graph)$name, tree$rowname)], digits = 2))
  
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18)) +
    ggtitle(paste("Decision Tree", tree_num))
  
  return(plot)
}

# Create and save plots for the first 200 trees
output_directory <- "F:\\Modeling\\CECI Update\\TreePlots"
dir.create(output_directory, showWarnings = FALSE)

for (i in 500:1000) {
  plot <- plot_tree(rf_model_cv$finalModel, i)
  jpeg_filename <- file.path(output_directory, paste0("tree_", i, ".jpg"))
  ggsave(filename = jpeg_filename, plot = plot, width = 30, height = 30, dpi = 300, units = "in")
}
