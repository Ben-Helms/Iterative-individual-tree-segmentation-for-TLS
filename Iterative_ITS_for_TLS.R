#==========================================================================================#
# Individual Tree Segmentation for TLS data via Spanner 
# Author - Benjamin Helms
# last edited 5-5-2026
#==========================================================================================#

# Load Packages 
source("package_setup.R")

# Set Directories  
inDir <- "scanInput" # Add input filepath with las files
outDir <- "treeData" # Add output filepath for tree data 

# Get List of LAS Files  
files <- list.files(pattern = ".las", path = inDir, full.names = TRUE)

# Build data frames 
# Data frame for final tree list
Tree_metrics <- data.frame()

# Build master list for all re-segmentation attempts
All_Intermediate_Metrics_List <- list()

# Build empirical cumulative distribution function for height-DBH ratio
FIA_predictions <- read.csv("./FIA_data/FIA_Prediction_Data.csv") # found in the RProject. Can use your own data if desired

# Calculate height/DBH allometry from FIA trees 
FIA_predictions <- FIA_predictions %>%
  mutate (DBH.m. = DIA*0.0254,
          HT_DBH_Ratio = HT.m./DBH.m.)

# Build height-DBH ECDF
ecdf_HT_DBH <- ecdf(FIA_predictions$HT_DBH_Ratio)
plot(ecdf_HT_DBH)

# Calculate CW/DBH allometries from Bechtold (2004)
FIA_predictions <- FIA_predictions %>%
  mutate(
    CW.m. = case_when(
      SPCD == 122 ~ (2.3089 + 1.1388 * DIA - 0.0089 * DIA^2)*0.3048,
      SPCD == 202 ~ (5.7753 + 1.0639 * DIA - 0.0109 * DIA^2)*0.3048
    )) %>%
  mutate(CW_DBH_Ratio = CW.m./DBH.m.)

# Build CW-DBH ECDF
ecdf_CW_DBH <- ecdf(FIA_predictions$CW_DBH_Ratio)
plot(ecdf_CW_DBH)

# Find the 90th and 10th percentiles
Tenth_percent_threshold_HT_DBH <- quantile(FIA_predictions$HT_DBH_Ratio, probs = 0.10, na.rm = TRUE)
Nintyth_percent_threshold_HT_DBH <- quantile(FIA_predictions$HT_DBH_Ratio, probs = 0.90, na.rm = TRUE)
Tenth_percent_threshold_CW_DBH <- quantile(FIA_predictions$CW_DBH_Ratio, probs = 0.10, na.rm = TRUE)
Nintyth_percent_threshold_CW_DBH <- quantile(FIA_predictions$CW_DBH_Ratio, probs = 0.90, na.rm = TRUE)

#======================================================================================================#

# Start tree extraction for loop 
for (file_path in files) {

  # Extract plot name, will need to change given your plot names
  plot_name <- tools::file_path_sans_ext(basename(file_path))
  plot_name <- sub("_C_\\d+.*$", "", plot_name) # subject to change depening on plot name

  # Read LAS file
  las <- readTLS(file_path)
  
  # Clip, Thin, Classify Ground, Normalize, and Filter Noise
  las <- clip_circle(las, x = 0, y = 0, radius = 15.00) # Clipped beyond the plot radius to avoid cutting tree crowns. 
  las_thin <- decimate_points(las, random_per_voxel(res = 1, n = 2500))
  las_thin <- classify_ground(las_thin, csf(sloop_smooth = FALSE, class_threshold = 0.1,
                                            cloth_resolution = 0.01, rigidness = 1L,
                                            iterations = 500L, time_step = 0.3))
  las_norm <- normalize_height(las_thin, tin())
  las_norm <- classify_noise(las_norm, ivf(0.25, 3))
  las_norm <- filter_poi(las_norm, Classification != LASNOISE)

  # Initial inventory/segmentation - get tree locations and radii via spanner
  inventory <- get_raster_eigen_treelocs(
    las = las_norm,
    res = 0.05,
    pt_spacing = 0.0254,
    dens_threshold = .2,
    neigh_sizes = c(0.333, 0.166, 0.5),
    eigen_threshold = 0.6666,
    grid_slice_min = 1.3,
    grid_slice_max = 1.8,
    minimum_polygon_area = 0.02,
    cylinder_fit_type = "ransac",
    max_dia = 1.0,
    SDvert = 0.25,
    n_best = 25,
    n_pts = 20,
    inliers = 0.9,
    conf = 0.99,
    max_angle = 20
  )

  # Segment individual trees via spanner
  segmentedLAS <- segment_graph(
    las = las_norm,
    tree.locations = inventory,
    k = 50,
    distance.threshold = .5,
    use.metabolic.scale = FALSE,
    ptcloud_slice_min = 1.3,
    ptcloud_slice_max = 1.8,
    subsample.graph = 0.1,
    return.dense = FALSE
  )

  # Filter trees within plot radius (11.34m)
  inventory <- inventory %>% 
  mutate(dist.m. = sqrt(st_coordinates(.)[, "X"]^2 + st_coordinates(.)[, "Y"]^2)) %>%
  filter(dist.m. <= 11.34)
  
  # Filter segmented LAS to include trees within the plot radius
  Individual_Tree_Seg <- filter_poi(segmentedLAS, treeID %in% inventory$TreeID)
  
  # Calculate plot-level initial inventory metrics via spanner 
  tree_metrics <- process_tree_data(treeData = inventory, segmentedLAS = Individual_Tree_Seg, return_sf = FALSE)

  # Standardize all metric names for consistency and remove saplings before flagging
  tree_metrics <- tree_metrics %>%
    mutate(
      CW.m. = 2 * sqrt(crown_area / pi),
      HT.m. = height, 
      DBH.m. = diameter, 
      CBH.m. = crown_base_height, 
      `Ht-DBH_ratio` = HT.m./DBH.m., 
      `CW-DBH_ratio` = CW.m./DBH.m.  
    ) %>%
    # Remove trees that don't meet inventory criteria
    filter(HT.m. >= 1.37) %>%
    filter(DBH.m. >= 0.127) %>%
    dplyr::select(TreeID, HT.m., DBH.m., CBH.m., CW.m., `Ht-DBH_ratio`, `CW-DBH_ratio`)

  # Talley the initial tree segmentation to compare improvements after resegmentation
  initial_segmented_count <- nrow(tree_metrics)

  # Create initial_tree_metrics data frame, calculate flags, and assign Segmentation_Attempt = 0
  initial_tree_metrics <- tree_metrics %>%
    mutate(
      `Ht-DBH_Quant` = ecdf_HT_DBH(`Ht-DBH_ratio`),
      `CW-DBH_Quant` = ecdf_CW_DBH(`CW-DBH_ratio`),
      flag = `Ht-DBH_ratio` < Tenth_percent_threshold_HT_DBH |
        `Ht-DBH_ratio` > Nintyth_percent_threshold_HT_DBH |
        `CW-DBH_ratio` < Tenth_percent_threshold_CW_DBH |
        `CW-DBH_ratio` > Nintyth_percent_threshold_CW_DBH,
      Segmentation_Attempt = 0
    ) %>%
    rename("treeID" = "TreeID")
  
  # Record Attempt 0 in the intermediate metrics list
  initial_list_entry <- initial_tree_metrics %>%
    mutate(
      "Plot Name" = plot_name,
      "Attempt" = 0,
      "Status" = "Initial Segmentation"
    ) %>%
    dplyr::select("Plot Name", "Attempt", "Status", "treeID", "HT.m.", "DBH.m.", "CBH.m.", "CW.m.",
           "Ht-DBH_ratio", "CW-DBH_ratio", "flag", "Ht-DBH_Quant", "CW-DBH_Quant", "Segmentation_Attempt")

  All_Intermediate_Metrics_List <- append(All_Intermediate_Metrics_List, list(initial_list_entry))

  # Split the initial trees for re-segmentation of flagged trees
  good_trees_initial <- initial_tree_metrics %>%
    filter(flag == FALSE)

  current_tree_metrics_to_process <- initial_tree_metrics %>%
    filter(flag == TRUE) # Only flagged trees proceed into the re-segmentation while loop

  # Add good trees to the final tree list
  final_tree_metrics <- good_trees_initial

  # Use the initial segmented LAS for re-segmentation.
  current_segmentedLAS <- Individual_Tree_Seg
  
  # Start while loop for re-segmentation 
  
  max_attempts <- 10
  attempt <- 1
  consecutive_failure_count <- 0

  while(TRUE) {

    # Separate Good and Flagged Trees from the current tree list
    good_trees_current_round <- current_tree_metrics_to_process %>%
      filter(flag == FALSE)

    flagged_trees_next_round <- current_tree_metrics_to_process %>%
      filter(flag == TRUE)
    
    # Add Good trees to the final results
    final_tree_metrics <- dplyr::bind_rows(
      final_tree_metrics, 
      good_trees_current_round
    )
    
    # Log/count flagged and good trees 
    num_flagged <- nrow(flagged_trees_next_round)
    num_good <- nrow(good_trees_current_round)

    # Update consecutive failure count
    if (num_good == 0) {
      consecutive_failure_count <- consecutive_failure_count + 1
    } else {
      consecutive_failure_count <- 0
    }

    # Initialize log status

    # Define while loop stopping conditions
    # Stop condition A: No trees flagged
    if (num_flagged == 0) {
      break
    }

    # Stop condition B: Tree list contains only flagged trees twice in a row. 
    if (consecutive_failure_count == 2) {
      cat(paste("Redundancy detected. >2 consecutive failures. Stopping.\n"))
      final_tree_metrics <- bind_rows(final_tree_metrics, flagged_trees_next_round)
      break
    }

    # Stop condition C: 10 attempts reached
    if (attempt >= max_attempts) {
      cat(paste("Maximum attempts (", max_attempts, ") reached. Stopping.\n"))
      final_tree_metrics <- bind_rows(final_tree_metrics, flagged_trees_next_round)
      break
    }

    # Filter las to only include flagged trees 
    las_to_resegment <- filter_poi(current_segmentedLAS, treeID %in% flagged_trees_next_round$treeID)

    # Inventory/Segmentation - tryCatch is used to keep the for loop from stopping if no trees are found or if the point cloud is too thin to resegment
    inventory_flagged <-  tryCatch(
      {
        get_raster_eigen_treelocs(
        las = las_to_resegment,
        res = 0.05,
        pt_spacing = 0.0254,
        dens_threshold = .02,
        neigh_sizes = c(0.333, 0.166, 0.5),
        eigen_threshold = 0.6666,
        grid_slice_min = 0.6666,
        grid_slice_max = 2,
        minimum_polygon_area = 0.025,
        cylinder_fit_type = "ransac",
        max_dia = 0.5,
        SDvert = 0.25,
        n_best = 25,
        n_pts = 25,
        inliers = 0.9,
        conf = 0.99,
        max_angle = 20)
      },
      error = function(e) {
        message("Segmentation fail for", plot_name)
        return(NULL)
      }
    )
    
    if (is.null(inventory_flagged) || nrow(inventory_flagged) == 0) {
      message(paste("Resegmentation complete:", plot_name, "— No additonal trees found."))
      final_tree_metrics <- bind_rows(final_tree_metrics, flagged_trees_next_round) 
      break 
    }

    # Segment Individual Trees
    segmentedLAS_flagged_new <- segment_graph(
      las = las_to_resegment,
      tree.locations = inventory_flagged,
      k = 50,
      distance.threshold = .5,
      use.metabolic.scale = FALSE,
      ptcloud_slice_min = 1.3,
      ptcloud_slice_max = 1.8,
      subsample.graph = 0.1,
      return.dense = FALSE
    )
    
    # Ensure treeID's match between inventory and segmentation 
    treeID_in_seg <- unique(segmentedLAS_flagged_new$treeID)
    inventory_flagged <- inventory_flagged %>%
      filter(TreeID %in% treeID_in_seg)

    # Process plot-level inventory metrics via spanner and clean/rename columns
    tree_metrics_flagged_new <- process_tree_data(treeData = inventory_flagged, segmentedLAS = segmentedLAS_flagged_new, return_sf = FALSE) %>%
      mutate(
        CW.m. = 2 * sqrt(crown_area / pi),
        HT.m. = height,
        DBH.m. = diameter,
        CBH.m. = crown_base_height,
        `Ht-DBH_ratio` = HT.m./DBH.m.,
        `CW-DBH_ratio` = CW.m./DBH.m.
      ) %>%
      filter(HT.m. >= 1.37) %>%
      filter(DBH.m. >= .127) %>%
      dplyr::select(TreeID, HT.m., DBH.m., CBH.m., CW.m., `Ht-DBH_ratio`, `CW-DBH_ratio`)

    # Re-flag the resegmented metrics for the next round's segmentation attempt
    current_tree_metrics_to_process <- tree_metrics_flagged_new %>%
      mutate(
        `Ht-DBH_Quant` = ecdf_HT_DBH(`Ht-DBH_ratio`),
        `CW-DBH_Quant` = ecdf_CW_DBH(`CW-DBH_ratio`),
        flag = `Ht-DBH_ratio` < Tenth_percent_threshold_HT_DBH |
          `Ht-DBH_ratio` > Nintyth_percent_threshold_HT_DBH |
          `CW-DBH_ratio` < Tenth_percent_threshold_CW_DBH |
          `CW-DBH_ratio` > Nintyth_percent_threshold_CW_DBH,
        Segmentation_Attempt = attempt
      ) %>%
      rename("treeID" = "TreeID")

    # Combine all trees evaluated in this attempt (both good & flagged) with the finalized trees
    full_tally_for_log <- bind_rows(final_tree_metrics, current_tree_metrics_to_process) %>%
      mutate(
        "Plot Name" = plot_name,
        "Attempt" = attempt,
        "Status" = paste0("Attempt ", Segmentation_Attempt, " Result")
      )

    # Create the log entry based on the full tally
    iteration_log_entry_df <- full_tally_for_log %>%
      dplyr::select(
        "Plot Name", "Attempt", "Status", "treeID", "HT.m.", "DBH.m.", "CBH.m.", "CW.m.",
        "Ht-DBH_ratio", "CW-DBH_ratio", "flag", "Ht-DBH_Quant", "CW-DBH_Quant", "Segmentation_Attempt"
      )

    # ----------------------------------------------------------------------

    # Append the iteration data frame to the master list
    All_Intermediate_Metrics_List <- append(All_Intermediate_Metrics_List, list(iteration_log_entry_df))

    # Update the segmented LAS data for the next round of segmentation 
    current_segmentedLAS <- segmentedLAS_flagged_new

    # Add to the attempt counter for the next loop iteration
    attempt <- attempt + 1

  } # End of while loop
  
  # Summarize and store individual tree metrics
  if (nrow(final_tree_metrics) > 0) {
    plot_summary <- final_tree_metrics %>%
      mutate("Plot Name" = plot_name,
             X = st_coordinates(.)[, "X"],
             Y = st_coordinates(.)[, "Y"]) %>%
      st_drop_geometry()
    
    Tree_metrics <- rbind(Tree_metrics, plot_summary) 
    
  } else {
    message(paste("No trees found for plot:", plot_name, "after filtering. Skipping metrics calculation."))
  }

  beepr::beep(sound = 1)

} # End of for the loop

# Create the final and intermediate metrics data frame after the loop is complete
if(length(All_Intermediate_Metrics_List) > 0) {
  Master_Intermediate_Metrics <- bind_rows(All_Intermediate_Metrics_List)
} else {
  Master_Intermediate_Metrics <- data.frame()
  warning("All_Intermediate_Metrics_List is empty. No intermediate metrics were recorded.")
}

Master_Intermediate_Metrics <- Master_Intermediate_Metrics %>%
  mutate(X = st_coordinates(.)[, "X"],
         Y = st_coordinates(.)[, "Y"]) %>%
  st_drop_geometry()

# Prepare to save the data

prep_data <- function(df) {
  df %>%
    group_by(`Plot Name`) %>%
    mutate(Unique_treeID = row_number()) %>%
    relocate(Unique_treeID, .after = `Plot Name`) %>%
    ungroup() %>%
    dplyr::select("Plot Name", "Unique_treeID", "flag", "Attempt", "Segmentation_Attempt", 
                  "HT.m.", "DBH.m.", "CBH.m.", "CW.m.", "Ht-DBH_ratio", 
                  "CW-DBH_ratio", "Ht-DBH_Quant", "CW-DBH_Quant", "X", "Y")
}

# --- USER SETTINGS ---
# print the max number of segmentation attempts
print(max(Master_Intermediate_Metrics$Segmentation_Attempt))

# define which attempts to download
attempts_to_save <- c(0,1) # 0 is Initial. Add any numbers from 1-10.
save_final <- TRUE             # Set to TRUE to save the final Tree_metrics
# ---------------------

# Save Final Metrics
if(save_final) {
  final_out <- Tree_metrics %>%
  group_by(`Plot Name`) %>% 
  mutate(Unique_treeID = row_number()) %>% 
  relocate(`Plot Name`, Unique_treeID) %>%
  dplyr::select(-c("treeID")) %>%
  ungroup()
  
  write.table(final_out, file.path(outDir, "Forest_Structure_Estimation_Output.csv"), row.names = FALSE, append = TRUE, sep = ",")
}

#Save Selected Attempts
for (i in attempts_to_save) {
  attempt_df <- Master_Intermediate_Metrics %>% filter(Attempt == i)
  
  if (nrow(attempt_df) > 0) {
    # Process and Save
    processed <- prep_data(attempt_df)
    file_label <- if(i == 0) "Initial" else paste0("Attempt_", i)
    
    write.table(processed, 
              file = file.path(outDir, paste0(file_label, "_Forest_Structure_Estimation.csv")), 
              row.names = FALSE, append = TRUE, sep = ",")
  }
}
