# Iterative-individual-tree-segmentation-for-TLS
This is an R-based workflow used to segment individual trees in point cloud data, then calculate forest structure parameters including diameter at breast height (DBH), tree height (height), canopy base height (CBH), crown width (CW), and tree density. 

This workflow is designed for LAS files from ground-based LiDAR scanners (i.e., terrestrial laser scanners and mobile laser scanners), but can process files from any scanning device. However, it uses a bottom-up segmentation approach based on tree bole identification that achieves the most accurate results when a high density of points is present in the understroy. 

In addition, this workflow is designed for plot-level forest inventories, which utilize a circular plot radius with the scan collected at plot center (0, 0). A co-registered multi-scan point cloud bundle can be used so long as the center of the co-registered bundle corresponds to the center of the plot. If mobile laser scanners are being used, the plot center should be georeferenced to the (0, 0) point of the scan. Any plot radius can be used; however, as the plot radius increases, so does occlusion, which can result in improper segmentation ([Gollob et al., 2019](https://doi.org/10.3390/rs11131602)). Thus far, this workflow has only been tested using 11.4 m radius plots (1/10th acre).

## Processing steps and description:

## **Inputs:** 

  1 - A list of LAS files.
  
  2 - A CSV containing manually collected tree diameter and height data used as the reference data to build an empirical cumulative distribution function (ECDF) for the iterative resegmentation process. 

  * ***Note***: The tree list should include the dominant tree species being inventoried and correspond to the geographic region of interest. A tree list containing *Pinus ponderosa* (Ponderosa Pine) and *Pseudotsuga menziesii* (Douglas-fir) data obtained from Colorado Front Range FIA plots is included in the repository.

## **Output:** 

### Tree_metrics - A tree list in the form of a CSV file containing the following columns:
  
  * **Plot name** - The exact name of the plot defined at the beginning of the processing loop

  * **Unique_treeID** - A tree ID for each tree in a plot

  * **HT.m.** - Individual tree height in meters 

  * **DBH.m.** - Individual tree diameter at breast height in meters 

  * **CBH.m.** - Individual tree crown base height in meters 

  * **CW.m.** - Individual tree crown width in meters 

  * **Ht-DBH_ratio** - Individual tree height to DBH ratio

  * **CW-DBH_ratio** - Individual tree crown width to DBH ratio

  * **Ht-DBH_Quant** - The quantile of the ECDF that the height-DBH ratio fell within 

  * **CW-DBH_Quant** - The quantile of the ECDF that the crown width-DBH ratio fell within

  * **flag** - True/false data denoting if the tree was flagged based on the quantile values. If true, the Ht-DBH or CW-DBH quantile fell outside of the central 80th percentile of the ECDF (i.e., quantile < 0.10 or quantile > 0.90). If false, the Ht-DBH or CW-DBH quantile fell inside the central 80th percentile of the ECDF (i.e., 0.10 ≤ quantile ≤ 0.90).

  * **Segmentation_Attempt** - The segmentation attempt where the tree was added to the final tree list. If flag = true, the tree was added after the final segmentation attempt.
   
  * ***Note*** - A tree list for each segmentation attempt (Initial_attempt - Attempt_10) is also calculated and can be saved as a CSV. Each attempt includes all trees segmented up to the most current attempt (i.e., the tree list for attempt 2 contains all trees segmented in the initial attempt, attempt 1, and attempt 2)  
    
## Required packages and justification:

## Considerations
