# code

```text
code/
├── 00_BBS Preparation/
│   ├── 00_BBS Initial Cleaning_updated.R
│   ├── 01_Creating spatial BBS routes.R
│   ├── 02_BBS Final Cleaning.R
│   ├── 03_Determining Sample Size.R
│   └── 04_BBS Conversion to site x species matrix.R
├── 01_Generating spatial polygons/
│   ├── 00_Generating_polygons_cookies.R
│   ├── 1_1 Estimating EH - pairs-4routes.R
│   ├── 1_2 Estimating alpha beta gamma - pairs.R
│   └── 1_3 Estimating alpha beta gamma - 4routes.R
├── 02_Landcover Harmonization and Extraction/
│   ├── 00_Projecting_Conusb.R
│   ├── 00_Projecting_Conusb.sh
│   ├── 01_Resample_NLCD_250m.R
│   ├── 02_Reclassify_NLCD_to_match_CONUSb.R
│   ├── 03_Crossvalidation of NLCD and Conusb.R
│   ├── 04_Evaluation of Crossvalidation.R
│   ├── 04_Visualization of Crossvalidation.R
│   ├── 05_Reclassify_based_on_accuracy.R
│   ├── 06_RE-Crossvalidation of NLCD and Conusb.R
│   ├── 07_RE-Evaluation of Crossvalidation.R
│   ├── 08_RE-Visualization of Crossvalidation.R
│   ├── 09_Relating_NLCD_and_CONUS_to_BBS_routes.R
│   ├── 09_Relating_NLCD_and_CONUS_to_BBS_routes.sh
│   ├── 10_Aggregating land cover_area_for BBS routes.R
│   ├── 10_Aggregating land cover_proportions_ for BBS routes.R
│   └── 11_Merging NLCD and CONUS output from cluster.R
├── 03_Calculating beta diversity and EH/
│   ├── 1_1 Estimating EH - pairs-4routes.R
│   ├── 1_2 Estimating alpha beta gamma - pairs.R
│   └── 1_3 Estimating alpha beta gamma - 4routes.R
├── 04_Hierarchical Modelling/
│   ├── 00_SamplingFunction.R
│   ├── Four Routes/
│   │   ├── 01_Main modelling_4routes.R
│   │   ├── 02_Model Diagnostics_4routes.R
│   │   └── 03_Main_modelling_4routes_without_outliers.R
│   └── Route Pairs/
│       ├── 01_Main modelling_pairs.R
│       ├── 02_Model Diagnostics.R
│       ├── 03_Main_modelling_pairs_without_outliers.R
│       ├── 04_Extracting Trend means and CIs_main model.R
│       └── 05_Converting Trend mean and CI to percent change.R
├── 05_Habitat Association Model/
│   ├── 00_Defining species and abundances.R
│   ├── 01_Wrangling abundances and habitat associations.R
│   ├── 02_Preparing data for modeling.R
│   ├── 03_Habitat association model.R
│   └── 04_Habitat association model_wo house sparrow.R
└── 06_Figures/
    ├── Figure 2a_S16a.R
    ├── Figure 2b_S16c.R
    ├── Figure 2c_S16e.R
    ├── Figure 2d.R
    ├── Figure 3 Post pred hab associations.R
    ├── Figure 4_Metacommunity and Hab Assoc Plot.R
    ├── Figure 5 - Land cover change by metacommunity.R
    ├── Figure 6 - Distance shifted.R
    ├── Figure S1 - Map of routes and cookies.R
    ├── Figure S10_Trend posterior means and CIs.R
    ├── Figure S11 - Pareto K_pairs.R
    ├── Figure S11-14 - Slopes v Polygon extent.R
    ├── Figure S12b_pairs_wo_outliers.R
    ├── Figure S12d_pairs_wo_outliers.R
    ├── Figure S12f_pairs_wo_outliers.R
    ├── Figure S13_dist_pairs_wo_outliers.R
    ├── Figure S14a_four_routes.R
    ├── Figure S14b_four_routes_wo_outliers.R
    ├── Figure S14c_four_routes.R
    ├── Figure S14d_four_routes_wo_outliers.R
    ├── Figure S14e_four_routes.R
    ├── Figure S14f_four_routes_wo_outliers.R
    ├── Figure S15 - Pareto K_4routes.R
    ├── Figure S16_dist_4routes.R
    ├── Figure S16_dist_4routes_wo_outliers.R
    ├── Figure S17_Rel Ab Changes for 6 hab assoc species.R
    ├── Figure S18 Post pred hab associations house sparrow.R
    └── Figure S19 - Rel Ab Change vs Hab Assoc by species.R
```
