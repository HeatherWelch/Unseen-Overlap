# Unseen-Overlap
Code repository to support "Unseen overlap between fishing vessels and top predators in the northeast Pacific"  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10553257.svg)](https://doi.org/10.5281/zenodo.10553257)  


## Relevant papers  
Welch, Heather, et al. "Hot spots of unseen fishing vessels." Science Advances 8.44 (2022): eabq2109.  
Welch, Heather, et al. "Impacts of marine heatwaves on top predator distributions are variable but predictable." Nature Communications 14.1 (2023): 5188.  
Kroodsma, David A., et al. "Tracking the global footprint of fisheries." Science 359.6378 (2018): 904-908.  

## Scripts
1. 1_spatial_extrapolation_extraction.R: spatially extrapolate hours and dates along gap event, extract species distributions to gap and fishing vessel activity data
2. 2_calculate_overlap.R: using outputs from 1. create master dataframes containing overlap metrics for observed fishing vessel activity, and the lower and upperbound estimates of unseen fishing vessel activity
3. 3_fig1_and_stats.R: code used to run analyses for Figure 1
4. 4_fig2.R: code used to run analyses for Figure 2
5. 5_fig3.R: code used to run analyses for Figure 3
6. 6_fig4.R: code used to run analyses for Figure 4  

## Datasets
1. AIS_metrics.csv: spatially allocated observed and unseen fishing vessel activity in the northeast Pacific
2. Species distribution datasets (n=5; [species].csv): distribution information for albatrosses, sooty shearwaters, tunas, sharks, mammals and leatherbacks     

