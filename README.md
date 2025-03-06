# 3D-RSD

This repository contains the implementation for **Fast Non-Line-of-Sight (NLOS) Imaging with Non-Planar Relay Surfaces**. For more details, please refer to the paper:  [https://ieeexplore.ieee.org/abstract/document/10233262](https://ieeexplore.ieee.org/abstract/document/10233262).

---

## Quick start
Clone the repo, select the `scene_index` from {0,2,6,10} in `RSD_main.m` and run. 

## Data description
The data in the 'data' folder contains the following fields:
- **ts**: the temporal resolution of the captured histogram.
- **total_rect_data**: the sum of the histograms captured by all the 16*16 SPAD array pixels.
- **camera_pos**: the scanned laser positions on the relay surface (noted that the name is flipped).
- **laser_pos**: the average SPAD position on the relay surface.

## Citation  
If you find this work useful, please cite:  
*C. Gu, T. Sultan, K. Masumnia-Bisheh, L. Waller, and A. Velten, "Fast Non-Line-of-Sight Imaging with Non-Planar Relay Surfaces," 2023 IEEE International Conference on Computational Photography (ICCP), Madison, WI, USA, 2023, pp. 1-12.*  
[DOI: 10.1109/ICCP56744.2023.10233262](https://doi.org/10.1109/ICCP56744.2023.10233262)  
