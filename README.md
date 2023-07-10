# 3D-RSD
Fast Non-line-of-sight Imaging with Non-planar Relay Surfaces

## Quick start
Clone the repo, change the scene index in 'RSD_main.m' and run. 
## Data description
The data in the 'data' folder contains the following fields:
- **ts**: the temporal resolution of the captured histogram.
- **total_rect_data**: the sum of the histograms captured by all the 16*16 SPAD array pixels.
- **camera_pos**: the scanned laser positions on the relay surface (noted that the name is flipped).
- **laser_pos**: the average SPAD position on the relay surface.