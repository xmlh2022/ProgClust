# ProgClust
ProgClust is a clustering algorithm for simultaneously clustering common cells and rare cells. It iteratively performs clustering and dynamically selects fano factor-based feature space or Gini index-based feature space to achieve clustering purpose according to the clustering purpose.

This contains three code files written in R and a folder containing data.

main.R: This file is used to run the ProgClust algorithm.

GiniClust2_packages.R: ProgClust is written based on the GiniClust2 algorithm. This file contains the libraries needed to run GiniClust2.

ProgClust_functions.R: This file contains all the functions needed to run the algorithm.

data: This folder contains 6 simulation data simulated from intestinal cells [1], and an inDrop data from mouse embryonic stem cells 4 days post-LIF[2].

>[1] Grün, D., Lyubimova, A., Kester, L., Wiebrands, K., Basak, O., Sasaki, N., et al. (2015). Single-cell messenger RNA sequencing reveals rare intestinal cell types. Nature, 525(7568), 251-255.

>[2] Klein, A. M., Mazutis, L., Akartuna, I., Tallapragada, N., Veres, A., Li, V., et al. (2015). Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells. Cell, 161(5), 1187-1201.


