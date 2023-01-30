# About scDecouple
scDecouple aims to decouple the true cellular response of the perturbation from the influence of infection bias. It models the distribution of perturbed cells and iteratively finds the maximum likelihood of cell cluster proportions as well as the real cellular response for each gRNA. 

# Background
scCRISPR-seq is an emerging high-throughput CRISPR screening technology that combines CIRPSR screening with single-cell sequencing technologies. It provides rich information on gene regulation. When performing scCRISPR-seq in a population of heterogeneous cells, the observed cellular response in perturbed cells may be caused not only by the perturbation, but also by the infection bias of guide RNAs (gRNAs) mainly contributed by intrinsic differences of cell clusters. The mixing of these effects poisons gene regulation studies. 

# Pipeline
scDecouple contains four steps: data preprocessing, PC selection, decoupling, and analysis. 
1. Cells are normalized and log-transformed. Variable genes are selected and transformed to PC space. 
2. PCs with high multimodality and explained variance are selected. The multimodality is defined by dip statistic. 
3. The decoupling process is performed on selected PCs to decouple observed changes after the establishment of the gaussian mixture model of control and perturbation groups (. 4. we calculated the cellular response on other PCs by observed FC and performed inverse transformation on all PCs to estimate cellular responses per gene, followed by pathway enrichment and perturbation-related gene ranking. 

# Author
Qiuchen Meng
mqc17@mails.tsinghua.edu.cn