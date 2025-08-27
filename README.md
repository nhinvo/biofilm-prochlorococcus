## *Prochlorococcus* Biofilm Dynamics: Taxonomic Classification and Variant Calling  
Repository for pipelines, scripts, and figures for variant calling and taxonomic classification analyses in the publication: 
1. Anjur-Dietrich, M.I., Jones, K.J., Mullet, J.I., **Vo, N.N.**, Castro, K.G., Parker, S.M,, Chisholm, S.W. (2025). Biofilm formation and dynamics in the marine cyanobacterium Prochlorococcus. [bioRxiv](https://doi.org/10.1101/2025.08.05.668435).  


## Microbial Abundance Analysis
All codes for microbial abundance analysis are located in `TaxonomicClassification/`.
- Raw reads are processed through a modified version of [ProSynTax-workflow](https://github.com/jamesm224/ProSynTax-workflow/tree/main) to obtain classification.    
    - Modified version of ProSynTax located in `TaxonomicClassification/Modified-ProSynTax-Workflow/`.  
- Scripts in `TaxonomicClassification/StandardizeData/` were used to perform limit of detection filtering and format data for plotting.  
- Scripts in `TaxonomicClassification/Plots/` were used to generate figures. 
- Quantitative (absolute) read counts were obtained using scripts in `TaxonomicClassification/AbsoluteGenomeEquivalents/`.  

**Figure 6:** 
![Figure 6](https://github.com/nhinvo/nhinvo.github.io/blob/master/assets/img/publications/biofilm-fig6.png)