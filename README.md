# sampling_random_trees
This is a python 3 script designed to randomly sample gene trees, produce a species tree from each subset and compare those phylogenies to a reference species tree. To assure a minimum level of phylogenetic consistency, all subsets had the information of each specimen in at least two gene trees. The species tree for each subsample is estimated by ASTRAL v. 5.7.7 (Zhang et al., 2018) and the gene concordance factor is calculated using IQ-TREE v. 2.1.2 (Minh, et al., 2020). Moreover, the topologies of the rooted species trees of the subsets are compared to a reference species tree using Robinson-Foulds distance. This program also check for the monophyly of predefined groups. 

######################################
Created By: Carlos Congrains and Scott Geib. Email: carloscongrains@gmail.com
######################################

REQUIREMENTS:

We developed and tested this script using Python version 3. It requires ETE3 module (http://etetoolkit.org/download/). The following programs must be installed in the system:

ASTRAL v. 5.7.7 (https://github.com/smirarab/ASTRAL)
IQ-TREE v. 2.1.2 (http://www.iqtree.org)

USAGE:








We performed this assessment using a custom python script (), which employed tools implemented in Environment for Tree Exploration (ETE) v. 3. (Huerta-Cepas et al., 2016)


trees_file,str(number_nodes_gCF),str(sum_gCF),str(avg_gCF),str(number_nodes_fixed_gCF),str(number_nodes_fixed_sum_gCF),str(number_nodes_fixed_avg_gCF),str(number_nodes_quartet_s),str(sum_quartet_s),str(avg_quartet_s),str(number_nodes_bootstrap_s),str(sum_bootstrap_s),str(avg_bootstrap_s),str(rf_result[0]),str(rf_result[1])]


References

Huerta-Cepas, J., Serra, F., & Bork, P. (2016). ETE 3: Reconstruction, analysis, and visualization of phylogenomic data. Molecular Biology and Evolution, 33(6), 1635–1638. https://doi.org/10.1093/molbev/msw046

Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5), 1530–1534. https://doi.org/10.1093/molbev/msaa015

Zhang, C., Rabiee, M., Sayyari, E., & Mirarab, S. (2018). ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinformatics, 19(6), 153. https://doi.org/10.1186/s12859-018-2129-y
