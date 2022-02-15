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

python sampling_trees.py --input_dir_trees example/gene_trees/ --input_speciestree example/reference_speciestree.tree --number_groups 10 --number_trees 20 --clades_monophyly_check example/clades_for_checking_monophyly  --astralpath full_path_to_astral.5.7.7.jar --outgroup A.bistrigata --output_DIR test_group_random_20_10 --output test_random_test_20_10.out

OUTPUT:

The main output is a tab separated csv file. The fields are:

1 groupID: file name of the subset of gene trees
2 number_branches_gcF: Number of branches
3 total_gCF: Sum of the gene concordance factor of each branch of the species tree inferred from a subset
4 avg_gCF: Average of the gene concordance factor of each branch of the species tree inferred from a subset
5 number_branches_f_gcF: Number of branches
6 total_f_gCF: Sum of the gene concordance factor of each branch of the reference species tree (--input_speciestree)
7 avg_f_gCF: Average of the gene concordance factor of each branch of the reference species tree (--input_speciestree)
8 number_branches_astral: Number of branches
9 total_PP: Sum of the local posterior probabilities of the species tree inferred from a subset
10 avg_PP: Average of the local posterior probabilities of the species tree inferred from a subset
11 number_branches_qs: Number of branches
12 total_qs: Sum of the quartet supports of the species tree inferred from a subset
13 avg_qs: Average of the quartet supports of the species tree inferred from a subset
14 RF: Robinson-Foulds distance of the species tree of a subset and the reference species tree (--input_speciestree)
15 max_RF: Maximum Robinson-Foulds distance of the species tree of a subset and the reference species tree (--input_speciestree)

The other fields would be results of checking the monophyly of the groups provided by the user (--clades_monophyly_check).


References

Huerta-Cepas, J., Serra, F., & Bork, P. (2016). ETE 3: Reconstruction, analysis, and visualization of phylogenomic data. Molecular Biology and Evolution, 33(6), 1635–1638. https://doi.org/10.1093/molbev/msw046

Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5), 1530–1534. https://doi.org/10.1093/molbev/msaa015

Zhang, C., Rabiee, M., Sayyari, E., & Mirarab, S. (2018). ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinformatics, 19(6), 153. https://doi.org/10.1186/s12859-018-2129-y
