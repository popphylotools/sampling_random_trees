#!/usr/bin/env python

"""
Retrieve fasta sequences from a list of ids and a multifasta file (database).
"""

import argparse,os,sys

import random
from ete3 import Tree
import subprocess
from collections import Counter
import csv


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--input_dir_trees', help='Path to dir containing only input unrooted trees in newick format. The file names must contain the extension .tre')
parser.add_argument('--input_speciestree', help='Path to a speciestree to compare with the subsampled trees.')
parser.add_argument('--number_groups', help='Number of groups of random trees.')
parser.add_argument('--number_trees', help='Number of trees per group.')
parser.add_argument('--clades_monophyly_check', help='File containing two fields sepated by tab. In the first field should be the name of the clade (unique and without spaces) to be checked as monophyletic. In the second fields should be the taxon names (without spaces or special characters) that belong to that are expected to form a monophyletic group. The taxon names must fully match names in the trees.')
parser.add_argument('--astralpath', help='Path to Astral program.')
parser.add_argument('--outgroup', help='Taxon name of the outgroup.')
parser.add_argument('--output_DIR', help='Path to a directory to save the trees.')
parser.add_argument('--output', help='A tab delimited file to store the results.')


if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

working_dir = os.getcwd()

if args.input_dir_trees:
    input_dir_trees = args.input_dir_trees
    if input_dir_trees[-1] != "/":
        input_dir_trees += "/"
    if input_dir_trees[0] != "/":
        input_dir_trees = os.path.join(working_dir,input_dir_trees)

if args.number_groups:
    number_groups = int(args.number_groups)

if args.number_trees:
    number_trees = int(args.number_trees)

if args.output_DIR:
    output_DIR = args.output_DIR
    if output_DIR[-1] != "/":
        output_DIR += "/"
    if not os.path.exists(output_DIR):
        os.makedirs(output_DIR)
    if output_DIR[0] != "/":
        output_DIR = os.path.join(working_dir,output_DIR)

if args.clades_monophyly_check:
    group_file = args.clades_monophyly_check
    selected_clades = {}
    with open(group_file, "r") as f:
        for line in f:
            line_list = line.rstrip().split('\t')
            selected_clades[line_list[0]] = []
            for taxon_name in line_list[1].rstrip().split(","):
                selected_clades[line_list[0]].append(taxon_name)

if args.outgroup:
    outgroup = args.outgroup

if args.astralpath:
    astralpath = args.astralpath

if args.input_speciestree:
    input_speciestree = args.input_speciestree
    if input_speciestree[0] != "/":
        input_speciestree = os.path.join(working_dir,input_speciestree)

if args.output:
    output = args.output
    if output[0] != "/":
        output = os.path.join(working_dir,output)

'''
FUNCTIONS
'''

def gen_random_idx(number_groups,number_elements,number_total_elements,list_ids,tree_dict,output_DIR,species_tree):

    paths_trees = []
    group_clusterIDs = []
    tmp_groups = 1    
    
    group_idxs =[]

    while len(group_idxs) < number_groups:
        group = []
        while len(group) < number_elements:
            random_index = random.randint(0,number_total_elements-1)
            if random_index not in group:
                group.append(random_index)

        group.sort()
        if group not in group_idxs:
            control = 0
            
            #List to get the frequency of the taxa in a group
            sample_info = []        
            
            group_trees_newick = []
            clusterIDs = []
            for idx in group:
                id = list_ids[idx]
                tree_ete3 = tree_dict[id]
                group_trees_newick.append(tree_ete3.write(format=1) + "\n")
                clusterIDs.append(id)
                for node in tree_ete3.traverse():
                    if node.is_leaf():
                        sample_info.append(node.name)

            #Count the taxa in a group
            samples_counter_dict = Counter(sample_info)
        
            for node in species_tree.traverse():
                if node.is_leaf():
                    if node.name in samples_counter_dict.keys():
                        if samples_counter_dict[node.name] < 2:
                            control = 1
                            break
                    else:
                        control = 1
                        break

            if control == 0:
                path_file_name = os.path.join(output_DIR,"random_{}_{}.tre".format(tmp_groups,number_groups))
                with open(path_file_name, "w") as output_file:
                    for tree_newick in group_trees_newick:
                        #Save the results in a new file
                        output_file.write(tree_newick)
        
                group_clusterIDs.append(clusterIDs)
                paths_trees.append(path_file_name)
                tmp_groups = tmp_groups + 1
                group_idxs.append(group)

    return paths_trees,group_clusterIDs

def run_astral(speciestree_dir,astralpath,path_trees):
    path_trees = os.path.normpath(path_trees)
    tree_basename = os.path.basename(path_trees)
    output_speciestree = os.path.join(speciestree_dir,tree_basename + "_spptree.tre")
    log_err_file = os.path.join(speciestree_dir,tree_basename + "_spptree.log")

    #Get the species_tree
    with open(log_err_file, "wb") as err_log:
        cmd = ["java","-jar",astralpath,"-i",path_trees,"-o",output_speciestree] 
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err_log , cwd=speciestree_dir)
        out = p.communicate()
    
    #Annotate the species tree (quartet support)
    output_speciestree_qs = os.path.join(speciestree_dir,tree_basename + "_spptree_qs.tre")
    log_err_file_qs = os.path.join(speciestree_dir,tree_basename + "_spptree_qs.log")

    with open(log_err_file_qs, "wb") as err_log:
        cmd = ["java","-jar",astralpath,"-q",output_speciestree,"-t","1","-i",path_trees,"-o",output_speciestree_qs] 
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err_log , cwd=speciestree_dir)
        out = p.communicate()

    return output_speciestree,output_speciestree_qs

def run_gene_concordance(speciestree_dir,astral_tree,path_trees):
    tree_basename = os.path.basename(path_trees)
    output_prefix = os.path.join(speciestree_dir,tree_basename + "_sptree_gfc")
    log_file_path = os.path.join(speciestree_dir,tree_basename + "_sptree_gfc.log")

    with open(log_file_path, "wb") as log_file:
        cmd = ["iqtree","-t",astral_tree,"--gcf",path_trees,"--prefix",output_prefix] 
        p = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.DEVNULL , cwd=speciestree_dir)
        out = p.communicate()

    return output_prefix

def run_gene_concordance_fixed(speciestree_dir,astral_tree_g,path_trees):
    tree_basename = os.path.basename(path_trees)
    output_prefix = os.path.join(speciestree_dir,tree_basename + "_sptree_gfc_fixed_sptree")
    log_file_path = os.path.join(speciestree_dir,tree_basename + "_sptree_gfc_fixed_sptree.log")

    with open(log_file_path, "wb") as log_file:
        cmd = ["iqtree","-t",astral_tree_g,"--gcf",path_trees,"--prefix",output_prefix] 
        p = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.DEVNULL , cwd=speciestree_dir)
        out = p.communicate()

    return output_prefix

def get_gfc_results(cGF_prefix):
    gCF = []
    with open(cGF_prefix + ".cf.stat", 'r') as file:
        reader = csv.reader(file,delimiter="\t")
        for line in reader:
            if len(line) == 12 and line[1] != 'NA' and line[1] != 'gCF':
                gCF.append(float(line[1]))
    sum_gCF = round(sum(gCF),2)
    number_nodes_gCF = len(gCF)
    avg_gCF = round(float(sum_gCF)/number_nodes_gCF,2)

    return number_nodes_gCF,sum_gCF,avg_gCF

def get_quartetsupport_results(astral_tree,astral_tree_qs):
    t_astral_ete3 = Tree(astral_tree)
    tqs_ete3 = Tree(astral_tree_qs)
    number_nodes = len(t_astral_ete3)
    quartet_s = []
    bootstrap_s = []
    for node in tqs_ete3.traverse():
        if len(node) > 1 and len(node) < number_nodes-1:
            quartet_s.append(node.support)

    for node in t_astral_ete3.traverse():
        if len(node) > 1 and len(node) < number_nodes-1:
            bootstrap_s.append(node.support)

    sum_quartet_s = round(sum(quartet_s),2)
    number_nodes_quartet_s = len(quartet_s)
    avg_quartet_s = round(float(sum_quartet_s)/number_nodes_quartet_s,2)

    sum_bootstrap_s = round(sum(bootstrap_s),2)
    number_nodes_bootstrap_s = len(bootstrap_s)
    avg_bootstrap_s = round(float(sum_bootstrap_s)/number_nodes_bootstrap_s,2)

    return number_nodes_quartet_s,sum_quartet_s,avg_quartet_s,number_nodes_bootstrap_s,sum_bootstrap_s,avg_bootstrap_s

def test_monophyly(tree_ete3,selected_clades):
    results=[]

    for clade in selected_clades.keys():
        if tree_ete3.check_monophyly(values=selected_clades[clade], target_attr="name")[0]:
            results.append("yes")
        else:
            results.append("no")
    return results


'''
Main
'''

tree_dict = dict()
list_ids = []
final_results = []
species_tree = Tree(input_speciestree)

#Create a header
header = ["groupID","number_branches_gcF","total_gCF","avg_gCF","number_branches_f_gcF","total_f_gCF","avg_f_gCF","number_branches_astral","total_PP","avg_PP","number_branches_qs","total_qs","avg_qs","RF","max_RF"]

for node_name in selected_clades.keys():
    header.append(node_name)
header.append("group")
final_results.append(header)


#Using the outgroup to convert species tree in a rooted tree. 
species_tree_1 = Tree(input_speciestree)
species_tree_1.set_outgroup( species_tree_1&outgroup )

#Create a dictionary of trees
for tree_file in os.listdir(input_dir_trees):
    path_tree = os.path.join(input_dir_trees, tree_file)
    id = tree_file.split(".")[0]
    t = Tree(path_tree)
    tree_dict[id] = t
    list_ids.append(id)

number_total_trees = len(list_ids)
paths_trees,group_clusterIDs = gen_random_idx(number_groups,number_trees,number_total_trees,list_ids,tree_dict,output_DIR,species_tree)

for path_trees_idx in range(len(paths_trees)):
    current_result = []
    path_trees = paths_trees[path_trees_idx]
    trees_file = os.path.basename(path_trees)
    #clusterIDs = ','.join(group_clusterIDs[path_trees_idx])
    astral_tree,astral_tree_qs = run_astral(output_DIR,astralpath,path_trees)
    cGF_prefix = run_gene_concordance(output_DIR,astral_tree,path_trees)
    cGF_fixed_prefix = run_gene_concordance_fixed(output_DIR,input_speciestree,path_trees)
    number_nodes_gCF,sum_gCF,avg_gCF = get_gfc_results(cGF_prefix)
    number_nodes_fixed_gCF,number_nodes_fixed_sum_gCF,number_nodes_fixed_avg_gCF = get_gfc_results(cGF_fixed_prefix)
    number_nodes_quartet_s,sum_quartet_s,avg_quartet_s,number_nodes_bootstrap_s,sum_bootstrap_s,avg_bootstrap_s = get_quartetsupport_results(astral_tree,astral_tree_qs)
    
    tree_random_astral = Tree(astral_tree)
    tree_random_astral.set_outgroup( tree_random_astral&outgroup )
    monophyly_test_results = test_monophyly(tree_random_astral,selected_clades)
    rf_result = tree_random_astral.robinson_foulds(species_tree_1)
    
    current_result = [trees_file,str(number_nodes_gCF),str(sum_gCF),str(avg_gCF),str(number_nodes_fixed_gCF),str(number_nodes_fixed_sum_gCF),str(number_nodes_fixed_avg_gCF),str(number_nodes_quartet_s),str(sum_quartet_s),str(avg_quartet_s),str(number_nodes_bootstrap_s),str(sum_bootstrap_s),str(avg_bootstrap_s),str(rf_result[0]),str(rf_result[1])]
    for result_monophyly in monophyly_test_results:
        current_result.append(result_monophyly)
    current_result.append(number_trees)

    final_results.append(current_result)

with open(output, 'w') as file:
    writer = csv.writer(file,delimiter ="\t")
    writer.writerows(final_results)

