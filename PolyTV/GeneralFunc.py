#!/usr/bin/env python3
# Author: CC

import numpy as np
import matplotlib.pyplot as plt
import math, random, time, shutil, glob
import sys,re,subprocess,os, argparse
import warnings

""" Cost Functions"""


def cost_epitopes(epitopes_list, vaccine, alleles, netMHCpan_path):

    fasta = ">current_vaccine\n" + vaccine

    output = subprocess.run([netMHCpan_path, "--", "-a", alleles, "-l", "9"], stdout=subprocess.PIPE, universal_newlines=True, input=fasta)
    """ Count the neo-epitopes from the output from NetMHCpan-4.1 and capture failed """
    neo_epitopes = 0
    for line in output.stdout.split("\n"):
        result = re.search(r'\s+\d{1,}',line)
        if (result is not None) and (line.endswith('SB')):
            line_peptide = line.strip().split()[2]
            result = next((True for epitope in epitopes_list if line_peptide in epitope), False)
            if result is False:
                neo_epitopes += 1
    return neo_epitopes 


def cost_linker(vaccine, epitopes_list, weight_linker, thre_linker):
    """ Cost function for the linkers  """
    
    # conserve the index [1st AA, last AA+1] of each epitope
    idx_epitopes= []
    
    for epitope in epitopes_list:
        idx = vaccine.index(epitope)
        idx_epitopes.append([idx, len(epitope)])
    idx_epitopes.sort(key=lambda x: x[0])
    
    # conserve the index [last AA of 1st epi+1, first AA of 2nd epi] for each pair
    # could add AA before the 1st AA of the linker until the 1st AA of the 2nd epi
    tail_head = list()
    for i in range(len(idx_epitopes)-1):
        tail = idx_epitopes[i][0] + idx_epitopes[i][1]
        head = idx_epitopes[i+1][0]
        tail_head.append([tail,head])
    
    len_linkers = []
    for head_tail in tail_head:
        len_linkers.append(head_tail[1]-head_tail[0])
    #print(len_linkers)
    
    cost_l = 0
    for i in len_linkers:
        if i > thre_linker:
            cost_l += (i-thre_linker)
    cost_l = weight_linker * cost_l
    
    return cost_l

def write_o(save_dir, epitopes_list, vaccine, alleles, netMHCpan_path, output_file, num_neo, running_time=None):
    
    fasta = ">current_vaccine\n" + vaccine
    output = subprocess.run([netMHCpan_path, "--", "-a", alleles, "-l", "9"], stdout=subprocess.PIPE, universal_newlines=True, input=fasta)

    neo_list = []
    candidate_list = []
    for line in output.stdout.split("\n"):
        need_info = []
        result = re.search(r'\s+\d{1,}',line)
        if (result is not None) and (line.endswith('B')):
            line_list = line.strip().split()
            line_peptide = line_list[2]
            result = next((True for epitope in epitopes_list if line_peptide in epitope), False)
            need_info.append(line_list[1])
            need_info.append(line_list[2])
            need_info.append(line_list[12])
            if result is False and line.endswith('SB'): 
                neo_list.append(' | '.join(need_info))
            if result is True:
                candidate_list.append(' | '.join(need_info))

    outfile = open(save_dir + '/' + output_file, 'a')
    outfile.write('Polytope vaccine with {} neo-epitopes: '.format(num_neo)+vaccine+'\n')
    outfile.write('Table1. Presented peptides in the selected-epitopes'+'\n')
    outfile.write('Alleles'+'  |  '+'Epi'+'  |  '+"%"+"Rank"+'\n')
    for i in range(len(candidate_list)):
        outfile.write(candidate_list[i]+'\n')
    outfile.write('\n')
    outfile.close()

    
    if len(neo_list) > 0:
        outfile = open(save_dir + '/' + output_file, 'a')
        outfile.write('Table2. Presentation of neo-epitopes'+'\n')
        outfile.write('Alleles'+'  |  '+'Epi'+'  |  '+"%"+"Rank"+'\n')
        for j in range(len(neo_list)):
            outfile.write(neo_list[j]+'\n')  
        outfile.write('\n')
        outfile.close()


    if running_time is not None:
        outfile = open(save_dir + '/' + output_file, 'a')
        outfile.write('Running time is: ' + running_time)
        outfile.write('\n')
        outfile.close()


def plot(save_dir, Big_steps, name_fig, picked_fitness, fitness_list_MC=None):
    if fitness_list_MC is not None:
        picked_fitness.extend(fitness_list_MC)
        plt.plot()
        plt.plot([i for i in range(Big_steps+1)], picked_fitness[:Big_steps+1],color='orange', label = 'Greedy search')
        plt.plot([i for i in range(Big_steps, len(picked_fitness))], picked_fitness[Big_steps:], color='green', label = 'Monte carlo')
        plt.legend(loc="lower left", bbox_to_anchor=(0, 1.02, 1, 0.2),ncol=2)
        plt.ylabel("Number of neo-epitopes")
        plt.xlabel("Iterations")
        plt.savefig(save_dir + '/' + name_fig + '_MC.png',dpi=1000)

    else:
        plt.plot([i for i in range(len(picked_fitness))], picked_fitness,color='orange', label = 'Greedy search')
        plt.legend(loc="lower left", bbox_to_anchor=(0, 1.02, 1, 0.2))
        plt.ylabel("Number of neo-epitopes")
        plt.xlabel("Iterations")
        plt.savefig(save_dir + '/' + name_fig + '_GS.png',dpi=1000)
