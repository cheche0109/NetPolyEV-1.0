#!/usr/bin/env python3
# Author: CC
import numpy as np
import matplotlib.pyplot as plt
import math, random, time, shutil, glob
import sys,re,subprocess,os, argparse
import warnings

sys.path.append(os.path.join(sys.path[0], 'PolyTV'))
import GeneralFunc as GF
######################################################
#### small move for creating linkers autonactically ##
######################################################


### Do for the epitopes ###
def reorder_epi(vaccine, epitopes_list):
    """ Pick random 2 selected epitopes to reorder """
    """ include the situation: the length of the 2 selected epitopes are different """

    random_2_epi = random.sample(epitopes_list,2) # random choose 2 epitopes
    idx_first = vaccine.find(random_2_epi[0]) # idx the pos in vaccine
    idx_second = vaccine.find(random_2_epi[1])

    # idx the 2 epitopes in vaccine sequentially.
    if idx_first < idx_second:
        idx_1st = idx_first
        epi1 = random_2_epi[0]
        idx_2nd = idx_second
        epi2 = random_2_epi[1]
    else:
        idx_1st = idx_second
        epi1 = random_2_epi[1]
        idx_2nd = idx_first
        epi2 = random_2_epi[0]

    # replace the 1st epitope with the 2nd & cut the vaccine
    vac_ = vaccine[: idx_1st] + epi2
    # replace the 2nd epitope with the 1st & cut the vaccine
    _cine = vaccine[idx_1st+len(epi1): idx_2nd] + epi1 + vaccine[idx_2nd + len(epi2): ]

    #first half + last half = whole vaccine
    vaccine = vac_ + _cine
    return vaccine


### Do for the linkers ###
### only the the pair epitopes which has linker between them could apply delete(vaccine) and change(vaccine)
## This pre-condition is considered in random_1_func(vaccine, epitopes_list)

#cannot if the vaccine doesn't have any spacers
def delete(vaccine):
    """ Delete random 1 amino acid within one of the linkers """
    idx_spacer= []

    for idx in range(len(vaccine)):
        base = vaccine[idx]
        if base.islower():
            idx_spacer.append(idx)
    delete_idx = random.choice(idx_spacer)

    vaccine = vaccine[:delete_idx] + vaccine[delete_idx+1:]
    return vaccine

def add(vaccine, epitopes_list):
    """ Add random 1 amino acid within one of the linkers """
    AAs = ['G', 'A', 'V', 'L','I','P','F','Y','W','S','T', 'C','M','N','Q','D','E','K','R','H']
    AAs = [x.lower() for x in AAs]
    
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

    random_pair = random.choice(tail_head) # select random pair
    ender = random_pair[0] # the idx for the first epitope in the selected pair
    starter = random_pair[1] # the idx for the second epitope in the selected pair
    if (starter-ender) == 0: # for the pair doesn't contain linker
        vaccine = vaccine[:ender]+random.choice(AAs)+vaccine[ender:]
    else:
        idx_insert = random.randint(0, starter-ender) + ender # range: from before first to after last
        vaccine = vaccine[:idx_insert]+random.choice(AAs)+vaccine[idx_insert:]
    return vaccine

# cannot if the vaccine doesn't have any spacers
def change(vaccine):
    """  Change random 1 amino acid within one of the linkers """

    AAs = ['G', 'A', 'V', 'L','I','P','F','Y','W','S','T', 'C','M','N','Q','D','E','K','R','H']
    AAs = [x.lower() for x in AAs]
    idx_spacer= []
    for idx in range(len(vaccine)):
        base = vaccine[idx]
        if base.islower():
            idx_spacer.append(idx)
    #print(idx_spacer)
    change_idx = random.choice(idx_spacer)
    vaccine = vaccine[:change_idx] + random.choice(AAs) + vaccine[change_idx+1:]    
    return vaccine


# totally 4 functions for small moves
# if the vaccine has no spacers, then can only call 2 functions: reorder or add

def random_1_func(vaccine, epitopes_list):

    vaccine_no_linker = 0
    for epitope in epitopes_list:
        vaccine_no_linker += len(epitope)

    vaccine_len = len(vaccine)
    
    idx = random.randint(1,4) # both included
    
    if idx == 1 and vaccine_len > vaccine_no_linker:
        #print('delete')
        vaccine = delete(vaccine)
    elif idx == 2:
        #print('add')
        vaccine = add(vaccine,epitopes_list)
    elif idx == 3 and vaccine_len > vaccine_no_linker:
        #print('change')
        vaccine = change(vaccine)
    elif idx == 4:
        #print('reorder')
        vaccine = reorder_epi(vaccine,epitopes_list)
    else:
        idx = random.randint(1,2)
        if idx == 1:
            #print('else-add')
            vaccine = add(vaccine,epitopes_list)
        else:
            #print('else-reorder')
            vaccine = reorder_epi(vaccine, epitopes_list)
    #print(vaccine)
    return vaccine



##########small moves######
def small_move(T_steps, eta, epitopes_list, alleles, netMHCpan_path, vaccine, weight_linker, thre_linker):

    start_small_move = time.time()

    initT = 10
    miniT = 0.1
    #itera = 30
    #eta = 0.11 # Cooling schedule
    t = initT

    #Generate the first random vaccine
    vaccine_old = vaccine
    #fsa_generator(vaccine_old, vaccine_path)
    neo_epitopes_old = GF.cost_epitopes(epitopes_list, vaccine_old, alleles, netMHCpan_path)
    #cost_of_linker = cost_linker(vaccine_old, epitopes_list, weight_linker, thre_linker)
    E_old = neo_epitopes_old #+ cost_of_linker

    #Conserve the number of neo-epitopes for every iterations from the Monte Carlo simulated anealing approach
    fitness_list2 = []
    fitness_list2.append(E_old)

    while t > miniT:
        
        for i in range(T_steps): # monte carlo accept probability
            vaccine_new = random_1_func(vaccine_old, epitopes_list)
            neo_epitopes_new = GF.cost_epitopes(epitopes_list, vaccine_new, alleles, netMHCpan_path)
            cost_of_linker = GF.cost_linker(vaccine_new, epitopes_list, weight_linker, thre_linker)
            E_new = neo_epitopes_new + cost_of_linker
            #shutil. rmtree('fsafile')
            E = E_new - E_old
            
            if E <= 0 or math.exp(-E/t) > random.random():
                """Because of the minimization optimization, so we have minus here"""
                # replace the old one with the new one
                vaccine_old = vaccine_new
                neo_epitopes_old = GF.cost_epitopes(epitopes_list, vaccine_old, alleles, netMHCpan_path)
                cost_of_linker = GF.cost_linker(vaccine_new, epitopes_list, weight_linker, thre_linker)
                E_old = neo_epitopes_old + cost_of_linker
            fitness_list2.append(E_old)
            if E_old == 0:
                break

        t = t - eta
        if E_old == 0:
            break

    end_small_move = time.time()

    running_time = "{:.0f}mins".format((end_small_move - start_small_move)/60)

    return running_time, vaccine_old, neo_epitopes_old, fitness_list2


