#!/usr/bin/env python3
# Author: CC

import numpy as np
import matplotlib.pyplot as plt
import math, random, time, shutil, glob
import sys,re,subprocess,os, argparse
import warnings

sys.path.append(os.path.join(sys.path[0], 'PolyTV'))
import GeneralFunc as GF
###############################################
#### big jumps for having the linker input ####
###############################################

def polytope_generator_linker(epitopes_list, linkers_list):
    """Shuffle the seleted epitopes to reoder and insert linkers from the function linker()"""

    counts = len(epitopes_list)
    random.shuffle(epitopes_list) # shuffle to get random order of the epitopes

    polytope_list = list() 
    for i in range(0, counts-1): 
        polytope_list.append(epitopes_list[i])
        polytope_list.append(random.choice(linkers_list)) # to insert a linker between epitopes
    polytope_list.append(epitopes_list[-1])
    polytope = ''.join(polytope_list)
    return polytope


######################################################
#### big jumps for creating linkers autonactically ###
######################################################

def linker():
    """Generate a linker by randomly picking 0-5 that every AA can be repetitive"""

    AAs = ['G', 'A', 'V', 'L','I','P','F','Y','W','S','T', 'C','M','N','Q','D','E','K','R','H']
    AAs = [x.lower() for x in AAs]
    num_AAs = random.randint(0,5) # the length of linker: [0,5]
    picked_AAs = random.choices(AAs,k=num_AAs) # randomly choose num_AAs AAs: repeatable
    linker = ''.join(picked_AAs)
    return linker


def polytope_generator(epitopes_list):
    """Shuffle the seleted epitopes to reoder and insert linkers from the function linker()"""

    counts = len(epitopes_list)
    random.shuffle(epitopes_list) # shuffle to get random order of the epitopes

    polytope_list = list() 
    for i in range(0, counts-1): 
        polytope_list.append(epitopes_list[i])
        polytope_list.append(linker()) # to insert a linker between epitopes
    polytope_list.append(epitopes_list[-1])
    polytope = ''.join(polytope_list)
    return polytope


def big_jump(Big_steps, epitopes_list, alleles, netMHCpan_path, linkers_list=None):


    initial_list = list()
    #Generate the first random vaccine
    if linkers_list==None:
        """ let the programming create linkers for you"""
        vaccine_old = polytope_generator(epitopes_list) # no linkers.txt
    else:
        """ have your own preferred linkers"""
        vaccine_old = polytope_generator_linker(epitopes_list, linkers_list) # big jump for have linkers.txt
    
    neo_epitopes_old = GF.cost_epitopes(epitopes_list, vaccine_old, alleles, netMHCpan_path)
    initial_list.append([neo_epitopes_old, vaccine_old])


    #Conserve the number of neo-epitopes for every iterations from the Monte Carlo simulated anealing approach
    fitness_list1 = list()
    fitness_list1.append(neo_epitopes_old)
    #itera = 100
    for i in range(Big_steps):
        if linkers_list==None:
            vaccine_new = polytope_generator(epitopes_list)
            """ have your own preferred linkers"""        
        else:
            """ let the programming create linkers for you"""
            vaccine_new = polytope_generator_linker(epitopes_list, linkers_list)# big jump for have linkers.txt
        
        #fsa_generator(vaccine_new, vaccine_path)
        neo_epitopes_new = GF.cost_epitopes(epitopes_list, vaccine_new, alleles, netMHCpan_path)
        E = neo_epitopes_new - neo_epitopes_old
        if E <= 0:
            # replace the old one with the new one
            vaccine_old = vaccine_new
            neo_epitopes_old = GF.cost_epitopes(epitopes_list, vaccine_old, alleles, netMHCpan_path)
        fitness_list1.append(neo_epitopes_old)
        if neo_epitopes_old == 0:
            break

    return vaccine_old, neo_epitopes_old, fitness_list1, initial_list


def repeat_GS(max_GS, Big_steps, epitopes_list, alleles, netMHCpan_path, linkers_list=None):
    Vaccine_Neo_epi = list()
    circulation = 1
    Fitness_List = dict()
    Initial_List = list()

    start_big_jump = time.time()
    
    while circulation <= max_GS:

        if linkers_list==None:
            vaccine_old, neo_epitopes_old, fitness_list, initial_list= big_jump(Big_steps, epitopes_list, alleles, netMHCpan_path, linkers_list=None)
        else:
            vaccine_old, neo_epitopes_old, fitness_list, initial_list= big_jump(Big_steps, epitopes_list, alleles, netMHCpan_path, linkers_list)
        Vaccine_Neo_epi.append([vaccine_old,neo_epitopes_old,circulation])
        Initial_List.append(initial_list[0])
        Fitness_List[circulation] = fitness_list
        circulation += 1
        if neo_epitopes_old == 0:
            break

    Vaccine_Neo_epi.sort(key=lambda x: x[1])
    vaccine = Vaccine_Neo_epi[0][0]
    num_neo = Vaccine_Neo_epi[0][1]
    picked_circulation = Vaccine_Neo_epi[0][2]
    initial_vaccine = Initial_List[picked_circulation-1][1]
    initial_neo = Initial_List[picked_circulation-1][0]

    for i in Fitness_List.keys():
        if i == picked_circulation:
            picked_fitness = Fitness_List[i]
            break

    end_big_jump = time.time()

    if end_big_jump-start_big_jump > 60:
        running_time = "{:.0f}mins".format((end_big_jump-start_big_jump)/60)
    else:
        running_time = "{:.0f}s".format(end_big_jump-start_big_jump)

    return running_time, initial_vaccine, initial_neo, vaccine, num_neo, picked_fitness




