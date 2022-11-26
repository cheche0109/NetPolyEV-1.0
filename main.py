#!/usr/bin/env python3
# Author: CC

import numpy as np
import matplotlib.pyplot as plt
import math, random, time, shutil, glob
import sys,re,subprocess,os, argparse
import warnings

from PolyTV.CheckIN import check_epitopes_file, check_alleles_file, check_alleles_command, check_linkers_file
from PolyTV.GreedySearch import repeat_GS
from PolyTV.MonteCarlo import small_move
from PolyTV.GeneralFunc import write_o, plot

########################################################

netMHCpan_path = "/Users/chenchen/Documents/DTU/PolyTope_Optimizer/netMHCpan-4.1/netMHCpan"

########################################################

parser = argparse.ArgumentParser(description='Polytope vaccine optimizer')
parser.add_argument('-e', dest='epitopesFile', type=str, required=True, help='Input selected epitopes file')
parser.add_argument('-a', dest='allelesFile', type=str, required=False, help='Input selected alleles file')
parser.add_argument('-as', '--Alleles', type=str, required=False, help='Type selected alleles')
parser.add_argument('-l', dest='LinkerFile',  type=str, required=False, help='Input selected linker file')
parser.add_argument('-fig', dest='nameFigure',  type=str, required=True, help='Name the plot')
parser.add_argument('-o', dest='OutputFile',type=str, required=True, help='Output file recording the vaccine and its presented epitopes')
parser.add_argument('-dir', action='store', dest='saveDirectory', type=str, required=True, help='Directory: plot(s) and output file(s)')

parser.add_argument("-se", action="store", dest="SEED", type=int, default=24, help="Random seed")
parser.add_argument("-bt", action="store", dest="Big_STEPS", type=int, default=100, help="Big steps in greedy search")
parser.add_argument("-gs", action="store", dest="Repeat_GS", type=int, default=5, help="Times for repeat greedy search")
parser.add_argument("-st", action="store", dest="T_STEPS", type=int, default=30, help="Finite steps at each temperature")
parser.add_argument("-cf", action="store", dest="CoolDown_Factor", type=float, default=0.11, help="Cooldown factor for the MC")

parser.add_argument("-wl", action="store", dest="Cost_Linker", type=float, default=0, help="Weight for the cost(linker)")
parser.add_argument("-tl", action="store", dest="Thre_Linker", type=int, default=6, help="The threshold of the length of the linker")
args = parser.parse_args()

epitopes_file = args.epitopesFile
alleles_file = args.allelesFile
linkers_file = args.LinkerFile
output_file = args.OutputFile
save_dir = args.saveDirectory
Big_steps = args.Big_STEPS
max_GS = args.Repeat_GS
T_steps = args.T_STEPS
eta = args.CoolDown_Factor
weight_linker = args.Cost_Linker
thre_linker = args.Thre_Linker
name_fig = args.nameFigure
num_seed = args.SEED

########################################################

if __name__ == '__main__':

    random.seed(num_seed)

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    """ Check inputs """
    epitopes_list = check_epitopes_file(epitopes_file)

    if args.allelesFile:
        alleles = check_alleles_file(alleles_file)
    elif args.Alleles:
        alleles = args.Alleles
        alleles = check_alleles_command(alleles)

    if args.LinkerFile:
        linkers_list = check_linkers_file(linkers_file)

    """ GreedySearch: Algorithm1 """
    if args.LinkerFile:
        running_time_GS, initial_vaccine, initial_neo, vaccine_o1, num_neo, picked_fitness = repeat_GS(max_GS, Big_steps, epitopes_list, alleles, netMHCpan_path, linkers_list)
    else:
        running_time_GS, initial_vaccine, initial_neo, vaccine_o1, num_neo, picked_fitness = repeat_GS(max_GS, Big_steps, epitopes_list, alleles, netMHCpan_path, linkers_list=None)

    outfile = open(save_dir + '/' + output_file, 'a')
    outfile.write('### Randomly produced the initial vaccine ###'+'\n')
    outfile.close()
    write_o(save_dir, epitopes_list, initial_vaccine, alleles, netMHCpan_path, output_file, initial_neo)

    outfile = open(save_dir + '/' + output_file, 'a')
    outfile.write('### Optimized vaccine after the greedy search ###'+'\n')
    outfile.close()
    write_o(save_dir, epitopes_list, vaccine_o1, alleles, netMHCpan_path, output_file, num_neo, running_time_GS)

    """ MonteCarlo: Algorithm2 & Plot the whole optimizing process """

    if num_neo > 0:
        running_time_MC, vaccine_o2, neo_epitopes, fitness_list_MC = small_move(T_steps, eta, epitopes_list, alleles, netMHCpan_path, vaccine_o1, weight_linker, thre_linker)
    
        outfile = open(save_dir + '/' + output_file, 'a')
        outfile.write('### Optimized vaccine after the monte carlo ###'+'\n')
        outfile.close()
        write_o(save_dir, epitopes_list, vaccine_o2, alleles, netMHCpan_path, output_file, neo_epitopes, running_time_MC)
        plot(save_dir, Big_steps, name_fig, picked_fitness, fitness_list_MC)
    else:
        plot(save_dir, Big_steps,name_fig, picked_fitness, fitness_list_MC=None)



