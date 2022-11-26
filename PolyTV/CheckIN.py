#!/usr/bin/env python3
# Author: CC

import numpy as np
import matplotlib.pyplot as plt
import math, random, time, shutil, glob
import sys,re,subprocess,os, argparse
import warnings


def check_epitopes_file(epitopes_file):
    AAs = ['G', 'A', 'V', 'L','I','P','F','Y','W','S','T', 'C','M','N','Q','D','E','K','R','H']

    try:
        infile = open(epitopes_file, 'r')
    except IOError as error:
        sys.stderr.write("Can't read file, reason: " + str(error) + "\n")
        sys.exit(1)

    " Check for the epitopes file: 1st: the lens need to be 8-15; 2nd: no newline between epitopes; 3rd: all need to be AA"
    epitopes_list = []
    for line in infile:
        line = line.strip()
        if line != "\n":
            epitope = ''
            if len(line) < 8 or len(line) > 15:
                raise ValueError("The length of epitope should between 8 and 15(both included\n")
            for base in line:
                if base not in AAs:
                    raise ValueError("Your <epitopes file> should only contain amino acids and should be one epitope per line\n")
                else:
                    epitope += base
            epitopes_list.append(epitope)
    infile.close()
    return epitopes_list

def check_alleles_file(alleles_file):
    """Generate alleles file for running netMHCpan"""
    classI = list()
    infile = open('./data/MHC_allele_names.txt','r')
    for line in infile:
        line = line.strip()
        classI.append(line)
    infile.close()

    with open(alleles_file,'r') as f:
        alleles = f.readline()
        your_alleles = alleles.split(',')
        #print(your_alleles)
        if (all(x in classI for x in your_alleles)):
            pass
        else:
            raise ValueError("The alleles you choose should from classI and all alleles should be one line\n")
    return alleles

def check_alleles_command(alleles):

    classI = list()
    infile = open('./data/MHC_allele_names.txt','r')
    for line in infile:
        line = line.strip()
        classI.append(line)
    infile.close()

    your_alleles = alleles.split(',')
    if (all(x in classI for x in your_alleles)):
        pass
    else:
        raise ValueError("The alleles you choose should from classI and all alleles should be one line\n")
    return alleles

def check_linkers_file(linkers_file):
    AAs = ['G', 'A', 'V', 'L','I','P','F','Y','W','S','T', 'C','M','N','Q','D','E','K','R','H']
    try:
        infile = open(linkers_file, 'r')
    except IOError as error:
        sys.stderr.write("Can't read file, reason: " + str(error) + "\n")
        sys.exit(1)  

    linkers_list = []
    for line in infile:
        line = line.strip()
        linker = ''
        for base in line:
            if base not in AAs:
                raise ValueError("Your <linkers file> should only contain amino acids and should be one linker per line\n")
            else:
                if 0 < len(line) < 10:
                    linker += base
                else:
                    raise ValueError("The length of linkers should between 1 and 9(both included\n")
        linkers_list.append(linker)
    linkers_list = [x.lower() for x in linkers_list]
    linkers_list.append('')
    infile.close()
    return linkers_list
