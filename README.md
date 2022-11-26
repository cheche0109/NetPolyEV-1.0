
                                   
				 PolyVaccine      |     INSTRUCTION

DESCRIPTION

The PolyVaccine is used for constructing the polytope vaccine, and our program is proposed for focusing on epitope assembly to design polytope vaccines by mainly implementing monte carlo simulated annealing approach, both the order of selected CD8+ T cell epitopes and the spacers with flexible length and sequence are optimized for the polytope vaccine, and the epitope-HLA class I binder is predicted by NetMHCpan-4.1. The method from the tool that is trained with data retrieved from mass spectrometry (MS) experiments known as eluted ligands (EL) is selected for our prediction work because it includes information both about the binding of peptide-MHC, the most selective step in the antigen presentation pathway, and about the prior steps of the antigen presentation pathway. 9 amino acids long, the most common length of the presented MHC class I ligands, is chosen for predicting the epitopes of the polytope vaccine for us. If the peptides are strong binders(SBs) or weak binders(WBs) will be informed based on the %Rank score. In our work, we choose the default threshold %Rank < 0.5 for SBs and %Rank < 2 for WBs for class I.

Webserver:https://services.healthtech.dtu.dk/service.php?PolyVaccine-1.0
We recommend you use our .py programming for the construction because we have 1-submission restriction for our webserver.

#################################
######## ! First step ! #########
#################################
First download the NetMHCpan-4.1 from https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1
and change the path for yourself in our .py file PolyVaccine.py
This should be the only thing you change in our codes.

#################################
####### Your input file #########
#################################
At least you need to have 2 input data: epitopes.txt and alleles.txt or typing in the alleles in the command line. 1) The selected epitopes should be one per line and the selected alleles should be one line and separated by ',', check the examples in our data directory  2) All alleles should be from Class I which are included in MHC_allele_names.txt in the data directory. 2) The selected epitopes should match the selected alleles, otherwise, our program will stop because the NetMHCpan-4.1 cannot run due to the mismatch problem. 3) Of course the epitope should only contain the amino acids, otherwise, our program cannot run.

Optional input file: linkers.txt, you can also have the input file for your preferred linkers, 1) The linkers should be one linker per line in the file. 2) the linker should only contain amino acids, otherwise, our program cannot run. 3) the length of the linker should be less than 10. 3) Only Algorithm 1 will be applied to your preferred linkers, the Algorithm 2 will modify your linkers. 4) Check the example file in the data directory.

#################################
####### Input parameters ########
#################################
mandatory inputs:
-e epitopes.txt
-a alleles.txt      (or: -as XXXXXXXXXXXX)
(XXXXXXXXXXXX should be replaced with your selected alleles and should in the format such as HLA-A02:01,HLA-A01:01,HLA-A03:01,HLA-A24:02)

Optinal inputs:
-l linkers.txt
-bt default value is 100, the iteration of the Algorithm 1 (Greedy search)
-st default value is 30, the finite steps at each temperature in the Algorithm 2 (Monte carlo)
-se default value is 24, the random seed
-gs default value is 6, the times you wanna repeat the Algorithm1, the default value is 6 means this Algorithm will be repeated for 5 times and the one produces the optimized vaccine with fewest neo-epitopes will be choosen for keeping optimizing.
-cf default value is 0.11, the cooldown factor for the Algorithm 2
-wl default value is 0, the weighting you give to the cost function for the length of linkers
-tl defualt value is 6, the longest linker you wanna have when -wl is given


###########################################
### Examples for running our program: #####
###########################################
No cost function for the length of the linkers and no linkers input:
nohup python PolyVaccine.py -e ./data/epitopes_covid.txt -a ./data/alleles.txt > ./test/covid.log &
nohup python PolyVaccine.py -e ./data/epitopes_covid_15mer.txt -a ./data/alleles.txt > ./test/covid_15mer.log &

nohup python PolyVaccine.py -e ./data/epitopes_hcv.txt -a ./data/alleles.txt > ./test/hcv.log &
nohup python PolyVaccine.py -e ./data/epitopes_flu.txt -a ./data/alleles.txt > ./test/flu.log &
nohup python PolyVaccine.py -e ./data/epitopes_smallpox.txt -a ./data/alleles.txt -st 55 > ./test/smallpox.log &
nohup python PolyVaccine.py -e ./data/epitopes_rsv.txt -as HLA-A02:01,HLA-A01:01,HLA-A03:01,HLA-A24:02,HLA-A11:01,HLA-B07:02,HLA-B51:01,HLA-B08:01,HLA-B35:01,HLA-B44:02 > ./test/rsv.log &

No cost function for the length of the linkers but gave linkers input:
nohup python PolyVaccine.py -e ./data/epitopes_hcv.txt -a ./data/alleles.txt -l ./data/linkers.txt> ./test/hcvL.log &

Cost function for the length of the linkers:
nohup python PolyVaccine.py -e ./data/epitopes_covid.txt -a ./data/alleles.txt -wl 1 > ./test/covid_LCost.log &

NOTE:
The result log file you obtained, should correspond with the one with the suffix '_test.log' in our test directory.

#############
### Note: ###
#############
1)
test directory: storing for the .log files got from the PolyVaccine.py

2)
In case of failing to read the selected alleles file (some computers need the newline after the alleles, some not), the user can type the alleles on the command line after '-as' and the example already given. test directory: storing for the .log files got from the PolyVaccine.py

3)
The generated iteration plot will just be under the 'PolyVaccine' folder, either named GS.png which means the monte carlo wasn't applied, or GS_MC.png which means both the greedy search and monte carlo were applied.

############
# PROBLEMS #
############

Contact Dr.Morten Nielsen, morni@dtu.dk or the student supervised by him Chen Chen, chenchen0109cc@gmail.com

DTU, summer of 2022
