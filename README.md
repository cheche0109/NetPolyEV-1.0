# NetPolyEV-1.0
Construction the polytope vaccines and filter the MHC-peptide binding by applying NetMHCpan-4.1
## Webserver
    https://services.healthtech.dtu.dk/services/NetPolyEV-1.0/
## Contributor
    Supervisor: Morten Nielsen (morni@dtu.dk)
    Student: Chen Chen (chenchenus2020@outlook.com)
## INSTALLATION
    Python > 3.5
    NetMHCpan-4.1 (https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1)
    Need to change netmhcpan path in main.py
## RUNNING
    python main.py
    optional arguments:
    -e EPITOPESFILE       Input selected epitopes file
    -a ALLELESFILE        Input selected alleles file
    -as ALLELES           Type selected alleles on commnad
    -l LINKERFILE         Input selected linker file
    -fig NAMEFIGURE       Name the optimizing process plot
    -o OUTPUTFILE         Output file recording the vaccine and its presented epitopes
    -dir SAVEDIRECTORY    Directory contains the plot(s) and output file(s)
    -se SEED              Random seed
    -bt BIG_STEPS         Iterations in greedy search, default is 100
    -gs REPEAT_GS         Times for repeat greedy search, default is 5
    -st T_STEPS           Finite steps at each temperature, default is 30
    -cf COOLDOWN_FACTOR   Cooldown factor for the MC (Initial temperature is 10, and the lowest temperature is 0.1)
    -wl COST_LINKER       Weight for the cost(linker), default is 0
    -tl THRE_LINKER       The threshold of the length of the linker, default is 9
    
## DISCLAIMER
    The executable software and the source code of NetPolyEV-1.0 is distributed free of charge as it is to any non-commercial users. 
