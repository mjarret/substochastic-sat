import subprocess
import os
import argparse
import csv

# Some parsers for command line options
parser = argparse.ArgumentParser(description='Autoselect parameters.')
parser.add_argument('pp',  metavar='pp',  nargs='?',  help='base directory',  default ='./')
args=parser.parse_args()

path = args.pp
outfile = "./chooser.out"
infofile = open(outfile, 'a+') # writing file
writer = csv.writer(infofile)

def toOptimize(weight,  runtime,  popsize,  starttime):
    nspecies = 1
    ntrials = 1
    hits=1
    for i in os.listdir(path):
        if i.endswith(".cnf"):
            filename = path + i
            p = subprocess.Popen(["./substochastic", filename , str(weight),  str(runtime),  str(popsize),  str(nspecies),  str(ntrials),  str(starttime)],  stdout = subprocess.PIPE)
            out, err = p.communicate()
            hits*=int(out.split('\n')[-2])
    if (hits == 1): return 0
    return hits

for weight in range(10, 25, 1):
    for runtime in range(6500, 10000, 250):
        for starttime in range(0, int(.40*runtime), int(0.05*runtime)):
            writer.writerow((weight, runtime, starttime, toOptimize(weight, runtime, 64, starttime)))
