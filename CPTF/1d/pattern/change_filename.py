import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat
from math import exp
import os
import argparse
import itertools
import sys
import copy
import argparse
import glob

parser = argparse.ArgumentParser(description = 'ens for drawing pattern')
parser.add_argument('--p',type=float,help='healing probability')
parser.add_argument('--Db',type=float,help='direction bias')
parser.add_argument('--seed',type=float,help='seed')
#parser.add_argument('--theta',type=float,help='exponent theta')
parser.add_argument('--theta',type=int,help='exponent theta')
#parser.add_argument('--ens',type=int,help='ENS NUMBER')
args = parser.parse_args()

t_min = 0; t_max = 1000;

'''
seed_half/theta1.5/pattern_L60000T_i10000T_f100000ens0p0.6theta1.5R0.65.txt
seed_half/theta1.5/pattern_L10000T_i10000T_f100000ens0p0.6theta1.5R0.65.txt
seed_half/theta1/pattern_L100T_i0T_f10000ens0p0.56theta1R0.5.txt
seed1/theta1/pattern_L100T_i0T_f10000ens4p0.48theta1R0.5Db0.5.txt
seed_half/theta1/pattern_L100T_i0T_f10000ens2p0.56theta1R0.5.txt
'''

txtdirname = "./seed_half/theta"+str(args.theta)
fname =  "/pattern_L100T_i0T_f10000ens*R0.5.txt" #if Db=0.5, Db is null in file name, file name is differnt in seed_half folder...

fname =  "./seed_half/*/pattern_L100T_i0T_f10000ens*R0.5.txt" #if Db=0.5, Db is null in file name, file name is differnt in seed_half folder...

fname_array = glob.glob(txtdirname+fname)
fname_array = glob.glob(fname)
print( fname_array )

for each_file in fname_array :
	print(each_file)
	print("changing...")
	os.rename(each_file, each_file[:-4]+"Db0.5.txt")
	print(each_file)

