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



fname =  "./seed_half/*/pattern_L100T_i0T_f10000*.txt" 

fname_array = glob.glob(fname)
print( fname_array )

for each_file in fname_array :
	print(each_file)
	try:
		pattern = np.loadtxt( each_file )
	except:
		print("error")
	#print(np.size(pattern))
	#print(each_file)
	#print("changing...")
	#os.rename(each_file, each_file[:-4]+"Db0.5.txt")
