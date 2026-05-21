import sys
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from multiprocessing import Pool
import seaborn as sns
import numpy as np

aim_file = sys.argv[1]
data_file = sys.argv[2]

data = open(data_file).read().replace(">","$>").split("$")[1:]
data_dict = {}
for gene in data:
	spa = gene.splitlines()
	data_dict[spa[0].strip().split()[0].strip(">")] = "".join([l.strip() for l in spa[1:]]).upper().replace("U","T")

aim = open(aim_file).read().replace(">","$>").split("$")[1:]
aim_dict = {}
for aim in aim:
	spa = aim.splitlines()
	aim_dict[spa[0].strip().split()[0].strip(">")] = "".join([l.strip() for l in spa[1:]]).upper().replace("U","T")

deal_list = []
for ag in aim_dict:
	for dg in data_dict:
		deal_list.append([ag, dg])

def align_func(d_list):
#	d_list[0] = a_g
#	d_list[1] = d_g
	alignments = pairwise2.align.globalms(aim_dict[d_list[0] ],data_dict[d_list[1] ],2, -1, -2, -2,penalize_extend_when_opening=True,penalize_end_gaps=False)
	align_info = format_alignment(*alignments[0],full_sequences=True).splitlines()[1]
#	print(align_info)
	align_len = len(align_info.strip())
	ident_len = align_info.count("|")
	ident_value = round(ident_len/align_len,4)
	return d_list[0] + "\t" + d_list[1] + "\t"  + str(ident_value)

pool_temp = []
if __name__ == '__main__':
	with Pool(processes=4) as pool:
		result = pool.map(align_func, deal_list)
		pool_temp = result

empty_df = pd.DataFrame(columns=list(aim_dict.keys()), index=list(data_dict.keys()))
for line in pool_temp:
	spa = line.split()
	empty_df.loc[spa[1],spa[0]] = float(spa[2])
	
empty_df.to_csv(aim_file + ".csv",index=True)
