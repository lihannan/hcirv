import sys
import pandas as pd
in_file = sys.argv[1]

df = pd.DataFrame()

data = open(in_file).read().replace(">","$>").split("$")[1:]


seq_dict = {}
for seq in data:
	spa = seq.splitlines()
	df[spa[0].split(".")[0].strip(">")] = list(spa[1].strip())

temp = ["position\tA\tT\tG\tC\t-\tmutate_ratio\tref_base\tmutate_base"]
for i in range(0,df.shape[0]):
	temp_list = []
	temp_dict = {}
	temp_list.append(i+1)
	temp_list.append(df.iloc[i].tolist().count("A"))
	temp_dict["A"] = df.iloc[i].tolist().count("A")
	temp_list.append(df.iloc[i].tolist().count("T"))
	temp_dict["T"] = df.iloc[i].tolist().count("T")
	temp_list.append(df.iloc[i].tolist().count("G"))
	temp_dict["G"] = df.iloc[i].tolist().count("G")
	temp_list.append(df.iloc[i].tolist().count("C"))
	temp_dict["C"] = df.iloc[i].tolist().count("C")
	temp_list.append(df.iloc[i].tolist().count("-"))
	temp_dict["-"] = df.iloc[i].tolist().count("-")
	if sum(temp_list[1:]) == 0:
		continue
	n = 1
	deal_list = sorted(temp_list[1:])
	ref_base = temp[0].split("\t")[1:][temp_list[1:].index(sorted(temp_list[1:])[-n] )]
	del temp_dict[ref_base]
	mutate_base = ""
	for base in temp_dict:
		if temp_dict[base] != 0:
			mutate_base += base

	mutate_ratio = (sum(temp_list[1:]) - max(temp_list[1:])) / sum(temp_list[1:])
	temp_list.append(mutate_ratio)
	temp_list = [str(t) for t in temp_list]
	temp_list.append(ref_base)
	temp_list.append(mutate_base)
	temp.append("\t".join(temp_list))

out = open(in_file + ".info","w")
out.write("\n".join(temp))
out.close()
