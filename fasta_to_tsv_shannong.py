import sys
in_file = sys.argv[1]

out_temp = []
for line in open(in_file):
	if not line.startswith(">"):
		temp_list  = [s for s in line.strip()]
		out_temp.append("\t".join(temp_list))
		
out = open(in_file + ".seq","w")
out.write("\n".join(out_temp))
out.close()
