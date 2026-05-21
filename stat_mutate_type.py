#Note: This script is not suitable for continuous base mutations. If continuous mutations are required, manual verification is required
import sys
in_file = sys.argv[1]

codon = { "GCT":"A",\
"GCC":"A",\
"GCA":"A",\
"GCG":"A",\
"TGT":"C",\
"TGC":"C",\
"GAT":"D",\
"GAC":"D",\
"GAA":"E",\
"GAG":"E",\
"TTT":"F",\
"TTC":"F",\
"GGT":"G",\
"GGC":"G",\
"GGA":"G",\
"GGG":"G",\
"CAT":"H",\
"CAC":"H",\
"ATT":"I",\
"ATC":"I",\
"ATA":"I",\
"AAA":"K",\
"AAG":"K",\
"TTA":"L",\
"TTG":"L",\
"CTT":"L",\
"CTC":"L",\
"CTA":"L",\
"CTG":"L",\
"ATG":"M",\
"AAT":"N",\
"AAC":"N",\
"CCT":"P",\
"CCC":"P",\
"CCA":"P",\
"CCG":"P",\
"CAA":"Q",\
"CAG":"Q",\
"CGT":"R",\
"CGC":"R",\
"CGA":"R",\
"CGG":"R",\
"AGA":"R",\
"AGG":"R",\
"TCT":"S",\
"TCC":"S",\
"TCA":"S",\
"TCG":"S",\
"AGT":"S",\
"AGC":"S",\
"ACT":"T",\
"ACC":"T",\
"ACA":"T",\
"ACG":"T",\
"GTT":"V",\
"GTC":"V",\
"GTA":"V",\
"GTG":"V",\
"TGG":"W",\
"TAT":"Y",\
"TAC":"Y"}

#type A
#orf1_region = [154,1134]	#981
#orf2_region = [1352,1996]	#645
#type B
orf1_region = [154,1134]	#981
orf2_region = [1380,2024]	#645

temp_dict = {}
data = open(in_file).readlines()[1:]
for line in data:
	spa = line.strip("\n").split("\t")
	temp_dict[int(spa[0])] = spa

def rev(seqs):
	seqs_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
	rev_seqs = ""
	rev_seqs += seqs_dict[ seqs[-1] ]
	rev_seqs += seqs_dict[ seqs[-2] ]
	rev_seqs += seqs_dict[ seqs[-3] ]
	return rev_seqs

temp = ["position\tA\tT	G\tC\t-\tmutate_ratio\tref_base\tmutate_base\ttype"]
for base in temp_dict:
	mutate_type = "black"
	if float(temp_dict[base][6]) != 0:
		if int(temp_dict[base][0]) >= orf1_region[0] and int(temp_dict[base][0]) <= orf1_region[1]:
			if (int(temp_dict[base][0]) - orf1_region[0]) % 3 == 2:
				base1 = "".join(temp_dict[base-2][7:])
				base2 = "".join(temp_dict[base-1][7:])
				b3_ref = "".join(temp_dict[base][7])
				b3_mutate = "".join(temp_dict[base][8])
				for b1 in base1:
					for b2 in base2:
						for b3 in b3_mutate:
							if codon[b1 + b2 + b3_ref] != codon[b1 + b2 + b3]:
								mutate_type = "red"
				
			if (int(temp_dict[base][0]) - orf1_region[0]) % 3 == 0:
				b1_ref = "".join(temp_dict[base][7])
				b1_mutate = "".join(temp_dict[base][8])
				base2 = "".join(temp_dict[base+1][7:])
				base3 = "".join(temp_dict[base+2][7:])
				for b2 in base2:
					for b3 in base3:
						for b1 in b1_mutate:
							if codon[b1_ref + b2 + b3] != codon[b1 + b2 + b3]:
								mutate_type = "red"
							
			if (int(temp_dict[base][0]) - orf1_region[0]) % 3 == 1:
				base1 = "".join(temp_dict[base-1][7:])
				b2_ref = "".join(temp_dict[base][7])
				b2_mutate = "".join(temp_dict[base][8])
				base3 = "".join(temp_dict[base+1][7:])
				for b1 in base1:
					for b3 in base3:
						for b2 in b2_mutate:
							if codon[b1 + b2_ref + b3] != codon[b1 + b2 + b3]:
								mutate_type = "red"
			temp.append("\t".join(temp_dict[base]) + "\t" +mutate_type )
		elif int(temp_dict[base][0]) >= orf2_region[0] and int(temp_dict[base][0]) <= orf2_region[1]:
			if (orf2_region[1] - int(temp_dict[base][0])) % 3 == 2:
				b1_ref = "".join(temp_dict[base][7])
				b1_mutate = "".join(temp_dict[base][8])
				base2 = "".join(temp_dict[base+1][7:])
				base3 = "".join(temp_dict[base+2][7:])
				for b2 in base2:
					for b3 in base3:
						for b1 in b1_mutate:
							
							if codon[rev(b1_ref + b2 + b3)] != codon[rev(b1 + b2 + b3)]:
								mutate_type = "red"

			if (orf2_region[1] - int(temp_dict[base][0])) % 3 == 1:
				base1 = "".join(temp_dict[base-1][7:])
				b2_ref = "".join(temp_dict[base][7])
				b2_mutate = "".join(temp_dict[base][8])
				base3 = "".join(temp_dict[base+1][7:])
				for b1 in base1:
					for b3 in base3:
						for b2 in b2_mutate:
							if codon[rev(b1 + b2_ref + b3)] != codon[rev(b1 + b2 + b3)]:
								mutate_type = "red"
							
			if (orf2_region[1] - int(temp_dict[base][0])) % 3 ==  0:
				base1 = "".join(temp_dict[base-2][7:])
				base2 = "".join(temp_dict[base-1][7:])
				b3_ref = "".join(temp_dict[base][7])
				b3_mutate = "".join(temp_dict[base][8])
				for b1 in base1:
					for b2 in base2:
						for b3 in b3_mutate:
							if codon[rev(b1_ref + b2 + b3)] != codon[rev(b1 + b2 + b3)]:
								mutate_type = "red"
			temp.append("\t".join(temp_dict[base]) + "\t" + mutate_type)
							
		else:
			temp.append("\t".join(temp_dict[base]) + "\t" + mutate_type)
	else:
		temp.append("\t".join(temp_dict[base]) + "\t" + mutate_type)
			
		
out = open(in_file + ".type","w")
out.write("\n".join(temp))
out.close()