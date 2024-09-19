#!/usr/bin/python3

import sys

strains_to_cov = {}
		
#report_file = "~{straingst_report}"
report_file = sys.argv[1]
coverage_cutoff = float(sys.argv[2])
outfile = "STRAIN_REF"
cov_file = "STRAIN_COV"

with open(report_file, 'r') as infile:
	next(infile) #skip headers 
	for line in infile:
		l=line.split('\t')
		strains_to_cov[l[1]] = l[5]


output_strain = max(strains_to_cov, key = strains_to_cov.get)
out_cov = strains_to_cov[max(strains_to_cov, key = strains_to_cov.get)]

# if float(out_cov) > 0.8:
# 	print(output_strain)
# else:
# 	print("Insufficient_COV")

with open (outfile, 'w') as f:
	if float(out_cov) > coverage_cutoff:
		f.write(output_strain)
	else:
		f.write("Insufficient_COV")
