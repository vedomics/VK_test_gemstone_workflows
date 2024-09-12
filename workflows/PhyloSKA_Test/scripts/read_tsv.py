#!/usr/bin/python3

import sys

strains_to_cov = {}
		
#report_file = "~{straingst_report}"
report_file = sys.argv[1]

with open(report_file) as infile:
	next(infile) #skip headers 
	for line in infile:
		l=line.split('\t')
		strains_to_cov[l[1]] = l[5]

output_strain = max(strains_to_cov, key = strains_to_cov.get)
out_cov = strains_to_cov[max(strains_to_cov, key = strains_to_cov.get)]

with open("strain.txt", "w") as f:
	f.write(output_strain)

with open("cov.txt", "w") as f:
	f.write(out_cov)
		