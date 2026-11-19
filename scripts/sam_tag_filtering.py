#!/usr/bin/python

import sys
import re


print("new_filter_1")
sam_file = sys.argv[1]
filtered_sam_file = sys.argv[2]
with open(sam_file, 'r') as f_in:
    with open(filtered_sam_file, 'w') as f_out:
        for line in f_in:
            if line[0] != '@':
                line_list = line.split('\t')
                NM = int([x for x in line_list if x.startswith('NM')][0].split(':')[2])
                MD = [x for x in line_list if x.startswith('MD')][0].split(':')[2].rstrip()
                MD = re.split('\d+',MD)[1:-1]
                if len(re.findall('fwd',filtered_sam_file)) == 1:
                    if NM >= 2 and MD.count('T') <= 1 and ((MD.count('A') + MD.count('C'))/NM) >= 0.7:
                        f_out.write(line)
                if len(re.findall('rev',filtered_sam_file)) == 1:
                    if NM >= 2 and MD.count('A') <= 1 and ((MD.count('T') + MD.count('G'))/NM) >= 0.7:
                        f_out.write(line)
            else:
                f_out.write(line)