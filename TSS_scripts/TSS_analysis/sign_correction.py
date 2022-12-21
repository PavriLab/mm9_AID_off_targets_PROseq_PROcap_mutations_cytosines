#!/usr/bin/env python3
'''
Description:
    Python script that corrects "multisignals", signals that expand to more than
    one position, in bedgraph files. The second argument is the input bedgraph
    and the third is the output file. Can be used with other files types as well.

Usage: ./sign_correction.py input output

'''
import sys

if len(sys.argv)!=3:
    print('Invalid number of arguments.\nUsage: ./sign_correction.py input output')
    sys.exit()



input_bg=sys.argv[1]
output_bg=sys.argv[2]

with open(input_bg,'r') as bg, \
open(output_bg, 'w') as out:
    size=0
    for i in bg:
        i=i.strip()
        i=i.split('\t')
        i[2]=int(i[2])
        i[1]=int(i[1])
        size=i[2]-i[1]
        if size==1:
            print(*i, sep='\t', file=out)
        elif size>1:
            for j in range(0,size):
                print(i[0],i[1]+j,i[1]+j+1,i[3], sep='\t', file=out)
