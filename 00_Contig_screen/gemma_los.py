import os
import re
from numpy import transpose
import pandas as pd
import numpy as np
import sys
i = 0
b = 0
c = []
c1 = []
#message1 = []
d = {}
f = {}
d1 = {}
d2 = {}
fi = 'GWAS.assoc_1.txt'
files = sys.argv[1]
a=0
c=0
with open(files) as file_object:
    for line in file_object:
        line1 = line.split("\n")
        line2 = line1[0].split("\t")
        if a>0:
            if line2[0]==config_name:
                d_value=d_value+int(line2[3]) - int(line2[2])
                if c == 0:
                    if int(line2[1])/2 < d_value:
                        print(config_name)
                        c=1
            else:
                config_name = line2[0]
                d_value = int(line2[3]) - int(line2[2])
                c=0
                if c == 0:
                    if int(line2[1]) / 2 < d_value:
                        print(config_name)
                        c=1
        else:
            a=a+1
            config_name=line2[0]
            d_value=int(line2[3])-int(line2[2])
            if int(line2[1])/2 < d_value:
                print(config_name)
                c=1

        # if '\tLachesis_group' in line:
        #     line1 = line.split("\n")
        #     line2 = line1[0].split("\t")
        #     line1 = line2[1].split(":")
        #     line1 = line1[0].split("achesis_group")
        #     #print(line1[1])
        #     message = line1[1]
        #     for i in range(11):
        #         message = message + '\t' + line2[i+1]
        #     with open(fi, 'a') as file_object:
        #         file_object.write(message + '\n')
        # else:
        #     with open(fi, 'a') as file_object:
        #         file_object.write(line)