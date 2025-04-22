import os
import re
import sys

i = 0
b = 0
c = []
c1 = []
d = {}
f = {}
d1 = {}
d2 = {}
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
