import os
import re
import sys

files = sys.argv[1]

a=0
c=0
diff=0
with open(files, 'r') as file:
    for line in file:
        line1 = line.split("\n")
        line2 = line1[0].split("\t")
        a=a+int(line2[3])
        c=c+1
f=3*(a/c)
d=0
with open(files, 'r') as file_object:
    for line in file_object:
        line1 = line.split("\n")
        line2 = line1[0].split("\t")
        if int(line2[3]) > f:
            if d==0:
                name=line2[0]
                start=line2[1]
                end=line2[2]
                d=d+1
            else:
                if name==line2[0]:
                    if 500000 >= int(line2[1])-int(end):
                        end=line2[2]
                    else:
                        if diff < int(end)-int(start):
                            diff=int(end)-int(start)
                            maxs=start
                            maxe=end
                        #print(name+"\t"+start+"\t"+end)
                        start=line2[1]
                        end=line2[2]
                else:
                    if diff < int(end) - int(start):
                        diff = int(end) - int(start)
                        maxs = start
                        maxe = end
                    print(name + "\t" + maxs + "\t" + maxe + "\t" + str(diff))
                    diff=0
                    name = line2[0]
                    start = line2[1]
                    end = line2[2]
    if diff < int(end) - int(start):
        diff = int(end) - int(start)
        maxs = start
        maxe = end
    print(name + "\t" + maxs + "\t" + maxe + "\t" + str(diff))