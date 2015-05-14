import math
import numpy as np
f = open("1A81H.pdb", "r")
w = open ("wrappers.txt", "r")

num = 0
res = []
a_prev = "blah"
b_prev = "blah"
for line in w:
    if (line.find("HB_") == 0):
        a = line[49:52].strip()
        b = line[72:75].strip()
        if(a.isdigit() and b.isdigit()):
            if(a != a_prev and b != b_prev):
                a_prev = a
                b_prev = b
                res.append(int(a))
                res.append(int(b))
print res
print len(res)
for i in xrange(0, len(res),2):
    print i
    for line in f:
        if (line.find("ATOM") == 0):
            line =  line.split()
            if int(line[5]) == res[i]:
                pass
                #get position
            elif int(line[5]) == res[i+1]:
                pass
                #get position
