import os
import matplotlib.pyplot as plt

qdict = {}
with open('test-mapq.txt') as f:
    lines = f.readlines()
    for x in lines:
        mapq = int(x.strip())
        if mapq in qdict.keys():
            qdict[mapq] = qdict[mapq]+1
        else: 
            qdict[mapq] = 0
    
    plt.hist(qdict)
    plt.show()