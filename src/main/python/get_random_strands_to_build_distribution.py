import os
from random import random
dir = 'my_dir'

for name in os.listdir(dir):
    count = 0
    ifile = open(os.path.join(dir,name))
    line = '\n'
    print name
    while line != '':
        count += 1 # one-based line numbers
        line = ifile.readline()
        if random() <= 0.000016442722119879474:
            ofile = open(name.replace('.fa','') + '_' + str(count) + '.fa','w')
            ofile.write(line)
            ofile.close()
    ifile.close()
