#!/usr/bin/python
import numpy
import os
from scipy import stats

DIST_DIR = '/titan/cancerregulome9/workspaces/users/trobinso/kde2/results'
TARGET_DIR = '/titan/cancerregulome9/workspaces/users/trobinso/xargs2'
OUTFILE = 'scored_regions.out'
FILTER_P = 0.01
INFINITY = float('inf')

chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y']
datapoints = {}

for filename in os.listdir(DIST_DIR):
    if filename[0:3] == 'chr':
        id = filename.split('_')[0].replace('chr','')
        if id in chrs:
            print filename
            file = open(os.path.join(DIST_DIR,filename), 'r')
            line = '\n'
            while line != '':
                line = file.readline()
                split = line.split('\t')
                if len(split) > 7:
                    if split[7] not in datapoints:
                        datapoints[split[7]] = []
                    datapoints[split[7]].append(float(split[5]))

for key in datapoints:
    try:
        datapoints[key] = stats.kde.gaussian_kde([datapoints[key]])
        print datapoints[key].dataset
    except:
        print key,'-> None'
        datapoints[key] = None

outfile = open(os.path.join(TARGET_DIR,OUTFILE),'w')
for filename in os.listdir(TARGET_DIR):
    if filename[0:3] == 'chr':
        id = filename.split('_')[0].replace('chr','')
        print filename
        file = open(os.path.join(TARGET_DIR,filename), 'r')
        line = '\n'
        while line != '':
            line = file.readline()
            split = line.split('\t')
            
            if id in chrs and len(split) > 7:
                score = float(split[5])
                if datapoints[split[7]] != None: 
                    pvalue = (datapoints[split[7]]).integrate_box_1d(score,INFINITY)
                    if pvalue <= FILTER_P:
                        outfile.write(line.rstrip('\r\n') + '\t' + str(pvalue) + '\r\n')

outfile.close()
