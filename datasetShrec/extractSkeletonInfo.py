import os
from os.path import isdir, join
import glob
from shutil import copyfile

pathData='C:\\Users\\artix\\Downloads\\HandGestureDataset_SHREC2017\\HandGestureDataset_SHREC2017\\'
curDir='C:\\Users\\artix\\Desktop\\datasetFrancese';

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result

ff=find_all('skeletons_world.txt',pathData)
for f in ff:
    gest=f.split("HandGestureDataset_SHREC2017\\HandGestureDataset_SHREC2017")[1]
    files=curDir+gest
    print(curDir+gest)
    os.makedirs(os.path.dirname(files), exist_ok=True)
    copyfile(f, files)
