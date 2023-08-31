#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import os
import numpy as np
import sys


# In[47]:


def copy_hits(loop_name):
    family_dir = 'align_result/'
    loop_dir = family_dir + loop_name
    target_folder = 'mimic_result/' + loop_name
    if os.path.isdir(target_folder) == False:
        os.mkdir(target_folder)
    csv_name = os.listdir(loop_dir)
    csv_name = [i for i in csv_name if 'combi' not in i]
    for i in range(len(csv_name)):
        folder_name = csv_name[i].replace('.csv','')
        store_folder = target_folder + '/' + folder_name
        if os.path.isdir(store_folder) == False:
            os.mkdir(store_folder)
        data = pd.read_csv(loop_dir + '/' + csv_name[i])
        x = data.iloc[:,1]
        rms = data.iloc[:,0]
        label = []
        for j in range(len(x)):
            split = []
            for k in range(len(x[j])):
                if x[j][k] == '_':
                    split.append(k)
            label.append(split)
        #print(len(label),label)
        for l in range(len(x)):
            folder_2 = str(x[l][:label[l][1]])
            folder_2 = folder_2.replace('_','R')
            #print(folder_2)
            file_name = str(x[l])
            target = store_folder + '/' + file_name
            cmd = 'cp ' + 'conformers_for_matching/' + folder_name + '/' + folder_2 + '/' + file_name + ' '+  target
            os.system(cmd)
            print(l,cmd)


# In[48]:


if __name__== "__main__":
    copy_hits(sys.argv[1])

