#!/usr/bin/env python
# coding: utf-8

# In[51]:


# the way to use it: python overlay_pymol.py input
import os
import pandas as pd
import numpy as np
import pymol.cmd as cmd
import sys


# In[52]:


# get csv files to pair mol2 with pdb
def get_paired(loop_name):
    family_dir = 'align_result/'
    loop_dir = family_dir + loop_name
    csv_files = os.listdir(loop_dir)
    paired = []
    for i in range(len(csv_files)):
        sub_paired = []
        path = loop_dir + '/' + csv_files[i]
        data = pd.read_csv(path)
        mol2 = data.iloc[:,1]
        mol2 = mol2.values.tolist()
        pdb = data.iloc[:,2:]
        pdb = pdb.values.tolist()
        new_pdb = []
        for j in range(len(pdb)):
            a = ''
            for k in range(len(pdb[j])):
                a += str(pdb[j][k])
            new_pdb.append(a)
        sub_paired.append(mol2)
        sub_paired.append(new_pdb)
        paired.append(sub_paired)
    return paired


# In[53]:


def overlay(loop_name):
    paired = get_paired(loop_name)
    paired_ = pd.DataFrame(paired)
    paired_.to_csv(loop_name + 'paired.csv')
    mol2_dir = 'mimic_result/' + loop_name
    pdb_dir = 'loop_target/' + loop_name
    #__main__.pymol_argv = [ 'pymol']
    for i in range(len(paired)):
        length = int(paired[i][0][0][:2].replace('_',''))
        loop_dir = pdb_dir + '/' + loop_name + '_' + str(length)
        conformer_dir = mol2_dir + '/' + str(length)
        loop_list = os.listdir(loop_dir)
        for j in range(len(loop_list)):
            loop = loop_list[j]
            mol2_list = []
            for k in range(len(paired[i][1])):
                paired[i][1][k] = paired[i][1][k].replace('nan','')
                if paired[i][1][k] == loop:
                    mol2_list.append(k)
            loop_path = loop_dir + '/' + loop
            print(loop_path)
            cmd.load(loop_path, loop)
            group_object = ''
            for l in range(len(mol2_list)):
                mol2_path = conformer_dir + '/' + paired[i][0][mol2_list[l]]
                print(mol2_path)
                cmd.load(mol2_path, paired[i][0][mol2_list[l]])
                cmd.align(paired[i][0][mol2_list[l]],loop)
                #print(rmsd)
                group_object += paired[i][0][mol2_list[l]] + ' '
            print(group_object)
            print(loop + ' ' + group_object)
            cmd.group(loop + '_', loop + ' ' + group_object)


# In[54]:


if __name__== "__main__":
    overlay(sys.argv[1])
    cmd.save(sys.argv[1] + '.pse')
    cmd.quit()


