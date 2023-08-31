#!/usr/bin/env python
# coding: utf-8

# In[15]:


import os
import pymol
from pymol import cmd, stored
import numpy as np
import __main__
import pandas as pd
import sys
import multiprocessing as mp
from multiprocessing import process
from multiprocessing import pool


# In[16]:


def alignment(loop_name, n, length):
    length = str(length)
    __main__.pymol_argv = [ 'pymol', '-qc'] 
    #pymol.finish_launching()
    dir_1 = os.listdir('conformers_for_matching')
    #print(dir_1)
    dir_1.sort()
    dir_path = '2014_hotloop_align_result/' + length + '/'
    if os.path.isdir(dir_path) == False:
        os.mkdir(dir_path)
    for i in range(len(dir_1)):
        #print(i)
        #print(dir_1[i])
        if dir_1[i] == length:
            print(length)
            name_1 = 'conformers_for_matching/' + dir_1[i]
            folder = dir_1[i]
            dir_2 = os.listdir(name_1)
            dir_2.sort()
            output_name = []
            output_rmsd = []
            loop_path = loop_name + '_' + str(folder)
            for j in range(len(dir_2)):
                print(dir_2[j])
                name_2 = name_1 + '/' + dir_2[j]
                dir_3 = os.listdir(name_2)
                dir_3.sort()
                file_name = []
                RMSD = []
                for k in range(len(dir_3)):
                    try:
                        #print(dir_3[k])
                        name_3 = name_2 + '/' + dir_3[k]
                        #print(name_3)
                        #loop_list = os.listdir('2014_hotloop_screen/' + length + '/'+ loop_name)
                        #print(loop_list[s])
                        cmd.load('2014_hotloop_screen/' + length + '/'+ loop_name,'mobile')
                        #print('2014_hotloop_screen/' + length + '/'+ loop_name,'mobile')
                        cmd.load(name_3,'target')
                        rmsd = cmd.align('mobile & i. 1-' + length +  '& ! h.', 'target & i. 1-' + length + ' &! h.', cycles = 0)
                        #print('mobile & i. 1-' + length +  '& ! h.', 'target & i. 1-' + length + ' &! h.')
                        complex_pair = []
                        complex_pair.append(dir_3[k])
                        complex_pair.extend(loop_name.replace('.pdb',''))
                        file_name.append(complex_pair)
                        #print(rmsd[0])
                        RMSD.append(rmsd[0])
                        cmd.delete('target')
                        cmd.delete('mobile')
                        cmd.reinitialize()
                    except:
                        print('{} is error'.format(dir_3[k]))
                index_min = np.argmin(RMSD)
                output_name.append(file_name[index_min])
                output_rmsd.append(RMSD[index_min])
            #print(output_rmsd)
            final_output_name = []
            final_output_rmsd = []
            for i in range(n):
                min_id = np.argmin(output_rmsd)
                final_output_name.append(output_name[min_id])
                final_output_rmsd.append(output_rmsd[min_id])
                output_name.remove(output_name[min_id])
                output_rmsd.remove(output_rmsd[min_id])
            df = pd.DataFrame(final_output_name, final_output_rmsd)
            df.to_csv(dir_path  + loop_name.replace('.pdb','') + '.csv')


# In[17]:


def copy_hits(loop_name, length):
    length = str(length)
    family_dir = '2014_hotloop_align_result/' + length + '/'
    loop_dir = family_dir
    target_folder = 'mimic_result/' + loop_name.replace('.pdb','')
    if os.path.isdir(target_folder) == False:
        os.mkdir(target_folder)
    csv_name = os.listdir(loop_dir)
    for i in range(len(csv_name)):
        if csv_name[i] == loop_name.replace('.pdb','.csv'):
            folder_name = csv_name[i].replace('.csv','')[-1]
            store_folder = target_folder + '/' 
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
                target_name = file_name.replace('.mol2','')
                target = store_folder + '/' + target_name + '.mol2'
                cmd = 'cp ' + 'conformers_for_matching/' + folder_name + '/' + folder_2 + '/' + file_name + ' '+  target
                os.system(cmd)
                print(l,cmd)


# In[18]:


# get csv files to pair mol2 with pdb
def get_paired(loop_name, length):
    length = str(length)
    family_dir = '2014_hotloop_align_result/'+ length + '/'
    loop_dir = family_dir
    csv_files = os.listdir(loop_dir)
    paired = []
    #print(csv_files)
    for i in range(len(csv_files)):
        #print(csv_files[i],loop_name.replace('.pdb','.csv'))
        if csv_files[i] == loop_name.replace('.pdb','.csv'): 
            sub_paired = []
            path = loop_dir + csv_files[i]
            #print(path)
            data = pd.read_csv(path)
            #print(data)
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


# In[19]:


def overlay(loop_name, length):
    cmd.reinitialize()
    #print(length)
    paired = get_paired(loop_name,length)
    #print(paired)
    length = str(length)
    mol2_dir = 'mimic_result/' + loop_name.replace('.pdb','') + '/'
    pdb_dir = '2014_hotloop_screen/' + length + '/' + loop_name
    #__main__.pymol_argv = [ 'pymol']
    for i in range(len(paired)):
        length = int(paired[i][0][0][:2].replace('_',''))
        loop_dir = pdb_dir
        conformer_dir = mol2_dir
        #loop_list = os.listdir(loop_dir)
        loop = loop_dir
        mol2_list = []
        for k in range(len(paired[i][1])):
            if paired[i][1][k] == loop_name.replace('.pdb',''):
                mol2_list.append(k)
        loop_path = loop_dir
        #print(loop_path)
        cmd.load(loop_path, loop_name.replace('.pdb',''))
        group_object = ''
        for l in range(len(mol2_list)):
            mol2_path = conformer_dir + '/' + paired[i][0][mol2_list[l]]
            print(mol2_path)
            cmd.load(mol2_path, paired[i][0][mol2_list[l]])
            cmd.align(paired[i][0][mol2_list[l]],loop_name.replace('.pdb',''))
            #print(rmsd)
            group_object += paired[i][0][mol2_list[l]] + ' '
        print(group_object)
        print(loop + ' ' + group_object)
        cmd.group(loop + '_', loop + ' ' + group_object)
    #print(length)
    cmd.save('2014_hotloop_overlays/' + str(length) + '/' + loop_name.replace('.pdb','_overlays') + '.pse')
    #cmd.quit()


# In[25]:


def one_loop_scan(y):
    key_loops = os.listdir('2014_hotloop_screen/8/')
    loop_name = key_loops[y]
    alignment(loop_name,50,8)
    copy_hits(loop_name,8)
    overlay(loop_name,8)


# In[26]:


def screen_hotloop():
    key_loops = os.listdir('2014_hotloop_screen/8/')
    length = len(key_loops)
    times = int(length/1)
    for i in range(times):
        pool = mp.Pool(processes=1)
        x = pool.map(one_loop_scan, range(i*1,(i+1)*1))
        pool.close()


# In[27]:


if __name__== "__main__":
    screen_hotloop()

