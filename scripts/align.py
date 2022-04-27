#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pymol
from pymol import cmd
import numpy as np
import __main__
import pandas as pd
import sys


# In[ ]:


def modify_pdb(path, length):
    name = path.replace('.pdb','')
    path = 'loop_target/' + path
    length = length
    dir_path = 'loop_target/' + name
    if os.path.isdir(dir_path) == False:
        os.mkdir(dir_path)
    dir_path = dir_path + '/' + name + '_' + str(length)
    if os.path.isdir(dir_path) == False:
        os.mkdir(dir_path)
    file = open(path, 'r')
    txt = file.readlines()
    res_start = int(txt[0][23:26])
    res_end = int(txt[-2][23:26])
    distance = res_end - res_start + 1
    change = distance - length + 1
    for i in range(change):
        start = res_start + i
        output_name = dir_path + '/' + str(start) + '_' + str((start+length)) + '.pdb'
        output = open(output_name, 'w')
        for j in range(len(txt)):
            number = txt[j][23:26].replace(' ','')
            if number.isnumeric():
                res_id = int(number)
                if txt[j][:4] == 'ATOM' and res_id >= start and res_id <= (start+length-1):
                    id = res_id - start + 1
                    if id < 10:
                        data = txt[j][:17] + 'UNK' + txt[j][20:23] + ' ' + str(id) + txt[j][26:]
                        output.write(data)
                    else:
                        data = txt[j][:17] + 'UNK' + txt[j][20:23] + str(id) + txt[j][26:]
                        output.write(data)
        output.write('END')
        output.close()


# In[ ]:


def alignment(loop_name, n):
    n = int(n)
    for i in range(4,11):
        pdb_name = loop_name + '.pdb'
        modify_pdb(pdb_name, i)
    __main__.pymol_argv = [ 'pymol', '-qc'] 
    #pymol.finish_launching()
    dir_1 = os.listdir('conformers_for_matching')
    dir_1.sort()
    dir_path = 'align_result/' + loop_name
    if os.path.isdir(dir_path) == False:
        os.mkdir(dir_path)
    for i in range(len(dir_1)):
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
                    print(dir_3[k])
                    name_3 = name_2 + '/' + dir_3[k]
                    loop_list = os.listdir('loop_target/'+ loop_name + '/' +loop_path)
                    for s in range(len(loop_list)):
                        print(loop_list[s])
                        cmd.load('loop_target/' + loop_name + '/' +loop_path + '/' + loop_list[s],'mobile')
                        cmd.load(name_3,'target')
                        rmsd = cmd.align('mobile', 'target', cycles = 0)
                        complex_pair = []
                        complex_pair.append(dir_3[k])
                        complex_pair.extend(loop_list[s])
                        file_name.append(complex_pair)
                        RMSD.append(rmsd[0])
                        cmd.delete('target')
                        cmd.delete('mobile')
                    cmd.reinitialize()
                except:
                    print('{} is error'.format(dir_3[k]))
            index_min = np.argmin(RMSD)
            output_name.append(file_name[index_min])
            output_rmsd.append(RMSD[index_min])
        final_output_name = []
        final_output_rmsd = []
        for i in range(n):
            min_id = np.argmin(output_rmsd)
            final_output_name.append(output_name[min_id])
            final_output_rmsd.append(output_rmsd[min_id])
            output_name.remove(output_name[min_id])
            output_rmsd.remove(output_rmsd[min_id])
        df = pd.DataFrame(final_output_name, final_output_rmsd)
        df.to_csv(dir_path + '/' + folder + '.csv')


# In[ ]:


if __name__== "__main__":
    alignment(sys.argv[1],sys.argv[2])

