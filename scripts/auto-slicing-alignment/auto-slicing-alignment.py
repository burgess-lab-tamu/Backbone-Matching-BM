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



#modify the loop PDB file to a unified format
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



# auto-slicing alignment
#loop name is the name of the loop PDB file. For example, if the loop structure file is uPA.pdb, then loop name is 'uPA'
# n is the number of hit cyclo-organopeptides to save'; for example, if n = 50, then top 50 cyclo-organopeptides are saved in the final output
def alignment(loop_name, n, k):
    n = int(n)
    k = int(k)
    w = k
    if k >= 10:
        for i in range(4,11):
            pdb_name = loop_name + '.pdb'
            modify_pdb(pdb_name, i)
    else:
        for i in range(4,k+1):
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
        if int(dir_1[i]) <= w:
            print(dir_1[i], w)
            name_1 = 'conformers_for_matching/' + dir_1[i]
            folder = dir_1[i]
            dir_2 = os.listdir(name_1)
            dir_2.sort()
            output_name = []
            output_rmsd = []
            loop_path = loop_name + '_' + str(folder)
            for j in range(len(dir_2)):
                #print(dir_2[j])
                name_2 = name_1 + '/' + dir_2[j]
                dir_3 = os.listdir(name_2)
                dir_3.sort()
                file_name = []
                RMSD = []
                for k in range(len(dir_3)):
                    try:
                        #print(dir_3[k])
                        name_3 = name_2 + '/' + dir_3[k]
                        loop_list = os.listdir('loop_target/'+ loop_name + '/' +loop_path)
                        for s in range(len(loop_list)):
                            #print(loop_list[s])
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


#copy best hits in a folder
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




#export pse files of loop fragments with their best overlayed cyclo-organopeptides
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
    cmd.save(loop_name + '.pse')



if __name__== "__main__":
    alignment(sys.argv[1],sys.argv[2], sys.argv[3])
    copy_hits(sys.argv[1])
    get_paired(sys.argv[1])
    overlay(sys.argv[1])

