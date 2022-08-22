# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:34:39 2021

@author: daisy
"""

import time

start_time = time.time()

folder_path = 'Z:/Daisy/2022/hbmut_experiment/2022-04-08/Series001/mrna/cyt'

import sassfunctions as sass
import os
import pandas as pd
import sys

folder_list = ['nucl', 'exonic', 'intronic']
file_list = ["_Position.", "_Intensity_", "Ellipsoid_Axis"]

files_dic = sass.get_imaris_files(folder_path, folder_list, file_list)

positions_dic = {}

for folder in files_dic:

        folder_pos = [x for x in files_dic[folder]['files'] if "_Position." in x]
        
        if len(folder_pos) != 1:

            sys.exit(f'wrong number of positions files detected: {folder} || {folder_pos}')

        positions_dic[folder] = sass.get_positions(os.path.join(folder_path, files_dic[folder]['path'], folder_pos[0]))
        
nuc_df = pd.DataFrame(positions_dic['nucl'][1.0]).T
spot_df = pd.DataFrame(positions_dic['spot'][1.0]).T

spot_to_nuc_dic = sass.assign_spots_to_nuc_fast_centroid(nuc_df, spot_df)
output = sass.make_output(spot_to_nuc_dic, spot_df, nuc_df)

#%%

out_dir = folder_path + '/data_output'

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
        
output.to_csv(out_dir + '/position_data.csv')

print("--- %s seconds ---" % (time.time() - start_time))
