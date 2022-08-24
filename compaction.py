import pandas as pd
import numpy as np
import math
import seaborn as sns
import os
import argparse

def get_dicts(directory):
	current_folder = directory
	found_files = []
	dir_list = os.listdir(directory)
	for i in dir_list:
		if 'suntag' in i:
			st_file_list = os.listdir(directory + '/' + i)
			for j in st_file_list:
				if 'Position' in j:
					st_pos_df = pd.read_csv(directory + '/'  + i + '/' + j, header = 2)
					found_files.append('suntag')
		if 'hb' in i:
			hb_file_list = os.listdir(directory + '/' + i)
			for j in hb_file_list:
				if 'Position' in j:
					hb_pos_df = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
					found_files.append('hb')
		if 'ng' in i:
			ng_file_list = os.listdir(directory + '/' + i)
			for j in ng_file_list:
				if 'Position' in j:
					ng_pos_df = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
					found_files.append('ng_pos')
				if 'Intensity_Sum_Ch'in j:
					if 'Ch=1' in j:
						ng_intensity_df_1 = pd.read_csv(directory + '/' +i + '/' + j, header = 2)
					if 'Ch=2' in j:
						ng_intensity_df_2 = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
					if 'Ch=3' in j:
						ng_intensity_df_3 = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
					if 'Ch=4' in j:
						ng_intensity_df_4 = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
					found_files.append('ng_int')

	if 'suntag' not in found_files:
		print ('Error: No Suntag File')
	if 'hb' not in found_files:
		print ('Error: No Hunchback File')
	if 'ng_pos' not in found_files:
		print ('Error: No NG Position File')
	if 'ng_int' not in found_files:
		print ('Error: No NG Intensity File')
	
	st_pos_dict = {}
	for r in st_pos_df.index:
		st_pos_dict[st_pos_df.loc[r, 'ID']] = [st_pos_df.loc[r, 'Position X'], 
		st_pos_df.loc[r, 'Position Y'], st_pos_df.loc[r, 'Position Z']]
		
	hb_pos_dict = {}
	for r in hb_pos_df.index:
		hb_pos_dict[hb_pos_df.loc[r, 'ID']] = [hb_pos_df.loc[r, 'Position X'], 
		hb_pos_df.loc[r, 'Position Y'], hb_pos_df.loc[r, 'Position Z']]
			
	ng_pos_dict = {}
	for r in ng_pos_df.index:
		ng_pos_dict[ng_pos_df.loc[r, 'ID']] = [ng_pos_df.loc[r, 'Position X'], 
		ng_pos_df.loc[r, 'Position Y'], ng_pos_df.loc[r, 'Position Z']]
	
	ng_intensity_dict = {}
	for r in ng_intensity_df_1.index:
		ng_intensity_dict[ng_intensity_df_1.loc[r, 'ID']] = [ng_intensity_df_1.loc[r, 'Intensity Sum']]
	for r in ng_intensity_df_2.index:
		for j in ng_intensity_dict.keys():
			if ng_intensity_df_2.loc[r, 'ID'] == j:
				ng_intensity_dict[j].append(ng_intensity_df_2.loc[r, 'Intensity Sum'])
	for r in ng_intensity_df_3.index:
		for j in ng_intensity_dict.keys():
			if ng_intensity_df_3.loc[r, 'ID'] == j:
				ng_intensity_dict[j].append(ng_intensity_df_3.loc[r, 'Intensity Sum'])
	for r in ng_intensity_df_4.index:
		for j in ng_intensity_dict.keys():
			if ng_intensity_df_4.loc[r, 'ID'] == j:
				ng_intensity_dict[j].append(ng_intensity_df_4.loc[r, 'Intensity Sum'])
		
			
	return st_pos_dict, hb_pos_dict, ng_pos_dict, ng_intensity_dict, current_folder
	
		
	

def get_dist(st_pos_dict, hb_pos_dict, ng_pos_dict, ng_intensity_dict):
	less_than_um = []
	closest_spot = {}

	for i in st_pos_dict:
		x = st_pos_dict[i][0]
		y = st_pos_dict[i][1]
		z = st_pos_dict[i][2]
		distance_dict = {}
		for j in hb_pos_dict:
			x1 = hb_pos_dict[j][0]
			y1 = hb_pos_dict[j][1]
			z1 = hb_pos_dict[j][2]
			distance = math.sqrt((x - x1)**2 + (y - y1)**2 + (z - z1)**2)
			if distance < 1:
				less_than_um.append([i, j, distance])
			distance_dict[distance] = j
			min_distance = min(distance_dict.keys())
			closest_spot[i] = [distance_dict[min_distance], min_distance]
	
	less_than_300 = []
	coloc_nG = {}

	for i in st_pos_dict:
		x = st_pos_dict[i][0]
		y = st_pos_dict[i][1]
		z = st_pos_dict[i][2]
		distance_dict = {}
		for j in ng_pos_dict:
			x1 = ng_pos_dict[j][0]
			y1 = ng_pos_dict[j][1]
			z1 = ng_pos_dict[j][2]
			distance = math.sqrt((x - x1)**2 + (y - y1)**2 + (z - z1)**2)
			if distance < 0.22:
				less_than_220.append([i, j, distance])
				coloc_nG[i] = [j, distance]
				
	all_data = []
	for i in closest_spot:
		if i in coloc_nG.keys():
			all_data.append([i, closest_spot[i][0], closest_spot[i][1], 'coloc', coloc_nG[i][0], coloc_nG[i][1]])
		else:
			all_data.append([i, closest_spot[i][0], closest_spot[i][1], 'non_coloc'])
	
	for i in ng_intensity_dict:
		for j in all_data:
			if len(j) == 6:
				if i == j[4]:
					j.append(ng_intensity_dict[i][0])
					j.append(ng_intensity_dict[i][1])
					j.append(ng_intensity_dict[i][2])
					j.append(ng_intensity_dict[i][3])	
	return all_data
			
def make_file(all_data):
	df = pd.DataFrame(all_data, columns = ['suntag_id', 'closest_hb_spot_ID', 'closest_hb_spot_dis', 'coloc_status', 'coloc_nG_ID', 'coloc_ng_dis', 'coloc_ng_intens_Ch=1', 'coloc_ng_intens_Ch=2', 'coloc_ng_intens_Ch=3', 'coloc_ng_intens_Ch=4'])
	return df
	
	
parser = argparse.ArgumentParser()
parser.add_argument('directory', help = 'folder with Imaris files in')
args = parser.parse_args()

st_pos_dict, hb_pos_dict, ng_pos_dict, ng_intensity_dict, current_folder = get_dicts(args.directory)
all_data = get_dist(st_pos_dict, hb_pos_dict, ng_pos_dict, ng_intensity_dict)
df = make_file(all_data)
df.to_csv(current_folder + '/compaction_data.csv')

