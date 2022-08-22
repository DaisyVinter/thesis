import pandas as pd
import os
import numpy as np
import argparse

test_directory = 'Work/test'

def get_dicts(directory):
    dir_list = os.listdir(directory)
    found_files = []
    for i in dir_list:
        if 'surf' in i:
            nuc_file_list = os.listdir(directory + '/' + i)
            for j in nuc_file_list:
                if 'Position' in j:
                    nuc_spots = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
                    found_files.append('nuc')
        if 'all' in i:
            all_file_list = os.listdir(directory + '/' + i)
            for j in all_file_list:
                if 'Position' in j:
                    all_spots = pd.read_csv(directory + '/' + i + '/' + j, header = 2)
                    found_files.append('all')
    
    if 'all' not in found_files:
        return('No total position file')
    if 'nuc' not in found_files:
        return('No nuclear position file')
    if len(found_files) != 2:
        return('Wrong number of position files')
    
    return nuc_spots, all_spots, all_file_list
     
     
def find_cyt_spots(nuc_spots, all_spots):
    
    cyt_spots = pd.concat([all_spots, nuc_spots]).drop_duplicates(keep = False, subset = ['Position X', 'Position Y', 'Position Z'])
    cyt_IDs = np.array(cyt_spots['ID'])
    all_IDs = np.array(all_spots['ID'])
    nuc_IDs = np.setdiff1d(all_IDs, cyt_IDs)
    
    return nuc_IDs, cyt_IDs

def output_cyt_files(directory, file_list, nuc_IDs, cyt_IDs):
	dir_list = os.listdir(directory)
	if 'cytoplasmic_spots' not in dir_list:
		os.mkdir(directory + '/cytoplasmic_spots')
	if 'n#clear_spots' not in dir_list:
		os.mkdir(directory + '/n#clear_spots')
	for i in dir_list:
		if 'all' in i:
			for file in all_file_list:
				try:
					stats = pd.read_csv(directory + '/' + i +  '/' + file, header = 2, index_col = 'ID')
					stats.drop(nuc_IDs, inplace = True)
					stats.drop(stats.filter(regex="Unname"), axis = 1, inplace = True)
					stats['ID'] = stats.index
					stats.to_csv(directory + '/cytoplasmic_spots/' + file, index = False)
				except KeyError:	
					print(file + ': No IDs')
					pass
				
			for file in all_file_list:
				try:
					stats = pd.read_csv(directory + '/' + i +  '/' + file, header = 2, index_col = 'ID')
					stats.drop(cyt_IDs, inplace = True)
					stats.drop(stats.filter(regex="Unname"), axis = 1, inplace = True)
					stats['ID'] = stats.index
					stats.to_csv(directory + '/n#clear_spots/' + file, index = False)
				except KeyError:
					print(file + ': No IDs')
					pass
                 
                    

parser = argparse.ArgumentParser()
parser.add_argument('directory', help = 'folder with Imaris files in')
args = parser.parse_args()

nuc_spots, all_spots, all_file_list = get_dicts(args.directory)
nuc_IDs, cyt_IDs = find_cyt_spots(nuc_spots, all_spots)
output_cyt_files(args.directory, all_file_list, nuc_IDs, cyt_IDs)

