

def add_spots_to_name(directory):
	import os
	file_list = os.listdir(directory)
	for file in file_list:
		os.rename(directory + '/' + file, directory + '/spots_' + file) 
	os.rename(directory, directory + ' spots')

import argparse
parser = argparse.ArgumentParser(description = 'Add spots to files')
parser.add_argument('directory', help = 'folder with Imaris files in')
args = parser.parse_args()

add_spots_to_name(args.directory)