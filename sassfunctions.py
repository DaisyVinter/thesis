# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 15:54:30 2021

@author: daisy
"""
def get_positions(positions_file):

    """
    ======
    inputs
    ======

    positions_file: posiiton file from imaris containing the x, y and z coordinates of an imaris object

    ======
    output
    ======

    pos_dic: dictionary of positions
    
    """

    pos_dic = {}

    with open(positions_file) as o_pos:

        data_found = False

        for line in o_pos:

            split_line = line.strip().split(',')

            if split_line[0].startswith('Position X') and not data_found:

                data_found = True

                column_index = build_index_dic(split_line)

                continue

            elif not data_found:

                continue

            pos_id = split_line[column_index['ID']]
            pos_time = float(split_line[column_index['Time']])
            pos_x = float(split_line[column_index['Position X']])
            pos_y = float(split_line[column_index['Position Y']])
            pos_z = float(split_line[column_index['Position Z']])

            try:

                pos_dic[pos_time][pos_id] = (pos_x, pos_y, pos_z)

            except KeyError:

                pos_dic[pos_time] = {pos_id: (pos_x, pos_y, pos_z)}
  
    return pos_dic
#%%
def get_imaris_files(folder_path, folder_list, file_list):
    """
    Searches through the embyro folder (folder_path) and gets all the relevant files to run the data analysis.

    Folder names are in the folder_list varaible and the files will be passed as the file_list variable.

    Data is returned as a dictionary with the top level folders as keys containing the te relevant internal files.
    """

    import os

    file_dic = {}

    spots_added = False

    for folder in os.listdir(folder_path):
        for x in folder_list:
            if x in folder.lower():
                if 'nuc' not in x:
                    spots_added = True

                file_dic[x] = {'path': folder, 'files': []}
                files_in_folder = os.listdir(os.path.join(folder_path, folder))

                for file in files_in_folder:
                    for y in file_list:
                        if y.lower() in file.lower():
                            file_dic[x]['files'].append(file)

    if not spots_added:
        print('no spots', spots_added)
        for folder in os.listdir(folder_path):
            if 'spot' in folder.lower():
                file_dic['spot'] = {'path': folder, 'files': []}

                files_in_folder = os.listdir(os.path.join(folder_path, folder))

                for file in files_in_folder:
                    
                    if 'spot' in file.lower():
                        file_dic['spot']['files'].append(file)

    # pp(file_dic)
    # exit()
    return file_dic
#%%

def import_axis(elipsoid_files, folder_path):
    import os
    files_to_get = ['A', 'B', 'C']

    axis_dic = {}

    for x in elipsoid_files:

        for filex in files_to_get:

            if 'Ellipsoid_Axis_'+filex+'.' in x:
        
                with open(os.path.join(folder_path, x)) as openX:

                    data_found = False

                    for line in openX:

                        split_line = line.strip().split(',')

                        if split_line[0].startswith('Ellipsoid Axis '+filex+' X') and not data_found:

                            data_found = True

                            column_index = build_index_dic(split_line)
                            # pp(column_index)
                            continue

                        elif not data_found:

                            continue

                        nuc_id = split_line[column_index['ID']]

                        nuc_time = float(split_line[column_index['Time']])

                        pos_x = float(split_line[column_index['Ellipsoid Axis '+filex+' X']])
                        pos_y = float(split_line[column_index['Ellipsoid Axis '+filex+' Y']])
                        pos_z = float(split_line[column_index['Ellipsoid Axis '+filex+' Z']])

                        try:

                            axis_dic[filex][nuc_time][nuc_id] = (pos_x, pos_y, pos_z)

                        except KeyError:

                            try:

                                axis_dic[filex][nuc_time] = {nuc_id: (pos_x, pos_y, pos_z)}

                            except KeyError:

                                axis_dic[filex] = {nuc_time: {nuc_id: (pos_x, pos_y, pos_z)}}

    return axis_dic
#%%
def build_index_dic(ls):

    dic = {}

    for x in range(0, len(ls)):

        dic[ls[x]] = x

    return dic
#%%
def assign_spots_to_nuc_fast_centroid(nuc_dic, spot_dic):
    import scipy.spatial
    import numpy as np

    '''
    Inputs:
    nuc_dic : positions of nuclei from position dictionary (items in dic)
    spot_dic : positions of spots from spot dictionary (items in dic)
    nuc_id : nucleus ID from position dictionary (key in dic)
    spot_id : spot ID from position dictionary (key in dic)
    
    Assigns each spot to closest nucleus based on nucleus centroid
    
    Returns dictionary of each spots (key) closest nucleus (value)
    '''

    nuc_id = np.array(nuc_dic.index)
    spot_id = np.array(spot_dic.index)
    
    matrix_of_closest = scipy.spatial.distance.cdist(nuc_dic, spot_dic, metric = 'euclidean')
    index_of_closest = np.where(matrix_of_closest == np.amin(matrix_of_closest, axis = 0))
    spot_to_nuc_dic = {}
    counter = 0
    lenspot = len(spot_id)
    for spot in range(lenspot):                                         ####can I vectorise this?
        counter += 1
        print('assigning_spots:  ' + str(counter/lenspot*100)[:4], end = '\r')
        nuc = nuc_id[index_of_closest[0][spot]]
        spot = spot_id[index_of_closest[1][spot]]
        spot_to_nuc_dic[spot] = nuc 
        
    return spot_to_nuc_dic
#%%
##new function makes output csv using pandas
def make_output(spot_to_nuc_dic, spot_df, nuc_df):
    
    import pandas as pd
    
    spot_df.rename(columns = {0: 'x', 1: 'y', 2: 'z'}, inplace = True)

    spot_to_nuc = pd.DataFrame(spot_to_nuc_dic, index = ['nuc']).T

    spot_and_nuc = spot_df.join(spot_to_nuc)
    
    num_per_nuclei = spot_and_nuc['nuc'].value_counts() 

    no_spots = {'nuc': [], 'num_spots': 0}
    
    nuc_id = nuc_df.index
    
    counter = 0
    lennuc = len(nuc_id)

    for nuc in nuc_id:

        counter += 1
        print('finding number of spots:  ' + str(counter/lennuc*100)[:4], end = '\r')
    
        try:

            spot_and_nuc.loc[spot_and_nuc['nuc'] == nuc, 'num_spots'] = num_per_nuclei.loc[nuc]

        except KeyError:

            no_spots['nuc'].append(nuc)
            
    with_no_spot_nucs = pd.concat([spot_and_nuc.reset_index(), pd.DataFrame(no_spots)], ignore_index = True).rename(columns = {'index': 'spot', 0: 'spot_x', 1:'spot_y', 2: 'spot_z'})
    
    counter1 = 0 
    
    for i in nuc_id:                            #can this be vectorised? slowest part
        counter1 += 1
        print('finding nuclear position:  ' + str(counter1/lennuc*100)[:4], end = '\r')
        
        with_no_spot_nucs.loc[with_no_spot_nucs['nuc'] == i, 'nuc_x'] =  nuc_df.loc[i][0]   
        with_no_spot_nucs.loc[with_no_spot_nucs['nuc'] == i, 'nuc_y'] =  nuc_df.loc[i][1]   
        with_no_spot_nucs.loc[with_no_spot_nucs['nuc'] == i, 'nuc_z'] =  nuc_df.loc[i][2]   
        
    print(with_no_spot_nucs.head(10))
    print(with_no_spot_nucs.tail(10))
    return with_no_spot_nucs
#%%

def assign_nuc_axis_dic(axis_dic, nuc_dic):

    """
    Uses the axis of the nucleus to generate many points along the axis for each nucleus.

    """

    pos_dic = {}
    # tup_ls = []

    for t in axis_dic['C']:

        for nuc in axis_dic['C'][t]:

            nuc_pos = nuc_dic[t][nuc]

            axis_info = axis_dic['C'][t][nuc]

            fx_ls, fy_ls, fz_ls = make_further_ls(0, 25, 40, nuc_pos[0], nuc_pos[1], nuc_pos[2], axis_info)

            for i in range(0, len(fx_ls)):

                tup_key = (fx_ls[i], fy_ls[i], fz_ls[i])
                # tup_ls.append((fx_ls[i], fy_ls[i], fz_ls[i]))
                # print(type(tup_key))
                try:

                    pos_dic[t][tup_key].append(nuc)

                except KeyError:

                    try:

                        pos_dic[t][tup_key] = nuc

                    except KeyError:

                        pos_dic[t] = {tup_key: nuc}

    print("Generation finished")

    for x in pos_dic:

        print(x, len(pos_dic[x]))

    return pos_dic
#%%

def make_further_ls(start, stop, num_pos, x, y, z, axis_info):

    import numpy as np

    fx_ls = []
    fy_ls = []
    fz_ls = []

    pos_to_plot = np.linspace(start, stop, num_pos)

    for ni in pos_to_plot:

        fx_ls.append(x+axis_info[0]*ni)
        fy_ls.append(y+axis_info[1]*ni)
        fz_ls.append(z+axis_info[2]*ni)

    return fx_ls, fy_ls, fz_ls
#%%

def assign_spots_to_nuc_fast_axis(nuc_dic, spot_dic):
    
    import pandas as pd 
    import scipy.spatial.distance
    import numpy as np
    
    spot_df = pd.DataFrame(spot_dic).T
    nuc_axis_df = pd.DataFrame(nuc_dic.keys())
    spot_to_nuc_dic = {}
    
    try:
        
        dis_matrix = scipy.spatial.distance.cdist(np.array(nuc_axis_df), np.array(spot_df), metric = 'euclidean')   #uses too much memory  
    
        index_of_closest = np.where(dis_matrix == np.amin(dis_matrix, axis = 0))       
        
        spot_id = np.array(spot_df.index)
        nuc_id = np.array(list(nuc_dic.values()))
        
        counter = 0
        lenspot = len(spot_id)
        
        for spot in range(lenspot):
            counter += 1
            print('assigning_spots:  ' + str(counter/lenspot*100)[:4], end = '\r')
            nuc = nuc_id[index_of_closest[0][spot]]
            spot = spot_id[index_of_closest[1][spot]]
            spot_to_nuc_dic[spot] = nuc
    except MemoryError:
        
        print('Not enough memory: split into 4')
          
        spot_df_split = np.array_split(spot_df, 4)
        
        for i in range(4):
            
            dis_matrix = scipy.spatial.distance.cdist(nuc_axis_df, spot_df_split[i], metric = 'euclidean')
            index_of_closest = np.where(dis_matrix == np.amin(dis_matrix, axis = 0))   
            
            spot_id = np.array(spot_df_split[i].index)
            nuc_id = np.array(list(nuc_dic.values()))
            
            counter = 0
            lenspot = len(spot_id)
            
            for spot in range(lenspot):
                counter += 1
                print('block' + str(i) +': assigning_spots:  ' + str(counter/lenspot*100)[:4], end = '\r')
                nuc = nuc_id[index_of_closest[0][spot]]
                spot = spot_id[index_of_closest[1][spot]]
                spot_to_nuc_dic[spot] = nuc
            
    
            
    return spot_to_nuc_dic
           

#%%
def get_angles_to_axis(axis_dic, nuc_dic):

    import math
    import numpy as np
    angle_dic = {}
    
    axes = ['A', 'B', 'C']
    
    for ax in axes:
        for t in axis_dic[ax]:
            for nuc in axis_dic[ax][t]:
                if nuc not in angle_dic:
                    angle_dic[nuc] = {}
                axis_info = axis_dic[ax][t][nuc]
                angle_to_plane = np.arcsin(abs(axis_info[2])/(np.sqrt(np.power(axis_info[0],2) + np.power(axis_info[1],2) + np.power(axis_info[2],2))))
                angle_deg = angle_to_plane*180/math.pi
                angle_dic[nuc][angle_deg] = ax

    return angle_dic

#%%
def get_max_axis(nuc_axis_angle):

    max_axis = {}

    for nuc in nuc_axis_angle:
        angles = []
        for angle in nuc_axis_angle[nuc]:
            angles.append(angle)
        max_ang = max(angles)
        max_axis[nuc] = nuc_axis_angle[nuc][max_ang]
        
    return max_axis

#%%

def assign_nuc_axis_dic_new(axis_dic, nuc_dic):

    """
    Uses the axis of the nucleus to generate many points along the axis for each nucleus.
only gives A for some reason
    """
    pos_dic = {}
        # tup_ls = []
    
    #n = ['1157', '1158']
    axes = ['A', 'B', 'C']
    
    for ax in axes:
        
        for t in axis_dic[ax]:
            
            for nuc in axis_dic[ax][t]:
                           
            #for nuc in n:
    
                nuc_pos = nuc_dic[t][nuc]            ##gets nuclear position
                
                axis_info = axis_dic[ax][t][nuc]          ##gets vectors for all axes
    
                fx_ls, fy_ls, fz_ls = make_further_ls(-20, 20, 40, nuc_pos[0], nuc_pos[1], nuc_pos[2], axis_info)
    
                for i in range(0, len(fx_ls)):
    
                    tup_key = (fx_ls[i], fy_ls[i], fz_ls[i])   
    
                    try:
                        pos_dic[nuc][tup_key].append(ax)
    
                    except KeyError:
    
                        try:
    
                            pos_dic[nuc][tup_key] = [ax]              
    
                        except KeyError:
    
                            pos_dic[nuc] = {tup_key : [ax]}
    
    return pos_dic

#%%
def max_axis_points(pos_dic, max_axis):

    max_axis_pos = {}
       
    for nuc in pos_dic:
        for tup_key in pos_dic[nuc]:
            if max_axis[nuc] in pos_dic[nuc][tup_key]:
                max_axis_pos[tup_key] = [nuc]
                
    return max_axis_pos

