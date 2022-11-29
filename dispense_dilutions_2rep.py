import time
import sys
from method_labware import *
from pace_util import * #TODO: import everything explicitly, at least from pace_util

import csv
import numpy as np



matrix_csv = 'drug_plate.csv'

def read_matrix_from_csv(csvfile):
    matrix = np.genfromtxt(csvfile, delimiter=',', dtype='unicode')
    matrix = np.rot90(matrix, 1)
    return matrix

matrix_a = read_matrix_from_csv('drug_a_plate.csv')
matrix_b = read_matrix_from_csv('drug_b_plate.csv')

drug_dict = {}
well_idx = 0
with open('drug_source.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row[1]:
            drug_dict.update({row[1]:well_idx})
        well_idx += 1



source_dict_a = {'x': {'Vessel':(dilution_plate, [0,1,2,3,4,5,6,7]),
                     'Volume':20},
               '10': {'Vessel':(dilution_plate, [8]),
                     'Volume':20},
               '8': {'Vessel':(dilution_plate, [9]),
                     'Volume':20},
               '6': {'Vessel':(dilution_plate, [10]),
                     'Volume':20},
               '4': {'Vessel':(dilution_plate, [11]),
                     'Volume':20},
               '2': {'Vessel':(dilution_plate, [12]),
                     'Volume':20},
               '9': {'Vessel':(dilution_plate, [16]),
                     'Volume':20},        
               '7': {'Vessel':(dilution_plate, [17]),
                     'Volume':20},
               '5': {'Vessel':(dilution_plate, [18]),
                     'Volume':20}, 
               '3': {'Vessel':(dilution_plate, [19]),
                     'Volume':20}, 
               '1': {'Vessel':(dilution_plate, [20]),
                     'Volume':20}, 
                     
               }

source_dict_b = {str(11-i):{'Vessel': (dilution_plate, list(range(8*i+24, 8*i+29))), 'Volume':20} for i in range(1,10)}
source_dict_b.update({'1':{'Vessel':(dilution_plate, [93,94]), 'Volume':20}})
#source_dict_b = {k:{'Vessel':(dilution_plate, [v['Vessel'][1][0]+16]), 'Volume':v['Volume']} if k!='x' else v for k, v in source_dict_a.items()}

def build_dispense_lists(source_dict, matrix):
    dispense_dict = {}
    for reagent in source_dict:
        dispense_dict.update({reagent: []})
        plate_index = 0
        for y in range(len(matrix)):
            for x in range(len(matrix[y])):    
                if matrix[y][x] == reagent:
                    dispense_dict[reagent].append(plate_index)
                plate_index += 1
    return dispense_dict

dispense_dict_a = build_dispense_lists(source_dict_a, matrix_a)
dispense_dict_b = build_dispense_lists(source_dict_b, matrix_b)


def columnwise_dispense_from_list(target_plate, source_plate, source_dict, dispense_dict):
    source_pos_map = {source_dict[key]['Vessel'][1][idx]:key for key in source_dict for idx in range(len(source_dict[key]['Vessel'][1]))}
    source_pos_list = [None]*96
    for key in source_pos_map:
        source_pos_list.insert(key, source_pos_map[key])
        
    
    dispense_pos_map = {}
    for reagent in dispense_dict:
        for well in dispense_dict[reagent]:
            dispense_pos_map.update({well:reagent})
    
    dispense_pos_list = [None]*384
    for key in dispense_pos_map:
        dispense_pos_list[key] = dispense_pos_map[key]
    
    
    list_of_dispense_series = []
    list_of_source_cols = []
    counter = 0
    while counter < 384:
        accumulated = [] #Maybe tag each accumulated well as a tuple with 2nd value equals size of dispense sequence
        for col in range(len(dispense_pos_list)//16):
            for source_col in range(len(source_pos_list)//8): # Grab a column from dispense list
                skip_next_well = False
                well_number_range = range(col*16, col*16+16)  # and a column from source list,
                reagent_list = []                             # and compare reagents to find columnwise
                col_dispense_poss = []                        # dispense opportunities
                col_source_poss = []
                source_col_list = source_pos_list[source_col*8:source_col*8+8]
                if all([well==None for well in source_col_list]):
                    continue
                for well_idx in well_number_range:
                    if skip_next_well:
                        skip_next_well = False
                        continue
                    if dispense_pos_list[well_idx] in source_col_list and dispense_pos_list[well_idx] is not None:
                        reagent = dispense_pos_list[well_idx]
                        col_dispense_poss.append(well_idx)
                        reagent_list.append(reagent)
                        col_source_poss.append(source_col_list.index(reagent)+8*source_col)
                        source_col_list[:source_col_list.index(reagent)+1] = [None]*len(source_col_list[:source_col_list.index(reagent)+1])
                        skip_next_well = True
                if len(col_dispense_poss) > 0:
                    list_of_source_cols.append([col_source_poss, reagent_list])
                    counter += 1

    list_of_source_cols = sorted(list_of_source_cols, key = lambda x: len(x[0]), reverse = True)    
    move_on_flag = False
    while not move_on_flag:
        for source_col in list_of_source_cols:
            aspirate_poss_list = [(source_plate, well) for well in source_col[0]]
            print(aspirate_poss_list)
            active_col_dispense_lists = []
            for active_col in range(len(dispense_pos_list)//16):
                active_source_list = list(source_col[1])
                well_number_range = range(active_col*16, active_col*16+16)
                active_col_dispense_poss = []
                for well_idx in well_number_range:
                    reagent = dispense_pos_list[well_idx]
                    if skip_next_well:
                        skip_next_well = False
                        continue
                    if reagent in active_source_list and reagent is not None and well_idx not in accumulated:
                        accumulated.append(well_idx)
                        active_col_dispense_poss.append((well_idx, active_source_list.index(reagent)))
                        active_source_list[:active_source_list.index(reagent)+1] = [None]*len(active_source_list[:active_source_list.index(reagent)+1])
                        skip_next_well = True
                dispense_poss_list = [None]*8
                            
                for well in active_col_dispense_poss:
                    dispense_poss_list[well[1]] = (target_plate, well[0])
                if all([well == None for well in dispense_poss_list]):
                    move_on_flag = True
                    continue
                active_col_dispense_lists.append(dispense_poss_list)
                dispense_series = [aspirate_poss_list, active_col_dispense_lists]
                        
            source_lists = [series[0] for series in list_of_dispense_series]
            if dispense_series[0] not in source_lists and not all([not poss for poss in dispense_series[1]]):
                list_of_dispense_series.append(dispense_series)
            elif dispense_series[0] in source_lists and not all([not poss for poss in dispense_series[1]]):
                matching_source_list = source_lists.index(dispense_series[0])
                list_of_dispense_series[matching_source_list][1] = list_of_dispense_series[matching_source_list][1] + dispense_series[1]
                        
        counter += 1       
    print("accumulated")             
    print(len(accumulated))
    for l in range(len(list_of_dispense_series)):
        reduced_dispense_series = []
        for series in list_of_dispense_series[l][1]:
            if series not in reduced_dispense_series:
                reduced_dispense_series.append(series)
        list_of_dispense_series[l][1] = reduced_dispense_series

    
    return list_of_dispense_series

def tips_list_iterator(asp_list, tips_list):
    asp_list =  asp_list + [None]*(8-len(asp_list))
    pickup_list = []
    for i in range(8):
        if asp_list[i]:
            first_tip = [tip for tip in tips_list if tip][0]
            first_tip_index = first_tip[1]
            tips_list[first_tip_index] = None
            pickup_list.append(first_tip)
        else:
            pickup_list.append(None)
    return pickup_list, tips_list



#deck_res_source
def single_dispense_from_list(deck_res_source, deck_res_target, wells_list, tip_vol, dispense_vol):
    asp_poss = [deck_res_source] + [None]*7
    tip_vols = [tip_vol] + [None]*7
    aspirate(ham_int, asp_poss, tip_vols)
    current_vol = tip_vol
    for well in wells_list:
        disp_poss = [(deck_res_target, well)] + [None]*7
        dispense_vols = [dispense_vol] + [None]*7
        dispense(ham_int, disp_poss, dispense_vols)
        current_vol = current_vol - dispense_vol
        if current_vol < dispense_vol:
            replace_vol = tip_vol - current_vol
            replace_vols = [replace_vol] + [None]*7
            aspirate(ham_int, asp_poss, replace_vols)
            current_vol = replace_vol


def matrix_dispense(source_dict, matrix, target, tip_vol, tips_iter):
    #matrix = read_matrix_from_csv(matrix_csv)
    dispense_lists = build_dispense_lists(source_dict, matrix)
    for reagent in source_dict:
        dispense_list = dispense_lists[reagent]
        deck_res_source = source_dict[reagent]['Vessel']
        dispense_vol = source_dict[reagent]['Volume']
        tips = next(tips_iter)
        tip_pick_up(ham_int, tips)
        single_dispense_from_list(deck_res_source, target, dispense_list, tip_vol, dispense_vol)
        tip_eject(ham_int, tips)


if __name__ == '__main__':

    with HamiltonInterface(simulate=True) as ham_int: #TODO: simulate=False
        normal_logging(ham_int)
        
        initialize(ham_int)
        
        
        if 'make' in sys.argv:
            arg_idx = sys.argv.index('make')
            drug_a_name = sys.argv[arg_idx+1]
            drug_b_name = sys.argv[arg_idx+2]
            
            tip_pos_a = [(drug_source_tips, 0)] + [None]*7
            tip_pos_b = [(drug_source_tips, 1)] + [None]*7
            tip_pick_up(ham_int, tip_pos_a)
            
            drug_a_pos = [(drug_source_plate, drug_dict[drug_a_name])] + [None]*7
            dispense_a_pos = [(dilution_plate, 8)] + [None]*7
            drug_vol = [1] + [None]*7
            
            drug_source_class = 'LowVolumeFilter_Water_DispenseSurface_Empty'
            aspirate(ham_int, drug_a_pos, drug_vol, liquidClass = drug_source_class)
            dispense(ham_int, dispense_a_pos, drug_vol, liquidClass = drug_source_class)
            
            tip_eject(ham_int, tip_pos_a)
            
            tip_pick_up(ham_int, tip_pos_b)
            drug_b_pos = [(drug_source_plate, drug_dict[drug_b_name])] + [None]*7
            dispense_b_pos = [(dilution_plate, 24)] + [None]*7
            drug_vol = [1] + [None]*7
            
            aspirate(ham_int, drug_b_pos, drug_vol, liquidClass = drug_source_class)
            dispense(ham_int, dispense_b_pos, drug_vol, liquidClass = drug_source_class)
            
        
        tip_pos_a = [(dilution_tips, 0)] + [None]*7
        tip_pos_b = [(dilution_tips, 1)] + [None]*7
        
        
        dilution_series_a = [8, 16, 9, 17, 10, 18, 11, 19, 12, 20]
        dilution_series_b = [well_idx + 16 for well_idx in dilution_series_a]
        serial_dilutions_a = [[(dilution_plate, dilution_well)] + [None]*7 for dilution_well in dilution_series_a]
        serial_dilutions_b = [[(dilution_plate, dilution_well)] + [None]*7 for dilution_well in dilution_series_b]
        dilution_vols = [500] + [None]*7
        
        dilution_class = 'HighVolumeFilter_Water_DispenseSurface_Empty'
        
        serial_dilutions_list = (serial_dilutions_a, serial_dilutions_b)
        dilution_tips_list = [tip_pos_a, tip_pos_b]
        
        for tip_poss, serial_dilution_poss in zip(dilution_tips_list, serial_dilutions_list):
            tip_pick_up(ham_int, tip_poss)
            for dilution_num in range(len(serial_dilution_poss)-1):
                aspirate(ham_int, serial_dilution_poss[dilution_num], dilution_vols, liquidClass = dilution_class)
                ## add mixing here ##
                dispense(ham_int, serial_dilution_poss[dilution_num+1], dilution_vols, liquidClass = dilution_class)
            tip_eject(ham_int, tip_poss)
        
        list_of_dispense_series_a = columnwise_dispense_from_list(drug_a_plate, dilution_plate, source_dict_a, dispense_dict_a)
        list_of_dispense_series_b = columnwise_dispense_from_list(drug_b_plate, dilution_plate, source_dict_b, dispense_dict_b)

        
        def execute_columnwise_dispense(list_of_dispense_series, tips_poss_list, vol, max_vol):
            for series in list_of_dispense_series:
                source_list = series[0]
                source_list =  source_list + [None]*(8-len(source_list))
                target_lists = series[1]
                print(target_lists)
                vols_lists = [[vol if y else None for y in target_lists[i]] for i in range(len(target_lists))]
                int_vols_lists = [[vol if y else 0 for y in target_lists[i]] for i in range(len(target_lists))] 
                asp_vols = [min(max_vol, sum([vols[i] for vols in int_vols_lists])) if source_list[i] else None for i in range(8)]
                current_vols = list(asp_vols)
                tips, tips_poss_list = tips_list_iterator(source_list, tips_poss_list)
                tip_pick_up(ham_int, tips)
                aspirate(ham_int, source_list, asp_vols)
                for target_poss in target_lists:
                    vols_list = [vol if pos else None for pos in target_poss]
                    dispense(ham_int, target_poss, vols_list)
                    current_vols = [current_vols[i]-vols_list[i] if vols_list[i] else current_vols[i] for i in range(8)]
                    if any ([v < vol for v in current_vols if v is not None]):
                        reaspirate_vols = [max_vol - v if v is not None else None for v in current_vols]
                        aspirate(ham_int, source_list, reaspirate_vols)
                        current_vols = [current_vols[i] + reaspirate_vols[i] if current_vols[i] is not None else None for i in range(len(current_vols))]
                tip_eject(ham_int, tips)
            return tips_poss_list
        
        tips_poss_list = [(drug_a_tips, _) for _ in range(96)]
        tips_poss_list = execute_columnwise_dispense(list_of_dispense_series_a, tips_poss_list, 20, 300)
        tips_poss_list = execute_columnwise_dispense(list_of_dispense_series_b, tips_poss_list, 20, 300)
        
        
        tip_pick_up_96(ham_int, pinning_tips)
        
        for quadrant in range(4):
            aspirate_384_quadrant(ham_int, drug_a_plate, quadrant, 20)
            dispense_384_quadrant(ham_int, drug_b_plate, quadrant, 20)
        tip_eject_96(ham_int, pinning_tips)
        
        