# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 11:48:15 2022

@author: stefa
"""
import sys
import csv
import numpy as np
import math
import os
from pyhamilton import (HamiltonInterface,  LayoutManager, 
 Plate96, Tip96, initialize, tip_pick_up, tip_eject, ResourceType, Plate384,
 aspirate, dispense,  oemerr, resource_list_with_prefix, normal_logging)

 
def extend_2d_list(lst2d, *extension_params):
    extension_params = (extension_params)
    converted_list = []
    for i in range(len(lst2d)):
        converted_sublist = []
        for j in range(len(lst2d[i])):
            converted_sublist.append((lst2d[i][j] + extension_params))
        converted_list.append(converted_sublist)
    return converted_list


def extend_2d_list_by_function_on_index(lst2d, f):
    converted_list = []
    for i in range(len(lst2d)):
        converted_sublist = []
        for j in range(len(lst2d[i])):
            converted_sublist.append((lst2d[i][j], f(i,j)))
        converted_list.append(converted_sublist)
    return converted_list        

def extend_2d_list_by_function_on_field(lst2d, f, k):
    converted_list = []
    for i in range(len(lst2d)):
        converted_sublist = []
        for j in range(len(lst2d[i])):
            converted_sublist.append((f(lst2d[i][j][k])))
        converted_list.append(converted_sublist)
    return converted_list        


class Well:
    def __init__(self, idx, contents, resource):
        self.idx = idx
        self.contents = contents
        self.resource = resource      
    
      

class PlateMap:
    """
    The PlateMap class represents a plate used in a laboratory setting. A 
    PlateMap object is initialized with a matrix representing the contents of the plate, 
    a resource name, and a target flag indicating whether or not the plate is a target plate. 
    """
    
    def __init__(self, matrix, resource, target = False):
        self.matrix = matrix
        self.resource = resource
        self.add_idxs_to_matrix()
        self.add_resource_name_to_matrix()
        self.matrix = self.convert_matrix_to_wells()
        self.convert_wells_to_none()
        self.target = target
        self.num_columns = len(self.matrix[0])
        self.spacing = len(self.matrix)//8 - 1
    
    
    def convert_wells_to_none(self):
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                if self.matrix[i][j].contents == 'o':
                    self.matrix[i][j].contents = None
    
    def print_map(self):
        for row in self.matrix:
            print(row)
    
    def map_well_idx(self, i,j): #Row i, column j
        return i+j*len(self.matrix)
    
    def add_idxs_to_matrix(self):
        self.matrix = extend_2d_list_by_function_on_index(self.matrix, self.map_well_idx)
    
    def add_resource_name_to_matrix(self):
        self.matrix = extend_2d_list(self.matrix, self.resource)
    
    def column_by_index(self, i):
        return [row[i] for row in self.matrix]
    
    def convert_matrix_to_wells(self):
        converted_list = []
        for i in range(len(self.matrix)):
            converted_sublist = []
            for j in range(len(self.matrix[i])):
                well = Well(idx = self.map_well_idx(i, j), 
                            contents = self.matrix[i][j][0],
                            resource = self.resource)
                converted_sublist.append(well)
            converted_list.append(converted_sublist)
        return converted_list     

def rep_well(well):
    if well:
        return (well.resource, well.idx)
    else:
        return None
    
def rep_seq(seq):
    rep_a = [rep_well(well) for well in seq.aspirate]
    rep_d = [rep_well(well) for well in seq.dispense]
    return [rep_a, rep_d]

class Sequence:
    """A class representing a sequence of liquid transfers.
    This class represents a sequence of liquid transfers that can be executed
    by a liquid handling robot. The sequence is defined by a list of tuples,
    each tuple containing the source and destination wells for each transfer.
    The `aspirate` and `dispense` attributes contain the source and destination
    wells, respectively. The `multiple_dispense` attribute contains a list of
    dispensing sequences that result from combining this sequence with other
    sequences using the `combine_dispense` method.
    """

    
    def __init__(self, sequence_list):
        self.aspirate = self.pad_to_eight([transfer[0] for transfer in sequence_list])
        self.dispense = self.pad_to_eight([transfer[1] for transfer in sequence_list])
        self.multiple_dispense = [self.dispense]
    
    
    def pad_to_eight(self, lst):
        while len(lst) < 8:
            lst.append(None)
        return lst

        
    
    def combine_dispense(self, seq_to_add):
        aspiration_contents = [well.contents if well else None for well in self.aspirate]
        for dispense in seq_to_add.multiple_dispense:
            channel_contents = aspiration_contents.copy()
            new_dispense = [None]*8
            for well in dispense:
                if not well:
                    continue
                idx = channel_contents.index(well.contents)
                new_dispense[idx] = well
                channel_contents[idx] = None
            self.multiple_dispense.append(new_dispense)
    

class TargetPlate(PlateMap):
    
    """
    The TargetPlate class is a child class of PlateMap that 
    represents a target plate. It defines several methods that match target 
    wells with corresponding wells on a source plate and optimizes the matching 
    sequences by combining them if possible.
    """
    
    def __init__(self, matrix, resource):
        super().__init__(matrix, resource)
        self.target_list = [element for innerList in self.matrix for element in innerList]
    
    
    def match_targets_sources(self, source_map):
        
        """
        Matches target columns with source columns and returns a list of sequences
       of matched target and source wells.

        """
        matches = []
        target_spacing = self.spacing
        source_spacing = source_map.spacing
        for column_idx in range(self.num_columns):
            target_column = self.column_by_index(column_idx)
            for source_column in source_map.all_columns():
                if not self.test_columnwise_contents(source_column, target_column):
                    continue
                match = self.columnwise_match(source_column, source_spacing, target_column, target_spacing)
                matches += match
        return matches
    
    def test_columnwise_contents(self, col1, col2):
        return any(well.contents in (src.contents for src in col1) and well.contents != None for well in col2)
    
    def test_add_to_sequence(self, source, target, source_spacing, target_spacing, sequence):
        """
        Tests if the given source and target wells can be added to the given sequence. 
        This is true if their contents match and the spacing between them is correct given 
        the source and target spacings.
        """
        match = source.contents == target.contents and target.contents != None
        empty_sequence = len(sequence) == 0
        return match and (empty_sequence
        or (source.idx > sequence[-1][0].idx + source_spacing \
        and target.idx > sequence[-1][1].idx + target_spacing))
    
    def columnwise_match(self, source_col, source_spacing, target_col, target_spacing):
        unmatched_targets = target_col.copy()
        sequence = []
        sequences = []
        while self.test_columnwise_contents(source_col, unmatched_targets):
            sequence = []
            for source in source_col:
                for target in unmatched_targets:
                    if self.test_add_to_sequence(source, target, source_spacing, target_spacing, sequence):
                        sequence.append([source, target])
                        unmatched_targets.remove(target)
                    if len(sequence) == 8:
                        sequences.append(Sequence(sequence))
                        sequence = []
            if len(sequence)>0:
                sequences.append(Sequence(sequence))
        return sequences
    
    def combine_dispense_test(self, seq1, seq2):
        if seq1 == seq2:
            return False
        return all(well in seq1.aspirate or well==None for well in seq2.aspirate)
        
        
    def test_subset_sequences(self, sequences):
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if i != j and all(well in sequences[j].aspirate or well==None for well in sequences[i].aspirate):
                    return True
        return False

    
    def combine_sequences(self, sequence_list):
        sequences = sequence_list.copy()
        while self.test_subset_sequences(sequences):
            for seq1 in sequences:
                for seq2 in sequences:
                    if self.combine_dispense_test(seq1, seq2):
                        sequences.remove(seq1)
                        seq1.combine_dispense(seq2)
                        sequences.append(seq1)
                        sequences.remove(seq2)
                        break
        return sequences
    
    def index_sequence(self, disp):
        for channel in disp:
            if isinstance(channel, Well):
                return channel.idx
        return None
    
    def reorder_sequences(self, sequences):
        reordered_sequences = []
        for sequence in sequences:
            sequence.multiple_dispense = sorted(sequence.multiple_dispense, key=lambda disp: self.index_sequence(disp)) 
            reordered_sequences.append(sequence)
        return reordered_sequences
    
    def match_and_optimize(self, source):
        sequences = self.match_targets_sources(sources)
        sequences = self.combine_sequences(sequences)
        sequences = self.reorder_sequences(sequences)
        return sequences
        
            

    
class SourceMap:
    """
    The SourceMap class is used to represent the mapping between a 
    collection of plates and the columns on those plates.
    """


    def __init__(self, plate_size):
        self.carrier_list = []
        self.col_nums = int(math.sqrt(plate_size/(2/3)))
        self.spacing = int(plate_size/self.col_nums)//8 - 1
    
    def add_carrier(self, plate_list):
        
        if not isinstance(plate_list, list):
            plate_list = [plate_list]
        
        col_nums = [len(plate.matrix[0]) for plate in plate_list]
        
        if not all([col_num == self.col_nums for col_num in col_nums]):
            raise Exception("Every plate in a carrier must have the same number of columns")
        
        self.carrier_list.append(plate_list)
    
    def column_by_carrier(self, carrier_idx, column_idx):
        column = []
        carrier = self.carrier_list[carrier_idx]
        for plate in carrier:
            column += plate.column_by_index(column_idx)
        return column
    
    def all_columns(self):
        all_columns = []
        for carrier_idx in range(len(self.carrier_list)):
            for column_idx in range(self.carrier_list[carrier_idx][0].num_columns):
                all_columns.append(self.column_by_carrier(carrier_idx, column_idx))
        return all_columns
    


matrix_csv = 'drug_plate.csv'

def plate_from_csv(csvfile, name):
    matrix = np.genfromtxt(csvfile, delimiter=',', dtype='unicode')
    matrix = matrix.tolist()
    plateMap = PlateMap(matrix, name)
    return plateMap

def target_from_csv(csvfile, name):
    matrix = np.genfromtxt(csvfile, delimiter=',', dtype='unicode')
    matrix = matrix.tolist()
    plateMap = TargetPlate(matrix, name)
    return plateMap

class TipRack:
    
    def __init__(self, rack):
        self.rack = rack
        self.starting_tips = rack._num_items
        self.remaining_tips = rack._num_items
    
    def get_tips(self, num_tips):
        current_tip = self.starting_tips - self.remaining_tips
        tips_list = [(self.rack, tip) for tip in range(current_tip, current_tip + num_tips)]
        self.remaining_tips -= num_tips
        return tips_list
    
    def get_tips_seq(self, seq):
        num_tips = len([ch for ch in seq.aspirate if ch])
        tips = self.get_tips(num_tips)
        return tips


lmgr = LayoutManager('assets//deck.lay')
drug_a_tips = lmgr.assign_unused_resource(ResourceType(Tip96, "drug_a_tips"))
drug_a_tips = TipRack(drug_a_tips)


drug_a_plate = lmgr.assign_unused_resource(ResourceType(Plate384, "drug_a_plate"))
drug_b_plate = lmgr.assign_unused_resource(ResourceType(Plate384, "drug_b_plate"))
dilution_plate = lmgr.assign_unused_resource(ResourceType(Plate96, "dilution_plate"))

plate_a = target_from_csv('drug_a_plate.csv', drug_a_plate)
plate_b = target_from_csv('drug_b_plate.csv', drug_b_plate)
dilutions = plate_from_csv('drug_dilutions.csv', dilution_plate)

sources = SourceMap(96)
sources.add_carrier([dilutions])

sequences_a = plate_a.match_and_optimize(sources)
sequences_b = plate_b.match_and_optimize(sources)


def sequence_to_poss(seq):
    return [(well.resource, well.idx) if well else None for well in seq]

def list_diff(list_a, list_b):
    return [el[0]-el[1] if el[1] else el[0] for el in zip(list_a, list_b)]

if __name__ == '__main__': 
    with HamiltonInterface(simulate=True) as ham_int:
        normal_logging(ham_int, os.getcwd())
        initialize(ham_int)

        for sequence in sequences_a:
            tips = drug_a_tips.get_tips_seq(sequence)
            tip_pick_up(ham_int, tips)
            
            aspiration_poss = sequence_to_poss(sequence.aspirate)
            aspiration_vols = [300 if pos else None for pos in aspiration_poss]
            channel_vols = aspiration_vols.copy()
            
            aspirate(ham_int, aspiration_poss, vols = aspiration_vols, liquidClass = 'StandardVolumeFilter_Water_DispenseJet_Empty')
           
            for disp in sequence.multiple_dispense: 
                dispense_poss = sequence_to_poss(disp)
                dispense_vols = [20 if pos else None for pos in dispense_poss]
                if any([ch<20 for ch in channel_vols if ch!=None]):
                    vols = list_diff(aspiration_vols, channel_vols)
                    channel_vols = aspiration_vols.copy()
                    aspirate(ham_int, aspiration_poss, vols = vols, liquidClass = 'StandardVolumeFilter_Water_DispenseJet_Empty')
                channel_vols = list_diff(channel_vols, dispense_vols)
                dispense(ham_int, dispense_poss, vols = dispense_vols, liquidClass = 'StandardVolumeFilter_Water_DispenseJet_Empty')