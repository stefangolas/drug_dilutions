# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 20:29:12 2022

@author: stefa
"""
from method_labware import *
from pace_util import * #TODO: import everything explicitly, at least from pace_util


if __name__ == '__main__':

    with HamiltonInterface(simulate=simulating) as ham_int, ClarioStar() as reader_int: #TODO: simulate=False
        normal_logging(ham_int)
        
        initialize(ham_int)

        tip_pick_up_96(ham_int, drug_b_tips)
        aspirate_384(ham_int, drug_b_plate, 0, 80)
        for quadrant in range(1,4):
            dispense_384_quadrant(ham_int, drug_b_plate, quadrant, 20)
        tip_eject(ham_int, drug_b_tips)
        
        tip_pick_up_96(ham_int, drug_a_tips)
        aspirate_384(ham_int, drug_a_plate, 0, 80)
        for quadrant in range(1,4):
            dispense_384_quadrant(ham_int, drug_b_plate, quadrant, 20)