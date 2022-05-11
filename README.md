# Drug dilutions
This is a PyHamilton library for dispensing serial 2-fold dilutions of 2 different therapeutic drugs onto 384-well plates and then combining these into one 384-well plate with quad-pinning for combinatorial dosage screening. </br>

The main script `py dispense_dilutions_2rep.py` reads csv files containing the dose patterns on each drug plate represented visually in a 2D matrix.</br>

The script then automatically determines a time-efficient pattern for dispensing to the target plates and autonomously performs this dispense routine independent of user specifications</br>
