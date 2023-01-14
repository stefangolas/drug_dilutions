# Drug dilutions
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/colors.png" width="320"/>
This is a PyHamilton (Hamilton Microlab STAR) library for dispensing 2-fold serial dilutions of 2 different therapeutic drugs onto 384-well plates and then combining these into one 384-well plate with quad-pinning for combinatorial drug dosage screening. </br>
</br>


Script created for "High-throughput Approaches to Uncover Synergistic Drug Combinations in Leukemia" </br>
Preprint here: https://www.biorxiv.org/content/10.1101/2022.11.29.518409v1

## Updated Script

`robot_method.py` is an updated and heavily refactored script that provides classes and functions for creating your own optimized workflows using plate maps.

## Original Process

1. The user runs the main script by running `py dispense_dilutions_2rep.py` with an optional `--make i j` flag where `i` and `j` represent source wells for concentrated drug reagent in a 96-well plate. If `--make i j` is passed then the robot will dispense from these two wells into the respective starting wells for each drug diluton series. If not, the script will assume the user has performed this step manually.</br>

2. Once the initial well for each dilution series is loaded with the specified drug, the robot performs serial 2-fold dilutions.
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/bunguloj3.gif" width="280"/>
</br>

3. The script then reads csv files `drug_a_plate.csv` and `drug_b_plate.csv` containing the dose patterns for each target plate represented visually in a 2D matrix and calculates a time-efficient set of dispense series to load each well as specified.</br> 
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/hamiltonstar5.png" width="420"/>
</br>


4. The script then runs a dispense routine which automatically determines implementation details such as tip pick up and refilling tips from source wells.</br>
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/plas23.gif" width="210"/>
</br>

5. Finally, the reagents in the first plate are added to the second plate using the 96-channel head for quad-pinning.
