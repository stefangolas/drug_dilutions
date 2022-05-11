# Drug dilutions
This is a PyHamilton library for dispensing 2-fold serial dilutions of 2 different therapeutic drugs onto 384-well plates and then combining these into one 384-well plate with quad-pinning for combinatorial dosage screening. </br>

1. The user runs the main script by running `py dispense_dilutions_2rep.py` with an optional `--make i j` flag where `i` and `j` represent source wells for concentrated drug reagent in a 96-well plate. If `--make i j` is passed then the robot will dispense from these two wells into the respective starting wells for each drug diluton series. If not, the script will assume the user has performed this step manually.</br>

2. Once the initial well for each dilution series is loaded with its drug, the robot performs serial 2-fold dilutions.
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/bunguloj2.gif" width="240"/>


3. The script then reads csv files containing the dose patterns for each target plate represented visually in a 2D matrix and calculates a time-efficient set of dispense series to load each well as specified.</br> 
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/hamiltonstar5.png" width="480"/>



4. The script then runs a dispense routine which automatically determines implementation details such as tip pick up and refilling tips from source wells.</br>
<img src="https://github.com/stefangolas/drug_dilutions/blob/master/imgs/plas23.gif" width="240"/>

