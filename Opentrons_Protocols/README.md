# Opentrons_Code
----------------
A repository containing code used to replicate reactions performed on the Opentron platform.

## Guide to Using Opentrons Code
--------------------------------
Solvent benchmarking scripts: `Solvent_benchmark_acetonitrile_automated_screen.py`, `Solvent_benchmark_acetonitrile_guest_assay.py` 
Scripts for the experimental protocols on the OT-2: `automated_synthesis_opentron_script.py`, `automated_host_guest_assay_opentron_script.py`

`substance_locations`, the locations of plate (1-9) and the well plate numbers are specified for each substance, in addition to the amount in each well and the substance name.
The amount value is specified so the Opentrons knows to move onto the next well with the same substance.

`move_commands` contains the commands to move the Opentrons to the specified locations. Each move must contain the following information:
- `substance`: The name of the substance to be moved.
- `amount`: The amount of the substance to be moved in µL.
- `plate`: The plate number to be moved to on the deck.
- `location`: The well location on the plate to be moved to. Mulptiple locations can be specified from the same plate.


## Opentrons Parameters
-----------------------

Parameters selected for the Opentrons were carefully chosen based on performing a benchmark aspirating and dispensing different volumes of solvent around the system.
For each new solvent system, a benchmark should be performed prior to running experiments to investigate the optimal parameter for each solvent system.
Otherwise, issues the users may face include pipette dripping and incorrect aspirate and dispensing amounts.

The main parameters to control the Opentrons, alongside their default values used, include:
- Gantry speed - the speed in which the robot moves. 
- Plunger flow rates speed - controls the speed at which liquuids are aspirated and dispensed at. 
- Air gap - Aspirating an air gap of 15 uL has been found to reduce chance of the robot dripping and is included by default for each movement in the example script.
- Max aspiration volume - Controls the maximum amount aspirating in the pipette, not including the air gap. Default values: 285 uL for the 300 pipette.
These parameters have been validated on acetonitrile, but a optimal settings may not be translatable to different solvent systems.
