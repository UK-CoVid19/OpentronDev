#!/usr/bin/env python
# coding: utf-8

# ## RNA extraction protocol using BOMB.bio kit.
# 
# 
# The following code commands the OT2 to extract COVID-19 RNA from a liquid samples using the Beckman Coulter RNAdvance Viral XP kit.
# 
# 
# ## Protocol
# [Refer to OpenCell Station B protocol instructions for operation guide:
# https://docs.google.com/document/d/11Kfc2KW56N5ggyGUDIIb-YVgKA0DsH8MjnlfujjUXmc/edit?usp=sharing]
# 

## Resources & information
#
# Beckman Coulter protocol [Beckman Coulter RNAdvance Viral XP]:
# (https://bec-techdocs-prod.s3.us-west-2.amazonaws.com/techdoc/files/ifu/en/C58637AB.pdf?AWSAccessKeyId=ASIA2KJI7HMZX3TYHYFZ&Expires=1596550598&Signature=PVVD2ItFQlqHTa2DBaEFx6YRWXw%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEN3%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCIHlkgpKrGPkBPZerRbom5WB3MWFZR2U8aQ7M8eFQd68lAiEAowcmbF%2F56bnDt5MPamADeaqLab3fSWLxcOwxQv50BKkq2wEIpv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgw3MDkyOTI0MDc2MDMiDHdWLXepk%2F%2FUk552lCqvAVclZnz7pC38DXxT1UKrtWvqp%2BFPIUqMmtjsuLRwY91Dr8SSGx3fJ94XGvZ9DihEe3RsxvxNZmQbRfrXrvzeZO%2Bb7njq7BUXcIyzbDrGRADErz7yjmykpARUwLbj75gtqU7rqxmSArQVRiag4VfJz6wtoUw07QIqUbOZzcKpotVdJmyJMb0lxJJhFvxsFEbK6KN0WxTgKmuuDMAUz%2F7DPd6LBJwPzWkAJF7l%2F12DqpEwv6yl%2BQU64AFsbkbYn9r91PIrGWbfAkc%2F6mgI5alxALdDSaWWSIfmLECeAtX57riPZl2MT%2FIImZoVEQYn2pfXyJeh%2FHZbayGgcae5gzn0B%2F4WxpRxmUqhpETBsZU30lXEKeyiSsETNKqhWMn9RI0U6LDuBTJjk59wgtUGkhFrCq%2BOG9Fc2%2Fax2x7foZoNuCOxWb7xy97G%2BVox%2Fl2HLeV23f%2BTQLWceXaJqqd%2FfpDFKH%2BzNAPMBsQX2yxb5m4oiOkY87jUROJ0Uo2bWELLwDHcBVlpazwe%2F%2B5vWQGJpk20uNEyAX0j5BKq1A%3D%3D)
# Opentrons OT2 API v2 [OpentronsPythonAPIV2.pdf]:
# (https://docs.opentrons.com/OpentronsPythonAPIV2.pdf)

metadata = {
    'protocolName': 'Beckman Coulter RNAdvance Viral XP V1',
    'author': 'Aubin Fleiss <afleiss@ic.ac.uk>, Neil MacKenzie, Eyal Kazin <eyalkazin@gmail.com>, Alex Perkins <a.perkins19@ic.ac.uk>, Mara Donora <mara@opencell.bio>', 
    'source': 'Testing' #'Custom Protocol Request'
}


#####################
#USER DEFINED VALUES#
#####################
number_of_sample_columns = 12
test_mode = False
#####################
#                   #
#####################

# import standard modules
from collections import OrderedDict
import time
import numpy as np
# import Opentrons modules
from opentrons import labware, instruments, modules, robot, types


# magnetic module
magdeck = modules.load('magdeck', '9')
magdeck.disengage()


# Deep well plate
plate_name = 'fischerbrand_96_wellplate_2000ul'
if plate_name not in labware.list():
    custom_plate = labware.create(
        plate_name,                    # name of you labware
        grid=(12, 8),                    # specify amount of (columns, rows)
        spacing=(8.96, 8.96),             # distances (mm) between each (column, row)
        diameter = 8.5,                     # diameter (mm) of each well on the plate
        depth=41,                       # depth (mm) of each well on the plate
        volume=2000)
    
    
# PCR plate
axy_plate = 'axygen_96_wellplate_400ul'
if axy_plate not in labware.list():
    custom_2_plate = labware.create(
        axy_plate,                    # name of you labware
        grid=(12, 8),                    # specify amount of (columns, rows)
        spacing=(9, 9),             # distances (mm) between each (column, row)
        diameter = 5.3,                     # diameter (mm) of each well on the plate
        depth=20,                       # depth (mm) of each well on the plate
        volume=400)
    

# reagents plate
#use deep well for now
trough = labware.load(plate_name, '8', 'trough')

# ethanol plate
ethanol_plate = labware.load('fischerbrand_96_wellplate_2000ul', '6', share=True)

# fresh plate
pcr_plate= labware.load(axy_plate, '1', 'fresh plate')

# sample plate
sample_plate = labware.load(plate_name, '9', share=True)


# instanciate tip rack in remaining slots
tip_rack_1 = labware.load('opentrons_96_filtertiprack_200ul', '2')
tip_rack_2 = labware.load('opentrons_96_filtertiprack_200ul','4')
tip_rack_3 = labware.load('opentrons_96_filtertiprack_200ul','5')
tip_rack_4 = labware.load('opentrons_96_filtertiprack_200ul', '7')
tip_rack_5 = labware.load('opentrons_96_filtertiprack_200ul', '10')
tip_rack_6 = labware.load('opentrons_96_filtertiprack_200ul', '11')

#these tips are mapped to the sample wells and are ONLY used for the wash steps
tip_rack_ethanol_wash = labware.load('opentrons_96_filtertiprack_200ul', 3)


tips = [tip_rack_1, tip_rack_2, tip_rack_3, tip_rack_4, tip_rack_5, tip_rack_6] 


# ## Instanciate pipette and set flow rate

# load pipette
m300 = instruments.P300_Multi(mount='right', tip_racks=tips)

m300.set_flow_rate(aspirate=150, dispense=150)


NEW_TIP_MODE = 'never'

if test_mode:
    MIX_REPETITIONS = 2
    MIX_REPETITIONS_WATER = 5
else:
    MIX_REPETITIONS = 10
    MIX_REPETITIONS_WATER = 180

reagents = OrderedDict()



#  350 μl of silica-coated magnetic beads
reagents['magnetic_beads'] = {'well': 'A12', 
                              'transfer_volume': 350, 
                              'mix_volume': 100, 
                              'mix_repetitions': MIX_REPETITIONS,
                              'new_tip': NEW_TIP_MODE}


# 40 µl of nuclease-free water 
reagents['nuclease_free_water'] = {'well': 'A5', 
                                   'transfer_volume': 40, 
                                   'mix_volume': 20, 
                                   'mix_repetitions': MIX_REPETITIONS_WATER,
                                   'new_tip': NEW_TIP_MODE}


# reagents setup

for reagent_name in reagents:
    reagents[reagent_name]['setup'] = trough.wells(reagents[reagent_name]["well"])


# Define custom functions

def mix_wells(mix_locations, mix_reps):
    """ Function to mix [mix_locations] thoroughly by aspirating/rejecting liquid at different heights in a well,
    performed [mix_reps] times """

    m300.set_flow_rate(aspirate=300, dispense=550)
    
    for well in mix_locations:
        
        if not m300.tip_attached:
            m300.pick_up_tip()

        m300.move_to(well.top(20), strategy='arc') 
        
        for position in range(20,2,-2):
            m300.aspirate(volume=200, location=well.bottom(position), rate=1.0)
            m300.blow_out(well.bottom(position))

            m300.aspirate(volume=200, location=well.bottom(position), rate=1.0)
            m300.blow_out(well.top(-2))

        m300.move_to(well.top(20), strategy='arc')
    m300.set_flow_rate(aspirate=150, dispense=150)

def resuspend(well_to_mix):
    """ Function to resuspend contents of [well_to_mix] by pipetting liquid up and down while gradually descending into the well """

    m300.set_flow_rate(aspirate=150, dispense=150)
        
    if not m300.tip_attached:
        m300.pick_up_tip()

    m300.move_to(well_to_mix.top(10), strategy='arc') # fist move to the well
        
        # then aspirate and reject

    for position in np.arange(1.2,0.4, -0.2):
        print(position)
        m300.aspirate(volume=150, location=well_to_mix.top(-5), rate=1.0)
        m300.mix(5, 50, location=well_to_mix.bottom(position))
        m300.dispense(volume=150, location=well_to_mix.top(-5), rate=1.0)

    m300.move_to(well_to_mix.top(20), strategy='arc')

def resuspendLITE(well_to_mix):
    """ Function to resuspend contents of [well_to_mix] by pipetting liquid up and down while gradually descending into the well (less) """

    m300.set_flow_rate(aspirate=150, dispense=150)
        
    if not m300.tip_attached:
        m300.pick_up_tip()

    m300.move_to(well_to_mix.top(10), strategy='arc') # fist move to the well
        
        # then aspirate and reject

    for position in np.arange(0.8,0.4, -0.2):
        print(position)
        m300.aspirate(volume=150, location=well_to_mix.top(-5), rate=1.0)
        m300.mix(5, 50, location=well_to_mix.bottom(position))
        m300.dispense(volume=150, location=well_to_mix.top(-5), rate=1.0)

    m300.move_to(well_to_mix.top(20), strategy='arc')
        

def transfer_and_mix(reagent, samples):
    """ Custom function to transfer [reagent] from correct source wells to [samples] & mix """
    
    for s in samples:

        if not m300.tip_attached:
            m300.pick_up_tip()
            
        #Air gap of 10ul to help avoid dripping
        m300.transfer(reagent['transfer_volume'], reagent['setup'].bottom(0.6), s.top(-10), new_tip=reagent['new_tip'], air_gap=10)
        m300.set_flow_rate(aspirate=200, dispense=200)
        aspirate_volume = 200-reagent['mix_volume']
        m300.aspirate(volume=aspirate_volume, location=s.top(10), rate=1.0)
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)
        m300.dispense(volume=aspirate_volume, location=s.top(10), rate=1.0)
        m300.blow_out()
        m300.set_flow_rate(aspirate=150, dispense=150)
        m300.drop_tip()

def transfer_and_mixBeads(reagent, samples):
    """ Custom function to transfer [Beads] from correct source wells to [samples] & mix 
    (where [Beads] = [reagent])"""
    
    for s in samples:

        well_code = str(s).split(" ")[-1][:-1]
        if well_code in ['A1','A2','A3','A4']:
            sourcewell = trough.wells('A12')
        elif well_code in ['A5','A6','A7','A8']:
            sourcewell = trough.wells('A11')
        elif well_code in ['A9','A10','A11','A12']:
            sourcewell = trough.wells('A10')

        if not m300.tip_attached:
            m300.pick_up_tip()

        #Resuspends the beads before each transfer; resuspends more thoroughly every 4 transfers.
        if well_code in ['A1', 'A5', 'A9']:
            resuspend(sourcewell)
        else:
            resuspendLITE(sourcewell)

        #Air gap of 10ul to help avoid dripping
        m300.transfer(reagent['transfer_volume'], sourcewell.bottom(0.6), s.top(-10), new_tip=reagent['new_tip'], air_gap=10)
        m300.set_flow_rate(aspirate=200, dispense=200)
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)
        m300.blow_out()
        m300.set_flow_rate(aspirate=150, dispense=150)
        m300.drop_tip()


        
def trash_supernatant(volume, height, samples):
    """ function to remove [volume in ul] of supernatant from [samples], pipetting [height] units from the bottom of the well"""
    
    for s in samples:
        m300.pick_up_tip()
        if volume <=190:
            m300.aspirate(volume=10, location=s.top(10), rate=1.0)
        m300.transfer(volume, s.bottom(height), m300.trash_container.top(5), new_tip='never', air_gap=10, blow_out = True)
        # extra blowout:
        m300.delay(seconds = 1)
        m300.dispense(10)
        m300.delay(seconds = 1)
        # transfer function tends to eject a small volume of air after all liquid is trashed
        # which forms bubbles and may lead to cross contaminations (does not happen with all liquids
        # Keep eyes peeled at this stage)
        m300.drop_tip()
                
        
def text_in_a_box(line,border_char="#"):
    """ function to print some text in a box of asterisks"""
    
    new_text=str(line)
    line_len = len(line)
    new_text = "\n"+border_char*(line_len+4)+"\n"+border_char+" "+line+" "+border_char+"\n"+border_char*(line_len+4)
    
    return(new_text)



# IMPORTANT REMARKS

# the API is not too robust yet as regards sanity checks
# Consequently the robot is still a danger to itself 
# and has pronounced taste for self-destruction

# never ever remove the block below
# unless you want the robot to pipette wells located beyond plate boundaries
# crushing all your labware

# also, when using a multi-channel pipette, make sure you are ALWAYS 
# using well coordinates from first row (A1 to A12) of your 96-well plate
# unless you want to spent countless hours re-calibrating your robot after
# its arm collided on external walls

if number_of_sample_columns > 12:
    raise Exception("Please specify a valid number of sample columns.")
    

samples = sample_plate.rows('A')[0:number_of_sample_columns]

# home
robot.home()


# Add 350 µl of silica-coated magnetic beads
transfer_and_mixBeads(reagents['magnetic_beads'], samples)

m300.delay(seconds=300)


# Settle the magnetic beads on a magnetic stand and discard the supernatant
robot.comment("Activating magdeck for 90 seconds")
magdeck.engage(height=12)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(seconds=600)


# trash supernatant
trash_supernatant(volume=650, height=0.4, samples=samples)


# In order to map the tips to the samples I have had to write the code outside of a function.
# It should make it more readable anyway

magdeck.disengage()



# ethanol wash (400 ul), repeated 2 times
robot.comment(text_in_a_box("Ethanol wash steps. Uses specific tips. Loops 2x"))

reps_ = 2
for _ in range(reps_):
    
    magdeck.disengage()
    
    for well in samples:
        
        #maps tips to sample well - uses specific tip box
        well_code = str(well).split(" ")[-1][:-1]
        
        #if not m300.tip_attached:
        m300.pick_up_tip(tip_rack_ethanol_wash[well_code])

        m300.transfer(400, 
                      ethanol_plate.wells(well_code).bottom(2), 
                      well.top(-10), new_tip='never', air_gap = 10)
        
        m300.set_flow_rate(aspirate=200, dispense=250)
        m300.aspirate(100, well.top(20))
        m300.mix(MIX_REPETITIONS, 100, well)
        m300.dispense(100, well.top(-20))
        m300.set_flow_rate(aspirate=150, dispense=150)

        m300.return_tip()


    robot.comment("Activating magdeck for 120 seconds")
    magdeck.engage(height=12)
    
    if test_mode:
        m300.delay(seconds=5)
    else:
        m300.delay(seconds=120)

    #trash_supernatant(volume=400, height=2, samples=samples, pipette = 'ethanol')
    for well in samples:
        
        #uses same tips
        well_code = str(well).split(" ")[-1][:-1]
        
        m300.pick_up_tip(tip_rack_ethanol_wash[well_code])
        # trashes supernatant from the bottom of the well (0.2mm) if last repetition
        # ensures maximal ethanol removal before drying stage
        if _ == (reps_-1):
            m300.transfer(400, well.bottom(0.2), m300.trash_container.top(10), new_tip='never' , air_gap = 10, blow_out = True)
        else:
            m300.transfer(400, well.bottom(0.6), m300.trash_container.top(10), new_tip='never' , air_gap = 10, blow_out = True)

        # to remove bubbles before returning tips to box:
        m300.set_flow_rate(aspirate=20, dispense=150)
        m300.aspirate(40)
        m300.dispense(20)
        m300.delay(seconds=2)
        m300.dispense(20)
        m300.set_flow_rate(aspirate=150, dispense=150)
        # transfer function tends to eject a small volume of air after all liquid is trashed
        # which forms bubbles and may lead to cross contaminations (does not happen with all liquids
        # Keep eyes peeled at this stage)#
        
        m300.return_tip()

magdeck.disengage()


# Bead drying stage
if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(seconds=60)



# Add 40 µl of nuclease-free water to elute RNA, mix, incubate for 5 mins
transfer_and_mix(reagents['nuclease_free_water'], samples)
m300.delay(seconds=300)

#turn on Magdeck to remove beads
robot.comment("Activating magdeck for 120 seconds")
magdeck.engage(height=12)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(seconds=120)

#transfer 40ul of eluted sample to PCR plate
# pcr plate mapped to samples.
# aspirate from (near) bottom of well
# air gap of 10ul to protect sample

for well in samples:
        
        if not m300.tip_attached:
            m300.pick_up_tip()

        well_code = str(well).split(" ")[-1][:-1]
        m300.set_flow_rate(aspirate=30, dispense=30)
        m300.transfer(40, well.bottom(0.3), pcr_plate.wells(well_code).bottom(0.5), new_tip='always', air_gap=10, blow_out = True)
        
magdeck.disengage()