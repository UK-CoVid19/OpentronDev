#!/usr/bin/env python
# coding: utf-8

# ## RNA extraction protocol using BOMB.bio kit.
# 
# 
# The following code commands the OT2 to extract COVID-19 RNA from a liquid samples using the BOMB.bio kit.
# 
# 
# ## Protocol
# 
# 1. Specify number of samples to run, in groups of eight. e.g. eight samples = 1 sample channel, 32 samples = 4 sample channels etc.
# 
# 2. Run Cell 4 to generate Python script called: rna_extraction_jupyter_exported.py
# 
# 3. Turn on OT2, connect and calibrate the deck.
# 4. Run Python script and perform labware calibration with dummy labware.
# 
# 5. Click 'return tip and proceed to run'. Before starting:
# 6. Clear deck and clean/RNA Zap deck, labware, pipettes and walls.
# 7. Install RNA free labware with regents. Remember to remove lids.
# 
# 8. Click "start".
# 
# 

# In[ ]:


## Resources & information


#* BOMB.BIO protocol: [BOMB total RNA extraction mammalian GITC v1.0.pdf](https://bomb.bio/wp-content/uploads/2018/09/8.2_BOMB_total_RNA_extraction_mammalian_GITC_V1.0.pdf)
#<p>
#* Opentrons OT2 API v2: [OpentronsPythonAPIV2.pdf](https://docs.opentrons.com/OpentronsPythonAPIV2.pdf)


# In[ ]:


metadata = {
    'protocolName': 'RNA Extraction v0.2',
    'author': 'Aubin Fleiss <afleiss@ic.ac.uk>, Neil MacKenzie, Eyal Kazin <eyalkazin@gmail.com>, Alex Perkins <a.perkins19@ic.ac.uk>, M Donora <matthew@opencell.bio>', 
    'source': 'Testing' #'Custom Protocol Request'
}


# ## To simulate/export protocol
# 
# The next cell exports the jupyter notebook to rna_extraction_jupyter_exported.py. The exported file can be either:
# - used directly in the Opentrons app
# - simulated in the command line for instance : $ opentrons_simulate rna_extraction.py
# - simulated within this notebook. To do so run the cell after the next

# ## User parameters: 
# - <b>number_of_sample_columns</b> (integer) defines the number of columns in the test plate that contain samples. This parameter is probably the only one the user will provide once the protocol is established
# <p>
# - <b>DNAse_incubation</b> (bool) switch on (True) or off (False) DNAse treatment
# <p>
# - <b>test_mode</b> (bool) in test mode (True) incubation times are reduced to 2 sec, mixing steps are reduced to 2 reps. In nornal mode (False), all incubation and mixing steps are restored to their normal durations

# In[ ]:


number_of_sample_columns = 3

film = True

test_mode = False

# switches on the first isopropanol wash using the same tips as the subsequent ethanol washes

# x, y, z
# need to tune on monday


# ## Installing, updating, loading modules
# 

# In[8]:


# intalling opentrons module (needed only once)

#import sys
#!{sys.executable} -m pip install opentrons 


# In[9]:


# updating opentrons module (do it once every few weeks or after an API update)

#import sys
#!{sys.executable} -m pip install --upgrade opentrons


# In[10]:


# import standard modules
from collections import OrderedDict
import time
import numpy as np
# import Opentrons modules
from opentrons import labware, instruments, modules, robot, types


# ## Connect to robot

# In[12]:


# make connected robot blink lights
"""
# some functions to play around with lights
robot.get_rail_lights_on()
robot.turn_off_rail_lights()
robot.turn_on_rail_lights()

robot.turn_off_button_light()
robot.turn_on_button_light()

robot.get_lights()
robot.set_lights(button=None, rails=None)
"""
robot.identify(5) #blink lights for 5 seconds


# ## Instanciate and initialise modules

# In[14]:


# magnetic module
magdeck = modules.load('magdeck', '9')
magdeck.disengage()


# 
# ## Instanciate labware
# 

# In[15]:


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
nunc_plate = 'nunc_96_wellplate_400ul'
if nunc_plate not in labware.list():
    custom_2_plate = labware.create(
        nunc_plate,                    # name of you labware
        grid=(12, 8),                    # specify amount of (columns, rows)
        spacing=(8.9, 8.9),             # distances (mm) between each (column, row)
        diameter = 6.8,                     # diameter (mm) of each well on the plate
        depth=12.24,                       # depth (mm) of each well on the plate
        volume=400)
    

#liquid_bin = labware.load('agilent_1_reservoir_290ml', 11)
    
# reagents plate
#use deep well for now
trough = labware.load(plate_name, '8', 'trough')

# ethanol plate
ethanol_plate = labware.load('fischerbrand_96_wellplate_2000ul', '6', share=True)

# fresh plate
pcr_plate= labware.load(nunc_plate, '1', 'fresh plate')

# sample plate
sample_plate = labware.load(plate_name, '9', share=True)


# instanciate tip rack in remaining slots
tip_rack_1 = labware.load('opentrons_96_filtertiprack_200ul', '2')
tip_rack_2 = labware.load('opentrons_96_filtertiprack_200ul','4' )
tip_rack_3 = labware.load('opentrons_96_filtertiprack_200ul','5' )
tip_rack_4 = labware.load('opentrons_96_filtertiprack_200ul', '7')
tip_rack_5 = labware.load('opentrons_96_filtertiprack_200ul', '10')
tip_rack_6 = labware.load('opentrons_96_filtertiprack_200ul', '11')

#these tips are mapped to the sample wells and are ONLY used for the wash steps
tip_rack_ethanol_wash = labware.load('opentrons_96_filtertiprack_200ul', 3)


tips = [tip_rack_1, tip_rack_2, tip_rack_3, tip_rack_4, tip_rack_5, tip_rack_6] 


# ## Instanciate pipette and set flow rate

# In[11]:



# load pipette
m300 = instruments.P300_Multi(mount='right', tip_racks=tips)

m300.set_flow_rate(aspirate=150, dispense=150)


# 
# ## Instanciate reagents
# 

# In[12]:


#lysis_buffer = trough.wells('A1')
# isopropanol = trough.wells('A2')
# magnetic_bead = trough.wells('A3')   # e.g, silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock)
# ethanol_80percent = trough.wells('A4')
# dnaseI_reaction_mix = trough.wells('A5')  # enzyme that removes DNA
# rna_binding_buffer = trough.wells('A6')
# nuclease_free_water = trough.wells('A7')
# liquid_waste = trough.wells('A12')  #  elution waste

NEW_TIP_MODE = 'never'

if test_mode:
    MIX_REPETITIONS = 2
    MIX_REPETITIONS_WATER = 5
else:
    MIX_REPETITIONS = 15
    MIX_REPETITIONS_WATER = 180

reagents = OrderedDict()

# Add 360 μl of isopropanol + magbeads (premixed)
# For 96 wells, need to add to A10, A11, A12
reagents['isopropanolBeads'] = {'well': 'A12', 
                               'transfer_volume': 360,
                               'mix_volume': 190, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}

# # Add 40 μl of silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock), seal and shake at RT at 1400 rpm for 5 min
# reagents['magnetic_beads'] = {'well': 'A10', 
#                               'transfer_volume': 40, 
#                               'mix_volume': 40, 
#                               'mix_repetitions': MIX_REPETITIONS,
#                               'new_tip': NEW_TIP_MODE}

# Remove the plate from the magnetic stand and add 400 μl isopropanol. Shake at RT at 1400 rpm for 2 min
# For 96 wells, need to add to A7, A8, A9
reagents['isopropanol_400'] = {'well': 'A9', 
                               'transfer_volume': 400, 
                               'mix_volume': 190, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}

#Add 40 µl of nuclease-free water to elute RNA, mix at 1300 rpm for 5 min
reagents['nuclease_free_water'] = {'well': 'A6', 
                                   'transfer_volume': 40, 
                                   'mix_volume': 20, 
                                   'mix_repetitions': MIX_REPETITIONS_WATER,
                                   'new_tip': NEW_TIP_MODE}


# In[13]:


# reagents setup

for reagent_name in reagents:
    reagents[reagent_name]['setup'] = trough.wells(reagents[reagent_name]["well"])


# ## Define custom functions

# In[32]:


def mix_wells(mix_locations, mix_reps):
    """mixes a well thoroughly by aspirating/rejecting liquid at different heights in a well"""

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
    """resuspends the contents of a well by pipetting liquid up and down while gradually descending into the well"""

    m300.set_flow_rate(aspirate=150, dispense=150)
        
    if not m300.tip_attached:
        m300.pick_up_tip()

    m300.move_to(well_to_mix.top(10), strategy='arc') # fist move to the well
        
        # then aspirate and reject

    for position in np.arange(1.4,0.4, -0.2):
        print(position)
        m300.mix(5, 50, location=well_to_mix.bottom(position))

    m300.move_to(well_to_mix.top(20), strategy='arc')
        

def transfer_and_mix(reagent, samples):
    
    for s in samples:

        if not m300.tip_attached:
            m300.pick_up_tip()
            
        #dispenses 30mm from bottom of the well to prevent contamination of reagents with split volumes. TO DO: weak height and speed
        #Air gap of 10ul to help avoid dripping
        m300.transfer(reagent['transfer_volume'], reagent['setup'].bottom(0.4), s.top(-10), new_tip=reagent['new_tip'], air_gap=10)
        m300.set_flow_rate(aspirate=200, dispense=200)
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)# note that according to nucleic_acid_extration.ot2.py .mix volume differs from .transfer volume
        m300.blow_out()
        m300.set_flow_rate(aspirate=150, dispense=150)
        m300.drop_tip()

def transfer_and_mixBEADS(reagent, samples):
    
    for s in samples:

        well_code = str(s).split(" ")[-1][:-1]
        if well_code in ['A1','A2','A3','A4','A5']:
            sourcewell = trough.wells('A12')
        elif well_code in ['A6','A7','A8','A9','A10']:
            sourcewell = trough.wells('A11')
        elif well_code in ['A11','A12']:
            sourcewell = trough.wells('A10')

        if not m300.tip_attached:
            m300.pick_up_tip()

        if well_code in ['A1','A6','A11']:
            resuspend(sourcewell)
        #dispenses 30mm from bottom of the well to prevent contamination of reagents with split volumes. TO DO: weak height and speed
        #Air gap of 10ul to help avoid dripping
        m300.transfer(reagent['transfer_volume'], sourcewell.bottom(0.4), s.top(-10), new_tip=reagent['new_tip'], air_gap=10)
        m300.set_flow_rate(aspirate=200, dispense=200)
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)# note that according to nucleic_acid_extration.ot2.py .mix volume differs from .transfer volume
        m300.blow_out()
        m300.set_flow_rate(aspirate=150, dispense=150)
        m300.drop_tip()

        
def trash_supernatant(volume, height, samples):
    """ function to remove [volume in ul] of supernatant from [samples], pipetting [height] units from the bottom of the well"""
    # height to be tested, more or less reliable depending on API version
    
    for s in samples:
        m300.pick_up_tip()
        m300.transfer(volume, s.bottom(height), m300.trash_container.top(5), new_tip='never', air_gap=10, blow_out = True)
        m300.dispense(200)
        m300.delay(seconds = 2)
        m300.dispense(50)
        # transfer function tends to eject a small volume of air after all liquid is trashed
        # which forms bubbles and may lead to cross contaminations (does not happen with all liquids
        # Keep eyes peeled at this stage)
        m300.drop_tip()
        
        
"""
def trash_supernatant_V2(volume, height, samples):
    #function to remove [volume in ul] of supernatant from [samples], pipetting [height] units from the bottom of the well
    # height to be tested, more or less reliable depending on API version
    

    
    for s in samples:
        m300.pick_up_tip()
        m300.transfer(volume, s.bottom(height), location=types.Location(point=types.Point(500, 500, 250), labware=None), new_tip='never', air_gap=10, blow_out = True)
        # transfer function tends to eject a small volume of air after all liquid is trashed
        # which forms bubbles and may lead to cross contaminations (does not happen with all liquids
        # Keep eyes peeled at this stage)
        m300.drop_tip()
"""       
        
        
def text_in_a_box(line,border_char="#"):
    """ function to print some text in a box of asterisks"""
    
    new_text=str(line)
    line_len = len(line)
    new_text = "\n"+border_char*(line_len+4)+"\n"+border_char+" "+line+" "+border_char+"\n"+border_char*(line_len+4)
    
    return(new_text)


def blow_air(mins, samples):
    
    #same tip
    m300.pick_up_tip()

    #MD EDIT - to potentially reduce simulation time
    #each cycle takes 14 seconds
    #UPDATE for multiple rows this takes less time in between each sample than 14s.
    repetitions = int((mins*60)/(9.6+(4.2*number_of_sample_columns)))
    for r in range(repetitions):
        for s in samples:
            #continously blows 190ul of air over beads 5 times at a height of 5mm
            #m300.mix(5, 190, s.bottom(15))
            #MD EDIT
            if number_of_sample_columns <= 10:
                aspirate_speed = number_of_sample_columns*19
            else:
                aspirate_speed = 190
            m300.set_flow_rate(aspirate=aspirate_speed, dispense=100)
            m300.aspirate(190, s.top(15))
            m300.dispense(190, s.bottom(15))
    #MD EDIT - deprecated below
    # for the duration of the time specified: 
    # t_end = time.time() + 60 * minutes
    # while time.time() < t_end:

    #     for s in samples:
    #         #continously blows 190ul of air over beads 5 times at a height of 5mm
    #         #m300.mix(5, 190, s.bottom(15))
    #         #MD EDIT
    #         m300.set_flow_rate(aspirate=19, dispense=100)
    #         m300.aspirate(190, s.top(15))
    #         m300.dispense(190, s.bottom(15))
    m300.drop_tip()
            
            


# In[ ]:



    


# ## Recapitulate setup

# In[15]:


print("================================ Setup recap ================================")


# ### Instruments

# In[16]:


# list of (mount, instrument)

print("\nRecap Instruments")
print(robot.get_instruments())


# ### Containers

# In[17]:


# List all containers on the deck

print("\nRecap containers")
for elt in robot.get_containers():
    print("\t",elt.parent, elt)


# ### Reagents

# In[18]:


print("\nRecap reagents")
for reagent in reagents:
    print("\t",trough.parent, trough, reagents[reagent]["setup"], reagent)

print("\t",ethanol_plate.parent, ethanol_plate, "<All wells>", "ethanol")


# ### Pipettes

# In[19]:


# recap attached pipettes
print("\nRecap pipette(s)")
for elt in robot.get_attached_pipettes():
    print("\t", elt,robot.get_attached_pipettes()[elt])


# ### Modules

# In[20]:


print("\nRecap modules")
print("\tmagdeck", magdeck.status)


# In[21]:


print("================================= Run Start =================================")


# In[22]:


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


# In[23]:


# home
robot.home()

#------MD----- WHAT FORM DO SAMPLES ARRIVE IN? PRE-SUPERNATANT?

# In[24]:


# steps 1-2

#m300.move_to(location=types.Location(point=types.Point(30, 40, 100), labware=None))

# sample collection, nothing to do here


# In[25]:


# step 3
#trash_supernatant(volume=900, height=2, samples=samples)

#robot.comment(text_in_a_box("step 3"))
# Add 240 µl of lysis buffer, seal and shake at RT at 1400 rpm for 5 min
#transfer_and_mix(reagents['lysis_buffer'], samples)


# In[26]:


#step 4
robot.comment(text_in_a_box("step 4"))

# Add 360 µl of isopropanol + beads, seal and shake at RT at 1400 rpm for 5 min
transfer_and_mixBEADS(reagents['isopropanolBeads'], samples)


# ### here, a function to resuspend the beads before dispensing them is needed

# In[27]:


# step 5
robot.comment(text_in_a_box("step 5"))

# resuspend the beads
#V8 - BEADS NOW PRE-MIXED WITH IPA
# well_to_mix = reagents['magnetic_beads']["setup"]
# resuspend(well_to_mix)


# # Add 40 µl of silica-coated magnetic beads
# transfer_and_mix(reagents['magnetic_beads'], samples)


# In[28]:


# step 6
robot.comment(text_in_a_box("step 6"))

# Settle the magnetic beads on a magnetic stand and discard the supernatant
# this block can probably be factorised as a function
# considering the number of times it is used throughout the protocol

robot.comment("Activating magdeck for 30 seconds")
magdeck.engage(height=12)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(seconds=30)


# In[33]:


# volume & height from bottom to be adjusted based on tests
trash_supernatant(volume=650, height=0.4, samples=samples)


# In[4]:


# step 7

# In order to map the tips to the samples I have had to write the code outside of a function.
# It should make it more readable anyway

# there are two versions of the same isopropanol wash step. One that uses the mapped ethanol tips and one that uses it's own tips.
# TBC which one to use.


robot.comment(text_in_a_box("step 7 -  Isopropanol First wash"))
#Remove the plate from the magnetic stand and add 400 µl isopropanol
# Shake at 1400 rpm for 2 min


magdeck.disengage()

for well in samples:

    #gets the well code of for the sample.
    well_code = str(well).split(" ")[-1][:-1]
    if well_code in ['A1','A2','A3','A4','A5']:
        sourcewell = trough.wells('A9')
    elif well_code in ['A6','A7','A8','A9','A10']:
        sourcewell = trough.wells('A8')
    elif well_code in ['A11','A12']:
        sourcewell = trough.wells('A7')

    #picks up the tip in the sample positon on the ethanol_wash tip rack.
    m300.pick_up_tip(tip_rack_ethanol_wash[well_code])

    m300.transfer(reagents['isopropanol_400']['transfer_volume'], sourcewell.bottom(0.4), well.top(-10), new_tip='never', air_gap = 10)

    m300.set_flow_rate(aspirate=100, dispense=100)
    m300.aspirate(100, well.top(20))
    m300.mix(MIX_REPETITIONS, 100, well)
    m300.dispense(100, well.top(-20))
    m300.set_flow_rate(aspirate=150, dispense=150)

    m300.return_tip()


robot.comment("Activating magdeck for 30 seconds")
magdeck.engage(height=12)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(seconds=30)


# volume & height from bottom to be adjusted based on tests
#trash_supernatant(volume=900, height=2, samples=samples, pipette = 'ethanol')
for well in samples:

    #uses the same tips to discard the supernatant.
    well_code = str(well).split(" ")[-1][:-1]
    m300.pick_up_tip(tip_rack_ethanol_wash[well_code])
    m300.transfer(200, well.bottom(0.4), m300.trash_container.top(5), new_tip='never', blow_out = True, air_gap=10)
    m300.dispense(40)
    m300.delay(seconds = 2)
    m300.dispense(40)
    m300.transfer(250, well.bottom(0.4), m300.trash_container.top(5), new_tip='never', blow_out = True, air_gap=10)
    m300.dispense(40)
    m300.delay(seconds = 2)
    m300.dispense(40)

    # transfer function tends to eject a small volume of air after all liquid is trashed
    # which forms bubbles and may lead to cross contaminations (does not happen with all liquids
    # Keep eyes peeled at this stage)#
    m300.return_tip()

# # In[3]:


# steps 9-10-11, repeated 4 times
robot.comment(text_in_a_box("Ethanol wash steps. Uses specific tips. Loops 4x"))

repss = 4
for rep in range(repss):
    
    magdeck.disengage()
    
    for well in samples:
        
        #maps tips to sample well - uses specific tip box
        well_code = str(well).split(" ")[-1][:-1]
        
        #if not m300.tip_attached:
        m300.pick_up_tip(tip_rack_ethanol_wash[well_code])

        m300.transfer(200, 
                      ethanol_plate.wells(well_code).bottom(2), 
                      well.top(-10), new_tip='never', air_gap = 10)
        
        m300.set_flow_rate(aspirate=200, dispense=250)
        m300.aspirate(100, well.top(20))
        m300.mix(MIX_REPETITIONS, 100, well)
        m300.dispense(100, well.top(-20))
        m300.set_flow_rate(aspirate=150, dispense=150)

        m300.return_tip()


    robot.comment("Activating magdeck for 30 seconds")
    magdeck.engage(height=12)
    
    if test_mode:
        m300.delay(seconds=5)
    else:
        m300.delay(seconds=30)
    
    # volume & height from bottom to be adjusted based on tests
    #trash_supernatant(volume=900, height=2, samples=samples, pipette = 'ethanol')
    for well in samples:
        
        #uses same tips
        well_code = str(well).split(" ")[-1][:-1]
        
        m300.pick_up_tip(tip_rack_ethanol_wash[well_code])
        if rep == (repss-1):
            m300.transfer(300, well.bottom(0.2), m300.trash_container.top(10), new_tip='never' , air_gap = 10, blow_out = True)
        else:
            m300.transfer(300, well.bottom(0.4), m300.trash_container.top(10), new_tip='never' , air_gap = 10, blow_out = True)

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


# In[ ]:


#samples


# In[ ]:


# magdeck.disengage()


# ## kindly request human to move the plate to the temperature module

# In[1]:


# step 12
#robot.comment(text_in_a_box("step 12"))

#robot.pause('Please add new tips to 5,6,9.')
#m300.reset_tipracks()

#robot.comment("Please place plate on tempdeck")
#robot.pause()


# In[18]:


# whilst beads are drying -  blow air over them for the duration specified

if test_mode:
    blow_air(1, samples)
else:
    blow_air(35, samples)
    #m300.delay(minutes=2)
    m300.set_flow_rate(aspirate=100, dispense=100)


# In[19]:


#robot.comment("Please place plate back on magdeck")
#robot.pause()


## kindly request human to move the plate to the temperature module

#In[5]:
##################
#COMMENTS BELOW MD#

# step 20
robot.comment(text_in_a_box("step 20"))

# Add 40 µl of nuclease-free water to elute RNA, mix at 1300 rpm for 5 min
transfer_and_mix(reagents['nuclease_free_water'], samples)


# In[3]:


# step 21
robot.comment(text_in_a_box("step 21"))

#turn on Magdeck to remove beads
robot.comment("Activating magdeck for 30 seconds")
magdeck.engage(height=12)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(seconds=30)

#transfer 40ul of eluted sample to PCR plate
# pcr plate mapped to samples.
# aspirate from bottom of well
# dispense 0.1mm from bottom of well
# air gap of 10ul to protect sample

for well in samples:
        
        if not m300.tip_attached:
            m300.pick_up_tip()

        well_code = str(well).split(" ")[-1][:-1]
        m300.set_flow_rate(aspirate=30, dispense=30)
        m300.transfer(40, well.bottom(0.2), pcr_plate.wells(well_code).bottom(0.5), new_tip='always', air_gap=10, blow_out = True)
        
magdeck.disengage()
robot.comment("Done, at last!")


# ### I think I would do the last transfer manually to maximise liquid recovery while minimising the amount of beads

# ## <p style="text-align: center;"> The end </p>
