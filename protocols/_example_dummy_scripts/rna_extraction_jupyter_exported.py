#!/usr/bin/env python
# coding: utf-8

# ## To do list
# 
# * add vortexing function => added the function, but as we do not know if we really need it, I do not call it anywhere (yet) in the actual protocol. Also, it probably needs some tweaking.
# * if some dripping is observed during real tests, it might be useful to use the touch_tip commands for instance pipette.transfer(bla, bla, bla, touch_tip=True)  
# 

# ## Done / ready to be tested:
# * All reagents defined in one trough ✔️
# * except ethanol : deep well plate with each well designated to one sample => (the code is quick and dirty, but it should work) ✔️
# * fix reagents definition (ATM all reagents are initialised in well A1 ✔️
# * add one switch for DNAse (if we want to do it or not) ✔️
# * Add recap of setup after initiation  ✔️
# * switch for waiting times & mix reps (simulation vs IRL) ✔️
# * Add setup recaps ✔️
# * add a function to "vortex" and to resuspend the beads before dispensing otherwise they form clumps ✔️
# * extensive code re-organisation and documentation  ✔️

# ## Resources & information
# 
# 
# * BOMB.BIO protocol: [BOMB total RNA extraction mammalian GITC v1.0.pdf](https://bomb.bio/wp-content/uploads/2018/09/8.2_BOMB_total_RNA_extraction_mammalian_GITC_V1.0.pdf)
# <p>
# * Opentrons OT2 API v2: [OpentronsPythonAPIV2.pdf](https://docs.opentrons.com/OpentronsPythonAPIV2.pdf)
# 
# 

# In[ ]:


metadata = {
    'protocolName': 'RNA Extraction v0.2',
    'author': 'Aubin Fleiss <afleiss@ic.ac.uk>, Neil MacKenzie, Eyal Kazin <eyalkazin@gmail.com>',
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


number_of_sample_columns = 1

DNAse_incubation = True

test_mode = True


# ## Installing, updating, loading modules
# 

# In[ ]:


# intalling opentrons module (needed only once)

#import sys
#!{sys.executable} -m pip install opentrons 


# In[ ]:


# updating opentrons module (do it once every few weeks or after an API update)

#import sys
#!{sys.executable} -m pip install --upgrade opentrons


# In[ ]:


# import standard modules
from collections import OrderedDict

# import Opentrons modules
from opentrons import labware, instruments, modules, robot


# ## Connect to robot

# In[ ]:


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

# In[ ]:


# magnetic module
magdeck = modules.load('magdeck', '1')
magdeck.disengage()

# temperature module
tempdeck = modules.load('tempdeck', '4')
tempdeck.set_temperature(25)


# 
# ## Instanciate labware
# 

# In[ ]:



# sample plate
sample_plate = labware.load("corning_96_wellplate_360ul_flat", '1', share=True)

# reagents plate
trough = labware.load('trough-12row', '2', 'trough')

# ethanol plate
ethanol_plate = labware.load('96-deep-well', '5', share=True)

# fresh plate
fresh_plate = labware.load("corning_96_wellplate_360ul_flat", '3', 'fresh plate')

# instanciate tip rack in remaining slots
tips = [labware.load('opentrons-tiprack-300ul', str(slot)) for slot in range(6, 12)] 


# ## Instanciate pipette and set flow rate

# In[ ]:



# load pipette
m300 = instruments.P300_Multi(mount='right', tip_racks=tips)

m300.set_flow_rate(aspirate=150, dispense=300)


# 
# ## Instanciate reagents
# 

# In[ ]:


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
else:
    MIX_REPETITIONS = 15

reagents = OrderedDict()
# Add 240 μl of lysis buffer, seal and shake at RT at 1400 rpm for 5 min
reagents['lysis_buffer'] = {'well': 'A1', 
                            'transfer_volume': 240,
                            'mix_volume': 240, 
                            'mix_repetitions': MIX_REPETITIONS,
                            'new_tip': NEW_TIP_MODE}

# Add 320 μl of isopropanol, seal and shake at RT at 1400 rpm for 5 min
reagents['isopropanol_320'] = {'well': 'A2', 
                               'transfer_volume': 320,
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}

# Add 40 μl of silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock), seal and shake at RT at 1400 rpm for 5 min
reagents['magnetic_beads'] = {'well': 'A3', 
                              'transfer_volume': 40, 
                              'mix_volume': 40, 
                              'mix_repetitions': MIX_REPETITIONS,
                              'new_tip': NEW_TIP_MODE}

# Remove the plate from the magnetic stand and add 400 μl isopropanol. Shake at RT at 1400 rpm for 2 min
reagents['isopropanol_400'] = {'well': 'A2', 
                               'transfer_volume': 400, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}


# Add 150 µl of DNase I reaction mix and mix at 1300 rpm for 5 min at RT, centrifuge shortly and shake at 350 rpm for 15-60 min at 37 °C
reagents['DNaseI_reaction_mix_150'] = {'well': 'A5', 
                               'transfer_volume': 150, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}

#Add 600 µl RNA binding buffer to the digest and mix at 1000 rpm for 10 min
reagents['RNA_binding_buffer'] = {'well': 'A6', 
                               'transfer_volume': 600, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}


#Add 40 µl of nuclease-free water to elute RNA, mix at 1300 rpm for 5 min
reagents['nuclease_free_water'] = {'well': 'A7', 
                               'transfer_volume': 40, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP_MODE}


# In[ ]:


# reagents setup

for reagent_name in reagents:
    reagents[reagent_name]['setup'] = trough.wells(reagents[reagent_name]["well"])


# ## Define custom functions

# In[90]:


def mix_wells(mix_locations, mix_reps):
    """mixes a well thoroughly by aspirating/rejecting liquid at different heights in a well"""

    m300.set_flow_rate(aspirate=300, dispense=550)
    
    for well in mix_locations:
        
        if not m300.tip_attached:
            m300.pick_up_tip()

        m300.move_to(well.top(20), strategy='arc') 
        
        for position in range(20,2,-2):
            m300.aspirate(volume=300, location=well.bottom(position), rate=1.0)
            m300.blow_out(well.bottom(position))

            m300.aspirate(volume=300, location=well.bottom(position), rate=1.0)
            m300.blow_out(well.top(-2))

        m300.move_to(well.top(20), strategy='arc')
    m300.set_flow_rate(aspirate=150, dispense=300)

def resuspend(well_to_mix):
    """resuspends the contents of a well by pipetting liquid up and down while gradually descending into the well"""

    m300.set_flow_rate(aspirate=300, dispense=550)
        
    if not m300.tip_attached:
        m300.pick_up_tip()

    m300.move_to(well_to_mix.top(20), strategy='arc') # fist move to the well
        
        # then aspirate and reject

    for position in range(20,2,-2):
        m300.aspirate(volume=300, location=well_to_mix.bottom(position), rate=1.0)
        m300.blow_out(well_to_mix.bottom(position))

        m300.aspirate(volume=300, location=well_to_mix.bottom(position), rate=1.0)
        m300.blow_out(well_to_mix.top(-2))

    m300.move_to(well_to_mix.top(20), strategy='arc')
        

def transfer_and_mix(reagent, samples):
    for s in samples:

        if not m300.tip_attached:
            m300.pick_up_tip()
        
        m300.transfer(reagent['transfer_volume'], reagent['setup'], s, new_tip=reagent['new_tip'])
        m300.set_flow_rate(aspirate=300, dispense=550)
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)  # note that according to nucleic_acid_extration.ot2.py .mix volume differs from .transfer volume
        m300.set_flow_rate(aspirate=150, dispense=300)
        m300.drop_tip()

        
def trash_supernatant(volume, height, samples):
    """ function to remove [volume in ul] of supernatant from [samples], pipetting [height] units from the bottom of the well"""
    # height to be tested, more or less reliable depending on API version
    for s in samples:
        m300.pick_up_tip()
        m300.transfer(volume, s.bottom(height), m300.trash_container.top(), new_tip='never')
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


# In[99]:



    


# ## Recapitulate setup

# In[ ]:


print("================================ Setup recap ================================")


# ### Instruments

# In[ ]:


# list of (mount, instrument)

print("\nRecap Instruments")
print(robot.get_instruments())


# ### Containers

# In[ ]:


# List all containers on the deck

print("\nRecap containers")
for elt in robot.get_containers():
    print("\t",elt.parent, elt)


# ### Reagents

# In[ ]:


print("\nRecap reagents")
for reagent in reagents:
    print("\t",trough.parent, trough, reagents[reagent]["setup"], reagent)

print("\t",ethanol_plate.parent, ethanol_plate, "<All wells>", "ethanol")


# ### Pipettes

# In[ ]:


# recap attached pipettes
print("\nRecap pipette(s)")
for elt in robot.get_attached_pipettes():
    print("\t", elt,robot.get_attached_pipettes()[elt])


# ### Modules

# In[ ]:


print("\nRecap modules")
print("\tmagdeck", magdeck.status)
print("\ttempdeck", tempdeck.status)


# In[ ]:


print("================================= Run Start =================================")


# In[ ]:


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


# In[ ]:


# home
robot.home()


# In[ ]:


# steps 1-2

# sample collection, nothing to do here


# In[ ]:


# step 3

robot.comment(text_in_a_box("step 3"))
# Add 240 µl of lysis buffer, seal and shake at RT at 1400 rpm for 5 min
transfer_and_mix(reagents['lysis_buffer'], samples)


# In[ ]:


# step 4
robot.comment(text_in_a_box("step 4"))

# Add 320 µl of isopropanol, seal and shake at RT at 1400 rpm for 5 min
transfer_and_mix(reagents['isopropanol_320'], samples)


# ### here, a function to resuspend the beads before dispensing them is needed

# In[91]:


# step 5
robot.comment(text_in_a_box("step 5"))

# resuspend the beads
well_to_mix = reagents['magnetic_beads']["setup"]
resuspend(well_to_mix)


# Add 40 µl of silica-coated magnetic beads
transfer_and_mix(reagents['magnetic_beads'], samples)


# In[ ]:


# step 6
robot.comment(text_in_a_box("step 6"))

# Settle the magnetic beads on a magnetic stand and discard the supernatant
# this block can probably be factorised as a function
# considering the number of times it is used throughout the protocol

robot.comment("Activating magdeck for 5 minutes")
magdeck.engage(height=15)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(minutes=5)


# In[ ]:


# volume & height from bottom to be adjusted based on tests
trash_supernatant(volume=900, height=2, samples=samples)


# In[ ]:


# step 7
robot.comment(text_in_a_box("step 7"))

#Remove the plate from the magnetic stand and add 400 µl isopropanol
# Shake at 1400 rpm for 2 min
magdeck.disengage()
transfer_and_mix(reagents['isopropanol_400'], samples)


# In[ ]:


# step 8
robot.comment(text_in_a_box("step 8"))

#Settle the magnetic beads on a magnetic stand and discard the supernatant
robot.comment("Activating magdeck for 5 minutes")
magdeck.engage(height=15)
if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(minutes=5)


# In[ ]:


# volume & height from bottom of the well are to be adjusted based on tests
trash_supernatant(volume=900, height=2, samples=samples)


# In[ ]:


# steps 9-10-11, repeated 4 times
robot.comment(text_in_a_box("step 9-10-11"))


for rep in range(3):
    
    magdeck.disengage()
    
    for well in samples:
        
        if not m300.tip_attached:
            m300.pick_up_tip()

        
        well_code = str(well).split(" ")[-1][:-1]
        
        m300.transfer(400, 
                      ethanol_plate.wells(well_code).bottom(2), 
                      well, new_tip='never')
        
        m300.set_flow_rate(aspirate=300, dispense=550)
        m300.mix(MIX_REPETITIONS, 300, well) 
        m300.set_flow_rate(aspirate=150, dispense=300)

        m300.drop_tip()


    robot.comment("Activating magdeck for 5 minutes")
    magdeck.engage(height=15)
    
    if test_mode:
        m300.delay(seconds=5)
    else:
        m300.delay(minutes=5)
    
    # volume & height from bottom to be adjusted based on tests
    trash_supernatant(volume=900, height=2, samples=samples)
    


# In[ ]:


samples


# In[ ]:


magdeck.disengage()


# ## kindly request human to move the plate to the temperature module

# In[ ]:


# step 12
robot.comment(text_in_a_box("step 12"))

robot.comment("Please place plate on tempdeck")
robot.pause()


# In[ ]:


#tempdeck.set_temperature(50)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(minutes=10)#tempdeck.set_temperature(25)


# In[ ]:


robot.comment("Please place plate back on magdeck")
robot.pause()


# In[ ]:


# step 13
robot.comment(text_in_a_box("step 13"))

# Remove the plate from the magnets and add 150 µl of DNase I reaction mix
# and mix at 1300 rpm for 5 min at RT, centrifuge shortly and shake 
# at 350 rpm for 15-60 min at 37 °C

if DNAse_incubation:
    transfer_and_mix(reagents['DNaseI_reaction_mix_150'], samples)


# In[ ]:


# step 14
robot.comment(text_in_a_box("step 14"))

#Add 600 µl RNA binding buffer to the digest and mix at 1000 rpm for 10 min
transfer_and_mix(reagents['RNA_binding_buffer'], samples)


# In[ ]:


# step 15
robot.comment(text_in_a_box("step 15"))

robot.comment("Activating magdeck for 5 minutes")
magdeck.engage(height=15)
if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(minutes=5)


# In[ ]:



# volume & height from bottom to be adjusted based on tests
trash_supernatant(volume=900, height=2, samples=samples)


# In[ ]:


# steps 16-17-18, repeated 4 times
robot.comment(text_in_a_box("step 16-17-18"))


for rep in range(3):
    
    magdeck.disengage()
    
    for well in samples:
        
        if not m300.tip_attached:
            m300.pick_up_tip()

        
        well_code = str(well).split(" ")[-1][:-1]
        
        m300.transfer(400, 
                      ethanol_plate.wells(well_code).bottom(2), 
                      well, new_tip='never')
        
        m300.set_flow_rate(aspirate=300, dispense=550)
        m300.mix(MIX_REPETITIONS, 300, well) 
        m300.set_flow_rate(aspirate=150, dispense=300)

        m300.drop_tip()


    robot.comment("Activating magdeck for 5 minutes")
    magdeck.engage(height=15)
    
    if test_mode:
        m300.delay(seconds=5)
    else:
        m300.delay(minutes=5)
    
    # volume & height from bottom to be adjusted based on tests
    trash_supernatant(volume=900, height=2, samples=samples)
    


# ## kindly request human to move the plate to the temperature module

# In[ ]:


# step 19
robot.comment(text_in_a_box("step 19"))

robot.comment("Please place plate on tempdeck")
robot.pause()


# In[ ]:



tempdeck.set_temperature(50)

if test_mode:
    m300.delay(seconds=5)
else:
    m300.delay(minutes=30)

tempdeck.set_temperature(25)


# In[ ]:



robot.comment("Please place plate back on magdeck")
robot.pause()


# In[ ]:


# step 20
robot.comment(text_in_a_box("step 20"))

# Add 40 µl of nuclease-free water to elute RNA, mix at 1300 rpm for 5 min
transfer_and_mix(reagents['nuclease_free_water'], samples)

robot.comment("Done, at last!")


# ### I think I would do the last transfer manually to maximise liquid recovery while minimising the amount of beads

# ## <p style="text-align: center;"> The end </p>
