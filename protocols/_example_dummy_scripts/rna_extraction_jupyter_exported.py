#!/usr/bin/env python
# coding: utf-8

# ## To do list
# 
# * get a list of the reagents and their concentrations and volumes
# * get a list of all the equipment

# ## To simulate protocol
# 
# export file to rna_extraction.py manually then in terminal:<br>
# $ opentrons_simulate rna_extraction.py
# 
# OR go to the last two cells that automatically export to *.py then simulate it and print the output

# ## Resources
# 
# 
# following protocol from [BOMB total RNA extraction mammalian GITC v1.0.pdf](https://bomb.bio/wp-content/uploads/2018/09/8.2_BOMB_total_RNA_extraction_mammalian_GITC_V1.0.pdf)
# 
# 
# Snippets taken from nucleic_acid_extraction.ot2.py [see here](https://protocol-delivery.protocols.opentrons.com/protocol/1584)
# 

# In[1]:


metadata = {
    'protocolName': 'RNA Extraction v0.2',
    'author': 'Neil MacKenzie, Eyal Kazin <eyalkazin@gmail.com>, Aubin Fleiss <afleiss@ic.ac.uk>',
    'source': 'Testing' #'Custom Protocol Request'
}


# In[2]:


# intalling opentrons module (only once)
#import sys
#!{sys.executable} -m pip install opentrons 


# In[3]:


# updating opentrons module (not really needed)
#import sys
#!{sys.executable} -m pip install --upgrade opentrons


# In[4]:


# import standard modules
from collections import OrderedDict
from opentrons import labware, instruments, modules, robot


# ## To do list
# 
# * **figure out how to get LYSIS_BUFFER_NEWTIP='never' to work!**
# I have this somewhere in my scripts, I will have a look
# 
# * **figure out how to tell machine to stop, remember where it is and to continue after command** robot.pause() brings the machine to a halt, it is the same as clicking "Pause" in the app. The protocol then has to be resumed manually via the app (click "Resume"). To just introduce a delay use \[pipette instance].delay(opt) so for example here that would be m300.delay(minutes=5)
# 

# 
# ## Comments
# * transfer_volume is in μl
# * new_tip options: 'always', 'never', 'once'
# 

# In[5]:



# create custom labware
# is there a particular reason for using this custom labware ?
# I feel safer using deep well plates
# (since they are deep there is less risk of cross contamination by projection in my experience)
plate_name = 'MidSci-96-Well'
if plate_name not in labware.list():
    labware.create(
        plate_name,
        grid=(12, 8),
        spacing=(9, 9),
        diameter=5,
        depth=21,
        volume=200
    )

# labware
trough = labware.load('trough-12row', '2', 'trough')
fresh_plate = labware.load(plate_name, '3', 'fresh plate')

tips = [labware.load('opentrons-tiprack-300ul', str(slot)) for slot in range(5, 10)]
print('-' * 50)
print('Make sure the `tips` make sense!')
print(tips)
print('-' * 50)

# modules
magdeck = modules.load('magdeck', '1')
sample_plate = labware.load(plate_name, '1', share=True)

tempdeck = modules.load('tempdeck', '4')
tempdeck.set_temperature(25)



# instruments
m300 = instruments.P300_Multi(mount='right', tip_racks=tips)


# In[6]:


# here some volumes seem to be "capped" to 300 while
# some volumes in the BOMB protocol are larger than 300
# Shall we use the BOMB protocol values ? When using m300.transfer()
# the robot (surprisingly) knows how to split the large volume to 
# avoid flooding the pipette. 
# It is not much smarter than that though :)

# I do not remember how much volume there is per well in a through-12
# plate. I can check in the lab how much it is
# (I fear we might need more wells for isopro and ethanol for instance)

#lysis_buffer = trough.wells('A1')
# isopropanol = trough.wells('A2')
# magnetic_bead = trough.wells('A3')   # e.g, silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock)
# ethanol_80percent = trough.wells('A4')
# dnaseI_reaction_mix = trough.wells('A5')  # enzyme that removes DNA
# rna_binding_buffer = trough.wells('A6')
# nuclease_free_water = trough.wells('A7')
# liquid_waste = trough.wells('A12')  #  elusion waste

# the last comment line above (liquid_waste=...)
# is a GREAT idea. maybe if we trash the supernatants in a spare deep 
# well plate we can better avoid cross-contamination
# I think I will implement this in all my protocols :)

NEW_TIP = 'never'
MIX_REPETITIONS = 15

reagents = OrderedDict()
# Add 240 μl of lysis buffer, seal and shake at RT at 1400 rpm for 5 min
reagents['lysis_buffer'] = {'well': 'A1', 
                            'transfer_volume': 240,
                            'mix_volume': 240, 
                            'mix_repetitions': MIX_REPETITIONS,
                            'new_tip': NEW_TIP}

# Add 320 μl of isopropanol, seal and shake at RT at 1400 rpm for 5 min
reagents['isopropanol_320'] = {'well': 'A2', 
                               'transfer_volume': 320,
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP}

# Add 40 μl of silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock), seal and shake at RT at 1400 rpm for 5 min
reagents['magnetic_beads'] = {'well': 'A3', 
                              'transfer_volume': 40, 
                              'mix_volume': 40, 
                              'mix_repetitions': MIX_REPETITIONS,
                              'new_tip': NEW_TIP}

# Remove the plate from the magnetic stand and add 400 μl isopropanol. Shake at RT at 1400 rpm for 2 min
reagents['isopropanol_400'] = {'well': 'A2', 
                               'transfer_volume': 400, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP}


# Add 150 µl of DNase I reaction mix and mix at 1300 rpm for 5 min at RT, centrifuge shortly and shake at 350 rpm for 15-60 min at 37 °C
reagents['DNaseI_reaction_mix_150'] = {'well': 'A5', 
                               'transfer_volume': 150, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP}

#Add 600 µl RNA binding buffer to the digest and mix at 1000 rpm for 10 min
reagents['RNA_binding_buffer'] = {'well': 'A6', 
                               'transfer_volume': 600, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP}

#Add 600 µl RNA binding buffer to the digest and mix at 1000 rpm for 10 min
reagents['ethanol_80percent_400'] = {'well': 'A4', 
                               'transfer_volume': 600, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP}


#Add 40 µl of nuclease-free water to elute RNA, mix at 1300 rpm for 5 min
reagents['nuclease_free_water'] = {'well': 'A7', 
                               'transfer_volume': 40, 
                               'mix_volume': 300, 
                               'mix_repetitions': MIX_REPETITIONS,
                               'new_tip': NEW_TIP}


# In[7]:


# reagent setup

for reagent_name in reagents:
    reagents[reagent_name]['setup'] = trough.wells('A1')


# In[8]:


# I have written a function for "vortexing" wells
# It pipettes up/down the liquid at various (or even random)
# heights in the well, ensuring everything is homogeneously mixed
# I will look for it and add it so maybe we can use it in the function
# transfer_and_mix what do you think?


def transfer_and_mix(reagent, samples):
    for s in samples:
        m300.pick_up_tip()
        m300.transfer(reagent['transfer_volume'], reagent['setup'], s, new_tip=reagent['new_tip'])
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)  # note that according to nucleic_acid_extration.ot2.py .mix volume differs from .transfer volume
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


# In[9]:


# def run_custom_protocol(number_of_sample_columns: int = 12):
#run_custom_protocol(**{'number_of_sample_columns': 2})

# I like having a main function running the whole protocol
# but it is much more difficult to test each instruction individually

# I will split each line of the protocol into different cells for now



# ### number_of_sample_columns is basically the only parameter the user will change once the protocol is set up

# In[10]:


number_of_sample_columns=1


# In[11]:


if number_of_sample_columns > 12:
    raise Exception("Please specify a valid number of sample columns.")
    

samples = sample_plate.rows('A')[0:number_of_sample_columns]


# In[12]:


# steps 1-2 : sample collection


# In[13]:


# step 3

# Add 240 µl of lysis buffer, seal and shake at RT at 1400 rpm for 5 min
transfer_and_mix(reagents['lysis_buffer'], samples)


# In[14]:



# --- OpenTron does not have Seal and Shake modules ---
# figure out how to tell machine to stop, remember where it is and to continue after command
# transfer and mix isopropanol


# In[15]:


# step 4

# Add 320 µl of isopropanol, seal and shake at RT at 1400 rpm for 5 min
transfer_and_mix(reagents['isopropanol_320'], samples)


# In[16]:


# step 5

# Add 40 µl of silica-coated magnetic beads
transfer_and_mix(reagents['magnetic_beads'], samples)


# In[17]:


# step 6

#Settle the magnetic beads on a magnetic stand and discard the supernatant
# this block can probably be factorised, 
# considering the number of times it is used throughout the protocol
robot.comment("Activating magdeck for 5 minutes")
magdeck.engage(height=15)
m300.delay(minutes=5)# may need to increase to let the beads settle


# In[18]:


# volume & height from bottom to be adjusted based on tests
trash_supernatant(volume=900, height=2, samples=samples)


# In[19]:


# step 7

#Remove the plate from the magnetic stand and add 400 µl isopropanol
# Shake at 1400 rpm for 2 min
magdeck.disengage()
transfer_and_mix(reagents['isopropanol_400'], samples)


# In[20]:


# step 8

#Settle the magnetic beads on a magnetic stand and discard the supernatant
robot.comment("Activating magdeck for 5 minutes")
magdeck.engage(height=15)
m300.delay(minutes=5)# may need to increase to let the beads settle


# In[21]:


# volume & height from bottom to be adjusted based on tests
trash_supernatant(volume=900, height=2, samples=samples)


# In[22]:


# steps 9-10-11, repeated 4 times

for rep in range(3):
    
    magdeck.disengage()
    
    transfer_and_mix(reagents['isopropanol_400'], samples)
    
    #Settle the magnetic beads on a magnetic stand and discard the supernatant
    robot.comment("Activating magdeck for 5 minutes")
    magdeck.engage(height=15)
    m300.delay(minutes=5)# may need to increase to let the beads settle
    
    # volume & height from bottom to be adjusted based on tests
    trash_supernatant(volume=900, height=2, samples=samples)
    
magdeck.disengage()


# ## kindly request human to move the plate to the temperature module

# In[23]:


# step 12

robot.comment("Please place plate on tempdeck")
robot.pause()


# In[24]:


tempdeck.set_temperature(50)
m300.delay(minutes=10) # may need to adjust to let the beads dry
tempdeck.set_temperature(25)


# In[25]:


robot.comment("Please place plate back on magdeck")
robot.pause()


# In[26]:


# step 13

# Remove the plate from the magnets and add 150 µl of DNase I reaction mix
# and mix at 1300 rpm for 5 min at RT, centrifuge shortly and shake 
# at 350 rpm for 15-60 min at 37 °C
transfer_and_mix(reagents['DNaseI_reaction_mix_150'], samples)


# In[27]:


# step 14

#Add 600 µl RNA binding buffer to the digest and mix at 1000 rpm for 10 min
transfer_and_mix(reagents['RNA_binding_buffer'], samples)


# In[28]:


# step 15

robot.comment("Activating magdeck for 5 minutes")
magdeck.engage(height=15)
m300.delay(minutes=5)# may need to increase to let the beads settle
    
# volume & height from bottom to be adjusted based on tests
trash_supernatant(volume=900, height=2, samples=samples)


# In[29]:


# steps 16-17-18, repeated 4 times

for rep in range(3):
    
    magdeck.disengage()
    
    transfer_and_mix(reagents['ethanol_80percent_400'], samples)
    
    #Settle the magnetic beads on a magnetic stand and discard the supernatant
    robot.comment("Activating magdeck for 5 minutes")
    magdeck.engage(height=15)
    m300.delay(minutes=5)# may need to increase to let the beads settle
    
    # volume & height from bottom to be adjusted based on tests
    trash_supernatant(volume=900, height=2, samples=samples)
    
magdeck.disengage()


# ## kindly request human to move the plate to the temperature module

# In[30]:


# step 19

robot.comment("Please place plate on tempdeck")
robot.pause()


# In[31]:



tempdeck.set_temperature(50)
m300.delay(minutes=30) # may need to adjust to let the beads dry
tempdeck.set_temperature(25)


# In[32]:



robot.comment("Please place plate back on magdeck")
robot.pause()


# In[33]:


# step 20

# Add 40 µl of nuclease-free water to elute RNA, mix at 1300 rpm for 5 min
transfer_and_mix(reagents['nuclease_free_water'], samples)


# ### I think I would do the last transfer manually to maximise liquid recovery while minimising the amount of beads

# ## <p style="text-align: center;"> The end </p>

# ### export and run the whole protocol in the opentrons simulator

# In[ ]:




