# get a list of the reagents and their concentrations and volumes
# get a list of all the equipment

# to run from bash command line:
# > opentrons_simulate rna_extraction.py

# based on nucleic_acid_extraction.ot2.py from https://protocol-delivery.protocols.opentrons.com/protocol/1584

from collections import OrderedDict
from opentrons import labware, instruments, modules, robot

# TODO
# figure out how to get LYSIS_BUFFER_NEWTIP='never' to work!
# figure out how to tell machine to stop, remember where it is and to continue after command

# comments
# transfer_volume is in μl
# new_tip options: 'always', 'never', 'once'

NEW_TIP = 'never'
MIX_REPETITIONS = 15

reagents = OrderedDict()
# Add 240 μl of lysis buffer, seal and shake at RT at 1400 rpm for 5 min
reagents['lysis_buffer'] = {'well': 'A1', 'transfer_volume': 240, 'mix_volume': 240, 'mix_repetitions': MIX_REPETITIONS,
                            'new_tip': NEW_TIP}
# Add 320 μl of isopropanol, seal and shake at RT at 1400 rpm for 5 min
reagents['isopropanol_320'] = {'well': 'A2', 'transfer_volume': 300, 'mix_volume': 300, 'mix_repetitions': MIX_REPETITIONS,
                            'new_tip': NEW_TIP}
# Add 40 μl of silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock), seal and shake at RT at 1400 rpm for 5 min
reagents['magnetic_beads'] = {'well': 'A3', 'transfer_volume': 40, 'mix_volume': 40, 'mix_repetitions': MIX_REPETITIONS,
                            'new_tip': NEW_TIP}
# Remove the plate from the magnetic stand and add 400 μl isopropanol. Shake at RT at 1400 rpm for 2 min
reagents['isopropanol_400'] = {'well': 'A2', 'transfer_volume': 300, 'mix_volume': 300, 'mix_repetitions': MIX_REPETITIONS,
                            'new_tip': NEW_TIP}

metadata = {
    'protocolName': 'RNA Extraction v0.1',
    'author': 'Neil MacKenzie, Eyal Kazin <eyalkazin@gmail.com>',
    'source': 'Testing' #'Custom Protocol Request'
}

# create custom labware
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

tips = [labware.load('opentrons-tiprack-300ul', str(slot)) for slot in range(4, 10)]
print('-' * 50)
print('Make sure the `tips` make sense!')
print(tips)
print('-' * 50)

# modules
magdeck = modules.load('magdeck', '1')
sample_plate = labware.load(plate_name, '1', share=True)


# instruments
m300 = instruments.P300_Multi(mount='right', tip_racks=tips)


# reagent setup

for reagent_name in reagents:
    reagents[reagent_name]['setup'] = trough.wells('A1')


#lysis_buffer = trough.wells('A1')
# isopropanol = trough.wells('A2')
# magnetic_bead = trough.wells('A3')   # e.g, silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock)
# ethanol_80percent = trough.wells('A4')
# dnaseI_reaction_mix = trough.wells('A5')  # enzyme that removes DNA
# rna_binding_buffer = trough.wells('A6')
# nuclease_free_water = trough.wells('A7')
# liquid_waste = trough.wells('A12')  #  elusion waste

def transfer_and_mix(reagent, samples):
    for s in samples:
        m300.pick_up_tip()
        m300.transfer(reagent['transfer_volume'], reagent['setup'], s, new_tip=reagent['new_tip'])
        m300.mix(reagent['mix_repetitions'], reagent['mix_volume'], s)  # note that according to nucleic_acid_extration.ot2.py .mix volume differs from .transfer volume
        m300.drop_tip()

def run_custom_protocol(number_of_sample_columns: int = 12):
    if number_of_sample_columns > 12:
        raise Exception("Please specify a valid number of sample columns.")

    samples = sample_plate.rows('A')[0:number_of_sample_columns]

    # transfer lysis and mix (bursting the cells open)
    transfer_and_mix(reagents['lysis_buffer'], samples)

    # --- OpenTron does not have Seal and Shake modules ---
    # figure out how to tell machine to stop, remember where it is and to continue after command
    # transfer and mix isopropanol

    transfer_and_mix(reagents['isopropanol_320'], samples)

    transfer_and_mix(reagents['magnetic_beads'], samples)

    transfer_and_mix(reagents['isopropanol_400'], samples)






run_custom_protocol(**{'number_of_sample_columns': 2})
