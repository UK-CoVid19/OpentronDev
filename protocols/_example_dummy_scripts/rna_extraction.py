# get a list of the reagents and their concentrations and volumes
# get a list of all the equipment

# to run from bash command line:
# > opentrons_simulate rna_extraction.py

from opentrons import labware, instruments, modules, robot

# TODO
# figure out how to get LYSIS_BUFFER_NEWTIP='never' to work!
# figure out how to tell machine to stop, remember where it is and to continue after command

LYSIS_BUFFER_VOLUME = 240  # μl
LYSIS_BUFFER_NEWTIP = 'never' # other options: 'never', 'once'
LYSIS_MIX_REPETITIONS = 15

metadata = {
    'protocolName': 'RNA Extraction v0.1',
    'author': 'Neil MacKenzie, Eyal Kazin <protocols@opentrons.com>',
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
lysis_buffer = trough.wells('A1')
isopropanol = trough.wells('A2')
magnetic_bead = trough.wells('A3')   # e.g, silica-coated magnetic beads (BOMB protocol #2.1, 1:10 diluted from stock)
ethanol_80percent = trough.wells('A4')
dnaseI_reaction_mix = trough.wells('A5')  # enzyme that removes DNA
rna_binding_buffer = trough.wells('A6')
nuclease_free_water = trough.wells('A7')
liquid_waste = trough.wells('A12')  #  elusion waste

#def transfer_and_mix():

def run_custom_protocol(number_of_sample_columns: int = 12):
    if number_of_sample_columns > 12:
        raise Exception("Please specify a valid number of sample columns.")

    samples = sample_plate.rows('A')[0:number_of_sample_columns]

    # transfer lysis and mix (bursting the cells open)
    for s in samples:
        m300.pick_up_tip()
        m300.transfer(LYSIS_BUFFER_VOLUME, lysis_buffer, s, new_tip=LYSIS_BUFFER_NEWTIP)
        m300.mix(LYSIS_MIX_REPETITIONS, LYSIS_BUFFER_VOLUME, s)  # note that according to nucleic_acid_extration.ot2.py .mix volume differs from .transfer volume
        m300.drop_tip()

    # --- OpenTron does not have Seal and Shake modules ---
    # figure out how to tell machine to stop, remember where it is and to continue after command

    # transfer and mix isopropanol




run_custom_protocol(**{'number_of_sample_columns': 12})
