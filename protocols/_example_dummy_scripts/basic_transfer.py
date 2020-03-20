# source: https://docs.opentrons.com/v2/new_examples.html
# Moving 100 µL from one well to another

from opentrons import protocol_api

metadata = {'apiLevel': '2.0'}

# Moving 100 µL from one well to another:
def run(protocol: protocol_api.ProtocolContext):
    plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 1)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 2)
    p300 = protocol.load_instrument('p300_single', 'right', tip_racks=[tiprack_1])

    p300.transfer(100, plate['A1'], plate['B1'])

    """
    # 
    
    p300.pick_up_tip()
    p300.aspirate(100, plate.wells('A1'))
    p300.dispense(100, plate.wells('A1'))
    p300.return_tip()
    """


# Error in: Traceback (most recent call last): File "opentrons/server/helpers.py", line 32, in run_protocol File "<string>", line 1, in <module> ImportError: cannot import name 'protocol_api'