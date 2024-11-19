from opentrons import simulate
from opentrons import protocol_api, types


# metadata
metadata = {
    'protocolName': 'NMR solution transfer',
    'author': 'User Name <user.name@email.ac.uk>',
    'description': 'Protocol for transfering NMR solutions from plate to tubes',
    'apiLevel': '2.9'
}

def run(protocol: protocol_api.ProtocolContext):
    
    # labware
    tiprack1 = protocol.load_labware('opentrons_96_tiprack_300ul', 4)
    plate = protocol.load_labware('2mlcollection_96_wellplate_2000ul', 5)
    rack = protocol.load_labware('nmr_96_tuberack_1000ul', 6)

    # pipettes
    multi = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks = [tiprack1])

    mecn_vol = 500 #ul
    mecn_height_t0 = mecn_vol/(plate['A1'].length*plate['A1'].width)
    mecn_depth_t0 = mecn_height_t0 - 3


    for well, tube in zip(plate.rows()[0], rack.rows()[0]):
        multi.pick_up_tip()
        #multi.mix(10, 150, location = well.bottom(2), rate = 2)
        
        side_tubes = tube.top().move(types.Point(x=3.8/2, y=0, z=-3))
        for i in range(1): 

            multi.aspirate(250, well.bottom(mecn_depth_t0-mecn_depth_t0*(i+1)/3))
            multi.move_to(tube.top())
            multi.dispense(250, side_tubes).blow_out()
            #multi.blow_out()
        
        multi.drop_tip()

