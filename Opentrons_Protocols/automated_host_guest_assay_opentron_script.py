from pprint import pprint
from opentrons import protocol_api
import json
from pathlib import Path
from collections import defaultdict
from opentrons.types import Point

# metadata
metadata = {
    "protocolName": "Automated MOC host:guest assay",
    "author": "User Name <user.name@email.ac.uk>",
    "description": "Automated MOC host:guest assay with cage 1 (host) and 32 guests, with subset of cage 2 and 3 (hosts) and 4 guests.",
    "apiLevel": "2.9",
}

substance_locations = {
    "3": {
        "name": "Stock",
        "type": "analyticalsales_24_wellplate_80000ul",
        "A1": {
            "substance": "G1",
            "amount": 5000
        },
        "A2": {
            "substance": "G2",
            "amount": 5000
        },
        "A3": {
            "substance": "G3",
            "amount": 5000
        },
        "A4": {
            "substance": "G4",
            "amount": 5000
        },
        "A5": {
            "substance": "G5",
            "amount": 5000
        },
        "A6": {
            "substance": "G6",
            "amount": 5000
        },
        "B1": {
            "substance": "G7",
            "amount": 5000
        },
        "B2": {
            "substance": "G8",
            "amount": 5000
        },
        "B3": {
            "substance": "G9",
            "amount": 5000
        },
        "B4": {
            "substance": "G10",
            "amount": 5000
        },
        "B5": {
            "substance": "G11",
            "amount": 5000
        },
        "B6": {
            "substance": "G12",
            "amount": 5000
        },
        "C1": {
            "substance": "G13",
            "amount": 5000
        },
        "C2": {
            "substance": "G14",
            "amount": 5000
        },
        "C3": {
            "substance": "G15",
            "amount": 5000
        },
        "C4": {
            "substance": "G16",
            "amount": 5000
        },
        "C5": {
            "substance": "G17",
            "amount": 5000
        },
        "C6": {
            "substance": "G18",
            "amount": 5000
        },
        "D1": {
            "substance": "G19",
            "amount": 5000
        },
        "D2": {
            "substance": "G20",
            "amount": 5000
        },
        "D3": {
            "substance": "G21",
            "amount": 5000
        },
        "D4": {
            "substance": "G22",
            "amount": 5000
        },
        "D5": {
            "substance": "G23",
            "amount": 5000
        },
        "D6": {
            "substance": "G24",
            "amount": 5000
        }
    },
    "5": {
        "name": "Stock2",
        "type": "analyticalsales_24_wellplate_80000ul",
        "A1": {
            "substance": "G25",
            "amount": 5000
        },
        "A2": {
            "substance":"G26",
            "amount": 5000
        },
        "A3": {
            "substance": "G27",
            "amount": 5000
        },
        "A4": {
            "substance": "G28",
            "amount": 5000
        },
        "A5": {
            "substance": "G29",
            "amount": 5000
        },
        "A6": {
            "substance": "G30",
            "amount": 5000
        },
        "B1": {
            "substance": "G31",
            "amount": 5000
        },
        "B2": {
            "substance": "G32",
            "amount": 5000
        },
        "B3": {
            "substance": "cage2",
            "amount": 5000 # 5 mL
        },
        "B4": {
            "substance": "cage3",
            "amount": 5000 # 5 mL
        },
    },
    "6": {
        "name": "Stock3",
        "type": "fisher_6_wellplate_25000ul",
        "A1": {
            "substance": "cage1",
            "amount": 20000 # 20 mL
        },
    },
}



move_commands = [
    {
        "substance": "cage1",
        "amount": 520,
        "plate": 4,
        "location": [
            "A1", "A2", "A3", "A4",
            "B1", "B2", "B3", "B4",
            "C1", "C2", "C3", "C4",
            "D1", "D2", "D3", "D4",
            "E1", "E2", "E3", "E4",
            "F1", "F2", "F3", "F4",
            "G1", "G2", "G3", "G4",
            "H1", "H2", "H3", "H4",
        ]
    },
        {
        "substance": "cage2",
        "amount": 520,
        "plate": 4,
        "location": [
            "A5", 
            "B5", 
            "C5", 
            "D5", 
            "E5",
        ]
    },
        {
        "substance": "cage3",
        "amount": 520,
        "plate": 4,
        "location": [
            "F5", 
            "G5", 
            "H5", 
            "A6", 
            "B6",
        ]
    },    
    {
        "substance": "G1",
        "amount": 85,
        "plate": 4,
        "location": [
            "A1",
        ]
    },
    {
        "substance": "G2",
        "amount": 85,
        "plate": 4,
        "location": [
            "A2", 
        ]
    },
    {
        "substance": "G3",
        "amount": 85,
        "plate": 4,
        "location": [
            "A3",
        ]
    },
    {
        "substance": "G4",
        "amount": 85,
        "plate": 4,
        "location": [
            "A4", 
        ]
    },
    {
        "substance": "G5",
        "amount": 85,
        "plate": 4,
        "location": [
            "B1",
        ]
    },
    {
        "substance": "G6",
        "amount": 85,
        "plate": 4,
        "location": [
            "B2",
        ]
    },
    {
        "substance": "G7",
        "amount": 85,
        "plate": 4,
        "location": [
            "B3",
        ]
    },
    {
        "substance": "G8",
        "amount": 85,
        "plate": 4,
        "location": [
            "B4",
        ]
    },
    {
        "substance": "G9",
        "amount": 85,
        "plate": 4,
        "location": [
            "C1", 
        ]
    },
    {
        "substance": "G10",
        "amount": 85,
        "plate": 4,
        "location": [
            "C2",
        ]
    },
    {
        "substance": "G11",
        "amount": 85,
        "plate": 4,
        "location": [
            "C3",
        ]
    },
    {
        "substance": "G12",
        "amount": 85,
        "plate": 4,
        "location": [
            "C4",
            "A5",
            "F5"
        ]
    },
    {
        "substance": "G13",
        "amount": 85,
        "plate": 4,
        "location": [
            "D1",
        ]
    },
    {
        "substance": "G14",
        "amount": 85,
        "plate": 4,
        "location": [
            "D2",
        ]
    },
    {
        "substance": "G15",
        "amount": 85,
        "plate": 4,
        "location": [
            "D3",
            "B5",
            "G5"
        ]
    },
    {
        "substance": "G16",
        "amount": 85,
        "plate": 4,
        "location": [
            "D4",
            "C5",
            "H5"
        ]
    },
    {
        "substance": "G17",
        "amount": 85,
        "plate": 4,
        "location": [
            "E1",
        ]
    },
    {
        "substance": "G18",
        "amount": 85,
        "plate": 4,
        "location": [
            "E2",
        ]
    },
    {
        "substance": "G19",
        "amount": 85,
        "plate": 4,
        "location": [
            "E3", 
        ]
    },
    {
        "substance": "G20",
        "amount": 85,
        "plate": 4,
        "location":
        [
            "E4",
            "D5",
            "A6",
        ]
    },
    {
        "substance": "G21",
        "amount": 85,
        "plate": 4,
        "location": [
            "F1",
        ]
    },
    {
        "substance": "G22",
        "amount": 85,
        "plate": 4,
        "location": [
            "F2",
        ]
    },
    {
        "substance": "G23",
        "amount": 85,
        "plate": 4,
        "location": [
            "F3",
        ]
    },
    {
        "substance": "G24",
        "amount": 85,
        "plate": 4,
        "location": [
            "F4", 
        ]
    },
    {
        "substance": "G25",
        "amount": 85,
        "plate": 4,
        "location": [
            "G1", 
        ]
    },
    {
        "substance": "G26",
        "amount": 85,
        "plate": 4,
        "location": [
            "G2",
        ]
    },
    {
        "substance": "G27",
        "amount": 85,
        "plate": 4,
        "location": [
            "G3", 
        ]
    },
    {
        "substance": "G28",
        "amount": 85,
        "plate": 4,
        "location": [
            "G4", 
        ]
    },
    {
        "substance": "G29",
        "amount": 85,
        "plate": 4,
        "location": [
            "H1", 
        ]
    },
    {
        "substance": "G30",
        "amount": 85,
        "plate": 4,
        "location": [
            "H2", 
        ]
    },
    {
        "substance": "G31",
        "amount": 85,
        "plate": 4,
        "location": [
            "H3", 
        ]
    },
    {
        "substance": "G32",
        "amount": 85,
        "plate": 4,
        "location": [
            "H4",
        ]
    },
]


class Opentrons:
    def __init__(
        self,
        protocol: protocol_api.ProtocolContext,
        deck_info: dict,
    ):
        """
        Initialise the Opentrons object.

        Parameters
        ----------
        protocol : protocol_api.ProtocolContext
            The protocol context for the current protocol.
        deck_info : dict
            Dictionary of deck information.
            Contains information of substances and their quantities on the deck.

        Returns
        -------
        None
        """
        # Set gantry speeds
        self.protocol = protocol
        self.protocol.max_speeds["X"] = 250
        self.protocol.max_speeds["Y"] = 250
        self.protocol.max_speeds["Z"] = 250

        # TODO: Add dynamic loading of labware
        self.tiprack300 = self.protocol.load_labware(
            "opentrons_96_tiprack_300ul", 2
        )
        self.plate = self.protocol.load_labware(
            "2mlcollection_96_wellplate_2000ul", 4 # 96 well plate
        )
        self.labware = {}
        self.labware[4] = self.plate
        self.labware[2] = self.tiprack300
        # Add substances to the deck
        for deck_number in deck_info:
            self.labware[int(deck_number)] = self.protocol.load_labware(
                deck_info[deck_number]["type"], int(deck_number)
            )
        # Delete name and type from deck_info
        for deck_number in deck_info:
            del deck_info[deck_number]["type"]
            del deck_info[deck_number]["name"]
        # Loading pipettes

        self.left_pipette = protocol.load_instrument(
            "p300_single_gen2", "left", tip_racks=[self.tiprack300]
        )
        self.left_pipette.flow_rate.aspirate = 70
        self.left_pipette.flow_rate.dispense = 70
        self.left_pipette.swelled = False
        
        # For tracking substance amounts
        self.substances = defaultdict(list)
        for position in deck_info:
            position_int = int(position)
            for well_plate in deck_info[position]:
                substance_position = self.labware[position_int][well_plate]
                substance = deck_info[position][well_plate]
                substance_name = substance["substance"]
                amount = substance["amount"]
                self.substances[substance_name].append(
                    {
                        "position": substance_position,
                        "amount": amount,
                    }
                )

    def swell_tip(self, pipette, position):
        """
        Swells the tip in the `stock` labware location at a specified location.

        Notes
        -----
        Perform before a pipette is used for transfer.

        Parameters
        ----------
        pipette : pipette
            The pipette to be used.

        position:
            Position to swell the tip in.

        Returns
        -------
        None

        """
        for i in range(1):
            pipette.aspirate(100, position)
            self.protocol.delay(10)
            pipette.move_to(position.top())
            pipette.dispense(100, location=position)
        pipette.swelled = True

    def move_without_drip(self, position_to, position_from, pipette, amount):
        """
        Transfers substance from one location to another without dripping (hopefully).

        Notes
        -----
        Ideally, the swell function will be used before this function is called to reduce the
        probability of drips.

        Parameters
        ----------
        position_to: position
            Location of the target well plates to move substance to.
        position_from: position
            Location of the source well plates to move substance from.
        amount: float or int
            Amount of substance to be moved. (in uL)


        Returns
        -------
        None

        """
        # Check if pipette swelled before movement
        pipette.aspirate(amount, position_from)
        pipette.air_gap(15)
        pipette.dispense(location=position_to.top(z=1)).blow_out(location=position_to.top(z=1))

    def move_substance(self, amount, substance_name, position_to, pipette):
        """
        Moves a specified amount of substance from one location to another.

        Parameters
        ----------
        amount : float or int
            Amount of substance to be moved. (in mL)
        substance_name : str
            Name of the substance to be moved.
        position_to: position
            Location of the target well plates to move substance to.
        pipette: pipette
            The pipette to be used.
        value:

        Returns
        -------
        None

        """
        # Find the deck location of the substance
        try:
            substance_position = self.substances[substance_name][0]["position"]
        except:
            pprint(self.substances)
            print(f"Substance {substance_name} not found.")
            raise RuntimeError(
                f"Substance not found in deck. This could be due to the deck being empty, or the amount of {substance_name} needed exceeding the amount placed on the deck."
            )

        # Perform swelling for left pipette
        if pipette.swelled != substance_name:
            self.swell_tip(
                pipette=pipette,
                position=substance_position,
            )
            pipette.swelled = substance_name
            # Check to see if amount of substance is greater than the amount on the deck
        minimum_volume = 250
        if amount > (
            self.substances[substance_name][0]["amount"] - minimum_volume
        ):
            print(
                f"Amount of {substance_name} needed is greater than the amount on the deck.\n"
                f"Trying again by moving to next well containing {substance_name} and aspirating from there."
            )
            # Change the well plate to move from
            self.substances[substance_name].pop(0)
            if len(self.substances[substance_name]) == 0:
                raise RuntimeError(
                    f"No more {substance_name} left on the deck. Please check the deck and try again."
                )
            self.move_substance(
                amount=amount,
                substance_name=substance_name,
                position_to=position_to,
                pipette=pipette,
            )

        assert pipette.swelled == substance_name
        # Move substance
        self.move_without_drip(
            pipette=pipette,
            position_to=position_to,
            position_from=substance_position,
            amount=amount,
        )
        self.substances[substance_name][0]["amount"] -= amount


def run(protocol: protocol_api.ProtocolContext):
    """
    Run the protocol.
    """
    amount_added = defaultdict(lambda: 0)
    deck_info = substance_locations
    move_info = move_commands

    ot = Opentrons(protocol=protocol, deck_info=deck_info)

    # Extract list of moves from the move information
    moves = []
    for substance in move_info:
        for location in substance["location"]:
            moves.append(
                {
                    "substance": substance["substance"],
                    "location": location,
                    "amount": substance["amount"],
                    "plate": substance["plate"],
                }
            )
    added_substances = []
    positions_added = defaultdict(str)
    pipette_max_volume = 250
    for move in moves:
        location = move["location"]
        amount = int(move["amount"])
        plate = int(move["plate"])
        substance = move["substance"]
        # Get target well location
        target_well = ot.labware[plate].wells(location)[0]
        # Add the first substance to the list
        num_moves = amount // pipette_max_volume
        added = 0
        if len(added_substances) == 0:
            ot.left_pipette.pick_up_tip()
            added_substances.append(substance)
        # Check if substance matches the last substance added
        elif substance != added_substances[-1]:
            # Get a new tip
            ot.left_pipette.drop_tip()
            ot.left_pipette.pick_up_tip()
            added_substances.append(substance)
        for i in range(num_moves + 1):
            # Check if move is the last move
            if i == num_moves:
                amount_to_add = amount - added
            else:
                amount_to_add = pipette_max_volume
            if amount_to_add == 0:
                continue
            ot.move_substance(
                amount=amount_to_add,
                substance_name=substance,
                position_to=target_well,
                pipette=ot.left_pipette,
            )
            added += amount_to_add
            amount_added[location] += amount_to_add
        positions_added[location] += substance + " "
    for line in protocol.commands():
        continue
    pprint(amount_added)
    pprint(positions_added)
