from pprint import pprint
from opentrons import protocol_api
import json
from pathlib import Path
from collections import defaultdict
from opentrons.types import Point

# metadata
metadata = {
    "protocolName": "Solvent Benchmarking Script",
    "author": "User Name <user.name@email.ac.uk>",
    "description": "Solvent Benchmarking",
    "apiLevel": "2.9",
}

substance_locations = {
    "5": {
        "name": "Solvent",
        "type": "fisher_6_wellplate_25000ul",
        "A1": {
            "substance": "Test_Solvent",
            "amount": 10000
        },
        "A2": {
            "substance":"Water",
            "amount": 5000
        }
    }
}

move_commands = [
    {
        "substance": "Test_Solvent",
        "amount": 75,
        "plate": 4,
        "location": [
            "A1",
            "A2",
            "A3"
        ]
    },
    {
        "substance": "Water",
        "amount": 10,
        "plate": 4,
        "location": [
            "F1"
        ]
    },  
    {
        "substance": "Test_Solvent",
        "amount": 85,
        "plate": 4,
        "location": [
            "A4",
            "A5",
            "A6"
        ]
    },
    {
        "substance": "Water",
        "amount": 10,
        "plate": 4,
        "location": [
            "F1"
        ]
    }, 
    {
        "substance": "Test_Solvent",
        "amount": 150,
        "plate": 4,
        "location": [
            "B1",
            "B2",
            "B3"
        ]
    },
    {
        "substance": "Water",
        "amount": 10,
        "plate": 4,
        "location": [
            "F1"
        ]
    }, 
    {
        "substance": "Test_Solvent",
        "amount": 160,
        "plate": 4,
        "location": [
            "B4",
            "B5",
            "B6"
        ]
    },
    {
        "substance": "Water",
        "amount": 10,
        "plate": 4,
        "location": [
            "F1"
        ]
    },
    {
        "substance": "Test_Solvent",
        "amount": 250,
        "plate": 4,
        "location": [
            "C1",
            "C2",
            "C3"
        ]
    },
    {
        "substance": "Water",
        "amount": 10,
        "plate": 4,
        "location": [
            "F1"
        ]
    },
    {
        "substance": "Test_Solvent",
        "amount": 260,
        "plate": 4,
        "location": [
            "C4",
            "C5",
            "C6"
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
        # Load labware
        self.tiprack300 = self.protocol.load_labware(
            "opentrons_96_tiprack_300ul", 2
        )
        self.plate = self.protocol.load_labware(
            "analyticalsales_48_wellplate_2000ul", 4
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
        pipette.dispense(location=position_to.top(z=2)).blow_out(location=position_to.top(z=2))

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

        # Perform swelling
        if pipette.swelled != substance_name:
            self.swell_tip(
                pipette=pipette,
                position=substance_position,
            )
            pipette.swelled = substance_name
            # Check to see if amount of substance is greater than the amount on the deck
        minimum_volume = 500
        if amount > (
            self.substances[substance_name][0]["amount"] - minimum_volume
        ):
            print(
                f"Amount of {substance_name} needed is greater than the amount on the deck.\n"
                f"Trying again to move after changing the well plate to movw from of {substance_name}."
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
        # Move the substance
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
    substance_path = substance_locations
    move_info = move_commands
    amount_added = defaultdict(lambda: 0)

    ot = Opentrons(protocol=protocol, deck_info=substance_path)

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
    pipette_max_volume = 285
    for move in moves:
        location = move["location"]
        amount = int(move["amount"])
        plate = int(move["plate"])
        substance = move["substance"]
        # Get target well location
        target_well = ot.labware[plate].wells(location)[0]
        # Add the first substance to the list
        num_moves = amount // pipette_max_volume
        #if amount >= 230:
        #    print(num_moves)
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
