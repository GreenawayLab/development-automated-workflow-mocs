from pprint import pprint
from opentrons import protocol_api
import json
from pathlib import Path
from collections import defaultdict
from opentrons.types import Point

# metadata
metadata = {
    "protocolName": "Automated MOC synthesis screen refined",
    "author": "User Name <user.name@email.ac.uk>",
    "description": "MOC Metal and Ligand Screen for ligands A1, A3 and metal salts Zn(NTf2)2 and Zn(BF4)2",
    "apiLevel": "2.9",
}

substance_locations = {
    "5": {
        "name": "Stock",
        "type": "analyticalsales_24_wellplate_80000ul",
        "A1": {
            "substance": "TriamineA",
            "amount": 5000
        },
        "A2": {
            "substance": "TriamineA",
            "amount": 5000
        },
        "A3": {
            "substance": "Ald1",
            "amount": 5000
        },
        "A4": {
            "substance": "Ald1",
            "amount": 5000
        },
        "A5": {
            "substance": "Ald3",
            "amount": 5000
        },
        "A6": {
            "substance": "Ald3",
            "amount": 5000
        },
        "B1": {
            "substance": "Zn(NTf2)2",
            "amount": 5000
        },
        "B2": {
            "substance": "Zn(NTf2)2",
            "amount": 5000
        },
        "B3": {
            "substance": "Zn(BF4)2",
            "amount": 5000
        },
        "B4": {
            "substance": "Zn(BF4)2",
            "amount": 5000
        },
    },
    "6": {
        "name": "Stock2",
        "type": "fisher_6_wellplate_25000ul",
        "A1": {
            "substance": "MeCN",
            "amount": 20000
        },
        "A2": {
            "substance":"MeCN",
            "amount": 20000
        },
    }    
}

move_commands = [
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "A1"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "A3"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "A5"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "A7"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "B1"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "B3"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "B5"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "B7"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "C1"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "C3"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "C5"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "C7"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "D1"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "D3"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "D5"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "D7"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "E1"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "E3"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "E5"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "E7"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "F1"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "F3"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "F5"
        ]
    },
    {
        "substance": "Ald1",
        "amount": 122,
        "plate": 4,
        "location": [
            "F7"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "A2"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "A4"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "A6"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "A8"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "B2"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "B4"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "B6"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "B8"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "C2"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "C4"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "C6"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "C8"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "D2"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "D4"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "D6"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "D8"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "E2"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "E4"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "E6"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "E8"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "F2"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "F4"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "F6"
        ]
    },
    {
        "substance": "Ald3",
        "amount": 138,
        "plate": 4,
        "location": [
            "F8"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 745,
        "plate": 4,
        "location": [
            "A1"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 729,
        "plate": 4,
        "location": [
            "A2"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 745,
        "plate": 4,
        "location": [
            "A3"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 729,
        "plate": 4,
        "location": [
            "A4"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 678,
        "plate": 4,
        "location": [
            "A5"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 662,
        "plate": 4,
        "location": [
            "A6"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 678,
        "plate": 4,
        "location": [
            "A7"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 662,
        "plate": 4,
        "location": [
            "A8"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 698,
        "plate": 4,
        "location": [
            "B1"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 682,
        "plate": 4,
        "location": [
            "B2"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 719,
        "plate": 4,
        "location": [
            "B3"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 703,
        "plate": 4,
        "location": [
            "B4"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 631,
        "plate": 4,
        "location": [
            "B5"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 615,
        "plate": 4,
        "location": [
            "B6"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 652,
        "plate": 4,
        "location": [
            "B7"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 636,
        "plate": 4,
        "location": [
            "B8"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 651,
        "plate": 4,
        "location": [
            "C1"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 635,
        "plate": 4,
        "location": [
            "C2"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 693,
        "plate": 4,
        "location": [
            "C3"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 677,
        "plate": 4,
        "location": [
            "C4"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 584,
        "plate": 4,
        "location": [
            "C5"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 568,
        "plate": 4,
        "location": [
            "C6"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 626,
        "plate": 4,
        "location": [
            "C7"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 610,
        "plate": 4,
        "location": [
            "C8"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 603,
        "plate": 4,
        "location": [
            "D1"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 587,
        "plate": 4,
        "location": [
            "D2"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 667,
        "plate": 4,
        "location": [
            "D3"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 651,
        "plate": 4,
        "location": [
            "D4"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 536,
        "plate": 4,
        "location": [
            "D5"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 520,
        "plate": 4,
        "location": [
            "D6"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 600,
        "plate": 4,
        "location": [
            "D7"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 584,
        "plate": 4,
        "location": [
            "D8"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 556,
        "plate": 4,
        "location": [
            "E1"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 540,
        "plate": 4,
        "location": [
            "E2"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 640,
        "plate": 4,
        "location": [
            "E3"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 624,
        "plate": 4,
        "location": [
            "E4"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 489,
        "plate": 4,
        "location": [
            "E5"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 473,
        "plate": 4,
        "location": [
            "E6"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 573,
        "plate": 4,
        "location": [
            "E7"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 557,
        "plate": 4,
        "location": [
            "E8"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 508,
        "plate": 4,
        "location": [
            "F1"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 492,
        "plate": 4,
        "location": [
            "F2"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 614,
        "plate": 4,
        "location": [
            "F3"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 598,
        "plate": 4,
        "location": [
            "F4"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 441,
        "plate": 4,
        "location": [
            "F5"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 425,
        "plate": 4,
        "location": [
            "F6"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 547,
        "plate": 4,
        "location": [
            "F7"
        ]
    },
    {
        "substance": "MeCN",
        "amount": 531,
        "plate": 4,
        "location": [
            "F8"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "A1"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "A2"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "A3"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "A4"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "A5"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "A6"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "A7"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "A8"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "B1"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "B2"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "B3"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "B4"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "B5"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "B6"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "B7"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "B8"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "C1"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "C2"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "C3"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "C4"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "C5"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "C6"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "C7"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "C8"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "D1"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "D2"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "D3"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "D4"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "D5"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "D6"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "D7"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "D8"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "E1"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "E2"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "E3"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "E4"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "E5"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "E6"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "E7"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "E8"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "F1"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "F2"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "F3"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 133,
        "plate": 4,
        "location": [
            "F4"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "F5"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "F6"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "F7"
        ]
    },
    {
        "substance": "TriamineA",
        "amount": 200,
        "plate": 4,
        "location": [
            "F8"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 26,
        "plate": 4,
        "location": [
            "B3"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 26,
        "plate": 4,
        "location": [
            "B4"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 26,
        "plate": 4,
        "location": [
            "B7"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 26,
        "plate": 4,
        "location": [
            "B8"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 52,
        "plate": 4,
        "location": [
            "C3"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 52,
        "plate": 4,
        "location": [
            "C4"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 52,
        "plate": 4,
        "location": [
            "C7"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 52,
        "plate": 4,
        "location": [
            "C8"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 78,
        "plate": 4,
        "location": [
            "D3"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 78,
        "plate": 4,
        "location": [
            "D4"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 78,
        "plate": 4,
        "location": [
            "D7"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 78,
        "plate": 4,
        "location": [
            "D8"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 105,
        "plate": 4,
        "location": [
            "E3"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 105,
        "plate": 4,
        "location": [
            "E4"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 105,
        "plate": 4,
        "location": [
            "E7"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 105,
        "plate": 4,
        "location": [
            "E8"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 131,
        "plate": 4,
        "location": [
            "F3"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 131,
        "plate": 4,
        "location": [
            "F4"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 131,
        "plate": 4,
        "location": [
            "F7"
        ]
    },
    {
        "substance": "Zn(BF4)2",
        "amount": 131,
        "plate": 4,
        "location": [
            "F8"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 47,
        "plate": 4,
        "location": [
            "B1"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 47,
        "plate": 4,
        "location": [
            "B2"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 47,
        "plate": 4,
        "location": [
            "B5"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 47,
        "plate": 4,
        "location": [
            "B6"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 94,
        "plate": 4,
        "location": [
            "C1"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 94,
        "plate": 4,
        "location": [
            "C2"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 94,
        "plate": 4,
        "location": [
            "C5"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 94,
        "plate": 4,
        "location": [
            "C6"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 142,
        "plate": 4,
        "location": [
            "D1"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 142,
        "plate": 4,
        "location": [
            "D2"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 142,
        "plate": 4,
        "location": [
            "D5"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 142,
        "plate": 4,
        "location": [
            "D6"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 189,
        "plate": 4,
        "location": [
            "E1"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 189,
        "plate": 4,
        "location": [
            "E2"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 189,
        "plate": 4,
        "location": [
            "E5"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 189,
        "plate": 4,
        "location": [
            "E6"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 237,
        "plate": 4,
        "location": [
            "F1"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 237,
        "plate": 4,
        "location": [
            "F2"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 237,
        "plate": 4,
        "location": [
            "F5"
        ]
    },
    {
        "substance": "Zn(NTf2)2",
        "amount": 237,
        "plate": 4,
        "location": [
            "F6"
        ]
    }
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
        self.tiprack20 = self.protocol.load_labware(
            "opentrons_96_tiprack_20ul", 1
        )
        self.tiprack300 = self.protocol.load_labware(
            "opentrons_96_tiprack_300ul", 2
        )
        self.plate = self.protocol.load_labware(
            "analyticalsales_48_wellplate_2000ul", 4
        )
        self.labware = {}
        self.labware[4] = self.plate
        self.labware[1] = self.tiprack20
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
        pipette.dispense(location=position_to.top(z=2))
        pipette.blow_out(location=position_to.top(z=2))

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
        minimum_volume = 600
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
    pipette_max_volume = 270
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
