#!/usr/bin/env python

# ATHAD Python driver.
#
# The inherited version drove BOTH an Atmosphere and a Hydrosphere, printed a
# "Paleo Ocean Circulations" banner and reported bathymetry_path — and, at the bottom,
# actually ran the ocean rather than the atmosphere. None of that applies: ATHAD has no
# hydrosphere, no bathymetry and no time slices.

from pycond import Atmosphere


class Model(object):
    """
    ATHAD Model Object — Atmosphere of the Earth in the Hadean Eon.
    """
    def __init__(self):
        self.atm = Atmosphere()
        self.config_cond_xml = "config_cond.xml"

    def print_config(self):
        print("\n\n\n\n   ATHAD — Atmosphere of the Earth in the Hadean Eon")
        print("\n\n   configuration file name is                    ", self.config_cond_xml)
        print("   config path is                                ", self.atm.config_xml_path.decode('utf-8'))
        print("   output path for results is                    ", self.atm.output_path.decode('utf-8'))

    def run(self):
        print("\n   running ATHAD (single epoch, ~4.4 Ga)")
        self.atm.run()
        print("\n   ATHAD terminated successfully")
        print("\n")


if __name__ == '__main__':
    m = Model()
    m.print_config()
    m.run()
