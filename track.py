#!/usr/bin/python
# ex:set fileencoding=utf-8:

"""
python track.py files ed
"""

from __future__ import unicode_literals

import argparse
import os

from paragliding.parsers import FlightLog
from paragliding.parsers import Flight


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a directory with igc files')
    parser.add_argument('dir', type=str, help='directory')
    parser.add_argument('output', type=str, help='output filename')
    args = parser.parse_args()

    flights = FlightLog()
    for path, dirs, files in os.walk(args.dir):
        for file in files:
            if file[-4:] == ".igc":
                f = Flight(os.path.join(path, file), file)
                flights.add_flight(f)
    flights.make_tree()
    flights.write(args.output + '.kml')
