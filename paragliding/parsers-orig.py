#!/usr/bin/python
# ex:set fileencoding=utf-8:

from __future__ import unicode_literals

import numpy as np
import re
import xml.etree.ElementTree as ET

from itertools import combinations

from datetime import datetime
from datetime import timedelta
from pytz import utc

from .utils import averages
from .utils import moving
from .utils import binom

import logging
logger = logging.getLogger(__name__)


class FlightLog(ET.ElementTree):

    def __init__(self):
        root = ET.Element('kml', xmlns="http://earth.google.com/kml/2.2")
        self.document = ET.SubElement(root, 'Document')
        self.flights = []
        name = ET.SubElement(self.document, 'name')
        name.text = "Flights"
        super(FlightLog, self).__init__(element=root)

    def add_flight(self, flight):
        self.flights.append(flight)

    def make_tree(self):
        for flight in self.flights:
            self.document = flight.make_tree(self.document)


class Flight(object):

    colors_info = {
        '< -5': {
            'values': (-999, -5),
            'color': (255, 99, 71) # Tomato
        },
        '-5 - -1': {
            'values': (-5, -1),
            'color': (255, 165, 0) # Orange
        },
        '-1 - -0.1': {
            'values': (-1, -0.1),
            'color': (216, 196, 31) # Light orange
        },
        '-0.1 - 0.1': {
            'values': (-0.1, 0.1),
            'color': (157, 214, 35) # Yellow
        },
        '0.1 - 0.5': {
            'values': (0.1, 0.5),
            'color': (157, 214, 35) # Green
        },
        '> 0.5': {
            'values': (0.5, 999),
            'color': (31, 168, 58) # Dark green
        },
        'Unknown': {
            'values': (-9999999, -9999999),
            'color': (144, 144, 144) # Gray
        }
    }

    def __init__(self, file_or_filename, name, *args, **kwargs):

        self.location = None
        self.pilot = None
        self.glider = None
        self.date = None

        if name.endswith(".igc"):
            self.name = name.split('.igc')[0]
        else:
            self.name = name

        self.time = []
        self.lat = []
        self.lon = []
        self.gpsheight = []
        self.barheight = []
        self.cart = []

        self.datapoints = 0
        self.distances = {}

        self.read_igc(file_or_filename)

    def __str__(self):
        return self.name or "Trajectory"

    def __repr__(self):
        return "<%s: '%s' at 0x%x>" % (self.__class__.__name__, str(self), id(self))

    def get_color_name(self, delta):
        for name, info in self.colors_info.items():
            values = info['values']
            if delta >= values[0] and delta <= values[1]:
                return name
        return 'Unknown'

    def get_color(self, color_name, alpha=255):
        rgb = self.colors_info[color_name]['color']
        r = rgb[0]
        g = rgb[1]
        b = rgb[2]
        return format(alpha, '02x') + format(b, '02x') + format(g, '02x') + format(r, '02x')

    def read_igc(self, file_or_filename):
        """
        reads igc data into the object
        """

        if hasattr(file_or_filename, "readlines"):
            file_obj = file_or_filename
        else:
            file_obj = open(file_or_filename, "r")

        search = r"B([0-9]{6})([0-9]{6})"
        timedelta_days = 0
        time_old = None
        for line in file_obj.readlines():
            line = line.strip().decode("latin1")
            # line.decode("utf-8").strip()
            coord = re.match(r'B([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{5})(N|S)([0-9]{3})([0-9]{5})(W|E)(A|V)([0-9-]{5})([0-9-]{5})', line)
            if coord:
                a = coord.groups()

                time = timedelta(days=timedelta_days, hours=int(a[0]), minutes=int(a[1]), seconds=int(a[2]))
                if time_old and time < time_old:
                    time += timedelta(1)
                    timedelta_days += 1
                time_old = time

                lat = float(a[3]) + float(a[4]) / 60000
                if a[5] == "S":
                    lat *= -1

                lon = float(a[6]) + float(a[7]) / 60000
                if a[8] == "W":
                    lon *= -1

                barheight = int(a[10])
                gpsheight = int(a[11])

                self.time.append(self.date + time)
                self.lat.append(lat)
                self.lon.append(lon)
                self.barheight.append(barheight)
                self.gpsheight.append(gpsheight)
                continue

            if "HPSITSITE" in line:
                self.location = line.split(':',1)[1].strip()
                continue

            if "HOPLTPILOT" in line:
                self.pilot = line.split(':',1)[1].strip()
                continue

            if "HOGTYGLIDERTYPE" in line:
                self.glider = line.split(':',1)[1].strip()
                continue

            date = re.match(r"HFDTE([0-9]{2})([0-9]{2})([0-9]{2})", line)
            if date:
                d,m,y = date.groups()
                self.date = datetime(2000+int(y), int(m), int(d), tzinfo=utc)
                continue

            # TODO only log unparsed entries from igc file
            logger.debug(line.strip())

        self.lon = np.array(self.lon)
        self.lat = np.array(self.lat)
        self.barheight = np.array(self.barheight)
        self.gpsheight = np.array(self.gpsheight)

        self.datapoints = len(self.lon)
        self.datarange = np.arange(self.datapoints)

        file_obj.close()

    def make_tree(self, root):

        root_folder = ET.SubElement(root, 'Folder')
        name = ET.SubElement(root_folder, 'name')
        name.text = str(self)

        colored_folder = ET.SubElement(root_folder, 'Folder')
        data = ET.SubElement(colored_folder, 'name')
        data.text = "Colored"

        shadow_folder = ET.SubElement(root_folder, 'Folder')
        data = ET.SubElement(shadow_folder, 'name')
        data.text = "Shadow"

        # SHADOW

        marker = ET.SubElement(shadow_folder, 'Placemark')
        data = ET.SubElement(marker, 'name')
        style = ET.SubElement(marker, 'Style', id="ShadowLine")
        style = ET.SubElement(style, 'LineStyle')
        data = ET.SubElement(style, 'color')
        data.text = '48000000'
        data = ET.SubElement(style, 'width')
        data.text = '2.0'
        coordinates = ET.SubElement(marker, 'LineString')
        data = ET.SubElement(coordinates, 'tessellate')
        data.text = '1'
        data = ET.SubElement(coordinates, 'coordinates')
        data.text = ' '.join([
            "%.8f,%.8f,%d" % (
                self.lon[d],
                self.lat[d],
                1,
            )
            for d in range(len(self.time))
        ])

        smooth_height = averages(self.gpsheight, 20, binom)
        # speed = np.gradient(t.cart[0])**2 + np.gradient(t.cart[1])**2 + np.gradient(t.cart[2])**2 - np.gradient(t.gpsheight)**2
        # speed *= speed > 0
        # speed = np.sqrt(speed)*3.6

        # Track (color)

        for d in range(1, len(self.time)):
            delta = smooth_height[d] - smooth_height[d-1]
            color_name = self.get_color_name(delta)

            flight = ET.SubElement(colored_folder, 'Placemark')

            data = ET.SubElement(flight, 'name')
            # data.text = "%s | %s m | %s km/h | %.1f m/s" % ("time", "alt", "speed", delta)
            data.text = "%s m/s" % color_name
            data = ET.SubElement(flight, 'Style')
            style = ET.SubElement(data, 'LineStyle')
            data = ET.SubElement(style, 'color')
            data.text = self.get_color(color_name)
            data = ET.SubElement(style, 'width')
            data.text = "2.5"
            coord = ET.SubElement(flight, 'LineString')
            data = ET.SubElement(coord, 'tessellate')
            data.text = "1"
            data = ET.SubElement(coord, 'altitudeMode')
            data.text = "absolute"
            data = ET.SubElement(coord, 'coordinates')
            data.text = "%.8f,%.8f,%d %.8f,%.8f,%d"% (
                self.lon[d-1],
                self.lat[d-1],
                self.gpsheight[d-1],
                self.lon[d],
                self.lat[d],
                self.gpsheight[d],
            )

        return root

    def max_turning_points(self, iterator, coords=None, distance=-1, count=0):
        changed = False
        keep = set()
        data = set()
        for test_coords in iterator:

            for i in test_coords:
                data.add(i)

            new_distance = self.calc_turning_point_distance(test_coords)
            if new_distance > distance:
                coords = test_coords
                changed = True
                distance = new_distance
            elif new_distance > (0.95 + 0.02 * count) * distance:
                for i in test_coords:
                    keep.add(i)

        return coords, distance, keep, data, changed

    def calc_turning_points(self, points, guess=0, maxiter=20):
        count = 0
        if guess < points + 2:
            guess = points + 2
        coords, distance, keep, data, changed = self.max_turning_points(
            combinations(
                np.uint16(np.arange(guess) * (self.datapoints - 1) / (guess - 1)), points + 2
            )
        )
        keep = set()
        data = set()

        while changed and count < maxiter:
            iter_coords = keep.copy()
            for n in range(len(coords)):

                j = coords[n]
                if n == 0:
                    i = 0
                    k = coords[n + 1]
                elif n == (len(coords) - 1):
                    i = coords[n - 1]
                    k = self.datapoints - 1
                else:
                    i = coords[n - 1]
                    k = coords[n + 1]

                dist_i = self.get_distances(i)
                dist_j = self.get_distances(j)
                dist_k = self.get_distances(k)

                iter_coords.add(j)

                n_ik = (np.argmax((dist_i[i:k+1] + dist_k[i:k+1])) + i)
                iter_coords.add(n_ik)

                if n == 0:
                    n_k = np.argmax(dist_k[0:k+1])
                    iter_coords.add(n_k)
                elif n == (len(coords) - 1):
                    n_i = np.argmax(dist_i[i:k + 1]) + i
                    iter_coords.add(n_i)

                if i < j:
                    n_ij = (np.argmax((dist_i[i:j+1] + dist_j[i:j+1])) + i)
                    iter_coords.add(n_ij)

                if j < k:
                    n_jk = (np.argmax((dist_j[j:k+1] + dist_k[j:k+1])) + j)
                    iter_coords.add(n_jk)

            if not bool(iter_coords - data):
                # no new items to process
                logger.info("no new items to process")
                break

            coords, distance, keep, data, changed = self.max_turning_points(
                combinations(sorted(iter_coords), points + 2), coords, distance, count
            )
            count += 1

        return self.calc_turning_point_distance(coords), coords

    def calc_turning_point_distance(self, coords):
        d = 0
        for n in range(len(coords) - 1):
            i = coords[n]
            j = coords[n+1]
            d += self.get_distances(i)[j]
        return d

    def get_distances(self, i):
        if i in self.distances.keys():
            return self.distances[i]
        self.distances[i] = np.zeros(self.datapoints, np.float64)
        for n in range(self.datapoints):
            self.distances[i][n] = self.calc_distance(i, n)
        return self.distances[i]

    def calc_FAI_distance(self, i, j):
        # FAI earth-radius in meter
        R = 6371000.0

        latx = np.radians(self.lat[i])
        lonx = np.radians(self.lon[i])

        laty = np.radians(self.lat[j])
        lony = np.radians(self.lon[j])

        sinlat = np.sin((latx-laty)/2)
        sinlon = np.sin((lonx-lony)/2)

        return 2 * R * np.arcsin(np.sqrt(
            sinlat * sinlat + sinlon * sinlon * np.cos(latx) * np.cos(laty)
        ))

    def calc_distance(self, i, j):
        return self.calc_FAI_distance(i, j)

        # # WGS84 - http://de.wikipedia.org/wiki/Erdellipsoid
        # a = 6378137.
        # n = 298.257223563 # = 1/f = a/(a-b)
        # b = a*(1-1/n)
        # e = np.sqrt(a**2 - b**2)/a # numerische Exzentrizität
        # N = a/np.sqrt(1-e**2*np.sin(np.radians(self.lat))**2) # Krümmungsradius des Ersten Vertikals

        # self.cart = (
        #     (N+self.gpsheight)*np.cos(np.radians(self.lat))*np.cos(np.radians(self.lon)),
        #     (N+self.gpsheight)*np.cos(np.radians(self.lat))*np.sin(np.radians(self.lon)),
        #     (N*(1-e**2)+self.gpsheight)*np.sin(np.radians(self.lat)),
        # )
