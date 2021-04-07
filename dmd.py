import csv
import math
import argparse
import re
import logging
import sys
import time

parser = argparse.ArgumentParser(
    description='Calculates double-meridian distance (DMD).')
parser.add_argument('input_filename', metavar='<input>', type=str,
                    help='input CSV file')
parser.add_argument('output_filename', metavar='<output>', type=str,
                    help='output file', default='out.csv')
parser.add_argument('-v', '--verbosity', metavar='<verbosity>', type=str,
                    help='limits program diagnostic output level', default='INFO',
                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
args = parser.parse_args()

logging.basicConfig(
    format='%(asctime)s UTC [%(levelname)s] %(name)s.%(funcName)s: %(message)s',
    stream=sys.stderr,
    level=getattr(logging, args.verbosity.upper(), None))
logging.Formatter.converter = time.gmtime
logger = logging.getLogger("dmd")


class Line:
    def __init__(self, azimuth, distance):
        self.azimuth = azimuth
        self.distance = distance


class DegMinSec:
    def __init__(self, degrees=0, minutes=0, seconds=0, sign="+"):
        self.degrees = degrees
        self.minutes = minutes
        self.seconds = seconds
        self.sign = "+"

    def __repr__(self):
        return "DegMinSec({0},{1},{2},'{3}')".format(
            self.degrees, self.minutes, self.seconds, self.sign)


class SimpleRadialCurve:
    def __init__(self, azimuth_in, radius, azimuth_out):
        self.azimuth_in = azimuth_in
        self.azimuth_out = azimuth_out
        self.radius = radius

# The format is 'DDD.MMSSssss', but everything except the first D is optional.
# We could do this using regex, but it's more work than it's worth here.


def str_to_azimuth(input_string):
    logger.debug("input: input_string")
    dms_matcher = re.compile(
        r'^([+-])?([0-9]{1,3})?\.?([0-9]{1,2})?([0-9]+)?$')
    matched_groups = re.match(dms_matcher, input_string)
    if matched_groups is None:
        raise ValueError(
            '"{}" is not a valid DMS azimuth (note: not decimal degrees) formatted string.'.format(input_string))
    (sign, d, m, s) = matched_groups.group(1, 2, 3, 4)
    logger.debug('matched: %s %s %s %s', sign, d, m, s)
    sign = '+' if sign is None else sign
    d = 0.0 if d is None else float(d)
    m = 0.0 if m is None else float(m)
    s = 0.0 if s is None else (float(s) if len(
        s) <= 2 else (float(s[0:2] + "." + s[2:])))
    logger.debug('post-conversion: %s %f %f %f', sign, d, m, s)
    return DegMinSec(d, m, s, sign)


logger.debug(str_to_azimuth("123."))
logger.debug(str_to_azimuth("1"))
logger.debug(str_to_azimuth("."))
logger.debug(str_to_azimuth(".0"))
logger.debug(str_to_azimuth(".0123"))
logger.debug(str_to_azimuth(".01234567"))
logger.debug(str_to_azimuth("-.01234567"))


# Get input pairs:
# azimuth (line-out/radial-in), distance/radius, azimuth2 (radial-out)
with open(args.input_filename, newline='') as input_file:
    data_reader = csv.reader(input_file, delimiter=',', quotechar='\"')
    # Read in the data. DMD closures should never be prohibitively large (e.g.
    # not over 1000 points) so holding it in memory should be okay.
    for row in data_reader:
        # If the row has 2 entries, it should be [azimuth, distance]
        if len(row) == 2:
            row.append('0')
        # If the length of the individual row is 3, this means that
        # that the optional radial-out azimuth for simple radial curves is present.
        #
        # Interpret the current row as a simple radial curve.
        elif len(row) == 3:
            pass
    # Once the data is read, do the DMD calcuation.
    # Then, output the report.


def perform_dmd_calculation(formatted_data):
    pass


def output_report(results):
    pass
