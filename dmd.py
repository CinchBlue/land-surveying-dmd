import csv
import math
import argparse
import re
import logging
import sys
import time

def classlogger(class_):
    attribute_name = '_{}__log'.format(class_.__name__)
    setattr(class_, attribute_name, logging.getLogger(class_.__module__ + '.' + class_.__name__))
    return class_

class Line:
    """Represents a basic Line object."""

    def __init__(self, azimuth, distance):
        self.azimuth = azimuth
        self.distance = distance

    def __repr__(self):
        return "Line({0},{1})".format(
            self.azimuth, self.distance)


@classlogger
class DegMinSec:
    """
    Plain-Old-Data (POD) object for the tuple of (Degrees, Minutes, Seconds).
    Azimuth is a vector value, so the DMS values should be all positive,
    where the sign captures the 
    """

    def __init__(self, degrees=0.0, minutes=0.0, seconds=0.0, sign="+"):
        self.degrees = degrees
        self.minutes = minutes
        self.seconds = seconds
        self.sign = "+"

    def __repr__(self):
        return 'DegMinSec({0},{1},{2},"{3}")'.format(
            self.degrees, self.minutes, self.seconds, self.sign)

    def __str__(self):
        return '{}{:.0f}°{:.0f}′{}″'.format(
            self.sign, self.degrees, self.minutes, '{:f}'.format(self.seconds).rstrip('0').rstrip('.'))

    @classmethod
    def from_str(cls, input_string):
        """
        The format is 'DDD.MMSSssss', but everything except the first D is optional.
        We could do this using regex, but it's more work than it's worth here.
        """
        cls.__log.debug("input: input_string")
        dms_matcher = re.compile(
            r'^([+-])?([0-9]{1,3})?\.?([0-9]{1,2})?([0-9]+)?$')
        matched_groups = re.match(dms_matcher, input_string)
        if matched_groups is None:
            raise ValueError(
                '`{}` is not a valid DMS azimuth (note: not decimal degrees) formatted string.'.format(input_string))
        (sign, d, m, s) = matched_groups.group(1, 2, 3, 4)
        cls.__log.debug('matched: %s %s %s %s', sign, d, m, s)
        sign = '+' if sign is None else sign
        d = 0.0 if d is None else float(d)
        m = 0.0 if m is None else float(m)
        s = 0.0 if s is None else (float(s) if len(
            s) <= 2 else (float(s[0:2] + "." + s[2:])))
        cls.__log.debug('post-conversion: %s %f %f %f', sign, d, m, s)
        return DegMinSec(d, m, s, sign)



class InnerRadialCurve:
    def __init__(self, azimuth_in, radius, azimuth_out):
        self.azimuth_in = azimuth_in
        self.azimuth_out = azimuth_out
        self.radius = radius

    def __repr__(self):
        return "InnerRadialCurve({0},{1},{2})".format(
            self.azimuth_in, self.radius, self.azimuth_out)

    def to_radial_line_tuple(self):
        return (
            Line(self.azimuth_in, self.radius),
            Line(self.azimuth_out. self.radius))



def perform_dmd_calculation(obj_list):
    tracking_dep, tracking_lat = 0.0, 0.0
    error_dep, error_lat = 0.0, 0.0
    for item in obj_list:
        if item is Line:
            dep = tracking_dep + math.cos
    pass

def output_report(results):
    pass

def construct_argument_parser():
    parser = argparse.ArgumentParser(
    description='Calculates double-meridian distance (DMD).')
    parser.add_argument('input_filename', metavar='<input>', type=str,
                        help='input CSV file')
    parser.add_argument('output_filename', metavar='<output>', type=str,
                        help='output file', default='out.csv')
    parser.add_argument('-v', '--verbosity', metavar='<verbosity>', type=str,
                        help='limits program diagnostic output level', default='INFO',
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser

def construct_application_logger(args):
    logging.basicConfig(
        format='%(asctime)s UTC [%(levelname)s] %(name)s.%(funcName)s: %(message)s',
        stream=sys.stderr,
        level=getattr(logging, args.verbosity.upper(), None))
    logging.Formatter.converter = time.gmtime
    log = logging.getLogger(__name__)
    return log

def read_input_as_object_list(args):
    object_list = []
    # Get input pairs:
    # azimuth (line-out/radial-in), distance/radius, azimuth2 (radial-out)
    with open(args.input_filename, newline='') as input_file:
        data_reader = csv.reader(input_file, delimiter=',', quotechar='\"')
        # Read in the data. DMD closures should never be prohibitively large (e.g.
        # not over 1000 points) so holding it in memory should be okay.
        for row in data_reader:
            # If the row has 2 entries, it should be [azimuth, distance]
            if len(row) == 2:
                azimuth = DegMinSec.from_str(row[0])
                distance = float(row[1])
                object_list.append(Line(azimuth, distance))
            # If the length of the individual row is 3, this means that
            # that the optional radial-out azimuth for simple radial curves is present.
            #
            # Interpret the current row as a simple radial curve.
            elif len(row) == 3:
                azimuth_in = DegMinSec.from_str(row[0])
                radius = float(row[1])
                azimuth_out = DegMinSec.from_str(row[2])
                object_list.append(InnerRadialCurve(azimuth_in, radius, azimuth_out))

    # Once the data is read, return it.
    return object_list


def main():
    parser = construct_argument_parser()
    args = parser.parse_args()

    log = construct_application_logger(args)

    log.debug(DegMinSec.from_str("123."))
    log.debug(DegMinSec.from_str("1"))
    log.debug(DegMinSec.from_str("."))
    log.debug(DegMinSec.from_str(".0"))
    log.debug(DegMinSec.from_str(".0123"))
    log.debug(DegMinSec.from_str(".01234567"))
    log.debug(DegMinSec.from_str("-.01234567"))

    obj_list = read_input_as_object_list(args)
    calculation_stats = perform_dmd_calculation(obj_list)

    log.debug(obj_list)


if __name__ == "__main__":
    main()



# input: begin->radial, radial->begin
# delta = radial->begin
# inner curve angle will never be over 180. It can be exactly 180.
# points are always done clockwise.