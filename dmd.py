"""Double-Meridian Distance calculator

NOTE: You need to do the chord area like this:

Do DMD area with the CHORD LINE derived from the radial in/out lines.
Then add the SEGMENT (from CHORD to ARC) AREA.

A = (0.5) * r * r * (theta - sin theta)
theta = theta [radians]
"""

import csv
import math
import argparse
import re
import logging
import sys
import time
import copy


def classlogger(class_):
    """This is a convenience function that can be applied to a class as a "decorator" in order to add
    a `logging` module logger to it. 
    """
    attribute_name = '_{}__log'.format(class_.__name__)
    setattr(class_, attribute_name, logging.getLogger(
        class_.__module__ + '.' + class_.__name__))
    return class_


@classlogger
class AzimuthDMS:
    """
    Plain-Old-Data (POD) object for the tuple of (Degrees, Minutes, Seconds) representing an azimuth.
    Azimuth is a vector value, so the DMS values should be all absolute value,
    where the `self.sign` field captures the sign if necessary.
    """

    def __init__(self, degrees=0.0, minutes=0.0, seconds=0.0, sign="+"):
        self.degrees = degrees
        self.minutes = minutes
        self.seconds = seconds
        self.sign = sign
        self.normalize()

    def __repr__(self):
        return 'AzimuthDMS({0},{1},{2},"{3}")'.format(
            self.degrees, self.minutes, self.seconds, self.sign)

    def __str__(self):
        return '{}{:.0f}°{:.0f}\'{}\"'.format(
            self.sign, self.degrees, self.minutes, '{:f}'.format(self.seconds).rstrip('0').rstrip('.'))

    def normalize(self):
        """Normalizes the azimuth to less than 360 degrees and positive sign."""

        self.degrees = self.degrees % 360
        if self.sign == '-':
            # Convert to decimal degrees, and subtract from 360.0
            new_degrees = 360.0 - self.to_degrees()
            # Run the conversion from decimal degrees
            new_azimuth = AzimuthDMS.from_degrees(new_degrees)
            # Adjust self.
            self.degrees = new_azimuth.degrees
            self.minutes = new_azimuth.minutes
            self.seconds = new_azimuth.seconds
            self.sign = new_azimuth.sign

        return self

    @classmethod
    def from_str(cls, input_string):
        """Converts a string into an AzimuthDMS.

        The input format is 'DDD.MMSSssss', but everything except the first D is optional.
        To make this much less difficult to code, we use regular expressions to match exactly as we want.
        This is also a @classmethod in order to allow us to use the cls.__log for logging when we want.
        """

        dms_matcher = re.compile(
            r'^([+-])?([0-9]{1,3})?\.?([0-9]{1,2})?([0-9]+)?$')
        matched_groups = re.match(dms_matcher, input_string)
        if matched_groups is None:
            raise ValueError(
                '`{}` is an invalid DMS azimuth formatted string (note: not decimal degrees).'.format(input_string))
        (sign, d, m, s) = matched_groups.group(1, 2, 3, 4)
        sign = '+' if sign is None else sign
        d = 0.0 if d is None else float(d)
        m = 0.0 if m is None else float(m)
        # This is needed in the case where you have 'DDD.M'. The M is assumed to be the shortened form of 'DDD.M0'.
        m = 10.0*m if m < 10 and s is None else m
        s = 0.0 if s is None else (float(s) if len(
            s) <= 2 else (float(s[0:2] + "." + s[2:])))

        return AzimuthDMS(d, m, s, sign)

    def to_degrees(self):
        return self.degrees + self.minutes/60.0 + self.seconds/3600.0

    def to_radians(self):
        return math.radians(self.to_degrees())

    @staticmethod
    def from_degrees(decimal_degrees):
        (degrees_fraction, whole) = math.modf(decimal_degrees)
        degrees = whole % 360.0
        (minutes_fraction, minutes) = math.modf(degrees_fraction*60.0)
        seconds = minutes_fraction * 60.0
        return AzimuthDMS(degrees, minutes, seconds)

    @staticmethod
    def from_radians(radians):
        return AzimuthDMS.from_degrees(math.degrees(radians))

    def __add__(self, other):
        return AzimuthDMS.from_degrees(self.to_degrees() + other.to_degrees())

    def __sub__(self, other):
        return AzimuthDMS.from_degrees(self.to_degrees() - other.to_degrees())

    def __neg__(self):
        temp = copy.deepcopy(self)
        if temp.sign == '+':
            temp.sign = '-'
        else:
            temp.sign = '+'
        return temp.normalize()

    def to_normalized_lat_dep(self):
        """Converts AzimuthDMS into a normalized vector tuple of (latitude, departure)."""

        latitude = math.cos(self.to_radians())
        departure = math.sin(self.to_radians())
        return (latitude, departure)

    def reverse(self):
        return AzimuthDMS(degrees=self.degrees+180.0, minutes=self.minutes, seconds=self.seconds, sign=self.sign)


@classlogger
class Line:
    """
    Represents a basic Line object, with a AzimuthDMS azimuth.
    """

    def __init__(self, azimuth, distance):
        self.azimuth = azimuth
        self.distance = distance

    def __repr__(self):
        return "Line({0},{1})".format(
            self.azimuth, self.distance)

    def to_lat_dep(self):
        (norm_lat, norm_dep) = self.azimuth.to_normalized_lat_dep()
        return (self.distance * norm_lat, self.distance * norm_dep)

    def to_lat(self):
        return to_lat_dep[0]

    def to_dep(self):
        return to_lat_dep[1]

    def to_dmd_area_terms(self, last_dmd, last_dep):
        # Get the current line's latitude and depature
        (curr_lat, curr_dep) = self.to_lat_dep()
        # Calculate the DMD of the current line
        curr_dmd = last_dmd + last_dep + curr_dep
        # Calculate the current DMD area term
        curr_area = (curr_dmd * curr_lat)/2.0

        return (curr_dmd, curr_area)

@classlogger
class InnerRadialCurve:
    """Represents a simple radial curve less than 180 degrees.

    There are many ways to canonically represent a radial curve, but we choose to
    represent it with two azimuths going into and out of the center point,
    as well as the constant radius of the curve.
    """

    def __init__(self, azimuth_in, azimuth_out, radius):
        self.azimuth_in = azimuth_in
        self.azimuth_out = azimuth_out
        self.radius = radius

    def __repr__(self):
        return "InnerRadialCurve(azimuth_in={},azimuth_out={},radius={})".format(
            self.azimuth_in, self.azimuth_out, self.radius)

    def __str__(self):
        return "InnerRadialCurve(azimuth_in={},azimuth_out={},radius={},delta:{})".format(
            self.azimuth_in, self.azimuth_out, self.radius, self.get_delta())

    def get_in_out_lines(self):
        return (
            Line(self.azimuth_in, self.radius),
            Line(self.azimuth_out, self.radius))

    def curves_right(self):
        """returns True if this curve curves right or is straight, and False otherwise."""

        diff_azimuth = self.azimuth_out.to_degrees() - self.azimuth_in.to_degrees()
        return (180.0 <= diff_azimuth and diff_azimuth <= 360.0) or (-180.0 <= diff_azimuth and diff_azimuth <= 0.0)

    def get_chord_triangle_area(self):
        delta_rad = self.get_delta().to_radians()
        return (self.radius*self.radius) * math.sin(delta_rad) / 2.0

    def get_sector_area(self):
        delta_rad = self.get_delta().to_radians()
        return (delta_rad / 2.0) * self.radius * self.radius

    def get_segment_area(self):
        return self.get_sector_area() - self.get_chord_triangle_area()

    def get_delta(self):
        diff_azimuth = self.azimuth_out.to_degrees() - self.azimuth_in.to_degrees()
        return AzimuthDMS.from_degrees(
            abs(180.0 - abs(diff_azimuth)))

    def get_curve_length(self):
        return self.get_delta().to_radians() * self.radius

    def get_chord_line(self):
        # To get the chord line azimuth, you can add/subtract half delta with the azimuth out based on curvature.
        half_delta = AzimuthDMS.from_degrees(self.get_delta().to_degrees()/2.0)
        chord_azimuth = (self.azimuth_out + half_delta) if self.curves_right() else (self.azimuth_out - half_delta)
        # Long Chord Length = 2 * radius * sin(delta/2)
        chord_length = 2.0 * self.radius * math.sin(self.get_delta().to_radians() / 2.0)
        return Line(chord_azimuth, chord_length)


class DMDCalculationResult:
    def __init__(self, line_area, line_length, lat, dep, segment_area=0.0, curve_length=0.0):
        self.line_area = line_area
        self.line_length = line_length
        self.lat = lat
        self.dep = dep
        self.segment_area = segment_area
        self.curve_length = curve_length
        self.total_area = line_area + segment_area
        self.total_length = line_length + curve_length

    def __repr__(self):
        return 'DMDCalculationResult(line_area={:.8f},line_length={:.8f},lat={:.8f},dep={:.8f},segment_area={:.8f},curve_length={:.8f},total_area={:.8f},total_length={:.8f})'.format(
            self.line_area,
            self.line_length,
            self.lat,
            self.dep,
            self.segment_area,
            self.curve_length,
            self.total_area,
            self.total_length)


def perform_dmd_calculation(obj_list):
    """Performs the DMD calculation over a list of objects.

    We use helper functions that contain the majority of the calculation update logic.
    These functions are technically "closures", as they close over the outer environment's
    variables to allow them to do the update.
    """
    log = logging.getLogger(__name__)

    sum_line_area, sum_line_length = 0.0, 0.0
    sum_segment_area, sum_curve_length = 0.0, 0.0
    sum_lat, sum_dep = 0.0, 0.0
    last_dmd, last_dep = 0.0, 0

    def update_calculation_for_line(line, should_perform_line_sum=True):
        log = logging.getLogger(__name__)

        # We use `nonlocal` to be able to refer the the outer variables above in
        # order to be able to write to them.
        #
        # Otherwise, we would just create new local variables that would vanish
        # while running this function each time.
        nonlocal sum_line_area, sum_line_length, sum_lat, sum_dep, last_dmd, last_dep
        # Get the current line's lat, dep, DMD, area
        (lat, dep) = line.to_lat_dep()
        (dmd, area) = line.to_dmd_area_terms(last_dmd, last_dep)
        log.debug('line={}, lat={}, dep={}, dmd={}, area={}'.format(
            line, lat, dep, dmd, area))
        # Update sums
        sum_lat += lat
        sum_dep += dep
        sum_line_area += area
        # Update line-specific sums only if needed -- we will use this as a sub-calculation for
        # inner radial curve calculations too, which won't need this part.
        if should_perform_line_sum:
            sum_line_length += line.distance
        # Update the last dmd/dep fields before next iteration
        (last_dmd, last_dep) = dmd, dep

    def update_calculation_for_inner_radial_curve(curve):
        log = logging.getLogger(__name__)

        nonlocal sum_line_area, sum_segment_area, sum_curve_length

        # For inner radial curves, we want to calculate the chord line DMD
        # and then account for the curve's segment area from the chord after.
        chord_line = curve.get_chord_line()

        # Calculate the DMD for the chord_line
        update_calculation_for_line(chord_line, should_perform_line_sum=False)

        # Get the curve's segment area, and add/subtract from sum if it curves right/left.
        # (This is because we are doing the traverse clockwise.)
        segment_area = curve.get_segment_area()
        log.debug("curve={}, segment_area={}, curves_right={}".format(curve, segment_area, obj.curves_right()))
        if obj.curves_right():
            sum_segment_area += segment_area
        else:
            sum_segment_area -= segment_area

        sum_curve_length += curve.get_curve_length()

    # Iterate over the list of objects for the sums
    for obj in obj_list:
        if isinstance(obj, Line):
            update_calculation_for_line(obj)
        elif isinstance(obj, InnerRadialCurve):
            update_calculation_for_inner_radial_curve(obj)
        else:
            raise ValueError(
                '`{}` is not a valid value that can be used in a Double-Meridian Distance calculation.'.format(obj))
        #log.debug('{}'.format(DMDCalculationResult(
        #    line_area=abs(sum_line_area),
        #    line_length=sum_line_length,
        #    lat=sum_lat,
        #    dep=sum_dep,
        #    segment_area=sum_segment_area,
        #    curve_length=sum_curve_length)))
        
    # Need to take the absolute value of the sum for area
    sum_line_area = abs(sum_line_area)

    # Return the results
    return DMDCalculationResult(
        line_area=sum_line_area,
        line_length=sum_line_length,
        lat=sum_lat,
        dep=sum_dep,
        segment_area=sum_segment_area,
        curve_length=sum_curve_length)


def output_report(results):
    pass


def construct_argument_parser():
    parser = argparse.ArgumentParser(
        description="""Calculates traverse calculations using double-meridian distance (DMD).

    Points for the traverse are assumed to be taken in a clockwise fashion, with azimuths
    according to U.S. Geological Survey standards (north-based, clockwise).
    """)
    parser.add_argument('input_filename', metavar='<input>', type=str,
                        help='input CSV file')
    parser.add_argument('output_filename', metavar='<output>', type=str,
                        help='output file', default='out.csv')
    parser.add_argument('-v', '--verbosity', metavar='<verbosity>', type=str,
                        help='limits program diagnostic output level', default='INFO',
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser


def config_application_logger(args):
    logging.basicConfig(
        format='%(asctime)s UTC [%(levelname)s] %(funcName)s: %(message)s',
        stream=sys.stderr,
        level=getattr(logging, args.verbosity.upper(), None))
    logging.Formatter.converter = time.gmtime


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
                azimuth = AzimuthDMS.from_str(row[0])
                distance = float(row[1])
                object_list.append(Line(azimuth, distance))
            # If the length of the individual row is 3, this means that
            # that the optional radial-out azimuth for simple radial curves is present.
            #
            # Interpret the current row as a simple radial curve.
            elif len(row) == 3:
                if row[2] == '' or row[2] == None:
                    azimuth = AzimuthDMS.from_str(row[0])
                    distance = float(row[1])
                    object_list.append(Line(azimuth, distance))
                else:
                    azimuth_in = AzimuthDMS.from_str(row[0])
                    radius = float(row[1])
                    azimuth_out = AzimuthDMS.from_str(row[2])
                    object_list.append(InnerRadialCurve(
                        azimuth_in=azimuth_in, azimuth_out=azimuth_out, radius=radius))

    # Once the data is read, return it.
    return object_list


def main():
    # Parse the command line arguments
    parser = construct_argument_parser()
    args = parser.parse_args()

    # Configure the application logger
    config_application_logger(args)
    log = logging.getLogger(__name__)

    # Perform the calculation
    obj_list = read_input_as_object_list(args)
    calculation_stats = perform_dmd_calculation(obj_list)

    log.info(obj_list)
    log.info(calculation_stats)

    output_report(calculation_stats)


if __name__ == "__main__":
    main()
