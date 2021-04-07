import csv
import math
import argparse

parser = argparse.ArgumentParser(
    description='Calculates double-meridian distance (DMD).')
parser.add_argument('input_filename', metavar='i', type=str,
                    help='input CSV file')
parser.add_argument('output_filename', metavar='o', type=str,
                    help='output CSV file', default='out.csv')
args = parser.parse_args()

print(args)

# Get input pairs:
# azimuth (line-out/radial-in), distance/radius, azimuth2 (radial-out)
with open(parser.input_filename, newline='') as input_file:
    data_reader = csv.reader(input_file, delimiter=',', quotechar='\"')
    for row in data_reader:
        print(row)