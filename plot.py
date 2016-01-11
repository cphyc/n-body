#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description='Plot data from output of the n-body simulation.')
parser.add_argument('--output', type=str, default='output_int.dat',
                    help='The output file (default %(default)s)')
args = parser.parse_args()


if __name__ == '__main__':
