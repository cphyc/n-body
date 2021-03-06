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
    # Read data as CSV
    data = pd.read_csv(args.output, delim_whitespace=True)

    # Plot
    plt.plot(data.t, data.Ep, label='$Ep$')
    plt.plot(data.t, data.Ec, label='$Ec$')
    plt.plot(data.t, data.Ep+data.Ec, label='$E_\mathrm{tot}$')

    plt.xlabel('Time (arbitrary unit)')
    plt.ylabel('Energy (arbitrary unit)')

    plt.grid('on')
    plt.legend()

    plt.figure()

    deltaEc = data.Ec.max() - data.Ec.min()
    deltaE = data.Ep+data.Ec - (data.Ec.mean() + data.Ep.mean())

    plt.plot(data.t, deltaE/deltaEc, label='$\\frac{E_\mathrm{tot} - \\overline{E_\mathrm{tot}}}{Ec_\mathrm{max} - Ec_\mathrm{min}}$')

    plt.xlabel('Time (arbitrary unit)')
    plt.ylabel('Total energy dissipation relative to max Ec variation')

    plt.grid('on')
    plt.legend()

    plt.show()
