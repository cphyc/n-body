#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Plot data from output of the n-body simulation.')
parser.add_argument('--output', type=str, default='output.dat',
                    help='The output file (default %(default)s)')
args = parser.parse_args()

class AnimatedScatter(object):
    def __init__(self, numpoints):
        time = numpoints/100
        dt = 1e-2
        nbiter = int(time/dt)
        self.data = pd.read_csv(args.output, skiprows=[i for i in range(1,1+nbiter*(numpoints+1),numpoints+1)], delim_whitespace=True)
        self.numpoints = numpoints

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111,projection = '3d')
        self.ani = animation.FuncAnimation(self.fig, self.animate, init_func=self.setup_plot, interval=1, repeat=False)

    def setup_plot(self):
        init = 0
        end = self.numpoints

        self.scat = self.ax.scatter(self.data.x[init:end], self.data.y[init:end], self.data.z[init:end], animated=True)

        self.ax.set_xlim3d(-2, 2)
        self.ax.set_ylim3d(-2, 2)
        self.ax.set_zlim3d(-2, 2)

        return self.scat,

    def animate(self, i):
        init = i * self.numpoints
        end = (i+1) * self.numpoints

        self.scat._offsets3d = (self.data.x[init:end], self.data.y[init:end], self.data.z[init:end])

        return self.scat,

    def show(self):
        plt.show()

    def save(self):
        self.ani.save('test.mp4', writer='ffmpeg', fps=10, bitrate=10000, dpi=180)

if __name__ == '__main__':
    a = AnimatedScatter(1024)
    a.save()
#    a.show()
