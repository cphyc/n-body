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
parser.add_argument('--particles', type=int, default=1024,
                    help='Number of particles to use (default %(default)s)')
args = parser.parse_args()

class AnimatedScatter(object):
    def __init__(self, numpoints):
        time = 8
        dt = 1e-3
        dump_freq = 10
        nbiter = int(time/dt)
        self.data = pd.read_csv(args.output, skiprows=[i for i in range(1,1+nbiter*(numpoints+1)//dump_freq,numpoints+1)], delim_whitespace=True)
        self.numpoints = numpoints

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111,projection = '3d')
        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=nbiter//dump_freq, init_func=self.setup_plot, interval=1, repeat=False)

    def setup_plot(self):
        init = 0
        end = self.numpoints

        self.scat = self.ax.scatter(self.data.x[init:end], self.data.y[init:end], self.data.z[init:end], s=5, animated=True)

        self.ax.set_xlim3d(-1.5, 1.5)
        self.ax.set_ylim3d(-1.5, 1.5)
        self.ax.set_zlim3d(-1.5, 1.5)

        return self.scat,

    def animate(self, i):
        init = i * self.numpoints
        end = (i+1) * self.numpoints
        print('Frame', i)
        self.scat._offsets3d = (self.data.x[init:end], self.data.y[init:end], self.data.z[init:end])

        return self.scat,

    def show(self):
        plt.show()

    def save(self):
        self.ani.save('test.mp4', writer='ffmpeg', fps=10, bitrate=10000, dpi=180)

if __name__ == '__main__':
    a = AnimatedScatter(args.particles)
    a.save()
#    a.show()
