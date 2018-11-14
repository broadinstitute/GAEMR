"""
CREDIT FOR CODE GOES TO:
Majority of following code courtesy of example at Matplotlib tutorial
http://matplotlib.sourceforge.net/mpl_examples/api/radar_chart.py

for comments and further details, please go to above URL. Many thanks
to developers who contributed code there!

"""

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
import math
import sys
from matplotlib.pyplot import legend


#a unified implementation for radar chart usage
class RadarChart:
    """a generic radar chart"""

    def __init__(self,num_vars):
        """initialize class"""
        self.num_vars=num_vars
        
    #a class to generate radar plots
    def radar_factory(self,frame='circle'):
        theta = 2*np.pi * np.linspace(0, 1-1./self.num_vars, self.num_vars)
        theta += np.pi/2

        def draw_poly_patch(self):
            verts = unit_poly_verts(theta)
            return plt.Polygon(verts, closed=True, edgecolor='k')

        def draw_circle_patch(self):
            return plt.Circle((0.5, 0.5), 0.5)

        patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
        if frame not in patch_dict:
            raise ValueError, 'unknown value for `frame`: %s' % frame

        #private class to represent a radar chart axes
        class RadarAxes(PolarAxes):
            name = 'radar'
            RESOLUTION = 1
            draw_patch = patch_dict[frame]

            def fill(self, *args, **kwargs):
                closed = kwargs.pop('closed', True)
                return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

            def plot(self, *args, **kwargs):
                lines = super(RadarAxes, self).plot(*args, **kwargs)
                for line in lines:
                    self._close_line(line)

            def _close_line(self, line):
                x, y = line.get_data()
                if x[0] != x[-1]:
                    x = np.concatenate((x, [x[0]]))
                    y = np.concatenate((y, [y[0]]))
                    line.set_data(x, y)

            def set_varlabels(self, labels):
                self.set_thetagrids(theta * 180/np.pi, labels)

            def _gen_axes_patch(self):
                return self.draw_patch()

            def _gen_axes_spines(self):
                if frame == 'circle':
                    return PolarAxes._gen_axes_spines(self)
                spine_type = 'circle'
                verts = unit_poly_verts(theta)
                verts.append(verts[0])
                path = Path(verts)

                spine = Spine(self, spine_type, path)
                spine.set_transform(self.transAxes)
                return {'polar': spine}

        register_projection(RadarAxes)
        return theta


    def unit_poly_verts(self,theta):
        x0, y0, r = [0.5] * 3
        verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
        return verts
