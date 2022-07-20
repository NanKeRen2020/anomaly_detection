#! /usr/bin/env python

#/*
# * Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
# * All rights reserved.
# *
# * This program is free software: you can use, modify and/or
# * redistribute it under the terms of the GNU General Public
# * License as published by the Free Software Foundation, either
# * version 3 of the License, or (at your option) any later
# * version. You should have received a copy of this license along
# * this program. If not, see <http://www.gnu.org/licenses/>.
# */

import cairo
from math import pi
import math
import sys
import sys,os
import copy
import argparse
import numpy as np
import iio

minNFAColor=(1,0,0) #RED
colors= [(1, 165./255., 0), (0,1,0), (0,1,1), (1,1,1)]
colors_name= ["Orange", "Green", "Cyan", "White"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_img', type=str, help=('path to source image (png)'))
    parser.add_argument('vgg_layer', type=int, help=('VGG layer num used (0 if NN not used)'))
    parser.add_argument('scales', type=int, help=('Number of image scales'))
    parser.add_argument('path', type=str, help=('Path to the pixelNFA and radiusNFA'))
    parser.add_argument('pixelNFA_prefix', type=str, help=('pixelNFA files\' prefix'))
    parser.add_argument('radiusNFA_prefix', type=str, help=('radiusNFA files\' prefix'))
    parser.add_argument('output_name', type=str, help=('path to output image'))
    parser.add_argument('threshold_nfa', type=int, help=('nfa threshold'))
    args = parser.parse_args()

    # Create a version of the input which has been toned down
    orig = iio.read(args.input_img)
    for x in range(orig.shape[0]):
        for y in range(orig.shape[1]):
            orig[x][y] = 0.33*orig[x][y]
    iio.write("input_0_rescaled.png", orig)


    # Read input and create the visu
    image = cairo.ImageSurface.create_from_png(os.path.join(args.path,"input_0_rescaled.png"));
    cr = cairo.Context(image)
    cr.save()
    cr.set_source_surface(image, 0, 0)
    cr.get_source().set_filter(cairo.FILTER_NEAREST)
    cr.paint()
    cr.restore()

    threshold_nfa = args.threshold_nfa

    cr.set_line_width(0.5)

    nn = [0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    couche = args.vgg_layer

    detections=[]
    real_scale = 0
    prev_shape = -1
    for scale in range(args.scales):
        nfa = iio.read(os.path.join(args.path, args.pixelNFA_prefix+str(scale)+".tiff")) 
        radius = iio.read(os.path.join(args.path, args.radiusNFA_prefix+str(scale)+".tiff")) 
        if nfa.shape[0] != prev_shape:
            real_scale=scale
            prev_shape=nfa.shape[0]
        # For each position of the image
        paddingx = round(max(0,(orig.shape[0] - (2**nn[couche])*(2**real_scale)*nfa.shape[0])/2))
        paddingy = round(max(0,(orig.shape[1] - (2**nn[couche])*(2**real_scale)*nfa.shape[1])/2))
        for x in range(nfa.shape[0]):
            for y in range(nfa.shape[1]):
                #if it is add it to the list
                if nfa[x][y] < threshold_nfa:
                    detections.append((nfa[x][y], paddingy+(2**nn[couche])*(2**real_scale)*(y+1.5), paddingx+(2**nn[couche])*(2**real_scale)*(x+1.5), (2**nn[couche])*(2**real_scale)*radius[x][y]+1)) 

    detections.sort()

    if len(detections) > 0:
        # Plot the detection with the lowest NFA and remove it form the list
        (nfa_best, x_best, y_best, radius_best) = detections.pop(0)

        # Plot the remaining detections
        for det in range(len(detections)-1, -1, -1):
            (nfa, x, y, radius) = detections[det]
            # If it a detection similar to the best one then plot it also in the color of the minNFA color
            if nfa == nfa_best or np.isnan(nfa):
                cr.set_source_rgb(minNFAColor[0], minNFAColor[1], minNFAColor[2])
                cr.arc(x, y, radius, 0, 2*pi)
                cr.stroke()
            else:
                color = len(colors) - 1 - int(min(math.floor(math.log(max(-nfa,1))), len(colors)-1))
                cr.set_source_rgb(colors[color][0], colors[color][1], colors[color][2])
                cr.arc(x, y, radius, 0, 2*pi)
                cr.stroke()

        cr.set_line_width(1)
        cr.set_source_rgb(minNFAColor[0], minNFAColor[1], minNFAColor[2])
        cr.arc(x_best, y_best, radius_best, 0, 2*pi)
        cr.stroke()

        legend = open(os.path.join(args.path, "legend.txt"), "w")
        legend.write("Red corresponds to the minimum logNFA " + str(nfa_best) + "\n")
        for c in range(len(colors)):
                if c == 0:
                    legend.write(colors_name[c] + " corresponds to a logNFA between -infinity  and " + str(-math.exp(len(colors)-1)) + "\n")
                elif c == (len(colors)-1):
                    legend.write(colors_name[c] + " corresponds to a logNFA between " + str(-math.exp(1)) + " and " + str(0) + "\n")
                else:
                    legend.write(colors_name[c] + " corresponds to a logNFA between " + str(-math.exp(len(colors)-c)) + " and " + str(-math.exp(len(colors)-c-1)) + "\n")
        legend.close()

    image.write_to_png(args.output_name)


