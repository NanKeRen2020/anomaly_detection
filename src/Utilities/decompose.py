#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import os
import sys

import numpy as np
sys.path.append(os.path.dirname(__file__))
import iio

def gblur(img, sig):
   import scipy.ndimage as ndimage
   return ndimage.gaussian_filter(img, sigma=(sig, sig, 0), order=0)

def down(image):
   sub = gblur(image, 1*1.39)
   return sub[::2,::2,:]

def decompose(image, nn, prefix, levels, suffix, cur=0):
   iio.write(prefix+str(cur)+suffix, image)

   if levels==1:
      return [image]
    
   sub = down(image)
   if ((nn == "True") and ((sub.shape[0] < 150) or (sub.shape[1] < 150))):
       sub = image

   decompose(sub, nn, prefix, levels-1, suffix, cur+1)
   return

def main():
   # verify input
   if len(sys.argv) > 4:
      image = iio.read(sys.argv[1])
      if len(image.shape) == 2:
        image = np.expand_dims(image, 2)
      nn = sys.argv[2]
      prefix = sys.argv[3]
      levels = int(sys.argv[4])
      suffix = sys.argv[5]
   else:
      print("Incorrect syntax, use:")
      print("  > " + sys.argv[0] + " input nn prefix levels suffix")
      print("        multiscale decomposition using gaussian pyramid")
      sys.exit(1)

   decompose(image, nn, prefix, levels, suffix)


if __name__ == '__main__':
    main()
