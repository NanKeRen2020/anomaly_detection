/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIB_IMAGES_H_INCLUDED
#define LIB_IMAGES_H_INCLUDED

#include <vector>
#include <string>
#include <fftw3.h>

/**
 * @brief Structure containing size informations of an image.
 *
 * @param width     : width of the image;
 * @param height    : height of the image;
 * @param nChannels : number of channels in the image;
 * @param wh        : equal to width * height. Provided for convenience;
 * @param whc       : equal to width * height * nChannels. Provided for convenience.
 **/
struct ImageSize
{
	unsigned width;
	unsigned height;
	unsigned nChannels;
	unsigned wh;
	unsigned whc;
};

/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
int loadImage(
    const char* p_name
,   std::vector<float> &o_im
,   ImageSize &o_imSize
,   const bool p_verbose = false
);

/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 * @param p_min, p_max : range of data (usually [0, 255]).
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int saveImage(
    const char* p_name
,   std::vector<float> const& i_im
,   const ImageSize &p_imSize
);

/**
 * @brief Compute a difference image between i_im1 and i_im2.
 *
 * @param i_im1: reference image;
 * @param i_im2: image to compare;
 * @param o_imDiff: will contain the difference;
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   std::vector<float> &o_imDiff
);

#endif // LIB_IMAGES_H_INCLUDED
