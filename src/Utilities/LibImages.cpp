/*
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work: Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"
#include "mt19937ar.h"

#include <unistd.h> // getpid
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <cmath>

extern "C" {
#include "iio.h"
}

using namespace std;

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
,	std::vector<float> &o_im
,	ImageSize &o_imSize
,	const bool p_verbose
){
	//! read input image
	if (p_verbose) cout << endl << "Read input image...";

	float *imTmp = NULL;
	int w, h, c;
	imTmp =  iio_read_image_float_vec(p_name, &w, &h, &c);

	if (!imTmp) {
		cout << "error :: " << p_name << " not found or not a correct png image" << endl;
		return EXIT_FAILURE;
	}

	if (p_verbose) cout << "done." << endl;


	//! Some image informations
	if (p_verbose) {
		cout << "image size :" << endl;
		cout << " - width          = " << w << endl;
		cout << " - height         = " << h << endl;
		cout << " - nb of channels = " << c << endl;
	}

	//! Initializations
	o_imSize.width      = w;
	o_imSize.height     = h;
	o_imSize.nChannels  = c;
	o_imSize.wh         = w * h;
	o_imSize.whc        = w * h * c;
	o_im.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
		o_im[k] = imTmp[k];

	free(imTmp);

	return EXIT_SUCCESS;
}

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
){
    //! Allocate Memory
    float* imTmp = new float[p_imSize.whc];

    unsigned c = p_imSize.nChannels;
    unsigned w = p_imSize.width;
    unsigned h = p_imSize.height;

    for (unsigned k = 0; k < w*h*c; k++) {
        imTmp[k] = i_im[k];
    }

    iio_save_image_float_vec(p_name, imTmp, w, h, c);
    //! Free Memory
    delete[] imTmp;

    return EXIT_SUCCESS;
}

/**
 * @brief Compute a difference image between i_im1 and i_im2.
 *
 * @param i_im1: reference image;
 * @param i_im2: image to compare;
 * @param o_imDiff: will contain the difference;
 * @param p_sigma : standard deviation of the noise;
 * @param p_min, p_max : range of data (usually [0, 255]);
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   std::vector<float> &o_imDiff
){
    if (i_im1.size() != i_im2.size()) {
        cout << "Can't compute difference, i_im1 and i_im2 don't have the same size" << endl;
        cout << "i_im1 : " << i_im1.size() << endl;
        cout << "i_im2 : " << i_im2.size() << endl;
        return EXIT_FAILURE;
    }

    const unsigned size = i_im1.size();
    if (o_imDiff.size() != size)
        o_imDiff.resize(size);

    for (unsigned k = 0; k < size; k++)
        o_imDiff[k] = (i_im1[k] - i_im2[k]);

    return EXIT_SUCCESS;
}
