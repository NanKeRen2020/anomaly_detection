/*
 * Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
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
 * @file imagePatchManager.cpp
 * @brief Patch manager using the patches from the image
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "imagePatchManager.h"

ImagePatchManager::ImagePatchManager(std::vector<float> const& image, ImageSize ims, int sizeP)
{
	im = &image;
	imSize = ims;
	sizePatch = sizeP;
}

void ImagePatchManager::getAllPatches(std::vector<unsigned>& allPatches)
{
	for(unsigned x = 0, k = 0; x < imSize.width - sizePatch + 1; ++x)
	for(unsigned y = 0; y < imSize.height - sizePatch + 1; ++y, ++k)
		allPatches[k] = imSize.nChannels*(x + y*imSize.width);
}

int ImagePatchManager::getNbPatches()
{
	return (imSize.width - sizePatch + 1) * (imSize.height - sizePatch + 1);
}

float ImagePatchManager::distance(unsigned patch1, unsigned patch2)
{
	float dist = 0.f, dif;
	for (unsigned hy = 0; hy < sizePatch; ++hy)
	for (unsigned hx = 0; hx < sizePatch; ++hx)
	for (unsigned hc = 0; hc < imSize.nChannels; ++hc)
		dist += (dif = (*im)[patch1 + hx*imSize.nChannels + hy*imSize.nChannels*imSize.width + hc] - (*im)[patch2 + hx*imSize.nChannels + hy*imSize.nChannels*imSize.width + hc]) * dif;
	return std::sqrt(dist / (sizePatch * sizePatch * imSize.nChannels)) / 255.f;
}

ImagePatchManager::~ImagePatchManager()
{

}
