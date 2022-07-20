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

#ifndef IMAGEPM_HEADER_H
#define IMAGEPM_HEADER_H

#include <algorithm>
#include <cstdlib>
#include "../LibImages.h"
#include "patchManager.h"

/**
 *
 * Most basic patch manager. Patches directly match the image.
 *
 */

class ImagePatchManager: public PatchManager {
	public:
		/**
		 * @brief Basic constructor
		 *
		 * @param image: The video on which the manager will be based on
		 * @param ims: Image size
		 * @param sizeP: Spatial size for the patches
		 */
		ImagePatchManager(std::vector<float> const& image, ImageSize ims, int sizeP);

		/**
		 * check patchManager.h for information on these functions
		 * @{
		 */
		void getAllPatches(std::vector<unsigned>& allPatches);
		int getNbPatches();
		float distance(unsigned id1, unsigned id2);
		const ImageSize* infoIm() {return &imSize;};
		/**
		 * @}
		 */

		~ImagePatchManager();
	private:
		/// The video on which the manager is based
		const std::vector<float> *im;
		/// Image size
		ImageSize imSize;
		/// Spatial size of a patch
		int sizePatch;
		/// Number of chanels from the video to use (0 to nbChannels=
		int nbChannels;
};
#endif
