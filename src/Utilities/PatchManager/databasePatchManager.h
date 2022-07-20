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

#ifndef DATABASEPM_HEADER_H
#define DATABASEPM_HEADER_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "../LibImages.h"
#include "patchManager.h"

/**
 *
 * Patch manager defined on an external database
 *
 */

class DatabasePatchManager: public PatchManager {
	public:
		/**
		 * @brief Basic constructor
		 *
		 * @param image: The database on which the manager will be based on
		 * @param sizeP: Spatial size for the patches
		 * @param ims: Image size
		 */
		DatabasePatchManager(std::vector<float> const& image, ImageSize ims, int sizeP);

		/**
		 * check patchManager.h for information on these functions
		 * @{
		 */
		void getAllPatches(std::vector<unsigned>& allPatches);
		int getNbPatches();
		float distance(unsigned id1, unsigned id2);
		const ImageSize* infoIm() {return &imSize;};

		void setCurrentImage(const std::vector<float>& cur, ImageSize curSize);
		/**
		 * @}
		 */

		~DatabasePatchManager();
	private:
		/// The database on which the manager is based (for now it is on the form of a video
		const std::vector<float>* im;
		/// The video that will be used alongside the database
		const std::vector<float>* current;
		/// Spatial size of a patch
		int sizePatch;
		/// Size of the database image
		ImageSize imSize;
		/// Size of the current image
		ImageSize curSize;
};
#endif
