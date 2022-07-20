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

#ifndef PATCHMANAGER_HEADER_H 
#define PATCHMANAGER_HEADER_H 

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "../LibImages.h"
#include <cmath>


/**
 *
 * Basis for different patch manager
 *
 */

class PatchManager 
{
	public:
		/**
		 * @brief Load all patches with the internal indexing system 
		 *
		 * @param allPatches: Used to return the indexing
		 */
		virtual void getAllPatches(std::vector<unsigned>& allPatches) = 0;

		/**
		 * @brief Return the total number of patches
		 *
		 * @return number of patches
		 */
		virtual int getNbPatches() = 0;

		/** 
		 * @brief Compute the distance between the first patch and the second patch. The first patch is indexed directly from the video itself. The second patch is indexed by the local system.
		 *
		 * @param id1: Index of the first patch (in the video)
		 * @param id2: index of the second patch (in the local indexing)
		 *
		 * @return distance between the two patches
		 */
		virtual float distance(unsigned id1, unsigned id2) = 0;

		/** 
		 * @brief Get the size of the video on which the patch manager is based
		 *
		 * @return the size of the video
		 */
		virtual const ImageSize* infoIm() = 0;

	private:
};
#endif
