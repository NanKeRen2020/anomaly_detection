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

#ifndef LOCALPATCHPROP_HEADER_H
#define LOCALPATCHPROP_HEADER_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "../LibImages.h"
#include "../PatchManager/patchManager.h"
#include "../PartitionTree/forestManager.h"

/**
 *
 * Basis for different patch search
 *
 */

class LocalPatchPropagation {
	public:
		/**
		 * @brief Compute the nearest neighbors and return their position
		 *
		 * @param index: Used to return the position of the patches corresponding to the nearest neighbors
		 * @param pidx: The query patch
		 * @param excludeYourself: Exclude patches that are too close (spatially) to the query patch
		 *
		 * @return the number of nearest neighbors computed
		 **/
		virtual int estimateSimilarPatches(std::vector<std::pair<float, unsigned> > &index, const unsigned pidx, bool excludeYourself = false) = 0;
};
#endif
