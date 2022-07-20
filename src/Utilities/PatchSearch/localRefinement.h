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


#ifndef LR_HEADER_H
#define LR_HEADER_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "../comparators.h"
#include "../LibImages.h"
#include "../PatchManager/patchManager.h"
#include "../PartitionTree/forestManager.h"
#include "localPatchPropagation.h"

/**
 *
 * Local refinement search
 *
 */

class LocalRefinement : public LocalPatchPropagation {
	public:
		/**
		 * @brief basic constructor
		 *
		 * @param pm_: Current patch manager
		 * @param fm_: Forest to be used during the computation
		 * @param kNN_: Number of nearest neighbors to be computed
		 * @param prms: Denoising parameters
		 */
		LocalRefinement(PatchManager& pm_, ForestManager& fm_, int kNN_, Params prms) {pm = &pm_; fm = &fm_; kNN = kNN_; params = prms;};

		/**
		 * Check localPatchPropagation.h for more informations on these functions
		 * @{
		 */
		int estimateSimilarPatches(std::vector<std::pair<float,unsigned> > &index, const unsigned pidx, bool excludeYourself = false);
		/**
		 * @}
		 */
	private:
		/// Number of nearest neighbors to be computed
		ForestManager* fm;
		/// Current patch manager
		PatchManager* pm;
		/// Number of nearest neighbors to be computed
		unsigned kNN;
		/// Denoising parameters
		Params params;
};
#endif
