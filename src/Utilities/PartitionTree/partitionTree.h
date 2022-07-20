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


#ifndef PARTITIONTREE_HEADER_H
#define PARTITIONTREE_HEADER_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <unordered_map>
#include "../PatchManager/patchManager.h"
#include "../LibImages.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 *
 * Basis for different partition trees
 *
 */

class PartitionTree
{
	public:
		/**
		 * @brief Compute the nearest neighbors and return their position and the distance to each of them
		 *
		 * @param localResult: Used to return the position of the patches corresponding to the nearest neighbors alongside their distance to the query patch
		 * @param pidx: The query patch
		 * @param reranking: If using a dimensionality reduction, rerank the results using the true distances
		 *
		 * @return the stop condition of the bin in which pidx is
		 **/
		virtual void retrieveFromTree(std::vector<std::pair<float, unsigned> >& localResult, const unsigned pidx) = 0;

		/**
		 * @brief Return the bin containing the query patch
		 *
		 * @param pm: The new patch manager
		 **/
		virtual void updatePM(PatchManager& pm) = 0;
};



#endif
