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


#ifndef FORESTMANAGER_HEADER_H
#define FORESTMANAGER_HEADER_H

#include <vector>
#include "partitionTree.h"
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include "../comparators.h"
#include "../../model/params.h"
#include "../PatchManager/patchManager.h"

class ComparePairsFast {
	public:
		bool operator()(std::pair<int,unsigned>& a, std::pair<int,unsigned>& b)
		{
			return a.first > b.first;
		}
};

/**
 *
 * Used to manipulate multiple partition trees
 *
 */

class ForestManager
{
	public:	

		/**
		 * @brief Construct a basic forest manager from a list of trees
		 *
		 * @param pm_: The patch manager to be used with the forest (it is not necessarely the same as the one used for the construction of the trees)
		 * @param forest_: The list of trees to be used to create the forest
		 * @param kNN: The number of nearest neighbors to be computed during each search
		 * @param param: Parameter used to construct the trees
		 *
		 **/
		ForestManager(PatchManager& pm_, std::vector<PartitionTree*>& forest_, int kNN, Params param);

		/**
		 * @brief Empty constructor
		 **/

		ForestManager() {};

		/**
		 * @brief Compute the nearest neighbors and return their position and the distance to each of them
		 *
		 * @param index: Used to return the position of the patches corresponding to the nearest neighbors
		 * @param pidx: The query patch
		 * @param excludeYourself: If you want to exclude the current patch itself from the search step. 
		 * 	                   It is extremly useful when using projections
		 *
		 * @return the number of nearest neighbors computed, (it should always be kNN if there wasn't any problem)
		 **/
		int retrieveFromForest(std::vector<std::pair<float, unsigned> >& index, const unsigned pidx, bool excludeYourself = false);

		/**
		 * @brief Get the number of trees in the forest
		 *
		 * @return number of trees
		 **/
		int getNbTrees() { return forest.size(); };

		/**
		 * @brief Update the patch manager to a new one
		 *
		 * @param pm_: the new patch manager
		 **/
		void updatePM(PatchManager& pm_);

		
		/**
		 * @brief Get the patch manager used by the forest
		 *
		 * @return the patch manager
		 **/
		PatchManager* getpm() {return  pm;};

		/**
		 * @brief Get the number of similar patches searched using the forest
		 *
		 * @return the number of nearest neighbors
		 **/
		int getkNN() {return kNN;};

		/**
		 * @brief Get the parameters used to create the forest manager
		 *
		 * @return the parameter
		 **/
		Params get_params() {return params;};

		/**
		 * @brief Access the trees contained in the forest
		 *
		 * @return list of pointer to the trees
		 **/
		std::vector<PartitionTree*> getforest() {return forest;};

	private:
		PatchManager *pm;
		std::vector<PartitionTree*> forest;
		int kNN;
		Params params;

};
#endif
