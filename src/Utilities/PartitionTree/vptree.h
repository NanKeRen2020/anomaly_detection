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

#ifndef VPTREE_HEADER_H
#define VPTREE_HEADER_H

#include <queue>
#include "partitionTree.h"
#include "../../model/params.h"
#include "prmTree.h"

/**
 *
 * VP-tree
 *
 */

class VPTree : public PartitionTree
{
	public:
		/**
		 * @brief Constructor
		 *
		 * @param pm_: Patch manager on which the tree is constructed
		 * @param paramTree: Contains the parameters for the tree
		 *
		 * @{
		 */
		VPTree(PatchManager& pm_, PrmTree& paramTree);
		/**
		 * @}
		 */

		/**
		 *
		 * Check partitionTree.h for informations on the following functions
		 * @{
		 */
		void retrieveFromTree(std::vector<std::pair<float, unsigned> >& localResult, const unsigned pidx);

		void updatePM(PatchManager& pm_) {pm = &pm_;};
		/**
		 * @}
		 */

		~VPTree()
		{
			delete root;
		}

	private:
		/**
		 * @brief Structures of the node defining a KD-tree
		 *
		 * @param id: Each node of a tree has a unique id
		 * @param left,right: Beginning and end of the bin (if the node is a leaf) in the list of patches
		 * @param vantagePoint: Index corresponding to the vantage point of the split
		 * @param splitDistance: Distance to the vantage point separating the left children to the right children 
		 * @param diameter: diameter of the ball containing the Node
		 * @param child1: Pointer to the left child of this node in the tree
		 * @param child2: Pointer to the right child of this node in the tree
		 */
		struct Node
		{
			unsigned id;
			unsigned left, right;

			unsigned vantagePoint;
			float splitDistance;

			float diameter;

			struct Node* child1;
			struct Node* child2;	

			~Node()
			{
				if (child1) delete child1;
				if (child2) delete child2;
			}
		};

		/// Number of elements to be returned for each search
		int kNN;
		/// Max nb of elem in a leaf
		int leaf_max_size;
		/// Nb of selection for the vp
		const int nbSelectionVP;
		/// Nb of points to consider during the creation of the tree
		const int subSampling;

		/// Root of the KDTree 
		Node* root;

		Params param;

		/// PatchManager
		PatchManager* pm;

		/// index of the elements on which the tree is based
		std::vector<unsigned> idxPatches;

		/// Move the elements according to their position relatively to the pivot
		unsigned partition(unsigned left, unsigned right, unsigned pivot, std::vector<float>& distances, unsigned n, unsigned l);

		/// Fast computation of the median along a specific dimension for a set of element delimited by left and right 
		float quickselect(unsigned left, unsigned right, std::vector<float>& distances);

		/// Split the elements recursively and return a node each time, main function used to construct the tree 
		Node* splitTree(unsigned left, unsigned right);

		/// Return the index to the vantage point that will be used for the construction of the node
		unsigned computeVP(unsigned left, unsigned right);
};

#endif
