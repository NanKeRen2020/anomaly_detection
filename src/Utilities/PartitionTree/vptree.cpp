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
 * @file vptree.cpp
 * @brief Specific type of partition tree: VP-tree
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "vptree.h"

VPTree::VPTree(PatchManager& pm_, PrmTree& paramTree) : nbSelectionVP(paramTree.nbVP), subSampling(paramTree.subS)
{
	param = paramTree.prms;
	int k = param.nSimilarPatches; 
#ifdef _OPENMP
    std::srand(unsigned(std::time(NULL)) ^ omp_get_thread_num());
#else
	std::srand(unsigned(std::time(NULL))); 
#endif
	kNN = k;
	// We insure that we have at least k elements in the bins with a median split, keep at least a couple of elements so to avoid any problems
	leaf_max_size = std::max(2*k + 1, 5);
	pm = &pm_;
	int nbPatches = pm->getNbPatches();
	idxPatches.resize(nbPatches);
	pm->getAllPatches(idxPatches);
	root = splitTree(0, idxPatches.size()-1); 
}

static int rand_gen(int i)
{
	return std::rand() % i;
}

unsigned VPTree::partition(unsigned left, unsigned right, unsigned pivot, std::vector<float>& distances, unsigned n, unsigned l)
{
	float valuePivot = distances[pivot - l];

	unsigned temp = idxPatches[pivot];
	idxPatches[pivot] = idxPatches[right];
	idxPatches[right] = temp;
	float tempF = distances[pivot - l];
	distances[pivot - l] = distances[right - l];
	distances[right - l] = tempF;


	unsigned store = left;
	int allEqual = 0;

	// move all elements between $left$ and $right$ so they are at the correct position according to the pivot
	for(unsigned i = left; i <= right; ++i)
	{
		float candidate = distances[i - l];

		// trick in case all elements are equal to speed up the computation
		if(std::abs(candidate - valuePivot) < 0.000001)
			allEqual++;

		if(candidate < valuePivot)  
		{
			temp = idxPatches[store];
			idxPatches[store] = idxPatches[i];
			idxPatches[i] = temp;
			tempF = distances[store - l];
			distances[store - l] = distances[i - l];
			distances[i - l] = tempF;
			++store;
		}
	}

	temp = idxPatches[store];
	idxPatches[store] = idxPatches[right];
	idxPatches[right] = temp;
	tempF = distances[store - l];
	distances[store - l] = distances[right - l];
	distances[right - l] = tempF;

	if(allEqual >= (right - left) / 2)
		return n;

	// return the index where the pivot is
	return store;
}

// Fast computation of the median along a specific dimension for a set of element delimited by left and right 
float VPTree::quickselect(unsigned l, unsigned r, std::vector<float>& distances)
{
	unsigned right = r;
	unsigned left = l;
	if(l == r)
	{
		return distances[0];
	}	

	unsigned n = l + (r - l) / 2;
	int compteur = 0;
	while(1==1)
	{
		compteur++;
		unsigned pivot = left + (right - left) / 2;
		pivot = partition(left, right, pivot, distances, n, l);

		if(n  == pivot)
		{
			return distances[n - l];
		}
		else if(n < pivot)
		{
			right = pivot - 1;
		}
		else
		{
			left = pivot + 1;
		}
	}
}

VPTree::Node* VPTree::splitTree(unsigned left, unsigned right)
{
	Node* tempNode = new Node();

	// Test if we are creating a leaf or not
	if( (right - left + 1) <= leaf_max_size)
	{
		// It is a leaf, we set all informations regarding this leaf and finish
		tempNode->child1 = NULL;
		tempNode->child2 = NULL;

		tempNode->left = left;
		tempNode->right = right;
		tempNode->vantagePoint = 0;
		tempNode->splitDistance = 0.;
	}
	else
	{
		// It is not a leaf, create the split and the two children nodes
		std::vector<float> distances(right - left + 1);
		unsigned vantage = computeVP(left, right);

		float maxDist = 0.;
		for(unsigned i = 0, j = left; j <= right; ++i,++j)
		{
			distances[i] = pm->distance(vantage, idxPatches[j]);
			if(distances[i] > maxDist)
				maxDist = distances[i];
		}
		tempNode->left = left;
		tempNode->right = right;

		tempNode->vantagePoint = vantage;

		tempNode->splitDistance = quickselect(left, right, distances);
		// Create the nodes for the children 
		tempNode->child1 = splitTree(left, left + (right - left) / 2);
		tempNode->child2 = splitTree(left + (right - left) / 2 + 1, right);
	}
	
	// Return the pointer to the constructed node
	return tempNode;
}


unsigned VPTree::computeVP(unsigned left, unsigned right)
{
	std::random_shuffle(idxPatches.begin() + left, idxPatches.begin() + right + 1, rand_gen);

	std::vector<unsigned> elementToConsider(nbSelectionVP);
	for(int i = 0; i < nbSelectionVP; ++i)
		elementToConsider[i] = idxPatches[left + (unsigned)i];

	float maxVar = 0.;
	int interestingPoint = 0;
	if(nbSelectionVP == 1)
	{
		// Construct a VP tree for a forest, in this case no other computation need to be done. We choose a random VP
		return elementToConsider[0];
	}
	else
	{
		// Otherwise select the best candidates from a subsample of possible vantage points
		std::vector<unsigned> elementToConsider(nbSelectionVP);
		std::vector<unsigned> subSample(std::min(right-left+1, (unsigned)subSampling));

		std::vector<float> distancesToSubs(std::min(right - left + 1, (unsigned)subSampling), 0.);

		std::random_shuffle(idxPatches.begin() + left, idxPatches.begin() + right + 1, rand_gen);
		for(int i = 0; i < subSample.size(); ++i)
			subSample[i] = idxPatches[left + (unsigned)i];

		// Constructing a normal VP tree
		for(int i = 0; i < nbSelectionVP; ++i)
		{
			//Compute the distances to the point
			for(unsigned j = 0; j < subSample.size(); ++j)
			{
				distancesToSubs[j] = pm->distance(elementToConsider[i], subSample[j]); 
			}

			// compute the median
			float medianTemp = quickselect(left, left + distancesToSubs.size() - 1, distancesToSubs);	

			// compute the variance 
			float var = 0.;
			for(int j = 0; j < distancesToSubs.size(); ++j)
				var += (distancesToSubs[j] - medianTemp) * (distancesToSubs[j] - medianTemp);
			if(var > maxVar)
			{
				maxVar = var;
				interestingPoint = i;
			}
		}
	}	
	// Return the best VP candidate
	return elementToConsider[interestingPoint];
}

void VPTree::retrieveFromTree(std::vector<std::pair<float, unsigned> >& localResult, const unsigned pidx)
{
	Node* currentNode = this->root; 
	// Find the leaf in which the patch resides
	while(currentNode->child1 != NULL && currentNode->child2 != NULL)
	{
		if(pm->distance(pidx, currentNode->vantagePoint) > currentNode->splitDistance)
			currentNode = currentNode->child2;
		else
			currentNode = currentNode->child1;
	}

	const int sPx = param.sizePatch;

	for(unsigned i = currentNode->left, j = 0; i <= currentNode->right; ++i,++j)
	{
		float dist = pm->distance(pidx, idxPatches[i]);
		localResult.push_back(std::make_pair(dist, idxPatches[i]));
	}
}
