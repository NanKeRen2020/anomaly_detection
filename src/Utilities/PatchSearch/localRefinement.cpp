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

#include "localRefinement.h"

#define LR_X 9
#define LR_Y 9

int LocalRefinement::estimateSimilarPatches(std::vector<std::pair<float, unsigned> > &index, const unsigned pidx, bool excludeYourself)
{

	int sWx = LR_X;
	int sWy = LR_Y;
	const int sPx   = params.sizePatch;

	const ImageSize* im_sz = pm->infoIm();
	// Get the best patches from the forest
	int finalValue = fm->retrieveFromForest(index, pidx, excludeYourself);

	std::unordered_map<unsigned, int> alreadySeen;
	std::vector<std::pair<float, unsigned> > listMatches(0);
	for(int i = 0; i < kNN; ++i)
	{
		++alreadySeen[index[i].second];
		listMatches.push_back(index[i]);	
	}

	// Update the list using a local search centered around each patches (local refinement step)
	for(int i = 0; i < kNN; ++i)
	{
		// Define the small search region around the current match
		unsigned pc = index[i].second % im_sz->nChannels;
		unsigned px = (index[i].second / im_sz->nChannels) % im_sz->width;
		unsigned py = index[i].second / (im_sz->nChannels * im_sz->width);

		unsigned rangex[2];
		unsigned rangey[2];

		rangex[0] = std::max(0, (int)px - (sWx-1)/2);
		rangey[0] = std::max(0, (int)py - (sWy-1)/2);

		rangex[1] = std::min((int)im_sz->width  - sPx, (int)px + (sWx-1)/2);
		rangey[1] = std::min((int)im_sz->height - sPx, (int)py + (sWy-1)/2);

		sWx = rangex[1] - rangex[0] + 1;
		sWy = rangey[1] - rangey[0] + 1;

		// Search for the nearest neighbors
		for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			unsigned currentPatch = pc + im_sz->nChannels*im_sz->width*qy + im_sz->nChannels*qx;
			int seen = alreadySeen[currentPatch]++;
			if(seen == 0)
				listMatches.push_back(std::make_pair(pm->distance(pidx, currentPatch), currentPatch));
		}
	}
    

	// Only keep at most kNN elements
	unsigned nSimP = std::min(kNN, (unsigned)listMatches.size());
	std::partial_sort(listMatches.begin(), listMatches.begin() + nSimP, listMatches.end(), comparaisonFirst);
	for (unsigned n = 0; n < kNN; n++)
		index[n] = listMatches[n];
	return nSimP;
}
