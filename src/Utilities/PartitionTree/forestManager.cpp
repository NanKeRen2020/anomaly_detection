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
 * @file forestManager.cpp
 * @brief Functions managing a partition tree forest
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "forestManager.h"

ForestManager::ForestManager(PatchManager& pm_, std::vector<PartitionTree*>& forest_, int kNN_, Params prm)
{
	pm = &pm_;
	forest = forest_;
	kNN = kNN_;
	params = prm;
}

int ForestManager::retrieveFromForest(std::vector<std::pair<float, unsigned> >& indexes, unsigned pidx, bool excludeYourself)
{
	std::vector<std::vector<std::pair<float, unsigned> > > localResults(forest.size());
	std::unordered_map<unsigned, int> alreadySeen;


	// Collect candidates from the different trees
	for(int tree = 0; tree < forest.size(); tree++)
	{
		forest[tree]->retrieveFromTree(localResults[tree], pidx);
	}
	int totalDiffPatches = 0;

	std::vector<std::pair<float, unsigned> > file;

	// Search for the best match, useful if we want to exclude the current patch
	unsigned bestMatch = 0;
	float bestDist = 10000000000;
	for(int i = 0; i < localResults.size(); ++i)
	{
		for(int j = 0; j < localResults[i].size(); ++j)
		{
			if(localResults[i][j].first < bestDist)
			{
				bestMatch =localResults[i][j].second;
				bestDist = localResults[i][j].first;
			}
		}
	}

	// If the current patch has been found and we don't want it, reject it and a region around it
	if(excludeYourself && bestDist < 0.00001)
	{
		
		const ImageSize* im_sz = pm->infoIm();
		const int sPx   = params.sizePatch;
		const int excR  = params.excR;

		unsigned lx,ly,lt,lc;
		lc = bestMatch % im_sz->nChannels;
		lx = (bestMatch / im_sz->nChannels) % im_sz->width;
		ly = bestMatch / (im_sz->nChannels * im_sz->width);

		for(int dx = std::max(0,(int)lx-(excR/2)); dx < std::min(lx+(excR/2), im_sz->width); ++dx)
		for(int dy = std::max(0,(int)ly-(excR/2)); dy < std::min(ly+(excR/2), im_sz->height); ++dy)
		{
			++alreadySeen[lc + dx*im_sz->nChannels + dy*im_sz->nChannels*im_sz->width];
		}
	}

	// Create a list of unique candidates without excluded candidates 
	for(int i = 0; i < localResults.size(); ++i)
	{
		for(int j = 0; j < localResults[i].size(); ++j)
		{
			int seen = (alreadySeen[localResults[i][j].second]++);
			if(seen == 0)
			{
				++totalDiffPatches;
				file.push_back(localResults[i][j]);
			}
		}
	}
	unsigned nSimP = std::min((unsigned)kNN, (unsigned)file.size());
	std::partial_sort(file.begin(), file.begin() + nSimP, file.end(), comparaisonFirst);

	// Only save the best ones and return them
	for(int i = 0; i < nSimP; ++i)
		indexes[i] = file[i];
	return nSimP;
}

void ForestManager::updatePM(PatchManager& pm_)
{
	pm = &pm_;
	for(int i = 0; i < forest.size(); ++i)
		forest[i]->updatePM(pm_);
}
