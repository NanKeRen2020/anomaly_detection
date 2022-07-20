/*
 * Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
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
 * @file 
 * @brief 
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "model.h"
#include "../Utilities/comparators.h"

// colors
#define ANSI_BLK  "\x1b[30m"
#define ANSI_RED  "\x1b[31m"
#define ANSI_GRN  "\x1b[32m"
#define ANSI_YLW  "\x1b[33m"
#define ANSI_BLU  "\x1b[34m"
#define ANSI_MAG  "\x1b[35m"
#define ANSI_CYN  "\x1b[36m"
#define ANSI_WHT  "\x1b[37m"
#define ANSI_BBLK "\x1b[30;01m"
#define ANSI_BRED "\x1b[31;01m"
#define ANSI_BGRN "\x1b[32;01m"
#define ANSI_BYLW "\x1b[33;01m"
#define ANSI_BBLU "\x1b[34;01m"
#define ANSI_BMAG "\x1b[35;01m"
#define ANSI_BCYN "\x1b[36;01m"
#define ANSI_BWHT "\x1b[37;01m"
#define ANSI_RST  "\x1b[0m"

#define SQRT2_INV 0.7071067811865475
#define NTHREAD 32

/**
 * @brief Compute the final weighted aggregation.
 *
 * result: contain the aggregation, will contain the result at the end.
 * weight: contains the weights used during the aggregation;
 *
 * @return : none.
 **/
inline void computeWeightedAggregation(
	std::vector<float>& result
,	const ImageSize imSize
,	std::vector<float> const& weight
){
	for (unsigned y = 0; y < imSize.height; ++y)
	for (unsigned x = 0; x < imSize.width; ++x)
	for (unsigned c = 0; c < imSize.nChannels; ++c)
		result[c + x*imSize.nChannels + y*imSize.nChannels*imSize.width] /= weight[x + y*imSize.width];
}

void computeModel(
	const std::vector<float>& im
,	ImageSize imSize
,	const std::vector<float>& ref
,	ImageSize refSize
,	std::vector<float> &model
,	Params const& params
){
	//! Parameters initialization
	const unsigned sPx = params.sizePatch;

	//! Weight sum per pixel
	std::vector<float> weight(imSize.width*imSize.height, 0.f);

	const unsigned patch_num = params.nSimilarPatches;
	std::vector<std::pair<float, unsigned> > index(patch_num);

	//! Create the VPtree forest
	ImagePatchManager pm(ref, refSize, sPx);
	std::vector<PartitionTree*> trees; 
	PrmTree paramTree;
	paramTree.epsilon = 0.; 
	paramTree.nbVP = 1;
	paramTree.subS = 1000;
	paramTree.rand = 5;
    paramTree.prms = params;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
	shared(trees) 
#endif
	for(int i = 0; i < params.nbTrees; ++i)
		trees.push_back(new VPTree(pm, paramTree));	

	ForestManager fm(pm, trees, params.nSimilarPatches, params);

	DatabasePatchManager dbpm(ref, refSize, sPx);
	dbpm.setCurrentImage(im, imSize);
    fm.updatePM(dbpm);

	LocalRefinement lpp(dbpm, fm, patch_num, params);


	// Iterate on patches this way so there is no overlap during the aggregation. This allow to compute all patches in parallel.
	for (unsigned shifty = 0; shifty < sPx; ++shifty)
	for (unsigned shiftx = 0; shiftx < sPx; ++shiftx)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(NTHREAD) schedule(dynamic) \
        shared(lpp, model, weight) \
        firstprivate(index)
#endif
        for (unsigned py = shifty; py < imSize.height - sPx + 1; py+=sPx)
            for (unsigned px = shiftx; px < imSize.width - sPx + 1; px+=sPx)
            {
                const unsigned ij  = px + imSize.width*py;
                const unsigned ij3 = (px + imSize.width*py)*imSize.nChannels;

                //! Search for similar patches around the reference one
                unsigned nSimP = lpp.estimateSimilarPatches(index, ij3, params.same);

                // Denormalize the distances
                for(unsigned k = 0; k < nSimP; ++k)
                {
                	index[k].first *= 255.f;
                	index[k].first *= index[k].first;
                }

                float norm_factor = 0.f;
                float h2 = params.h*params.h;
                for (unsigned k = 0; k < nSimP; ++k)
                {
                    index[k].first = exp(-index[k].first/h2);
                    norm_factor += index[k].first;
                }
                if(norm_factor > 0)
                    norm_factor = 1.f / norm_factor;
                else
                {
                    for(int i = 0; i < nSimP; ++i)
                        index[i].first = 1;
                    norm_factor = 1.f / nSimP;
                }

                // Compute the average patch while aggregating it.
                for (unsigned y = 0; y < sPx; y++)
                for (unsigned x = 0; x < sPx; x++)
                {
                    for (unsigned c = 0; c < imSize.nChannels; c++)
                    {
                        for (unsigned k = 0; k < nSimP; ++k)
                            model[ij3 + y*imSize.width*imSize.nChannels + x*imSize.nChannels + c] +=
                                index[k].first * ref[index[k].second  + y*refSize.width*refSize.nChannels + x*refSize.nChannels + c] * norm_factor;
                    }

                    weight[ij + y*imSize.width + x] += 1;
                }
            }
    }

	//! Weighted aggregation
	computeWeightedAggregation(model, imSize, weight);

	//! Delete trees
	for(int i = 0; i < params.nbTrees; ++i)
		delete trees[i];
}
