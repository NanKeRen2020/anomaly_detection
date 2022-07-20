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

#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

#include "params.h"
#include "../Utilities/LibImages.h"
#include "../Utilities/PatchManager/patchManager.h"
#include "../Utilities/PartitionTree/prmTree.h"
#include "../Utilities/PartitionTree/forestManager.h"
#include "../Utilities/PartitionTree/vptree.h"
#include "../Utilities/PatchManager/databasePatchManager.h"
#include "../Utilities/PatchManager/imagePatchManager.h"
#include "../Utilities/PatchSearch/localRefinement.h"

#include <cstdlib>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy video;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: parameters of the method, contains:
 *			- sigma: standard deviation of the noise;
 *			- sizePatch: size of patches (sizePatch x sizePatch);
 *			- nSimilarPatches: number of similar patches;
 *			- sizeSearchWindow: size of the neighbourhood searching window;
 *			- useHomogeneousArea: if true, the trick of using homogeneous area will be used;
 *			- gamma: parameter used to determine if we are in an homogeneous area;
 *			- maxAvoid: parameter used to stop the paste trick;
 *			- beta: parameter used during the estimate of the denoised patch;
 *			- coefBaricenter: parameter to determine if the covariance matrix inversion is correct;
 *			- isFirstStep: true if it's the first step of the algorithm which is needed;
 *			- verbose: if true, print some informations, do nothing otherwise.
 *
 * @return Percentage of processed groups over number of pixels.
 **/
void computeModel(
	const std::vector<float>& im 
,	ImageSize imSize
,	const std::vector<float>& ref 
,	ImageSize refSize
,	std::vector<float> &model
,	Params const& params
);

/**
 * @brief Compute the final weighted aggregation.
 *
 * result: contain the aggregation, will contain the result at the end.
 * weight: contains the weights used during the aggregation;
 *
 * @return none.
 **/
inline void computeWeightedAggregation(
	std::vector<float>& result
,	const ImageSize imSize
,	std::vector<float> const& weight
);

#endif // MODEL_H_INCLUDED
