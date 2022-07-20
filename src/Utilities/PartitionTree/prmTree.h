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

#ifndef PRMTREE_HEADER_H
#define PRMTREE_HEADER_H

#include "../../model/params.h"

/**
 * @brief Structures of parameters useful to create a partition tree
 *
 * @param epsilon: Threshold to which the computation of a tree is stop (set to 0 for best use for now)
 * @param prms: Parameters from the denoising algorithms
 * @param nbVP: Number of candidate vantage points used during the computation of the VP-tree
 * @param subS: Number of distance to compute to estimate the median distance during the computation of the VP-tree
 * @param rand: Number of dimension to use for a randomized KD-tree
 **/
struct PrmTree
{
	float epsilon;
	Params prms;
	int nbVP;
	int subS;	
	int rand;
};
#endif
