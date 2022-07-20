/*
 * Original work Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */


#include "comparators.h"

bool comparaisonFirst(
	const std::pair<float, unsigned> &i_pair1
,	const std::pair<float, unsigned> &i_pair2
){
	return i_pair1.first < i_pair2.first;
}

bool comparaisonInverseFirst(
	const std::pair<float, unsigned> &i_pair1
,	const std::pair<float, unsigned> &i_pair2
){
	return i_pair1.first > i_pair2.first;
}
