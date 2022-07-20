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
 * @file nfa.cpp
 * @brief NFA functions.
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include "nfa.h"

using namespace std;

void coloredNoiseStatistic(std::vector<float>& residual, ImageSize& imSize, int R, float l, std::vector<float>& pixelNFA, std::vector<float>& radiusNFA, int M, int N, int HALFPATCHSIZE)
{
	int n0 = imSize.width;
	int n1 = imSize.height;

	float n02 = n0*n0;
	float n12 = n1*n1;

	for(int c = 0; c < imSize.nChannels; ++c)
	{
		std::vector<float> noise_gs(n0*n1);
		for(int x = 0; x < n0; ++x)
			for(int y = 0; y < n1; ++y)
				noise_gs[y + x*n1] = residual[x*imSize.nChannels + y*imSize.width*imSize.nChannels + c];

		int n1_dft = n1/2+1;

		std::vector<complex<float> > dft_rec_noise(n0*n1_dft);
		fftwf_plan plan = fftwf_plan_dft_r2c_2d(n0, n1,
				noise_gs.data(),
				reinterpret_cast<fftwf_complex*>(dft_rec_noise.data()),
				FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);

		/// Compute the measure kernel $s$ for each radius $r$
		std::vector<std::vector<float> > s(R, std::vector<float>(n0*n1));
		std::vector<std::vector<complex<float> > > ffts(R, std::vector<complex<float> >(n0*n1_dft));
		for(int r = 0; r < R; ++r)
		{
			float norms = 0;
			for(int i = 0, x = -n0/2+1; i < n0; ++i, ++x)
				for(int j = 0, y = -n1/2+1; j < n1; ++j, ++y)
				{
					if(x*x+y*y <= (r+1)*(r+1)*l*l)
					{
						s[r][j + i * n1] = 1;
						norms++;
					}
					else
						s[r][j + i * n1] = 0;
				}
			// normalize s in l_1
			for(int x = 0; x < n0; ++x)
				for(int y = 0; y < n1; ++y)
				{
					s[r][y + x*n1] = s[r][y + x*n1] / norms;
				}

			plan = fftwf_plan_dft_r2c_2d(n0, n1,
					s[r].data(),
					reinterpret_cast<fftwf_complex*>(ffts[r].data()),
					FFTW_ESTIMATE);
			fftwf_execute(plan);
			fftwf_destroy_plan(plan);
		}

		std::vector<std::vector<complex<float> > > dft_conv_noise(R, std::vector<complex<float> >(n0*n1_dft));
		for(int x = 0; x < n0; ++x)
			for(int y = 0; y < n1_dft; ++y)
			{
				/// Compute m(.,.)
				for(int r = 0; r < R; ++r)
					dft_conv_noise[r][y + x * n1_dft] = ffts[r][y + x * n1_dft] * dft_rec_noise[y + x * n1_dft];
			}

		std::vector<std::vector<float> > filtered_by_dft(R, std::vector<float>(n0*n1));
		/// Inverse back the dft
		for(int r = 0; r < R; ++r)
		{	
			fftwf_plan plan = fftwf_plan_dft_c2r_2d(n0, n1,
					reinterpret_cast<fftwf_complex*>(dft_conv_noise[r].data()),
					filtered_by_dft[r].data(),
					FFTW_ESTIMATE);
			fftwf_execute(plan);
			fftwf_destroy_plan(plan);
		}

		for(int r = 0; r < R; ++r)
		{
			std::vector<float> backup(n0*n1);
			for(int x = 0; x < n0; ++x) 
				for(int y = 0; y < n1; ++y) 
					backup[y + x*n1] = filtered_by_dft[r][y + x * n1];

			for(int x = 0; x < n0; ++x) 
			{
				int xs = (x + (n0)/2) % n0;
				for(int y = 0; y < n1; ++y) 
				{
					int ys = (y + (n1)/2) % n1;

					filtered_by_dft[r][ys + xs*n1] = backup[y + x * n1];
				}
			}
		}

		std::vector<double> sigmaphi(R);
		for(int r = 0; r < R; ++r)
		{
			// the residual is supposed to be centered. This doesn't change after the application of a Gaussian 
			sigmaphi[r] = 0.;
			for(int x = 0, ii = 0; x < n0; ++x)
				for(int y = 0; y < n1; ++y, ++ii)
					sigmaphi[r] += (filtered_by_dft[r][y + n1*x] * filtered_by_dft[r][y + n1*x] - sigmaphi[r]) / (float)(ii + 1);
			sigmaphi[r] /= (n02*n12);
			sigmaphi[r] = sqrt(std::max(sigmaphi[r], 0.));
		}

		float tolog10 = log(10);

		for(int r = 0; r < R; ++r)
		for(int x = HALFPATCHSIZE; x < (imSize.width-HALFPATCHSIZE); ++x)
		for(int y = HALFPATCHSIZE; y < (imSize.height-HALFPATCHSIZE); ++y)
		{
			float temp;
			temp = (sigmaphi[r] < 1e-8f) ? 1. : filtered_by_dft[r][y + n1*x] / (SQRT2*sigmaphi[r]*n0*n1);
			//temp = (abs(temp) > 26) ? -100000000 : std::log(imSize.nChannels * 4./3.*R*M*N*std::erfc(std::abs(temp)))/tolog10;
			temp = (abs(temp) > 26) ? -100000000 : std::log(M*N*std::erfc(std::abs(temp)))/tolog10;
            if(temp < pixelNFA[x + y*imSize.width])
            {
                pixelNFA[x + y*imSize.width] = temp;
                radiusNFA[x + y*imSize.width] = (r+1)*l;
            }
        }
	}

	float logc = std::log(imSize.nChannels);
	// Add the complement coefficient to take into account the multichannel affect
	for(int x = HALFPATCHSIZE; x < (imSize.width-HALFPATCHSIZE); ++x)
	for(int y = HALFPATCHSIZE; y < (imSize.height-HALFPATCHSIZE); ++y)
	{
		pixelNFA[x + y*imSize.width] += logc;
	}
}

