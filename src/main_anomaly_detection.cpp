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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <sstream>
#include <sys/stat.h>

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/videoio.hpp"
#include <iostream>
#include <opencv2/dnn.hpp>
#include <vector>
#include <string>
#include <limits>
#include <numeric>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/normal.hpp>
#include <random>


#include "Utilities/LibImages.h"
#include "Utilities/cmd_option.h"
#include "Utilities/PatchManager/patchManager.h"
#include "Utilities/nfa.h"
#include "Utilities/PartitionTree/prmTree.h"
#include "Utilities/PartitionTree/vptree.h"
#include "model/model.h"

#define R 2
#define l 2

// Define the fixed parameters here
#define NBTREES 4
#define hsigma 10.

using namespace std;
using namespace boost::math;

/**
 * @file   main.cpp
 * @brief  Main executable
 *
 * @author THIBAUD EHRET  <ehret.thibaud@gmail.com>
 **/

template<typename _Tp>
vector<_Tp> convertMat2Vector(const cv::Mat &mat)
{
	return (vector<_Tp>)(mat.reshape(1, 1));//通道数不变，按行转为一行
}
 
/****************** vector转Mat *********************/
template<typename _Tp>
cv::Mat convertVector2Mat(vector<_Tp> v, int channels, int rows)
{
	cv::Mat mat = cv::Mat(v);//将vector变成单列的mat
	cv::Mat dest = mat.reshape(channels, rows).clone();//PS：必须clone()一份，否则返回出错
	return dest;
}


cv::Mat norm_0_255(const cv::Mat& src)
{
    cv::Mat dst;
    switch (src.channels())
    {
    case 1:
        cv::normalize(src, dst, 1, 0, cv::NORM_MINMAX);
        break;
    case 3:
        cv::normalize(src, dst, 1, 0, cv::NORM_MINMAX);
        break;
    default:
        src.copyTo(dst);
        break;
    }
    return dst;
}

void cumsum_hist(const cv::Mat& hist, cv::Mat& cum_sum_hist)
{
        cum_sum_hist.at<float>(0) = hist.at<float>(0);
        for(size_t k = 1; k < hist.rows; ++k)
            cum_sum_hist.at<float>(k) = hist.at<float>(k) + cum_sum_hist.at<float>(k-1);
}

float get_score(float alpha, const cv::Mat residual, bool laplace, cv::Mat& data_mat, cv::Mat& cdf_data, float& scale)
{
        data_mat = cv::Mat::zeros(residual.size(), CV_32FC1);
        std::vector<double> data_vec;
        for (auto it = residual.begin<uchar>(); it != residual.end<uchar>(); ++it)
        {
            *(data_mat.begin<float>() + (it - residual.begin<uchar>())) = std::pow(std::abs(*it), alpha);

        }

        //std::cout<< "zero mat: " << zero_mat.size() << std::endl << zero_mat << std::endl;
        cv::Mat mat_mean, mat_stddev;
        cv::meanStdDev(data_mat, mat_mean, mat_stddev); 
        scale  = mat_stddev.at<double>(0, 0)/2;

        std::cout << "scale value0 : " << scale << std::endl;

    	double sum = std::accumulate(data_vec.begin(), data_vec.end(), 0.0);
	    double mean =  sum / data_vec.size(); //均值
	    double accum  = 0.0;
	    std::for_each (data_vec.begin(), data_vec.end(), [&](const double d) {
		       accum  += (d-mean)*(d-mean);
	    });
        
        if (laplace)
	    scale = sqrt(accum/(data_vec.size())/2); //方差
        else
        scale = sqrt(accum/(data_vec.size())); //方差
        std::cout << "scale value1: " << scale << std::endl;
        if (scale == 0 ) scale = 1;
        
        boost::math::laplace_distribution<double> lap_dis(0, scale);
        boost::math::normal_distribution<double> norm_dis(0, scale);

        cdf_data = cv::Mat::zeros(data_mat.size(), data_mat.type());
        for (auto it = data_mat.begin<float>(); it != data_mat.end<float>(); ++it)
        {
            if (laplace)
            *(cdf_data.begin<float>() + (it - data_mat.begin<float>())) = boost::math::cdf(lap_dis, *it);
            else
            *(cdf_data.begin<float>() + (it - data_mat.begin<float>())) = boost::math::cdf(norm_dis, *it);
     
        }

        cv::Mat hist, cumsum;
        float range[] = { 0.01, 1.01 };
        const float* histRange[] = { range };
        const int channels[] = { 0 };
        const int histSize[] = { 20 };
        cv::calcHist(&cdf_data, 1, channels, cv::Mat(), hist, 1, histSize, histRange, true, false);
        hist.at<float>(9) = 0;
        hist = hist / cv::sum(hist)[0];
        cv::Mat culsum_hist(hist.size(), hist.type());
        culsum_hist.at<float>(0) = hist.at<float>(0);
        for(size_t k = 1; k < hist.rows; ++k)
        {
           culsum_hist.at<float>(k) = hist.at<float>(k) + culsum_hist.at<float>(k-1);
        }
        
        
        cv::Mat one_hist = cv::Mat::ones(hist.size(), hist.type());
        cv::Mat cumsum_one(hist.size(), hist.type());
        cumsum_hist(one_hist*0.05, cumsum_one);
        cv::subtract(culsum_hist, cumsum_one, culsum_hist);
        for (auto it = culsum_hist.begin<float>(); it != culsum_hist.end<float>(); ++it)
        {

            *it = (*it)*(*it);
        }
        return cv::sum(culsum_hist)[0];
}



cv::Mat convert(const cv::Mat& residual_mat, int component)
{
    std::vector<float> alphas{0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4};
    float best_score = std::numeric_limits<float>::max();
    //std::cout<< "residual: " << residual.size() << std::endl << residual << std::endl;
    boost::math::normal_distribution<double> norm_dis;
    cv::Mat data_mat, cdf_data, best_data, residual, residual_gassian;
    float best_alpha, scale, score;
    std::string best_model;
    std::vector<cv::Mat> residuals;
    cv::split(residual_mat, residuals);

    for (int i = 0; i < residual_mat.channels(); ++i)
    {
       residual = residuals[component - 1 - i].clone();

    for (auto alpha: alphas)
    {


        score = get_score(alpha, residual, 1, data_mat, cdf_data, scale);
        
        if (score < best_score)
        {
            
            cv::Mat one_data = cv::Mat::ones(data_mat.size(), data_mat.type());
            cv::Mat sf_data = cv::Mat::zeros(data_mat.size(), data_mat.type());;
            //cv::subtract(one_data, cdf_data, sf_data);
            for (auto it = data_mat.begin<float>(); it != data_mat.end<float>(); ++it)
            {
                
                //*(it) = boost::math::quantile(norm_dis_1, 0.5);
                if (*it > 0) *( sf_data.begin<float>() + (it - data_mat.begin<float>()) )
                = 0.5*exp(-1*(*it)/scale);
                else *( sf_data.begin<float>() + (it - data_mat.begin<float>()) )
                = 1 - 0.5*exp((*it)/scale);   


                *( sf_data.begin<float>() + (it - data_mat.begin<float>()) ) 
                = boost::math::quantile( boost::math::complement(norm_dis, 
                  *( sf_data.begin<float>() + (it - data_mat.begin<float>())) ) );

                
            }

            cv::Mat mask_sf = cv::Mat::zeros(data_mat.size(), data_mat.type() );
            for (auto it = sf_data.begin<float>(); it != sf_data.end<float>(); ++it)
            {
                if (*it > 0)
                *(mask_sf.begin<float>() + (it - sf_data.begin<float>())) = 1;                
            }
            best_data = sf_data.mul(mask_sf);      


            cv::Mat ppf_data = cv::Mat::zeros(data_mat.size(), data_mat.type() ); 

            ppf_data = ppf_data.mul(~mask_sf);
            cv::add(best_data, ppf_data, best_data);
 
            best_alpha = alpha;
            best_score = score;
            
            best_model = "laplace";
        }


        score = get_score(alpha, residual, 0, data_mat, cdf_data, scale);
        if (score < best_score)
        {
            best_score = score;
            best_alpha = alpha;

            best_data = data_mat/scale;
  
            best_model = "gaussian";
        }

        
    }

        std::cout << "the best is " << best_model << ", " << best_alpha << ", " << best_score << std::endl;
        residuals[i] = best_data;
        cv::merge(residuals, residual_gassian);
    }

    return residual_gassian;


}


static cv::Mat formatImagesForPCA(const std::vector<cv::Mat> &data)
{
    cv::Mat dst(static_cast<int>(data.size()), data[0].rows*data[0].cols, CV_32F);
    for(unsigned int i = 0; i < data.size(); i++)
    {
        cv::Mat image_row = data[i].clone().reshape(1,1);
        cv::Mat row_i = dst.row(i);
        image_row.convertTo(row_i, CV_32F);
    }
    return dst;
}

cv::Mat get_nn_output(const cv::Mat& image, const std::string& proto_path, const std::string& model_path, int layer_num, int& row)
{
 
    std::cout << "image: " << image.size() << ", " << image.channels() << ", "
              << image.size[0] << ", " << image.size[1] << ", " << image.size[2]
              << ", " << image.size[3] << ", " << image.depth() << std::endl;
    //cv::normalize(image, image, 1, 0);
    //image.convertTo(image, CV_16SC3);
    std::ifstream proto_file(proto_path);
    std::string str;
    std::vector<std::string> strs;
    while (std::getline(proto_file, str))
    {
        strs.push_back(str);
    }
    std::cout << "proto file: " << strs.size() << ", " << strs[3] << ", "
              << strs[4] << ", " << strs[5] << std::endl;
    strs[3] = "input_dim: " + std::to_string(image.channels());
    strs[4] = "input_dim: " + std::to_string(image.cols);
    strs[5] = "input_dim: " + std::to_string(image.rows);
    std::ofstream out_proto_file(proto_path);
    for (auto str: strs)
    {
        out_proto_file << str << std::endl;
    }
    out_proto_file.close();

    cv::dnn::Net net = cv::dnn::readNetFromCaffe(proto_path, model_path);
    net.setPreferableBackend(0);
    net.setPreferableTarget(0);
    cv::Mat blob;
    cv::dnn::blobFromImage(image, blob, 255.0, cv::Size(475, 213), cv::Scalar( 0.40760392,  0.45795686,  0.48501961), false, false);
    net.setInput(blob);
    std::vector<cv::String> layer_names = net.getLayerNames();
    layer_names = std::vector<cv::String>{"conv1_1", "conv1_2", "pool1","conv2_1", "conv2_2", "pool2",
                    "conv3_1", "conv3_2", "conv3_3", "conv3_4", "pool3",
                    "conv4_1", "conv4_2", "conv4_3", "conv4_4"};
    cv::Mat result = net.forward(layer_names[layer_num]);
    std::cout << "result: " << result.size() << ", " << result.channels() << ", "
              << result.size[0] << ", " << result.size[1] << ", " << result.size[2]
              << ", " << result.size[3] << ", " << result.depth() << std::endl;
  
    std::vector<cv::Mat> images;
    cv::dnn::imagesFromBlob(result, images);
    std::cout << "blob image: " << images.size() << ", " << images[0].size() << ", " 
              << images[0].channels() << std::endl;
    row = images[0].rows;
    std::vector<cv::Mat> channels;
    cv::split(images[0], channels);
    std::cout << "channels: " << channels.size() << ", " << channels[0].size() << ", "
              << channels[0].channels() << ", "  << channels[0].size[0] << ", " 
              << channels[0].size[1] << ", " << channels[0].size[2] << ", " << channels[0].size[3] << std::endl;  
    

    cv::Mat data = formatImagesForPCA(channels);
    std::cout << "data info: " << data.size() << ", " << data.channels() << ", "
              << data.size[0] << ", " << data.size[1] << ", " << data.size[2]
              << ", " << data.size[3] << ", " << data.depth() << std::endl; 
    return data;
}

cv::Mat get_pca_componet(const cv::Mat& data, int row, int component, cv::Mat& reshape_result)
{
             
    // perform PCA
    cv::PCA pca(data, cv::Mat(), cv::PCA::DATA_AS_COL, component); // trackbar is initially set here, also this is a common value for retainedVariance
    cv::Mat mean = pca.mean.clone();
    cv::Mat eigen_values = pca.eigenvalues.clone();
    cv::Mat eigenv_ectors = pca.eigenvectors.clone();

    for (int i = 0; i < pca.eigenvectors.rows; ++i)
    {
        pca.eigenvectors.rowRange(i, i+1) = pca.eigenvectors.rowRange(i, i+1)/std::sqrt(*( pca.eigenvalues.begin<float>() + i));
        if (i == 1 || i == 4)
        pca.eigenvectors.rowRange(i, i+1) = -1*pca.eigenvectors.rowRange(i, i+1);       
    }

    // Demonstration of the effect of retainedVariance on the first image
    cv::Mat point = pca.project(data); // project into the eigenspace, thus the image becomes a "point"
    std::cout << "pca point: " << point.size() << ", " << point.channels() << ", " 
            << point.size[0] << ", " << point.size[1] << ", " << point.size[2] << ", " 
            << point.size[3] << std::endl;      
    

    return point;

}



int main(int argc, char **argv)
{
   
	const bool same         = true;
	//! Model parameters
	const int patch_size  = 8;
	const int num_patches  = 16;
    int layer_num = atoi(argv[4]);


	cv::Mat image = cv::imread(argv[1]);
    cv::Mat image_gau, scale_image, reshape_result;
    int width, height, row, component;
    
    for (int n = 0; n < 4; ++n)
    {
        if (image.cols < 150 || image.rows < 150) continue;        
        if (layer_num < 0)
        {
            row = image.rows;
            std::vector<cv::Mat> channels;
            cv::split(image, channels);
            scale_image = formatImagesForPCA(channels);  
        }
        else
        {
	        scale_image = get_nn_output(image, argv[5], argv[6], layer_num, row);

        }


    component = std::min(5, scale_image.size[0]/2);  
    cv::Mat pca_proj = get_pca_componet(scale_image, row, component, reshape_result);
    std::cout << "pca_proj info: " << reshape_result.size() << ", " << reshape_result.channels() << ", " 
            << reshape_result.size[0] << ", " << reshape_result.size[1] << ", " << reshape_result.size[2] << ", " 
            << reshape_result.size[3] << std::endl;  
    
    
    ImageSize imSize, refSize;
    imSize.width = reshape_result.cols;
	imSize.height = reshape_result.rows;
	imSize.nChannels = reshape_result.channels();
    imSize.wh = imSize.width*imSize.height;
	imSize.whc = imSize.width*imSize.height*imSize.nChannels;

    refSize.width = reshape_result.cols;
	refSize.height = reshape_result.rows;
	refSize.nChannels = reshape_result.channels();
    refSize.wh = refSize.width*refSize.height;
	refSize.whc = refSize.width*refSize.height*refSize.nChannels;

	std::vector<float> candidate;
	std::vector<float> reference;


    if(imSize.nChannels != refSize.nChannels)
    {
		fprintf(stderr, "%s: the number of channels needs to be the same for the image and the reference.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;

    }

	//! Compute denoising default parameters
	Params prms;
	prms.sizePatch = patch_size;
	prms.nbTrees = NBTREES;
	prms.nSimilarPatches = num_patches;
	prms.same = same;
	prms.excR = 2*patch_size;
	prms.h = hsigma;

	std::vector<float> model(imSize.whc);
	std::vector<float> residual(imSize.whc);

	//! Run denoising algorithm
	computeModel(candidate, imSize, reference, refSize, model, prms);

	//! Compute the residual
	computeDiff(candidate, model, residual);

    cv::Mat residual_mat = convertVector2Mat(residual, imSize.nChannels, imSize.height);
    std::cout << "residual_mat info: " << residual_mat.size() << ", " << residual_mat.channels() << ", " 
             << residual_mat.size[0] << ", " << residual_mat.size[1] << ", " << residual_mat.size[2] << ", " 
             << residual_mat.size[3] << std::endl; 
    
    cv::Mat residual_gaussian = convert(residual_mat, component);

	//! Declarations
	std::vector<float> residual_gau(residual_gaussian.reshape(1, 1));
    std::vector<float> original(image.reshape(1, 1));
	ImageSize origSize;
    imSize.width = residual_gaussian.cols;
	imSize.height = residual_gaussian.rows;
	imSize.nChannels = residual_gaussian.channels();
    imSize.wh = imSize.width*imSize.height;
	imSize.whc = imSize.width*imSize.height*imSize.nChannels;

    origSize.width = image.cols;
	origSize.height = image.rows;
	origSize.nChannels = image.channels();
    origSize.wh = origSize.width*origSize.height;
	origSize.whc = origSize.width*origSize.height*origSize.nChannels;


	//! Detection in the residual
	std::vector<float> pixelNFA(imSize.width*imSize.height);
	std::vector<float> radiusNFA(imSize.width*imSize.height);
	ImageSize nfaSize;
	nfaSize.width = imSize.width;
	nfaSize.height = imSize.height;
	nfaSize.nChannels = 1;
	nfaSize.wh = nfaSize.width * nfaSize.height;
	nfaSize.whc = nfaSize.width * nfaSize.height * nfaSize.nChannels;
	coloredNoiseStatistic(residual, imSize, R, l, pixelNFA, radiusNFA, origSize.width, origSize.height, patch_size/2);
    std::string pixelNFA_image = argv[2] + std::to_string(n) + "_cpp.tiff";
    std::string radiusNFA_image = argv[3] + std::to_string(n) + "_cpp.tiff";
	saveImage(pixelNFA_image.c_str(), pixelNFA, nfaSize);
	saveImage(radiusNFA_image.c_str(), radiusNFA, nfaSize);

        cv::GaussianBlur(image, image_gau, cv::Size(5, 5), 1.39, 1.39);
        cv::pyrDown(image_gau, image);
    }
    


    
	return EXIT_SUCCESS;
}
