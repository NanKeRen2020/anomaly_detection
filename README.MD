Reducing Anomaly Detection in Images to Detection in Noise.

# ABOUT

* Author    : shengqiang8814@qq.com
* Licence   : GPL v3+, see gpl.txt

# OVERVIEW

This source code provides an implementation of the anomaly detection method developped in 


[1], which has been studied in details in http://www.ipol.im

detection result image
![image](https://github.com/NanKeRen2020/anomaly_detection/blob/main/detections.png)


# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux. 

- Compilation. 
Automated compilation requires the cmake and make programs, gcc/g++7.5 and OpenCV3.4.

- Library. 
This code requires the libpng, libtiff, libjpeg and fftw3 library. 

The visualization tool also requires the python cairo library. The scipy python library is also needed.


Usage:
==========

1. Download the code package and extract it. Go to that directory. 

2. Configure and compile the source code using cmake and make. 
It is recommended that you create a folder for building:

UNIX/LINUX/MAC:
$ mkdir build; cd build
$ cmake ..
$ make

Binaries will be created in build/bin folder.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). 
The number of threads used by the code is defined in model/model.cpp.

3. Usage instruction:
Running " ./main_anomaly_detection src1.png pixelNFA_ pixelNFA_ deploy.prototxt VGG_normalised.caffemodel 3" 

computes the result using the method and parameters from [1], the caffe model file can download from [2].


4. The "python visu.py src1.png 0 4 ./ pixelNFA_ radiusNFA_ detections.pdf -2 " 

produces a comprehensive visual result of the detections in a image 

"detections.png" merging all scales and where the circle radius represents 

the detection scale and the color the strength of the detection (NFA). 

White corresponds to a weak detection (NFA test value between 

$2\times 10^{-3}$ and $10^{-2}$), cyan to a good detection (NFA between $4\times 10^{-8}$ and $2\times 10^{-3}$), 

green to a very strong detection (NFA between $8.10^{-21}$ and $4\times 10^{-8}$) 

and orange to an extremely salient detection (NFA smaller than $8\times 10^{-21}$). 

Red corresponds to the detection with lowest NFA..





References
==========

[1] "Davy, Axel, et al. "Reducing Anomaly Detection in Images to Detection in Noise." 

    2018 25th IEEE International Conference on Image Processing (ICIP). IEEE, 2018.", 

[2] http://www.ipol.im/.
