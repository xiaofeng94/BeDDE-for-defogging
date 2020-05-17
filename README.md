# BeDDE & exBeDDE for dehazing evaluation
## BeDDE
BeDDE is a real-world benchmark dataset for evaluations of dehazing methods.
It consists of 208 image paris of hazy images and clear refernece images. 

Those images were collected from 34 provincial capital cities of China. 
Images in each city were collect between 9:00~10:00 over 40 days and only one image were taken in each day.

For each image pair, a manually labelled mask is provided to delineate regions with the same contents.
We evaluate dehazing results on those regions.

## exBeDDE
exBeDDE is an extension of BeDDE, designed to measure the performance of dehazing evaluation metrics. It contains 167 hazy images and 1670 dehazed images. All hazy images with their clear reference are from 12 cities of BeDDE with most images. All dehazed images are generated from the 167 hazy images using 10 dehazing methods described below.

## Visibility Index (VI) and Realness Index (RI)
We find it is more reasonable to evaluate dehazed reuslts from two separate aspects, visibility and realness, and acorrdingly propose to two criterion,  visibility index and realness index, to evaluate dehazing methods. Details of the criterion may be find in the [paper]() (coming soon) titled as _Dehazing Evaluation: Real-World Benchmark Datasets, New Criteria and Baselines_.


# Download
You can download [BeDDE](https://drive.google.com/file/d/12p-MY2ZygT5Tl8q0oFxDIUg9B5Jn042-/view?usp=sharing) and [exBeDDE](https://drive.google.com/file/d/1swAyQS-j9QNTvLwsCJgbFXnjscB86CeL/view?usp=sharing) on Google drive.

You may also get access to [BeDDE]() or [exBeDDE]() on BaiduYun disk (coming soon).

# Testing

## Our environment

- Matlab 2017b

## Test on BeDDE
1. Download BeDDE.rar and unzip it to `./Defogging_eval`

2. Set Matlab work folder to `./Defogging_eval`

3. Run `eval_defog_method.m`

Then, you will see the VI score of hazy images. You may modify the variables to try other dehazing methods or metrics.

## Test on exBeDDE

1. Download exBeDDE.rar and unzip it to `./Defogging_eval`

2. Set Matlab work folder to `./Defogging_eval`

3. Run `assess_IQA_metric.m`

It will take a while. After that, you will see the performance of our VI on hazy groups. You may modify the variables to test on dehazed groups or assess other metrics.

## Test your own dehazing method on BeDDE

1. Create a folder named after your method in each city folder of BeDDE. Put the images results of each city to the corresponding folder you just created. Name all the dehazed images after their original hazy images or make sure the name of a dehazed image starts with the name of its original hazy image.

2. Set the variable `method_name` in `eval_defog_method.m` to the name of your method folder in each city folder.

3. Set the variable `eval_method` to `VI`, `RI`, `VSI` or other metrics and run the script to get the corresponding score of your method on BeDDE.


# Dehazing methods
The selected 10 dehazing methods are 
Fast Visibility Restoration (FVR), 
Dark Channel Prior (DCP), 
Bayesian Defogging (BayD), 
Color Attenuation Prior (CAP), 
Non-Local image Dehazing (NLD), 
MSCNN, 
DehazeNet (DeN), 
AOD-Net, 
DCPDN, 
and GFN,

We test them on BeDDE as the dehazing benchmarks. Code of those methods in our experiments will be available in this repo as well. Details of each method can be found in the folder named after the method.
