# BeDDE & exBeDDE for defogging evaluation
BeDDE is a real-world benchmark dataset for evaluations of defogging methods.
It consists of 208 image paris of foggy images and clear refernece images. 

Those images were collected from 34 provincial capital cities of China. 
Images in each city were collect between 9:00~10:00 over 40 days and only one image were taken in each day.

For each image pair, a manually labelled mask is provided to delineate regions with the same contents.
We evaluate defogging results on those regions.

What's more, in order to measure the performance of defogging evaluation metrics, we build an extension of BeDDE, exBeDDE, using 167 foggy images and 1670 defogged images. All foggy images with their clear reference are from 12 cities of BeDDE with most images. All defogged images are generated from the 167 foggy images using 10 defogging methods described below.

We find it is more reasonable to evaluate defogged reuslts from two separate aspects, visibility and realness, and acorrdingly propose to two criterion,  visibility index and realness index, to evaluate defogging methods. Details of the criterion may be find in the [paper]() (coming soon) titled as _Defogging Evaluation: Real-World Benchmark Datasets, New Criteria and Baselines_.


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

Then, you will see the VI score of foggy images. You may modify the variables to try other defogging methods or metrics.

## Test on exBeDDE

1. Download exBeDDE.rar and unzip it to `./Defogging_eval`

2. Set Matlab work folder to `./Defogging_eval`

3. Run `assess_IQA_metric.m`

It will take a while. After that, you will see the performance of our VI on foggy groups. You may modify the variables to test on defogged groups or assess other metrics.

## Test your own defogging method on BeDDE

1. Create a folder named after your method in each city folder of BeDDE. Put the images results of each city to the corresponding folder you just created. Name all the defogged images after their original foggy images or make sure the name of a defogged image starts with the name of its original foggy image.

2. Set the variable `method_name` in `eval_defog_method.m` to the name of your method folder in each city folder.

3. Set the variable `eval_method` to `VI` or `RI` and run the script to get the VI or RI score of your method on BeDDE.


# Defogging methods
The selected 10 defogging methods are 
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

We also test them on BeDDE as the defogging benchmarks. Code of those methods in our experiments will be available here too.


## DCPDN

### Our environment

- Ubuntu 16.04
- Python v3.5
- PyTorch v0.3.1 & Torchvision v0.2.0

### Testing on sample images
1. Download the [pretrained model](https://github.com/hezhangsprinter/DCPDN#demo-using-pre-trained-model) from the original author. Place the model in `./demo_model`

2. Generate .h5 files for test images by 

`python generate_testsample.py`

If it successes, there should be .h5 files in the folder `demo_image`.

3. Test sample images by

```python demo.py --valDataroot ./demo_image --netG ./demo_model/netG_epoch_8.pth```

Afterwards, the restored images will be saved in the folder `demo_image`
