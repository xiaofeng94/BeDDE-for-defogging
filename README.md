# BeDDE & exBeDDE for dehazing evaluation
## BeDDE
|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/chengdu_clear_rs.jpg"  width=""/>|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/chengdu_3_rs.jpg" width=""/>|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/chengdu_2_rs.jpg" width=""/>|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/chengdu_6_rs.jpg" width=""/>|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/chengdu_21_rs.jpg" width=""/>|
|:---:|:---:|:---:|:---:|:---:|
|Reference image|Light|Light|Medium|Heavy|

BeDDE (read as /ˈbedi/) is a real-world benchmark dataset for evaluations of dehazing methods.
It consists of 208 pairs of hazy images and clear refernece images. 
For each pair, a manually labelled mask is provided to delineate regions with the same contents.
We evaluate dehazing results on those regions.

## exBeDDE
exBeDDE is an extension of BeDDE, designed to **measure the performance of dehazing evaluation metrics**. It contains 167 hazy images and 1670 dehazed images with mean opinion scores labeled by people. Its hazy images come from BeDDE, and the dehazed images are generated by 10 dehazing methods.

## Visibility Index (VI) and Realness Index (RI)
We find it is more reasonable to evaluate dehazed reuslts from two separate aspects, i.e., visibility and realness, and acorrdingly propose two criteria, i.e., visibility index (VI) and realness index (RI). More details can be found in the [paper](https://ieeexplore.ieee.org/document/9099036) (or [this link](https://sse.tongji.edu.cn/linzhang/files/TIP-BeDDE.pdf)) titled as _Dehazing Evaluation: Real-World Benchmark Datasets, New Criteria and Baselines_. (Presented in Trans. Image Process.)

# Download
You can download [BeDDE](https://drive.google.com/file/d/12p-MY2ZygT5Tl8q0oFxDIUg9B5Jn042-/view?usp=sharing) and [exBeDDE](https://drive.google.com/file/d/1swAyQS-j9QNTvLwsCJgbFXnjscB86CeL/view?usp=sharing) on Google drive.

You may also get access to [the two datasets](https://pan.baidu.com/s/1lUVtdhyrvFBwl8tmwxqUXg) (key: m53p) on BaiduYun disk.

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
1. Create a folder named after your method in each city folder of BeDDE. Put the images results of each city to the corresponding folder you just created. Name all the dehazed images after their original hazy images or make sure the name of a dehazed image starts with the name of its original hazy image. Your directory tree may look like, 
```
BeDDE
├── beijing
│   ├── fog
│   │   ├── beijing_1.png
│   │   ├── beijing_2.png
│   │   ...
│   │   └── beijing_15.png
│   ├── gt
│   ├── mask
│   ├── <your_method_name>
│   │   ├── beijing_1_<your_method_name>.png
│   │   ├── beijing_2_<your_method_name>.png
│   │   ...
│   │   └── beijing_15_<your_method_name>.png
├── changsha
│   ...
├── chengdu
│   ...
...
```

2. Set the variable `method_name` in `eval_defog_method.m` to the name of your method.

3. Set the variable `eval_method` to `VI`, `RI`, `VSI` or other metrics and then run the script to get the score for your method on BeDDE.

# Benchmarks
## Dehazing methods
All methods were evaluated on BeDDE. VI shows their abilities to restore visibility, RI and LPIPS refer to realness of the results.
![dehazing benchmarks](https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/dehazing_bm.jpg)

## Dehazing metrics
All metrics were assessed on exBeDDE's hazy groups for the visibility evaluation and on dehazing groups for the realness evaluation.
|Visibility evaluation|Realness evaluation|
|:---:|:---:|
|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/metric_bm_vi.jpg" width="300"/>|<img src="https://github.com/xiaofeng94/BeDDE-for-defogging/blob/master/Defogging_eval/figures/metric_bm_ri.jpg" width="300"/>|

# Dehazing methods
The selected 14 dehazing methods adopted by exBeDDE are 
Fast Visibility Restoration (FVR), 
Dark Channel Prior (DCP), 
Bayesian Defogging (BayD), 
Color Attenuation Prior (CAP), 
Non-Local image Dehazing (NLD), 
MSCNN, 
DehazeNet, 
AOD-Net, 
DCPDN, 
GFN, 
DisentGAN,
Learning a patch quality comparator  (PQC), 
EPDN,
and GridDehazeNet.

We used offical code for all dehazing methods (except DCP). Please try the code link for those methods listed as follows: [FVR](http://perso.lcpc.fr/tarel.jean-philippe/visibility/), [DCP](https://github.com/sjtrny/Dark-Channel-Haze-Removal), [BayD](https://www.cs.drexel.edu/~kon/codeanddata/defog/), [CAP](https://github.com/JiamingMai/Color-Attenuation-Prior-Dehazing), [NLD](https://github.com/danaberman/non-local-dehazing), [MSCNN](https://github.com/rwenqi/Multi-scale-CNN-Dehazing), [DehazeNet](https://github.com/caibolun/DehazeNet), [AOD-Net](https://github.com/Boyiliee/AOD-Net), [DCPDN](https://github.com/hezhangsprinter/DCPDN) (Also try [this code](https://github.com/xiaofeng94/BeDDE-for-defogging/tree/master/DCPDN) with less issues to fix), [GFN](https://github.com/rwenqi/GFN-dehazing), [DisentGAN](https://github.com/xyang35/Disentangled-GAN), [PQC](https://github.com/san-santra/dehaze_t_comparator), [EPDN](https://github.com/ErinChen1/EPDN), [GridDehazeNet](https://github.com/proteus1991/GridDehazeNet)

<!--Code of those methods used in our experiments is available upon email requests. If one is frequently required, we will release it here as well. Current released methods: [DCPDN](https://github.com/xiaofeng94/BeDDE-for-defogging/tree/master/DCPDN)-->

# References
If our datasets and criteria are helpful, please consider citing the following papers. [1] and [2] for BeDDE. [2] for exBeDDE, VI, and RI.

[1] S. Zhao, L. Zhang, _et al._ Evaluation of defogging: A real-world benchmark dataset, a new criterion and baselines. In _ICME_, pp.1840-1845, 2019.

[2] S. Zhao, L. Zhang, _et al._ Dehazing Evaluation: Real-world Benchmark Datasets, Criteria and Baselines. _IEEE Trans. Image Process._, 28(1):6947-6962, 2020.
