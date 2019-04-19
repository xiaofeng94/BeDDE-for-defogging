# BeDDE-for-defogging 
BeDDE is a real-world benchmark dataset for evaluations of defogging methods.
It consists of 208 image paris of foggy images and clear refernece images. 

Those images were collected from 34 provincial capital cities of China. 
Images in each city were collect between 9:00~10:00 over 40 days and only one image were taken in each day.

For each image pair, a manually labelled mask is provided to delineate regions with the same contents.
We evaluate defogging results on those regions.

# Download
It will be available at least by the end of May.

You may have a preview of BeDDE in this [Paper](http://sse.tongji.edu.cn/linzhang/ICME2019/BeDDE.pdf)

# Defogging methods
10 representative defogging methods, i.e.,
Fast Vis- ibility Restoration (FVR) [24], Dark Channel Prior (DCP), 
Bayesian Defogging (BayD), 
Color Attenuation Prior (CAP), 
Non-Local image Dehazing (NLD), 
MSCNN, 
DehazeNet (DeN), 
AOD-Net, 
DCPDN, 
and GFN,
are evaluated on BeDDE. 

Code of those methods in our experiments will be available here too.
