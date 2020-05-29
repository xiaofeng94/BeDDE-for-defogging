## Our environment

- Ubuntu 16.04
- Python v3.5
- PyTorch v0.3.1 & Torchvision v0.2.0

## Testing on sample images
1. Download the [pretrained model](https://github.com/hezhangsprinter/DCPDN#demo-using-pre-trained-model) from the original author. Place the model in `./demo_model`

2. Generate .h5 files for test images by 

`python generate_testsample.py`

If it successes, there should be .h5 files in the folder `demo_image`.

3. Test sample images by

```python demo.py --valDataroot ./demo_image --netG ./demo_model/netG_epoch_8.pth```

Afterwards, the restored images will be saved in the folder `demo_image`

## References
Original implementation: https://github.com/hezhangsprinter/DCPDN
