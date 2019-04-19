import os
import random
import time

import torch
import torchvision.utils as vutils
from torch.autograd import Variable

from misc import *
import dehaze22  as net

# ----------------------
test_image_root = '../my_testset/for_DCPDN'
save_root = '../my_testset/results_temp'

method_name = 'DCPDN'
useGPU = True
imageSize = 1024
# ----------------------


dataset_type = 'pix2pix_val'


valBatchSize = 1
workerNum = 1
manualSeed = random.randint(1, 10000)

ngf = 64
ndf = 64
inputChannelSize = 3
outputChannelSize= 3

netG = net.dehaze(inputChannelSize, outputChannelSize, ngf)

model_file = './demo_model/netG_epoch_8.pth'
netG.load_state_dict(torch.load(model_file))
print(netG)

netG.train()
if useGPU:
  netG.cuda()

inputHolder = torch.FloatTensor(valBatchSize, inputChannelSize, imageSize, imageSize)
if useGPU:
  inputData = Variable(inputHolder.cuda(),volatile=True)
else:
  inputData = Variable(inputHolder,volatile=True)

cityFolderNames = list()
allRootFiles = os.listdir(test_image_root)

for folderName in allRootFiles:
  if '.' != folderName and '..' != folderName and ('.' not in folderName):
    cityFolderNames.append(folderName)

totalTime = 0
totalImgCount = 0
for cityName in cityFolderNames:
  valDataroot = os.path.join(test_image_root, cityName, 'fog')

  valDataloader = getLoader(dataset_type,
                          valDataroot,
                          imageSize, # not used
                          imageSize, # not used
                          valBatchSize,
                          workerNum,
                          mean=(0.5, 0.5, 0.5), std=(0.5, 0.5, 0.5),
                          split='Train',
                          shuffle=False,
                          seed=manualSeed)

  curSaveFolder = os.path.join(save_root, cityName, method_name)
  if not os.path.exists(curSaveFolder):
    os.makedirs(curSaveFolder, mode=0o777)

  print('----- %s'%cityName)
  for i, data in enumerate(valDataloader):
    input_cpu, _, __, ___, file_name = data

    if useGPU:
      input_cpu = input_cpu.float().cuda()
    else:
      input_cpu = input_cpu.float()
    # get paired data
    inputData.data.resize_as_(input_cpu).copy_(input_cpu)

    t0 = time.time()
    x_hat, tran_hat, atp_hat, dehaze2= netG(inputData)
    elapsedT = time.time() - t0
    totalTime = totalTime + elapsedT
    print('elapsed time: %f'%elapsedT)
    totalImgCount = totalImgCount + 1

    zz=x_hat.data
    
    zz1=zz[0,:,:,:]
    vutils.save_image(zz1, os.path.join(curSaveFolder, '%s_%s.png'%(file_name[0].split('/')[-1][:-3],method_name)),
                         normalize=True, scale_each=False)

print('total time: %.4f'%totalTime)
print('total count: %d'%totalImgCount)
print('average time: %.4f'%(totalTime/totalImgCount))