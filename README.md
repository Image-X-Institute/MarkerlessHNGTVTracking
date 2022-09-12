# MarkerlessHNGTVTracking
Code used to create the results described in the paper "A conditional generative adversarial network approach for segmenting head and neck tumors in kV images acquired during radiation therapy"

For a more general version of the CT Head Deformation code see the the repository click [here](https://github.com/ACRF-Image-X-Institute/CTHeadDeformation) 

For any questions email mark.gardner@sydney.edu.au

## Setup/installation

- Install Matlab (tested with version >= 2018b).
- Install [elastix](https://elastix.lumc.nl/index.php) and add the path to the elastix.exe command to the system environment path.
- Install the python module [platipy](https://pyplati.github.io/platipy/) using the command:
```
pip install platipy
```
- Clone the [pix2pix](https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix) repository.
- Install [visdom](https://github.com/fossasia/visdom) using the command:
```
pip install visdom
```

## Acquiring the data
The data for this project was downloaded from the cancer imaging archive [(TCIA)](https://www.cancerimagingarchive.net/). The database used was the HNSCC database which requires permission to download the data from. The list of patient data that was acquired is found in the [PatientList.txt](https://github.com/ACRF-Image-X-Institute/MarkerlessHNGTVTracking/PatientList.txt) file.

## How to generate training/testing data

- Run 'GenerateMultipleVolumes.m' with the variable 'BaseFile' as the directory where the HNSCC data is saved. 

- Run 'GenerateDRRs.m', with the input file as the directory for each patient. 

## How to train and test the network

In a command prompt run visdom using the command 
```
python -m visdom
```

In a seperate command prompt navigate to the pix2pix repository and train the network using:

```
python train.py --dataroot <Location of Training data just created> --load_size 550 --crop_size 512 --input_nc 1 --output_nc 1 --netG unet_256 --dataset_mode aligned --batch_size 8 --no_flip --save_epoch_freq 1 --n_epochs 10
python test.py --dataroot <Location of Training data just created> --load_size 512 --crop_size 512 --input_nc 1 --output_nc 1 --netG unet_256 --dataset_mode aligned --batch_size 8 --no_flip --save_epoch_freq 1 --n_epochs 10
```
More info on the pix2pix repository can be found [here](https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix)

## How to analyse and generate test results

- Run 'DeepLearningRandEvaluation.m' for each patient with the inputvariable the directory where the images from test.py were saved
- Run 'DeepLearningSummaryMultiVol.m' to combine the results for all patients. 
- Plot the results by running the jupyter notebooks worksheet 'TumoutMotionPlots.ipynb'
