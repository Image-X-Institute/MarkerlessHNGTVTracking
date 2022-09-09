# MarkerlessHNGTVTracking
Code used to create the results described in the paper "A conditional generative adversarial network approach for segmenting head and neck tumors in kV images acquired during radiation therapy"

For a more general version of the CT Head Deformation code see the the repository click [here](https://github.com/ACRF-Image-X-Institute/CTHeadDeformation) 

For any questions email mark.gardner@sydney.edu.au

## Setup/installation

-Install Matlab (tested with version >= 2018b).
-Install [elastix](https://elastix.lumc.nl/index.php) and add the path to the elastix.exe command to the system environment path.
-Install the python module [platipy](https://pyplati.github.io/platipy/) usign the command:
```
pip install platipy
```
-Install the [pix2pix](https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix) module.

## Acquiring the data
The data for this project was downloaded from the cancer imaging archive [(TCIA)](https://www.cancerimagingarchive.net/). The database used was the HNSCC database which requires permission to download the data from. The list of patient data that was acquired is found in the [PatientList.txt](https://github.com/ACRF-Image-X-Institute/MarkerlessHNGTVTracking/PatientList.txt) file.

## How to generate training/testing data

## How to train and test the network

## How to analyse and generate test results
