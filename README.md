# MSLNet
Code for our paper MSLNet and Perceptual Grouping for Guidewire Segmentation and Localization:
https://www.mdpi.com/1424-8220/25/20/6426/pdf?version=1760878225

Files:
- train_MSLNetNN.ipynb trains the MSLNet. The data shoudl be organized according to the NNUNet dataset folder structure (https://github.com/MIC-DKFZ/nnUNet).
- test_MSLNetNN.ipynb uses a trained model to predict segmentations from data.
- evaluation.ipynb evaluates one or more directories containing segmentations.
- localization.ipynb performs perceptual grouping and evaluates the result. Uses the package pycvs to extract curves
