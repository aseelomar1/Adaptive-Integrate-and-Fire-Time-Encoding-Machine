# Adaptive Integrate-and-Fire Time Encoding Machine (AIF-TEM)

## Description
This repository contains MATLAB code for the Adaptive Integrate-and-Fire Time Encoding Machine (AIF-TEM), as presented in the EUSIPCO24 conference and subsequent journal publication. The code aims to improve the performance of Integrate-and-Fire time encoding machine through adaptive techniques.

## Table of Contents
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Function Overview](#function-overview)
- [License](#license)
- [Attribution](#attribution)
- [Acknowledgments](#acknowledgments)
- [Publications](#publications)

### Usage
In your MATLAB command window:
```matlab
addpath('mainFunc');  % Add mainFunc to path
samplingMain;         % Call the main function samplingMain or quantizationMain
```
### Directory Structure
```bash
project_root/
│
├── mainFunc/
│   ├── samplingMain.m          # Main function for sampling and recovery of a signal
│   └── quantizationMain.m      # Main function for sampling, quantization and recovery of a signal
│
├── HelperFunc/
│   ├── bias_quant.m            # Helper function for bias quantization
│   ├── welford_update.m        # Helper function for Welford's algorithm
|   ├── calcEnergyMaxCoeff.m    # Calculate maximom energy and amplitude
|   ├── compute_G_matrix.m      # Compute G matrix that used to recover the signal from TEM output
|   ├── genSinc.m               # Generate a sum of sinc signal
|   └── quantizationTEM.m       # Perform Quantization
│
└── TEM/
    ├── AIF_TEM.m                # Adaptive Integrate-and-Fire TEM
    ├── IF_TEM.m                 # Classical Integrate-and-Fire TEM
    └── recover_TEM              # Recover the signal from A/IF-TEM output
```
### Function Overview
- **`samplingMain.m`**: This function generates a signal that is the sum of N sinc functions multiplied by random amplitudes. The signal is then sampled using both the Integrate-and-Fire Time Encoding Machine (IF-TEM) and the Adaptive Integrate-and-Fire Time Encoding Machine (AIF-TEM). The output from both samplers is used to recover the signal from its time samples. The recovery distortion and oversampling are calculated.
- **`quantizationMain.m`**: This function generates a signal that is the sum of N sinc functions multiplied by random amplitudes. The signal is sampled using both IF-TEM and AIF-TEM, and the output from both samplers is quantized. The signal is then recovered from the quantized time samples of both samplers. The recovery distortion and the total number of bits used to encode the signal are calculated.

### License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

### Attribution
If you use this code in your work, please cite the original author:

**Author**: Aseel Omar  
**Emails**: [aseel.omar@campus.technion.ac.il](mailto:aseel.omar@campus.technion.ac.il), [aseel.to16@gmail.com](mailto:aseel.to16@gmail.com)

Your contributions and usage of this code are appreciated!

### Acknowledgments
- Special thanks to my advisor, Alejandro Cohen, for his support.

### Publications
- Aseel Omar, Alejandro Cohen, "Adaptive Integrate-and-Fire Time Encoding Machine," presented at the EUSIPCO conference, 2024. [Conference Version](https://eurasip.org/Proceedings/Eusipco/Eusipco2024/pdfs/0002442.pdf)
- Aseel Omar, Alejandro Cohen, "Adaptive Integrate-and-Fire Time Encoding Machine with Quantization," arXiv preprint, 2024. [arXiv Version](https://arxiv.org/abs/2403.02992)

