# Period Amplitude Analysis (PAA) EEGlab Plugin

Method for detecting slow waves using GUI in EEGlab for MATLAB

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Mac, PC and Linux compatible.  
Designed for use with EEGlab 2019_1 on Matlab R2019a.  
For use on continuous eeglab datasets (*.set).  
Requires sleep scoring and movement artifacts to be included in EEG.event structure.

### Installing

Simply unzip and copy 'PAA' to '~/eeglab/plugins/'

## Usage

using the eeglab interface:

* at the matlab command prompt, add eeglab root directory to the path and launch eeglab
* load an eeglab dataset
* navigate to Tools>Period Amplitude Analysis>

using the batch script:

* use the included PAA_standalone.m script to run multiple files as a batch. 
* manually specify user-defined parameters prior to running script

## Authors

Sleep Well

## Contact 

https://www.sleepwellpsg.com

## License

Copyright (C) Sleep Well, 2020.  
See the GNU General Public License for more details.

## Relevant Citations

Feinberg et al, 1978. https://doi.org/10.1016/0013-4694(78)90266-3
Geering et al, 1993. https://doi.org/10.1111/j.1365-2869.1993.tb00074.x
Bersagliere and Achermann, 2010. https://doi.org/10.1111/j.1365-2869.2009.00775.x
