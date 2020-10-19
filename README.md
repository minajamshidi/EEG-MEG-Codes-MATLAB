# EEG-MEG-Codes-MATLAB
This repository includes useful MATLAB functions for EEG analysis.

* ssd: is a frequency-domain extension of ssd (spatio-spectral decomposition). Please see the function heading.

* select_component: is a function which helps manual ICA artifact rejection. The EEGLab ICA ispection window is slow and it does not give the opportunity of interaction with data. With this function, one can interactively see the spectrum of the EEG, if the selected components are rejected, as well as the amount of alpha power that will be rejected, and other options (see the heading of the function). Additionally, with left and right click can select "active" and "passive" components, meaning that active components should be rejected anyway and passive components are kept for possible future usage. For example, in many cases researchers prefer to remove eye and heart artifact, but also save the muscle artifact components, in case they want to remove them. The former will be selected by left click as active and the latter with right click as passive.

* plot_spec: is a function for plotting spectrum of a multichannel signal. Note: The EEGlab function has some bugs. 
