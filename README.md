<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->


<!-- ABOUT THE PROJECT -->
## About The Project

Widespread distributed fracturing during earthquakes threatens infrastructure and lifelines. We combine high-resolution rupture maps from the four major surface-rupturing strike-slip earthquakes in southern California and northern Mexico since 1992 to incorporate the displacements produced by secondary faulting and fracturing into a probabilistic displacement hazard analysis framework. Through analysis of the length and spatial distribution of mapped fractures for each event, we develop an expression for the probability of finding a fracture that accommodates a displacement that exceeds a given threshold Do at a given distance away from the fault. In our model, Do is a user-input variable that may be adjusted for different engineering applications. To facilitate general use, we accompany this study with a Python Jupyter notebook that produces a hazard curve with robust uncertainties for end-user application. Our model is best applied to estimating expected secondary displacements for strike-slip earthquakes, similar to those analyzed, with widespread fracturing across immature fault zones.

<!-- GETTING STARTED -->
## Getting Started

This repository contains the scripts and data required to create the general model published in Rodriguez Padilla and Oskin 2022N. This information is stored in the General end-user model - Python directory. 

Additionally, all of the scripts to develop this model based on analysis of data from the Landers, Hector Mine, El Mayor Cucapah, and Ridgecrest earthquakes is available in the All scripts analysis (Matlab + Python) directory. The data for running this analysis are stored in a separate repository that is also open-access (). All figures in the manuscript can be generated using these scripts and data. 

### Prerequisites

The general model described runs on Python Jupyter Notebooks. 

Some of the scripts for running the models for each individual earthquake require Matlab as well as the Matlab Mapping Toolbox. Some of the scripts rely on functions downloable from Mathworks. The specific dependencies for each Matlab script to run are listed at the beginning of the corresponding script. 


<!-- ROADMAP -->
## Creating hazard curves
### Running the general end-user model

- [ ] In the "General end-user model - Python" directory
    - [ ] Run the "PFDHA_secondary_fracturing.ipynb" Jupiter Notebook
    - [ ] Input your desired value for Do (in meters) when the dynamic prompt comes up
    - [ ] The script exports two pdf files: the hazard curve and the uncertainty distribution for each parameter in the  model

### Running the complete analysis in Rodriguez Padilla and Oskin (202N)
The equation numbers correspond to those in the manuscript. 

- [ ] In the "All scripts analysis (Matlab + Python)" directory:
    - [ ] To generate the fracture density decays (second term of eq.10)
    	- [ ] Run the "fracture_density_decay.m" script
    	- [ ] The script exports a text file with the fracture density at each x position with distance away from the fault
    	- [ ] Run the "MCMC_density_decay_Poisson.ipynb" script
    	- [ ] The script outputs the best-fit and posterior distribution of the fit of equation 8 to the fracture density decay as well as the Markov chain
    - [ ] To estimate the scaling ratio of fracture length and displacement for each event in CA from the FDHI database, including the four events in our study: 
    	- [ ] Run the "max_slip_length_scaling_CA.m" script
    - [ ] To estimate P(D>Do) for the EMC, Landers, Hector Mine, and Ridgecrest data: 
    	- [ ] Run the "generate_hazard_curve.m" script
    	- [ ] This script outputs figures 2 and 4 in the manuscript
    - [ ] To generate uncertainties for the hazard curve for each event
	- [ ] Run the "generate_hazard_curve.m" 
	- [ ] This script outputs figures A2-A5 in the appendix

<!-- CONTACT -->
## Contact

Please report suggestions and issues:

Your Name - [@_absrp](https://twitter.com/_absrp) - arodriguezpadilla@ucdavis.edu

Project Link: [https://github.com/absrp/PFDHA_strikeslip](https://github.com/absrp/PFDHA_strikeslip)

<p align="right">(<a href="#readme-top">back to top</a>)</p>




