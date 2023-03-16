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

Widespread distributed fracturing during earthquakes threatens infrastructure and lifelines. We combine high-resolution rupture maps from the five major surface-rupturing strike-slip earthquakes in southern California and northern Mexico since 1992 to incorporate the displacements produced by secondary ruptures into a probabilistic displacement hazard analysis framework. Through analysis of the spatial distribution of mapped ruptures and displacements for each of these events, we develop a magnitude-dependent expression for the probability of finding a secondary rupture that accommodates a displacement that exceeds a displacement threshold at a given distance away from the principal fault. Our model is best applied to estimating expected secondary displacements for strike-slip earthquakes, similar to those analyzed, with widespread ruptures across immature fault zones. 

<!-- GETTING STARTED -->
## Getting Started

This repository contains the scripts and data required to create the general model published in Rodriguez Padilla and Oskin 202N. The general model runs in script PFDHA_secondary_ruptures.ipynb (see instructions for use below). 

Additionally, all of the scripts to develop this model based on analysis of data from the Landers, Hector Mine, El Mayor Cucapah, and Ridgecrest (foreshock and mainshock) earthquakes is available in the PFHDA_model_allcode directory. The data for running this analysis are stored in a separate repository that is also open-access (https://figshare.com/projects/A_probabilistic_displacement_hazard_assessment_framework_for_distributed_ruptures_from_strike-slip_earthquakes/162349). All figures in the manuscript can be generated using these scripts and data. 

### Prerequisites

The general model described runs on Python Jupyter Notebooks. 

Some of the scripts for running the models for each individual earthquake require Matlab as well as the Matlab Mapping Toolbox. Some of the scripts rely on functions downloable from Mathworks. The specific dependencies for each Matlab script to run are listed at the beginning of the corresponding script. 


<!-- ROADMAP -->
### Running the general end-user model

- [ ] In the "General end-user model - Python" directory
    - [ ] Run the "PFDHA_secondary_ruptures.ipynb" Jupyter Notebook
    - [ ] Input your desired value for S_0 (in meters) and Mw when the dynamic prompt comes up
    - [ ] The script exports a pdf file with the hazard curve and the uncertainty distribution for each parameter in the  model


<!-- CONTACT -->
## Contact

Please report suggestions and issues:

[@_absrp](https://twitter.com/_absrp) - arodriguezpadilla@ucdavis.edu

Project Link: [https://github.com/absrp/PFDHA_strikeslip](https://github.com/absrp/PFDHA_strikeslip)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

Manuscript Link: * in review, stay tuned *




