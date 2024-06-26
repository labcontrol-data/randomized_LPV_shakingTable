# randomized_LPV_shakingTable

<p align="center"><a href="http://www.labcontrol.xyz/dokuwiki" target="_blank" rel="noopener"><img src="https://github.com/labcontrol-data/images/blob/main/logo.png"></a></p>

# [Labcontrol, Brazil](http://www.labcontrol.xyz/dokuwiki)

[**Labcontrol**](http://www.labcontrol.xyz/dokuwiki)  is a research laboratory located at Universidade Tecnológica Federal do Paraná (UTFPR), Brazil. UTFPR is a Brazilian public university located in the Paraná state, Brazil. [**Labcontrol**](http://www.labcontrol.xyz/dokuwiki)  develops research on Control Systems and Automation. The Scientific Director of Labcontrol is [Prof. Dr. Alessandro N. Vargas](http://www.anvargas.com). The projects hosted in [**Labcontrol**](http://www.labcontrol.xyz/dokuwiki)  are [described in this link (click here).](http://www.anvargas.com/blog)

About
============

This page provides information about the project developed in Labcontrol called ["Randomization to the control design of a shaking table with MR-damper."](http://www.anvargas.com/blog/mrdamper.html)  Experiments were carried out in practice in a laboratory testbed, and the data contained in this Github repository were collected in those experiments. 

[![DOI](https://zenodo.org/badge/330236633.svg)](https://zenodo.org/badge/latestdoi/330236633)

**Please check more details about this project in the page detailing the ["Randomization to the control design of a shaking table with MR-damper."](http://www.anvargas.com/blog/mrdamper.html)**


`maincode.m` is a MATLAB(R) script that calls Simulink-Matlab and generates simulation data. The script also generates figures. The figures contain both simulation and real-time data collected in a laboratory testbed.

For more details about the experimental data, as long as the corresponding academic publications, please visit [the project page](http://www.anvargas.com/blog).


Installation
============

1. The code in this repo solves linear matrix inequalities (LMIs). You will need to install two tools in your Matlab.
    - Yalmip: [yalmip.github.io](https://yalmip.github.io/)
    - Mosek: [mosek.com](https://www.mosek.com/)
2. After you ensure the two tools from "Step 1" above are working properly in your Matlab, then extract the ZIP file (or clone the git repository) in your computer.
3. Add the folders `matlab-code/` and `data/` to your path in MATLAB/Octave: e.g. 
    - using the "Set Path" dialog in MATLAB, or 
    - by running the `addpath` function from your command window or `startup` script.

Make sure that you are running Matlab 2017a (or a newer version). Older versions may work, but it is uncertain.

Usage
=====

Typical usage of `maincode.m` consists of running it in your MATLAB. The code generates ten figures.

MATLAB
------
  1. Run in this order:
     - `Identification_LPV.m` This file processes the experimental data and computes the identification matrices of a time-varying linear system. The code runs Mayne's estimator [1], a particular case of Kalman Filter, to perform the identification.
     - `cleaner_matrices_LPV.m` This file processes a large database of matrices and prunes them. It prepares the files for the next step.
     - `finalVertices.m` Computes the vertices of the polytope (i.e., it generates A_1,...,A_n and B_1,...,B_n).
     - `maincode.m` Main code that ensures the closed-loop system is stable under the randomized approach.

More information
================

* For more information about `maincode.m`, visit the author's page: [Prof. Alessandro N. Vargas](http://www.anvargas.com). You are welcome to help improve the code.
* You are free to use the data in your research. When doing so, please contact the author [Prof. Alessandro N. Vargas](http://www.anvargas.com) 
and let him know about your project. Depending on your research area, the author can help you interpret the data according to your application. The author can also suggest papers and books that can be helpful in your research.

[![DOI](https://zenodo.org/badge/330236633.svg)](https://zenodo.org/badge/latestdoi/330236633)

Citation
------
How to cite the data of this repository:

```
@misc{vargasGithub2024,
    author       = {A. N. Vargas},
    title        = {Data, source code, and documents for the shaking table with {MR}-damper}},
    month        = {Jan},
    year         = 2021,
    doi          = {10.5281/zenodo.4445334},
    version      = {1.0.3},
    publisher    = {Zenodo},
    url          = {https://zenodo.org/badge/latestdoi/330236633}
};
```
[1] MAYNE, D. Q. Optimal non-stationary estimation of the parameters of a linear system
with Gaussian inputs. Journal of Electronics and Control, Taylor & Francis, v. 14, n. 1, p.
101–112, 1963.
