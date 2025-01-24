<h1>A minimal model for multigroup adaptive SIS epidemics</h1>

![](https://komarev.com/ghpvc/?username=SaraSottile)

<h2>Description</h2>

This repository contains the supplementary material related to the paper "A minimal model for multigroup adaptive SIS epidemics", co-authored by Massimo A. Achterberg, <a href="https://mattiasensi.github.io/" target="_blank">Mattia Sensi</a>  and <a href="https://sarasottile.github.io/" target="_blank">Sara Sottile</a>.

The preprint version of the paper can be found <a href="https://arxiv.org/abs/2407.17639" target="_blank">here</a>. Feedback is welcome! If you encouter a bug or have problems, or want to give advice, please contact us:

<p>M.A.Achterberg@tudelft.nl</p>
<p>mattia.sensi@polito.it</p>
<p>sara.sottile4@unibo.it</p>

<h2>Overview</h2>
This repository contains the MATLAB codes used to generate the figures in the paper. The provided scripts include the necessary parameters and methods to reproduce the figures. In particular, the random seeds provided will allow user to obtain the exact same values we used.

<h2>Structure</h2>
The repository is organized into three folders, each corresponding to a specific set of figures from the paper:

<li>Folder 1: MATLAB codes for Figure 2.</li>
<li>Folder 2: MATLAB codes for Figures 3-4.</li>
<li>Folder 3: MATLAB codes for Figure 5.</li>

Each folder includes:

<li>MATLAB scripts for figure generation.</li>
<li>Network files required as input, provided directly within the respective folder.</li>

<h2>Reproducibility</h2>
The random seed is fixed for reproducibility. At the beginning of each script, you will find the following command:

<br>rng(1); % Setting seed=1 for reproducibility  </br>

<h2>Usage</h2>
Navigate to the folder corresponding to the figure you want to generate.
Run the provided MATLAB scripts to reproduce the figures.
The required network files will be automatically read as inputs.

<h2>Notes</h2>
Ensure that MATLAB is installed and properly configured on your system before running the codes.


