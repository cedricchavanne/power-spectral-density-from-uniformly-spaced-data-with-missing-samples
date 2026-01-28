# power-spectral-density-from-uniformly-spaced-data-with-missing-samples
Matlab files for producing figures of paper :  
Chavanne, C. Asymptotically-unbiased nonparametric estimation of the power spectral density from uniformly-spaced data with missing samples.  

Matlab files are organized in three folders :  

makedata : Matlab scripts to generate synthetic data and sampling, and to import turbulence data from Kang et al. (2003).  
makeplots : Matlab scripts to analyze data and generate figures of paper.  
tools : Matlab functions to compute cross-correlation and cross-spectral density of two time series.  

Content of each folder is described below in more details.  

makedata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
do_synthetic_data_fbm.m : script to run synthetic_data_fbm.m with specified parameters. Generate data files fbm_1over2.mat and fbm_1over3.mat.  
do_synthetic_holes.m : script to run synthetic_holes.m with specified parameters. Generate data files fbm_holes_33percent.mat and Kang2003_holes_50percent.mat.  
import_Kang2003.m : import data from Kang et al. (2003) downloaded from https://engineering.jhu.edu/meneveau/hot-wire-data.html. Generate data file Kang2003.mat.  
synthetic_data_fbm.m : function to generate random synthetic data for fractional Brownian motion.  
synthetic_holes.m : function to generate random synthetic holes with Bernoulli and batch-Bernoulli distributions.  

To generate synthetic fractional Brownian motion data and synthetic holes with Bernoulli and batch-Bernoulli distributions, follow these steps :  
1) create a directory named 'data' at the same level as, but separate from, the directory where you put the above-mentioned Matlab folders.  
2) Run Matlab, go to the makedata folder, and execute 'do_synthetic_data_fbm'.  
3) execute 'do_synthetic_holes'. Remark : this script also generates synthetic holes for the turbulence data from Kang et al. (2003).  

To import data from Kang et al. (2003), follow these steps :  
1) if not already done, create a directory named 'data' at the same level as, but separate from, the directory where you put the above-mentioned Matlab folders.  
2) in the directory 'data', create a directory named 'Kang2003'.   
3) In a web browser, go to https://engineering.jhu.edu/meneveau/hot-wire-data.html.  
4) Click on hyperlink named 'M20H1', and download zip archive M20H1.zip.  
5) Extract zip archive content in directory 'Kang2003', keeping zip archive directory structure. This creates a directory named 'M20H1', which contains the data files.  
6) Run Matlab, go to the makedata folder, and execute 'import_Kang2003'.  
7) if not already done, execute 'do_synthetic_holes' to generate synthetic holes for the turbulence data from Kang et al. (2003).  

makeplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
plot_corr.m : plot Figure 2.  
plot_fbm.m : plot Figures 1, 4, 5, 6, 7.  
plot_Kang2003.m : plot Figure 8.  
test_fbm.m : compute power spectra for fractional Brownian motion data with missing samples. Generate data files fbm_1over2_missing_33percent_spectra.mat and fbm_1over3_missing_33percent_spectra.mat.  
test_Kang2003.m : compute power spectra for data from Kang et al. (2003) with 50% synthetic missing data. Generate data file Kang2003_spectra.mat.  
test_WK_fbm.m : plot Figure 3.  

To generate figures of paper, follow these steps :  
1) if not already done, create a directory named 'figs' at the same level as, but separate from, the directory where you put the above-mentioned Matlab folders.  
2) Run Matlab, go to the makeplots folder, and execute 'test_fbm' to generate data files fbm_1over2_missing_33percent_spectra.mat and fbm_1over3_missing_33percent_spectra.mat.  
3) execute 'test_Kang2003' to generate data file Kang2003_spectra.mat. This may take some time due to the plomb Matlab function.  
4) execute 'plot_fbm' to plot Figures 1, 4, 5, 6, 7.  
5) execute 'plot_corr' to plot Figure 2.  
6) execute 'test_WK_fbm' to plot Figure 3.  
7) execute 'plot_Kang2003' to plot Figure 8.  

tools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
binaverage.m : bin-average power spectral density using a specified number of frequency bins per decade.  
psd.m : compute the power cross-spectral density of two random processes with no missing samples (periodogram).  
psd_WKbias.m : compute the power cross-spectral density of two random processes with missing samples using Wiener-Khinchin's theorem with the biased estimator of the standard cross-correlation.  
psd_WKcirc.m : compute the power cross-spectral density of two random processes with missing samples using Wiener-Khinchin's theorem with the unbiased estimator of the circular cross-correlation.  
psd_WKunbi.m : compute the power cross-spectral density of two random processes with missing samples using Wiener-Khinchin's theorem with the unbiased estimator of the standard cross-correlation.  
xcorrbias.m : compute the biased standard cross-correlation of two random processes.  
xcorrcirc.m : compute the unbiased circular cross-correlation of two random processes.  
xcorrunbi.m : compute the unbiased standard cross-correlation of two random processes.  

Author:  
Cedric Chavanne  
Professor  
Institut des sciences de la mer  
Universite du Quebec a Rimouski  
Rimouski, QC, Canada  
cedric.chavanne@ensta.org
