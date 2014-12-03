matlab_SPoC
===========

#### This package contains Matlab code for 
* Source Power Correlation analysis (SPoC, Dähne et al. 2014a)
* multimodal Source Power Correlation analysis (mSPoC, Dähne et al., 2013)
* canonical Source Power Correlation analysis (cSPoC, Dähne et al., 2014b)
* Spatio-Spectral Decomposition (SSD) for dimensionality reduction (Nikulin et al., 2011, Haufe et al. 2014b)


#### Important Notes:

0. Code for mSPoC to be added soon!
1. Please make sure the util folder (and all of its subfolders) are on the Matlab path. Otherwise the optimization required for (m/c)SPoC will not work! Run the `startup_spoc.m` script to add folders to the path. 
2. Please read the documentation of the matlab functions `ssd.m`, `spoc.m`, `mspoc.m`, `cspoc.m` and run / look at the respective examples. I have tried to explain everything that you need to know to use the functions. If there is unclarity, please let me know and I will try to improve the documentation. 
3. It is highly recommened to use dimensionality reduction via SSD before applying (m/c)SPoC. Dimensionality reduction greatly increases the computational speed and improves the quality of the results!
Below you find a snippet of matlab code that shows an example of how to use SSD for preprocessing.
4. EEGLAB plugins are on the way!


#### References

S. Dähne, F. Biessman, F. C. Meinecke, J. Mehnert, S. Fazli, K. R. Müller, "Integration of Multivariate Data Streams With Bandpower Signals", *IEEE Transactions on Multimedia*, 15(5):1001-1013, 2013

S. Dähne, F. C. Meinecke, S. Haufe, J. Höhne, M. Tangermann, K. R. Müller, V. V. Nikulin, "SPoC: a novel framework for relating the amplitude of neuronal oscillations to behaviorally relevant parameters", *NeuroImage*, 86:111-122, 2014

S. Dähne, V. V. Nikulin, D. Ramirez, P. J. Schreier, K. R. Müller, S. Haufe, "Finding brain oscillations with power dependencies in neuroimaging data", *NeuroImage*, 96:334-348, 2014 


S. Haufe, F. Meinecke, K. Görgen, S. Dähne, J. Haynes, B. Blankertz, F. Biessmann, "On the interpretation of weight vectors of linear models in multivariate neuroimaging", *NeuroImage*, 87:96-110, 2014

S. Haufe, S. Dähne, V. V. Nikulin, "Dimensionality reduction for the analysis of brain oscillations", *NeuroImage*, 101:583-597, 2014 



#### SSD preprocessing:

```matlab

% We assume X contains the continous EEG/MEG data with 
% columns of X corresponding to channels and rows to time points, 
% i.e. size(X) = [n_samples, n_channels]. 
% If you want to apply cSPoC, X is a cell array with X{n} being
% the n'th dataset. The preprocessing has to be applied to each 
% dataset individually


% First we call the SSD function which returns the 
% original data bandpassed and projected onto the SSD components
% Also, in order to compute the spatial patterns of extracted 
% (m/c)SPoC sources,  we need the covariance matrix of the 
% bandpassed sensorspace data.
% See SSD code for more info on the paramters

[W_ssd, ~, ~, C, X_ssd]  = ssd(X, bands, ...); 


% Dimensionality reduction by choosing only the first few 
% SSD components

n_ssd_components = 15; 
X_ssd = X_ssd(:,1:n_ssd_components);
W_ssd = W_ssd(:,1:n_ssd_components);


% Apply any of the (m/c)SPoC methods to get the SPoC spatial filters.

W_spoc = spoc(X_ssd, z);

% Note that the spatial filters are now defined in SSD space and have 
% to be mapped to sensor space to be applicable to the original data.

W = W_ssd * W_spoc;

% The spatial patterns are computed using the covariance matrix of the
% bandpass filtered sensor space data

A = W * C;

```
