
% Add subdirectories to path
addpath(genpath(pwd));
 
% Add subset of required SPM12 code provided by MCI distribution
% (this may not provide full MCI functionality)
w=pwd;
addpath([w(1:length(w)-3),'spm-code']);



