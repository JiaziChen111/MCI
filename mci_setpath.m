
% Are you using the standalone version which does not require SPM
% on your search path ?
standalone=0;

% Set spm_dir if you're not using the standalone version
% SPM can be downloaded from http://www.fil.ion.ucl.ac.uk/spm/software/
spm_dir='C:\Users\bzv17fbu\ucl\spm12';

% Are you goint to be using Sundials for integration of differential equations
add_sundials=1; 

% Set sun_dir if you're going to be using Sundials 
% The four major components of the Sundials package (CVODE,
% CVODES,IDA,IDAS) can be downloaded from http://computation.llnl.gov/casc/sundials/
sun_dir='C:\Users\bzv17fbu\ucl\fil\sampling\';

% Add subdirectories to path
addpath(genpath(pwd));

if standalone
    % Add subset of required SPM12 code provided by MCI distribution
    % (this may not provide full MCI functionality)
    w=pwd;
    addpath([w(1:length(w)-3),'spm-code']);
else
    % Put SPM12 on path
    disp('Adding SPM12 to path ');
    disp(' ');
    add_dir{1}=[spm_dir];
    add_dir{2}=[spm_dir,'\toolbox\dcm_meeg'];
    add_dir{3}=[spm_dir,'\toolbox\spectral'];
    for i=1:length(add_dir)
        addstr=['addpath ',add_dir{i}];
        disp(addstr);
        eval(addstr);
    end
end

if add_sundials
    disp('Adding Sundials to path ');
    disp(' ');
    clear add_dir
    add_dir{1}=[sun_dir,'sundials-2.5.0\sundialsTB\cvodes'];
    add_dir{2}=[sun_dir,'sundials-2.5.0\sundialsTB\cvodes\cvm'];
    for i=1:length(add_dir)
        addstr=['addpath ',add_dir{i}];
        disp(addstr);
        eval(addstr);
    end
end



