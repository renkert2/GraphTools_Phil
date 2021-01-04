%DAEMOT_AUTO_INSTALL installs the DAEMOT toolset
%   DAEMOT_AUTO_INSTALL should be executed after extracting the
%   DAEMOT toolset to the desired installation directory.  When called,
%   the script will install the necessary directories for operation of the
%   DAEMOT toolset.
%
%   See also addpath genpath

%   Developed by Matthew Williams - UIUC
%   Modified by Herschel Pangborn - UIUC in 7/2019 to remove unnecessary adding of subfolders

spath  = fileparts(mfilename('fullpath'));    % path of executed script

% installing icons directory
dir_path{1} = fullfile(spath,'Icons');
if exist(dir_path{1},'dir')
    addpath(genpath(dir_path{1}))
    err(1) = true;
else
    err(1) = false;
end

% installing support functions directory
dir_path{2} = fullfile(spath,'Support Functions');
if exist(dir_path{2},'dir')
    addpath(genpath(dir_path{2}))
    err(2) = true;
else
    err(2) = false;
end

dir_path{3} = spath;
if exist(dir_path{3},'dir')
    addpath(dir_path{3})
    err(3) = true;
else
    err(3) = false;
end

% installing electrical maps directory
dir_path{4} = fullfile(spath,'Graph Generation');
if exist(dir_path{4},'dir')
    addpath(genpath(dir_path{4}))
    err(4) = true;
else
    err(4) = false;
end

% Error checking
if ( err(1) == false && err(2) == false && err(3) ==  false && err(4) == false)
    warning(['DAEMOT_auto_install.m failed to install the necessary '...
        'directories. Administrator privileges may be required. Check', ...
        ' necessary privileges of the installation directory.'])
elseif ( err(1) == false || err(2) == false || err(3) ==  false || err(4) == false)
    warning(['One or more directories failed to install.  Manual ',...
        'installation required for the following:'])
    for i=1:4
        if err(i) == false
            disp(dir_path{i})
        end
    end
else
    disp(['Successfully installed the DAEMOT toolset. '...
        'MATLAB/Simulink may need to be restarted.'])
end

savepath      % saves the path     

clear spath dir_path err
