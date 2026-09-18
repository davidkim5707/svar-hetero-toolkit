%==========================================================================
% setup.m -- add core/ and third_party/ to the MATLAB path.
%
% Run once per session from anywhere:  run('<path to repo>/setup.m')
% The replication drivers under replications/ add the same two folders
% themselves, so they do not need this script.
%==========================================================================
hars_root__ = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(hars_root__, 'core')), ...
        genpath(fullfile(hars_root__, 'third_party')));

if verLessThan('matlab', '9.8')
    warning('hars:setup:oldMatlab', ...
        'MATLAB R2020a or newer is required for exportgraphics; figures will not be written.');
end
fprintf('[setup] added core/ and third_party/ under %s\n', hars_root__);
clear hars_root__
