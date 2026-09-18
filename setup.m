%==========================================================================
% setup.m -- add the library to the MATLAB path.
%
% Adds code/core and code/third_party with their subfolders. Run once per
% session from anywhere:  run('<path to repository>/setup.m')
% The examples and the tests call this script themselves. The replication
% drivers add the same two folders without it.
%==========================================================================
repo_root__ = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repo_root__, 'code', 'core')), ...
        genpath(fullfile(repo_root__, 'code', 'third_party')));

if verLessThan('matlab', '9.8')
    warning('setup:oldMatlab', ...
        'MATLAB R2020a or newer is required for exportgraphics; figures will not be written.');
end
if ~license('test', 'Statistics_Toolbox')
    warning('setup:noStatistics', ...
        'The Statistics and Machine Learning Toolbox is required (mvnrnd, quantile, gamrnd).');
end
fprintf('[setup] added code/core and code/third_party under %s\n', repo_root__);
clear repo_root__
