%==========================================================================
% download_baselines.m -- fetch the one baseline file that is too large for
% the repository.
%
% Writes  baselines/irfs_ramirez_monetary.mat   (348 MB)
% Source  release "baselines-v1" of this repository on GitHub.
% Run from this folder. main_monetary.m runs without the file and then leaves
% the AR (2018) comparison out of Figure 1 and out of the timing table.
%==========================================================================
BASELINE_URL  = ['https://github.com/davidkim5707/svar-hetero-toolkit/releases/' ...
                 'download/baselines-v1/irfs_ramirez_monetary.mat'];
BASELINE_FILE = fullfile('baselines', 'irfs_ramirez_monetary.mat');
BASELINE_BYTES = 347867945;

if exist(BASELINE_FILE, 'file') == 2
    fprintf('%s is already present.\n', BASELINE_FILE);
else
    fprintf('Downloading %s (348 MB) ...\n', BASELINE_FILE);
    websave(BASELINE_FILE, BASELINE_URL, weboptions('Timeout', 3600));
end

file_info = dir(BASELINE_FILE);
if file_info.bytes == BASELINE_BYTES
    fprintf('CHECK file size %d bytes: PASS\n', file_info.bytes);
else
    error('download_baselines:size', 'CHECK file size %d bytes, expected %d: FAIL', ...
          file_info.bytes, BASELINE_BYTES);
end
