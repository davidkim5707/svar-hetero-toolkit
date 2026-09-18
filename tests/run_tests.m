%==========================================================================
% run_tests.m -- run the self-tests of the zero-restriction routines.
%
% Each test builds a synthetic cache, runs one routine of code/core/zero, and
% prints PASS or FAIL per check. A failing check raises an error.
% Runtime is under a minute.
%==========================================================================
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup.m'));

test_build_zero_constraint();
test_draw_Q_zero_columnwise();
test_ess_draw_Q_zero_columnwise();
test_log_volume_element_zero();
fprintf('All zero-restriction self-tests finished.\n');
