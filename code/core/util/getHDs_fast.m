function [HD] = getHDs_fast(IRFs, shocks, whichVariable)
% IRFs: [ny x horizon x nshock]
% shocks: [T x nshock]
% whichVariable: scalar (index of variable whose decomposition is desired)

    hmax = size(shocks,1) - 1;
    nshock = size(shocks, 2);
    HD = zeros(nshock, 1);  % One entry for each structural shock

    for j = 1:nshock
        sum_contrib = 0;
        for h = 0:hmax
            shock_t = shocks(end - h, j);            % shock j at time t−h
            irf_t = IRFs(whichVariable, h+1, j);     % response of var to shock j at horizon h
            sum_contrib = sum_contrib + irf_t * shock_t;
        end
        HD(j) = sum_contrib;
    end
end