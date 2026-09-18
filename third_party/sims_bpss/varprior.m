function [ydum,xdum,breaks]=varprior(nv, nx, lags, mnprior, vprior, urprior, ybar, xbar, nstat, mnstart)
%function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
% ydum, xdum:   dummy observation data that implement the prior
% breaks:       vector of points in the dummy data after which new dummy obs's start
%                   Set breaks=T+[0;breaks], ydata=[ydata;ydum], xdum=[xdata;xdum], where 
%                   actual data matrix has T rows, in preparing input for rfvar3
% nv,nx,lags: VAR dimensions
% mnprior.tight:Overall tightness of Minnesota prior
% mnprior.decay:Standard deviations of lags shrink as lag^(-decay)
% vprior.sig:   Vector of prior modes for diagonal elements of r.f. covariance matrix
% vprior.w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
%                   Should be an integer, and will be rounded if not.  vprior.sig is needed
%                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
%                   Set vprior.w=0 to achieve this.
% Note:         The original Minnesota prior treats own lags asymmetrically, and therefore
%                   cannot be implemented entirely with dummy observations.  It is also usually
%                   taken to include the sum-of-coefficients and co-persistence components
%                   that are implemented directly in rfvar3.m.  The diagonal prior on v, combined
%                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
%                   prior mean generates larger prior variances for own than for cross-effects even in 
%                   this formulation, but here there is no way to shrink toward a set of unconstrained 
%                   univariate AR's.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m
if nargin < 10
    mnstart = 1;
end
if nargin < 9
    nstat = ones(nv,1);
end
if nargin < 8
    xbar = ones(1,nx);
end
if nargin < 7
    ybar = zeros(1,nv);
end

if isstruct(mnprior) && ...
   isfield(mnprior, 'tight') && ~isempty(mnprior.tight) && ...
   isfield(mnprior, 'decay') && ~isempty(mnprior.decay)
    if nx > 0
        xdum = zeros(lags + 1, nx, lags, nv);
    else
        xdum = [];
    end
    ydum = zeros(lags+1, nv, lags, nv);

    for il = 1:lags
        if il <= mnstart
            decayfactor = 1;
        else
            decayfactor = (il - mnstart + 1)^mnprior.decay;
        end
        ydum(il+1,:,il,:) = decayfactor * diag(vprior.sig);
    end

    % Initial condition for stationarity
    ydum(1,:,1,:) = diag(vprior.sig .* nstat);
    % Apply tightness and reshape
    ydum = mnprior.tight * reshape(ydum, [lags+1, nv, lags*nv]);
    ydum = flip(ydum, 1);

    xdum = mnprior.tight * reshape(xdum, [lags+1, nx, lags*nv]);
    xdum = flip(xdum, 1);

    % Optional: define breaks if needed later
    breaks = (lags+1) * (1:(nv*lags))';
    lbreak = breaks(end);
else
    ydum = [];
    xdum = [];
    breaks = [];
    lbreak = 0;
end
% elle = 0;
% for j1 = 1 : nv
%     for j2  =  1 : lags   
%         elle = 1 + elle;
%         ydum(:,:,elle) = ydum(:,:,elle) * mnprior.unit_root_(j1);
%     end
% end

% Unit roots and co-persistence prior (urprior)
if ~isempty(urprior)
    % Lambda part (co-persistence prior)
    if isfield(urprior, 'lambda') && ~isempty(urprior.lambda)
        lambda = urprior.lambda;
        ydum_lambda = repmat(ybar, lags+1, 1) * abs(lambda);
        if lambda > 0
            xdum_lambda = repmat(xbar, lags+1, 1) * lambda;
        else
            xdum_lambda = zeros(lags+1, nx);
        end
        ydum = cat(3, ydum, ydum_lambda);
        xdum = cat(3, xdum, xdum_lambda);
        breaks = [breaks; lbreak + (lags+1)];
        lbreak = breaks(end);
    end
    
    % Mu part (own-persistence prior)
    if isfield(urprior, 'mu') && ~isempty(urprior.mu)
        mu = urprior.mu;
        ydum_mu = zeros(lags+1, nv, nv);
        for iv = 1:nv
            if nstat(iv)
                ydum_mu(:, iv, iv) = ybar(iv);
            end
        end
        ydum_mu = mu * ydum_mu;
        xdum_mu = zeros(lags+1, nx, nv);
        ydum = cat(3, ydum, ydum_mu);
        xdum = cat(3, xdum, xdum_mu);
        breaks = [breaks; lbreak + (lags+1)*(1:nv)'];
        lbreak = breaks(end);
    end
    
end

if ~isempty(vprior) && ~isempty(vprior.w) && vprior.w > 0
    ydum2 = zeros(lags+1,nv,nv);
    xdum2 = zeros(lags+1,nx,nv);
    ydum2(end,:,:) = diag(vprior.sig) * vprior.w;  
    ydum = cat(3,ydum,ydum2);
    xdum = cat(3,xdum,xdum2);
    breaks = [breaks;(lags+1)*(1:nv)'+lbreak];
    lbreak = breaks(end);
end

dimy = size(ydum);
ydum = reshape(permute(ydum,[1 3 2]),dimy(1)*dimy(3),nv);
xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx);
breaks = breaks(1:(end-1));



