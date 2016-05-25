%% Blocking Matrix LMS Filter
%
% Apply an LMS filter to the blocking matrix of a GSC beamformer.  The
% input to each filter is the output of a Delay and Sum beamformer,
% output is subtracted from the delayed original track, and the error
% signal is the rerouted output signal, causing a decrease in target
% signal power.  
%
% *Syntax*
%
% |[z, wSave, mcAdapt, snrSave] = bmlms(x, b, mu, order, beta, phi, 
%      psi, bInit, wInit, wForce, snrThresh, nSnr, snrInit)|
%
% *Required Inputs*
%
% * |x| - Matrix of aligned microphone signals, each column a track
% * |b| - Column vector of Delay and Sum beamformer output
% * |mu| - Step size of adaptation.  This algorithm uses the
%     normalized LMS variant so mu must be in [0,2] to ensure
%     stability.
% * |order| - Order of the adaptive filter.
% * |beta| - Leakage coefficient (forgetting factor).  Beta must be in
%     [0,1], where 0 implies zero memory and 1 implies no forgetting
%     (neither extreme is likely to be very helpful).
%
% *Optional Inputs*  
% 
% (Supply [] if not needed)
%
% * |phi| and |psi| - Vectors of upper and lower bounds,
%     respectively, for a set of Coefficient-Constained Adaptive
%     Filters (CCAF's).
% * |bInit| - Vector of values to be initially fed into the adaptive
%     filter inputs ahead of b.  For a windowed application this can
%     allow for seamless transition between windows without data at
%     window boundaries losing time in the filter.  If empty, zeros
%     are fed as the initial values.
% * |wInit| - Matrix of tap weights for the first iteration.  This is
%     good for ensuring that runs over audio windows flow smoothly, To
%     avoid forcing this initial value, supply []
% * |wForce| - Matrix of tap weights to force the filters into using,
%     rather than adapting.  This is useful for recreating a previous
%     run with the same filter behavior.  To avoid forcing weights,
%     supply []
% * |snrThresh| - Threshold of adaptation.  When not empty, this stage
%     of the beamformer computes an estimate of the SNR at each
%     iteration as the ratio of the FBF output power to the BM output
%     power on each track.  If this ratio is greater than |snrThresh|
%     the BM taps are updated; otherwise the MC taps are updated (see
%     outputs)
% * |nSnr| - Estimate SNR every nSnr samples
%
% *Outputs*
%
% * |z| - Matrix of Blocking Matrix output, ready to be passed to the
%     Mutiple-Input Canceller (MC).
% * |wSave| - Matrix of all tap weights generated by the adaptive
%     filters.  These are useful for debugging and for recreating the
%     same filtering effects with different inputs, and the last set
%     can be used as input to a subsequent period if time windowing.
% * |mcAdapt| - Binary matrix instructing the corresponding filters in
%     the MC to adapt.  One column is written each nSnr samples.
% * |snrSave| - Matrix of estimated SNR's in dB at every iteration
%
% *References*
%
% * Hoshuyama, Osamu, and Akihiko Sugiyama. "Robust Adaptive
% Beamforming."  Microphone Arrays : Signal Processing Techniques and
% Applications. Ed. Michael Brandstein and Darren Ward. New York:
% Springer, 2001.
%
% *Todo*
%
% * Play with limiting (use .1 threshold on norm(bwin))
% * Plot the SNR to see what a good snrThresh would be
%
% Written by Phil Townsend (jptown0@engr.uky.edu) June 10, 2008

%% Function Declaration
function [z, wSave, mcAdapt, snrSave] = bmlms(x, b, mu, order, ...
    beta, phi, psi, bInit, wInit, wForce, snrThresh, nSnr, snrInit)

%% Argument Error Checking
error(nargchk(13, 13, nargin));
if ~isreal(x) || length(size(x)) ~= 2 || ~all(all(isfinite(x)))
    error('x must be a real matrix');
elseif ~isvector(b) || ~isreal(b) || ~all(isfinite(b))
    error('b must be a real vector');
elseif ~isvector(bInit) || ~isreal(bInit) || ...
	~all(all(isfinite(bInit)))
    error('bInit must be a real vector');
elseif isempty(mu) || ~isscalar(mu) || ~isreal(mu) || ~isfinite(mu)
    error('mu must be a real scalar');
elseif isempty(order) || ~isscalar(order) || ~isreal(order) || ...
        ~isfinite(order) || order < 0 || ...
	abs(mod(order,floor(order))) > 0
    error('order must be a positive integer');
elseif isempty(beta) || ~isscalar(beta) || ~isreal(beta) || ...
        ~isfinite(beta) || beta < 0 || beta > 1
    error('beta must be a real scalar in [0,1]');
elseif ~isempty(wInit) && (length(size(wInit)) ~= 2 || ...
        ~isreal(wInit) || ~all(all(isfinite(wInit))))
	error('wInit must be a real matrix');
elseif ~isempty(wForce) && ...
	(~isreal(wForce) || length(size(wForce)) ~= 3 || ...
	 size(wForce,1) ~= order || size(wForce,2) ~= size(x,2) || ...
	 ~all(all(all(isfinite(wForce)))))
    error('wForce must be a real cubic matrix');
elseif ~isempty(snrThresh) && (~isscalar(snrThresh) || ...
        ~isreal(snrThresh) || ~isfinite(snrThresh))
    error('snrThresh must be a real scalar')
elseif ~isempty(nSnr) && ...
	(~isscalar(nSnr) || ~isreal(nSnr) || ...
	 ~isfinite(nSnr) || nSnr < 0 || ...
	 abs(mod(nSnr,floor(nSnr))) > 0)
    error('nSnr must be a positive integer')
elseif xor(isempty(snrThresh), isempty(nSnr))
    error(['snrThresh and nSnr must both be specified or' ...
        ' both empty'])
elseif ~isempty(snrInit) && (~isvector(snrInit) || ~isreal(snrInit)) 
    error('snrInit must be a real vector');
elseif ~isempty(phi) && (length(size(phi))~=2 || ~isreal(phi) || ...
        ~all(all(isfinite(phi))))
    error('phi must be a real matrix')
elseif ~isempty(psi) && (length(size(psi))~=2 || ~isreal(psi) || ...
        ~all(all(isfinite(psi))))
    error('psi must be a real matrix')   
end

%% Adaptive Filtering
% Apply the LMS algorithm.  If phi and psi weren't specified, we get
% the normalized LMS adaptation over each column
%
% $${\bf w}_i(n+1) = \beta {\bf w}_i(n) +
% \mu\frac{z_i(n)}{{\bf b}^T(n){\bf b}(n)}{\bf b}(n),
% \quad \beta \in [0, 1] , \quad \mu \in [0, 2]$$
%
% where the algorithm is "leaky" if beta < 1.
%
% If phi and psi are specified, we have a coefficient-constrainted
% adaptive filter (CCAF), where vectors phi and psi are used as upper
% and lower bounds on the taps, respectively.  This method gives more
% direct control over the power limit of the filter:
%
% $$w_k(n+1) = \phi_k, \quad w_k(n+1) > \phi_k$$
%
% $$w_k(n+1) = \psi_k, \quad w_k(n+1) < \psi_k$$
%
% Note that computations are carried out using matrices for better
% performance.  See the m-file ccafbounds.m for a way to calculate the
% bounds phi and psi.
%

% Initialize variables.  Each column of w holds the taps for an
% adaptive filter on a column of the input microphone tracks x.  wSave
% saves a 3-D matrix of all values of w.  Use wInit for initial value
% of w if supplied; zeros otherwise.
[N, M] = size(x);  % length of tracks, number of mics
if ~isempty(wInit), w = wInit; else w = zeros(order, size(x,2)); end
wSave = zeros(order, size(x,2), length(b));
z = zeros(length(b), size(x,2));
if nargin == 3,  bInit = zeros(order, 1);  end

% AMC Thresholding initial conditions
snrSave = zeros(N, M);
mcAdapt = zeros(length(b), size(x,2));
if ~isempty(snrThresh)
    % Use snrInit if supplied; otherwise set to infinity so that we
    % always adapt over the first power window.
    if isempty(snrInit)
        snr = Inf*ones(1, size(z,2)); %#ok
    else
        snr = snrInit; %#ok
    end
    % Matrix of columns to adapt in MC
    mcAdapt = zeros(size(z)); 
    % Initial value of amcInds
    amcInds = snr > snrThresh;
end

% This highly-vectorized code implements many LMS filters of order
% ord--one for each track of input x
for n = 1:length(b)
    % Select the data window for this iteration
    if n <= order  % feed in supplied past values ahead of b
        bWin = [bInit(end-order+n+1:end); b(1:n)];
    else  % need only data from given b
        bWin = b(n-order+1:n,:);
    end

    % Calculate output for this iteration and ensure stability
    z(n,:) = x(n,:) - (w'*bWin)';  % Subtract filtered b
    if any(~isfinite(z(n,:))), error('Output has blown up'), end

    % If wForce supplied, use those supplied weights.  Otherwise,
    % figure out taps for the next iteration.
    if ~isempty(wForce), w = wForce(:,:,n); continue, end

    % If enough samples have passed and we're doing AMC
    % thresholding, find the SNR for each track.
    if ~isempty(snrThresh) && mod(n, nSnr) == 0
        snr = 10*log10(sum(b(n-nSnr+1:n).^2) ./ ...
            sum(z(n-nSnr+1:n, :).^2));
        amcInds = snr > snrThresh;  % find which tracks adapt
    end
    snrSave(n,:) = snr; % save SNR values for debug output
    
    % Perform LMS tap update for all columns.
    % 1. If no thresholding, always adapt
    % 2. If thresholding, update if snr is above threshold,
    %    otherwise continue
    if isempty(snrThresh) % adapt everyone
        w = beta*w + mu*bWin*ones(1,size(z,2))*diag(z(n,:)) / ...
            norm(bWin)^2;
    else  % threshold certain tracks
        mcAdapt(n, :) = ~amcInds;
        w(:, amcInds) = beta*w(:, amcInds) + mu*bWin*...
            ones(1,length(find(amcInds))) * diag(z(n,amcInds)) / ...
            norm(bWin)^2;
    end

    % CCAF (Coefficient Constrained) Tap Constraints
    if ~isempty(phi)
        % Expand bound vectors and run weight comparison
        uBoundInds = w > phi;  lBoundInds = w < psi;
        % Include AMC thresholding if specified
        if ~isempty(snrThresh)
            amcIndMat = ones(order,1) * amcInds;
            uBoundInds = uBoundInds & amcIndMat;
            lBoundInds = lBoundInds & amcIndMat;
        end
        % Apply bounds to tap matrix
        w(uBoundInds) = phi(uBoundInds);
        w(lBoundInds) = psi(lBoundInds);
    end

    % Save the taps used on this iteration
    wSave(:,:,n) = w;

end  % for n = 1:length(b)

end  % function bmlms