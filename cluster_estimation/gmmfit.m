%GMMFIT
% This function fits a GMM with k cluster centers to the data specified in
% data. If vary is set to 1, six GMMs will be fit with varying initial
% conditions and the one with maximum likelihood will be returned.
%
% Inputs:
% data      The dataset to fit GMM to
% k         The number of clusters
% vary      1 to vary initial conditions, 0 to just fit 1 GMM
% verbose   1 to output debug info
function [g, t] = gmmfit(data, k, vary, verbose)
    if nargin < 4, verbose = 1; end
    if nargin < 3, vary = 1; end
    
    sigma = (max(data(:, 1)) - min(data(:, 1)))/(2*k);
    
    if verbose == 1, disp(['Fitting with ' num2str(k) ' components']); end
    
    % Build one GMM with uniformly spaced mixtures
    startstruc = struct;
    ncols = ceil(sqrt(k));
    mu_x = linspace(min(data(:, 1)), max(data(:, 1)), ncols)';
    mu_y = linspace(min(data(:, 1)), max(data(:, 1)), ncols)';
    
    mu = [reshape(repmat(mu_x, 1, ncols), [], 1), reshape(repmat(mu_y', ncols, 1), [], 1)];
    
    startstruc.mu = mu(1:k, :);
    startstruc.Sigma = diag([sigma sigma]);
    startstruc.PComponents = 1/k*ones(1, k);
    
    opts = statset('MaxIter', 500);
    
    tic;
    gs{1} = gmdistribution.fit(data, k, 'start', startstruc, 'Regularize', 2*eps, 'Options', opts, 'SharedCov', true);
    
    if vary == 1
        disp('Testing Random Start Conditions');
        % Build five GMMs using random instances from the set
        for i = 2:6
            gs{i} = gmdistribution.fit(data, k, 'Regularize', .01, 'Options', opts);
        end
    end
    t = toc;
    
    NlogLs = cellfun(@(gm)gm.NlogL, gs);
    
    [~, ind] = min(NlogLs);
    
    g = gs{ind};
    
    if verbose == 1, disp(['Fit took ' num2str(t) ' seconds.']); end

end