%GMM_STATS
% This function takes a cell vector g of gmdistribution objects and N, the
% size of the data, and  calculates various statistics based on their log
% likelihood
%
% Inputs:
% g     Cell array of gmdistributions
% data  Data on which elements in g were trained
% 
% Output:
% AICs  Akaike information criterion calculated from -log(L)
% AICcs Corrected AIC for small datasets
% BICs  Bayes Information Criterion
% gap   Gap statistic
% gAICs AIC as output from gmdistribution object (not sure how they
%       penalize k, becuase their outputs are not the same
function [AICs, AICcs, BICs, gAICs, NlogLs] = gmm_stats(g, N)
    NlogLs = cellfun(@(gm)gm.NlogL, g);         % -log(L)
    ks = cellfun(@(gm)gm.NComponents, g);       % K
    BICs = cellfun(@(gm)gm.BIC, g);             % BIC
    gAICs = cellfun(@(gm)gm.AIC, g);            % AIC
    
    AICs = 2*ks + 2*(NlogLs);                   % Calculated AIC
    AICcs = (2*ks.*(ks+1))./(N-ks-1);           % Corrected AIC

end