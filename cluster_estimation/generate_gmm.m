%GENERATE_GMM
% This function generates a Gaussian Mixture Model (GMM) with k modes, and
% samples it n times into data. mu_type, cov_type, and prior_type can be
% either 'random' or 'uniform' (default) 
function [data, g] = generate_gmm(k, n, mu_type, cov_type, prior_type)
    if nargin < 5, prior_type = 'uniform'; end
    if nargin < 4, cov_type = 'uniform'; end
    if nargin < 3, mu_type = 'uniform'; end
    if nargin < 2, n = 50; end
    if nargin < 1, k = 25; end
    
    if strcmp(mu_type, 'uniform') 
        mu = 3.5*(1:ceil(sqrt(k)));
        [mux, muy] = meshgrid(mu, mu);
        mus = [reshape(mux, [], 1), reshape(muy, [], 1)];
        mus = mus(1:k, :);
    else
        mus = sqrt(k)*10*rand(k, 2);
        rows = 1;
        while sum(rows) > 0
            d = pdist(mus);
            
            rows = find(d < 4);
            
            if sum(rows) > 0
                mus = sqrt(k)*10*rand(k, 2);
            end
        end
    end
    
    if strcmp(cov_type, 'uniform')
        covs = repmat(.2*eye(2), [1 1 k]);
    else
        covs = [];
        for i = 1:k
            covs = cat(3, covs, .25+.5*rand(1,2)*eye(2));
        end
    end
    
    if strcmp(prior_type, 'uniform')
        priors = 1/k*ones(k, 1);
    else
        priors = 4+randn(k, 1);
        priors = priors / sum(priors);
    end
    
    g = gmdistribution(mus, covs, priors);
    
    data = random(g, n);

end