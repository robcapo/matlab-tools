%GEN_CLUSTER
% This is a recursive function that implements the DBSCAN method of
% clustering. Algorithm asks for a dataset and a starting observation, as
% well as some search parameters, and will output a cluster as well as the
% path taken to create that cluster. For more information on DBSCAN, see:
% http://en.wikipedia.org/wiki/DBSCAN
%
% Inputs:
% x         All data
% test_ind  Index of initial point to start cluster
% r         Search (neighborhood) radius
% m         Minimum number of observations needed per cluster
% ax_h      Axis handle (for plotting)
% color     Color of current cluster
function [ind, path] = gen_cluster(x, test_ind, r, m, inds, ax_h, color)
    if nargin < 7, color = 'b'; end
    if nargin < 6, ax_h = []; end
    if nargin < 5
        ind = test_ind;
    else
        ind = [inds; test_ind];
    end
    
    xT = x(test_ind, :);
    
    % Find distances of all data to test observation
    d = pdist2(x, xT);
    
    % Find observations that are in the neighborhood defined by r
    neighborhood_inds = find(d <= r);
    
    
    % Plotting
    if ~isempty(ax_h)
        cla;
        hold on;
        scatter(ax_h, x(:, 1), x(:, 2), 'ok');
        scatter(ax_h, x(ind, 1), x(ind, 2), ['+' color]);
        circle(xT(:, 1), xT(:, 2), r);
        hold off;
        
        xlim([-2 5]);
        ylim([-2 5]);
        drawnow;
        pause(1);
    end
    
    
    if length(neighborhood_inds) > m        % must be greater because xT is included in x
        % Create path matrix as [t1, n1; t2, n2; ...] where t specifies a
        % test observation, and n specifies a neighboring observation. All
        % observations are inputted as indices of x
        path = [repmat(test_ind, size(neighborhood_inds)), neighborhood_inds];
        for i = 1:length(neighborhood_inds)
            if (sum(ind == neighborhood_inds(i)) == 0)
                % Run next level of recursion
                [new_inds, paths]= gen_cluster(x, neighborhood_inds(i), r, m, ind, ax_h, color);
                
                % Append results of next level of recursion
                path = [path; paths];
                ind = [ind; new_inds];
            end
            
        end
    else 
        path = [];
    end
    
    % If an ind belongs to a cluster, it only belongs once.
    ind = unique(ind);
    
end