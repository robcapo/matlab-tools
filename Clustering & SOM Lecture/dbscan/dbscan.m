% This script implements the DBSCAN algorithm for clustering arbitrary
% data. It generates three separated gaussians and attempts to separate
% them into clusters
clear all;close all;

% Generate data
x = [mvnrnd([0 0], [.1 0; 0 .1], 10);
     mvnrnd([1.5 3], [.7 0; 0 .1], 10);
     mvnrnd([3, 0], [.1 0; 0 .1], 10)];


% Keep track of indices that we haven't clustered yet
unchecked_inds = 1:length(x);

% Color map for different clusters
c = 'rgbcmyko';
i = 1;
clusters = {};

figure;
while ~isempty(unchecked_inds) % while there's still instances we haven't clustered ...
    % start with a random observation
    test_ind = round(1+(length(unchecked_inds)-1)*rand());
    
    % Generate a cluster from that observation and return all indices
    % belonging to the cluster, and the path that DBSCAN took to find it
    % (used for plotting)
    [ind, path] = gen_cluster(x, unchecked_inds(test_ind), 1, 5, [], gca, c(i));
    
    % Determine a threshold of 5 as the minimum number of observations to
    % create a cluster
    if length(ind) >= 5
        clusters{i, 1} = ind;
        clusters{i, 2} = path;
        clusters{i, 3} = c(i);

        % Plotting 
        hold on;
        cla;
        title(['Cluster ' num2str(i) ' Created']);
        for k = 1:size(clusters, 1)
            scatter(x(:, 1), x(:, 2), 'ok');
            scatter(x(clusters{k, 1}, 1), x(clusters{k, 1}, 2), ['*' clusters{k, 3}]);

            tpath = clusters{k, 2};

            for j = 1:length(tpath)
                plot([x(tpath(j, 1), 1); x(tpath(j, 2), 1)], [x(tpath(j, 1), 2); x(tpath(j, 2), 2)]);
            end
        end
        pause;
        % End plotting
        
        i = i + 1;
    end
    
    % Remove all instances in current cluster from unchecked instances
    unchecked_inds(ismember(unchecked_inds, ind)) = [];
    
end