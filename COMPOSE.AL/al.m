% This Active Learning Class contains functionality for actively learning
% in non-stationary environments. It is designed to be owned by an instance
% of the compose class defined in compose.m. The algorithm is described in
% detail in:
% "Active Learning in Non-Stationary Environments," International Joint Conference on Neural Networks (IJCNN 2013)
%
% Created by Rob Capo
% Rowan ECE
classdef al < handle
	% random comment
    properties (SetAccess = protected)
        x_l                 % [MATRIX] NxD labeled data
        labels              % [VECTOR] Labels corresponding to x_l
        
        x_ul                % [MATRIX] NxD unlabeled data
        x_query             % [VECTOR] Indices corresponding to x_ul for which instances to query
        
        density_scores      %*[VECTOR] Scores of instances wrt density
        overlap_scores      %*[VECTOR] Scores of instances wrt overlap
        uncertainty_scores  %*[VECTOR] Scores of instances wrt uncertainty
        scores              %*[VECTOR] Overall scores of instances
        
        default_opts = struct % [STRUCT] Default AL options. Specified in constructor
        
        opts = struct       % [STRUCT] Containing the options for COMPOSE AL:
                            %         Cd: Weight of density metric in score
                            %         Co: Weight of overlap metric in score
                            %         Cu: Weight of uncertainty metric in score
                            %        win: Parzen hypercube dimensions
                            %  noise_thr: Minimum pct of instances for >0 density
                            %        thr: Threshold of score for instance to be queried
                            % uncertainty_class: Type of classifier to use for uncertainty posteriors ('nb' or 'nnet')
    end
    
    methods
        % Constructor
        function obj = al(opts)
            %%%%% Set default options %%%%%
            obj.default_opts.Co = 1;
            obj.default_opts.Cu = 1;
            obj.default_opts.noise_thr = 0;
            obj.default_opts.max_query = 0;
            obj.default_opts.thr = 1;
            obj.default_opts.verbose = 0;
            
            obj.default_opts.uncertainty_class = 'pdf';
            obj.default_opts.h = 1; % hidden layer size for neural net
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.set_opts(opts);
        end
        
        % Find scores for AL object's unlabeled data
        function [scores queries] = query(obj, x)
            if ~isfield(obj.opts, 'win'), error('You must include a Parzen window size parameter ''win'' for density estimation. This parameter should be a row vector with the same number of elements as the dimensionality of the data. You can set this parameter using obj.set_opt(''win'', val);'); end
            if nargin < 2, x = obj.x_ul; end
            
            if obj.opts.verbose >= 1, tic; end
            obj.density_scores = obj.rank_density(x);
            if obj.opts.verbose >= 1, disp(['Density estimation took ' num2str(toc) ' seconds.']); end
            
            if obj.opts.verbose >= 1, tic; end
            obj.overlap_scores = obj.rank_overlap(x);
            if obj.opts.verbose >= 1, disp(['Overlap estimation took ' num2str(toc) ' seconds.']); end
            
            if obj.opts.verbose >= 1, tic; end
            obj.uncertainty_scores = obj.rank_uncertainty(x);
            if obj.opts.verbose >= 1, disp(['Uncertainty estimation took ' num2str(toc) ' seconds.']); end
            
            
            scores = (obj.density_scores - obj.opts.Co * obj.overlap_scores) .* exp(obj.opts.Cu*obj.uncertainty_scores);
            obj.scores = scores;
            
            queries = find(scores >= obj.opts.thr);
            
            if obj.opts.max_query > 0 && length(queries) > obj.opts.max_query
                [~, ind] = sort(scores(queries), 'descend');
                queries = queries(ind(1:obj.opts.max_query));
            end
            
            if obj.opts.verbose >= 1
                disp(['Range on density scores: ' num2str(min(obj.density_scores)) ' to ' num2str(max(obj.density_scores)) '.']);
                disp(['Range on overlap scores: ' num2str(min(obj.overlap_scores)) ' to ' num2str(max(obj.overlap_scores)) '.']);
                disp(['Range on uncertainty scores: ' num2str(min(obj.uncertainty_scores)) ' to ' num2str(max(obj.uncertainty_scores)) '.']);
                disp(['Range on instance scores: ' num2str(min(scores)) ' to ' num2str(max(scores)) '.']);
            end
        end
        
        % Set the data for the AL object
        function set_data(obj, x_ul, x_l, labels)
            if size(x_l, 1) ~= size(labels, 1), error('Please make sure the labeled data has the same number of rows as the labels'); end
            
            obj.x_ul = x_ul;
            obj.x_l = x_l;
            obj.labels = labels;
        end
        
        % Set options
        function o = set_opts(obj, opts)
            if ~isfield(opts, 'win') && ~isfield(obj.opts, 'win'), warning('AL:MissingWin', 'Please specify a Parzen window size parameter ''win'' for density estimation. This parameter should be a row vector with the same number of elements as the dimensionality of the data.'); end
            
            obj.opts = obj.merge_structs(obj.default_opts, obj.opts, opts);
            o = obj.opts;
        end
        
        % Set the value of a specific option
        function o = set_opt(obj, key, val)
            % Set default options if nothing has been set yet
            if numel(fieldnames(obj.opts)) == 0
                obj.opts = obj.default_opts;
            end
            
            obj.opts.(key) = val;
            o = obj.opts;
        end
    end
    
    % Internal ranking methods
    methods
        % Method to rank density value of instances in x
        function scores = rank_density(obj, x)
            r = size(x, 1);
            [ur, uc] = size(obj.x_ul);

            scores = zeros(r, 1);

            for i = 1:r
                x_center = x(i, :);

                % create window boundaries
                box_min = repmat(x_center - obj.opts.win/2, ur, 1);
                box_max = repmat(x_center + obj.opts.win/2, ur, 1);

                % find unlabeled instances inside the window
                x_in = obj.x_ul(sum(and(obj.x_ul >= box_min, obj.x_ul <= box_max), 2) / uc == 1, :);
                n_in = size(x_in, 1);

                if n_in > obj.opts.noise_thr * ur
                    % Calculate normalized Euclidean distances of all
                    % instances in the box (range from 0 to sqrt(2))
                    z = tic;
                    norm_euc = obj.mahal_dist(x_in, x_center, diag((obj.opts.win/2).^2));
                    
                    % remove instances that fall in the rectangular window
                    % but not the elipsoid window
                    norm_euc(norm_euc > 1) = [];
                    
                    
                    ul_dist_sum = mean(exp(-4*norm_euc));
                else
                    ul_dist_sum = 0;
                end

                scores(i) = ul_dist_sum;
            end
        end
        
        % Method to rank overlap value of instances in x
        function scores = rank_overlap(obj, x)
            [r, c] = size(x);
            lr = size(obj.x_l, 1);

            scores = zeros(r, 1);

            for i = 1:r
                x_center = x(i, :);

                % create window boundaries
                box_min = repmat(x_center - obj.opts.win/2, lr, 1);
                box_max = repmat(x_center + obj.opts.win/2, lr, 1);

                % find unlabeled instances inside the window
                x_in = obj.x_l(sum(and(obj.x_l >= box_min, obj.x_l <= box_max), 2) / c == 1, :);
                n_in = size(x_in, 1);

                if n_in > 0
                    % Calculate normalized Euclidean distances of all
                    % instances in the box (range from 0 to sqrt(2))
                    norm_euc = obj.mahal_dist(x_in, x_center, diag(obj.opts.win/2));
                    % remove instances that fall in the rectangular window
                    % but not the elipsoid window
                    norm_euc(norm_euc > 1) = [];
                    
                    l_dist_sum = mean(exp(-4*norm_euc));
                else
                    l_dist_sum = 0;
                end

                scores(i) = l_dist_sum;
            end
        end
        
        % Method to rank entropy values of instances in x
        function scores = rank_entropy(obj, x)
            win = max(x) - min(x);

            r = size(x, 1);
            ulabels = unique(obj.labels);

            scores = zeros(r,1);

            for i = 1:r
                x_center = x(i, :);
                probs = zeros(length(ulabels), 1);
                for j = 1:length(ulabels)
                    
                    probs(j) = mean(mvnpdf(obj.x_l(obj.labels == ulabels(j), :), x_center, diag(win)));
                    
                end
                
                scores(i) = sum(probs.*log2(1./probs));
            end
        end
        
        % Method to rank uncertainty values of instances in x
        function scores = rank_uncertainty(obj, x)
            ulabels = unique(obj.labels);
            nC = length(ulabels);
            switch (obj.opts.uncertainty_class)
                case 'nb' % Naive Bayes
                    nb = NaiveBayes.fit(obj.x_l, obj.labels, 'Prior', 'uniform', 'Distribution', 'kernel');
                    scores = posterior(nb, x);
                case 'nnet' % Neural Network
                    obj.opts.h
                    net = feedforwardnet(obj.opts.h);
                    net = train(net, obj.x_l', full(ind2vec(obj.labels')));
                    scores = net(x')';
                case 'pdf' % Gaussian PDF
                    [r, c] = size(x);
                    scores = zeros(r, nC);
                    
                    for i = 1:r
                        for j = 1:nC
                            scores(i, j) = mean(mvnpdf(obj.x_l(obj.labels == ulabels(j), :), x(i, :), diag(obj.opts.win.^2/4)));
                        end
                    end
            end
            scores = (sum(scores ./ repmat(max(scores, [], 2), 1, nC), 2) - 1) / (nC - 1);
        end
    end
    
    % Bootstrap methods:
    % merge_structs(a,b,c,...) Merge structures with b overwriting a, c
    %    overwriting b, etc.
    %
    % mahal_dist(x,mu,sig) Determines mahalanobis distance of all instances
    %    in x from the distribution specified by mu and sig.
    methods (Static)
        % Merges fields from structure arguments. Each struct overwrites
        % the struct before it
        function s = merge_structs(varargin)
            s = struct;
            
            % Iterate through each struct argument
            for i = 1:nargin
                tStruc = varargin{i};
                
                if isa(tStruc, 'struct')
                    names = fieldnames(tStruc);

                    % Iterate through all nonempty fields in current struct
                    for j = 1:numel(names)
                        if ~isempty(tStruc.(names{j}))

                            s.(names{j}) = tStruc.(names{j});

                        end
                    end
                else
                    warning('merge_structs:BadArg', ['Argument ' num2str(i) ' is not a struct. It has been ignored']);
                end
            end
        end
        
        % Determines mahalanobis distance of all instances in x from the
        % distribution specified by mu and sig.
        function d = mahal_dist(x, mu, sig)
            mu = repmat(mu, size(x, 1), 1);
            d = sqrt(sum(((x - mu) * sig^-1).^2, 2));
        end
    end
end