% This class contains functionality for learning in initially labeled
% non-stationary environments with batches of data. The algorithm is
% described in detail in:
% "COMPOSE: A Semi-Supervised Learning Framework for Initially Labeled Non-Stationary Streaming Data" IEEE Transactions on Neural Networks and Learning Systems, Special issue on Learning in Nonstationary and Dynamic Environments
%
% Created by Karl Dyer and Rob Capo
% Rowan ECE
classdef compose < handle
    
    %% Class Properties - Protected %%
    properties (SetAccess = protected)
        timestep    = 1     % [INTEGER] The current timestep of the dataset
        synthetic   = 0     % [INTEGER] 1 Allows synthetic data during cse and {0} does not allow synthetic data
        n_cores     = 1     % [INTEGER] The number of cores to use for parallel processing segments of code
        
        verbose     = 1     % [INTEGER] Level of feedback displayed during run {default}
                            %    0  : No Information Displayed
                            %   {1} : Command line progress updates
                            %    2  : Plots when possible and Command line progress updates
        
        data                % [CELL] cell array of timesteps each containing a matrix N instances x D features
        labels              % [CELL] cell array of timesteps each containing a vector N instances x 1 - Correct labels
        hypothesis          % [CELL] cell array of timesteps each containing a N instances x 1 - Classifier hypothesis
        core_support        % [CELL] cell array of timesteps each containing a N instances x 1 - binary vector indicating if instance is a core support (1) or not (0)
        
        classifier_func     % [STRING] string corresponding to classifier in ssl class
        classifier_opts     % [STRUCT] struct of options for the selected classifer in ssl class
        
        cse_func            % [STRING] string corresponding to function in cse class
        cse_opts            % [STRUCT] struct of options for the selected cse function in cse class
        
        performance         % [VECTOR] column vector of classifier performances at each timestep
        comp_time           % [MATRIX] matrix of computation time for column 1 : ssl classification, column 2 : cse extraction
        
        plot_mov = {}       % [VECTOR] containing all frames of the compose movie
        plot_callback       % [CALLBACK] containing a user defined function that will be called with:
                            % plot_callback(data, labels, obj);
                            
        dataset
        
        axes_h = 0          % [INTEGER] handle to axis
        figure_xlim         % [VECTOR] Containing the xlim for the plot
        figure_ylim         % [VECTOR] Containing the ylim for the plot
        colormap = [.7 .7 .7; colormap('lines')];
        
        OBJ_learner         % [OBJECT] Learner object from ssl class (see ssl.m)
        OBJ_cse             % [OBJECT] CSE object from cse class (see cse.m)
    end
    
    properties
        step                % [INTEGER] Containing what step the algorithm is on 1: pre-classification 2: pre-cse 3: post-cse 
    end
    
    %% Methods - Visible %%
    methods
        %% COMPOSE Constructor %%
        function obj = compose(dataset, verbose)
            if isempty(strfind(version,'7.9'))
                narginchk(1,2);                                    % Validate number of arguments passed
            end
         
            validateattributes(dataset,{'struct'},{'nonempty'})            % Validate dataset passed in
            if nargin == 2                                                 % If the user passes a verbosity option
                validateattributes(verbose,{'numeric'},{'scalar','integer','>=',0,'<=',2}); % validate the value and
                obj.verbose = verbose;                                     % load it into the object
            end
            
            obj.colormap(2, :) = []; % remove blue so first two colors are red and green
            
            obj.data = dataset.data;         % Initialize cell array with same number of time steps as dataset and add 4th column for classifier hypothesis
            obj.labels = dataset.labels;     % Load the dataset into the labels property
            obj.core_support = dataset.use;  % Load the dataset into the core support property
            obj.hypothesis = deal(cell(size(dataset.data,1),1)); % Initialize cell arrays for core supports and hypothesis
            obj.performance = zeros(1,size(dataset.data, 1)); % Initialize performance matrix with a spot for each timestep
            obj.comp_time = zeros(2,numel(size(dataset.data, 1)));   % Initialize computation time matrix with a spot for each time
            
            obj.dataset = dataset;
            
            % Find overall window in which data will drift
            all_data = cell2mat(dataset.data);
            obj.figure_xlim = [min(all_data(:, 1)) max(all_data(:, 1))];
            obj.figure_ylim = [min(all_data(:, 2)) max(all_data(:, 2))];
            
            if obj.verbose > 1 && size(obj.data,2) > 3      % If verbosity is set to plot (2) and the data has more than 3 features 
                obj.verbose = 1;                            % lower the verbosity to just command window updates (1) and 
                disp('Warning: Your data has more than three dimensions - verbosity level set to 1'); % display warning 
            end
        end
        
        %% Set Number of Cores for Parallel Processing %%
        function set_cores(obj, cores)
            if license('test','Parallel Computing Toolbox')                         % Determine if parallel computing toolbox is available
                if cores > feature('numcores')                                      % if the requested number of cores exceeds available cores
                    disp(['You do not have that many cores on this machine',...     % issue a warning that the cores has been reduced
                        ' - cores has be set to ', num2str(feature('numcores'))]);
                    obj.n_cores = feature('numcores');                              % set n_cores to the maximum number allowed
                else                                                                % otherwise
                    obj.n_cores = cores;                                            % set n_cores to request number
                end
            else
                obj.n_cores = 1;                                                    % If parallel computing toolbox is not available limit n_cores to 1
            end
        end
        
        %% Set Classifier & Options %%
        function set_classifier(obj, user_selection, user_opts)
            if isempty(strfind(version,'7.9'))
                narginchk(2,3);
            end
            
            if isempty(obj.OBJ_learner)                                 % If there is no learner object
                obj.OBJ_learner = ssl(0);                     % Create one
                obj.OBJ_learner.set_data(obj.data{obj.timestep},obj.labels{obj.timestep}); % Load the first batch of data into the learner object
            end
            
            if nargin < 3                                               % If user has not specified options
                obj.OBJ_learner.set_classifier(user_selection);         % load in the selected classifier with default options
            else                                                        % otherwise
                obj.OBJ_learner.set_classifier(user_selection,user_opts); % load in the selected classifier and options
            end
                                  
            obj.classifier_func = user_selection;                       % load record the selected classifier into the object
            obj.classifier_opts = obj.OBJ_learner.classifier_opts;      % load the classifier options into the object
        end
        
        %% Set plot movie frame %%
        function set_plot_mov(obj, frame)
            obj.plot_mov{end+1} = frame;
        end
        
        %% Set Density Estimation Type & Options %%
        function set_cse(obj, user_selection, user_opts)
            if isempty(strfind(version,'7.9'))
                narginchk(2,3)                              % Validate number of arguments
            end    
                
            if isempty(obj.OBJ_cse)                                 % If there is no learner object
                obj.OBJ_cse = cse(0,obj.synthetic,obj);       % Create one
            end
            
            if nargin < 3                                           % If user has specified options
                obj.OBJ_cse.set_boundary(user_selection);           % load in the selected cse and default settings
            else                                                    % OTHERWISE
                obj.OBJ_cse.set_boundary(user_selection,user_opts); % load in the selected cse with user settings
            end
                        
            obj.cse_func = user_selection;                          % load record the selected classifier into the object
            obj.cse_opts = obj.OBJ_cse.boundary_opts;               % load the classifier options into the object
            
        end
        
        function set_plotcallback(obj, func)
            if isa(func, 'function_handle')
                obj.plot_callback = func;
            end
        end
        
        function set_axes(obj, axes)
            if isnumeric(axes)
                obj.axes_h = axes;
            else
                disp('Warning: Axes handles must be numeric Using default axes handle');
            end
        end
        
        %% Run COMPOSE %%
        function filename = run(obj)
            start = obj.timestep;                               % Identify starting time from value in timestep
                                    
            for ts = start:size(obj.data, 1)                    % Iterate over all timesteps from the start to the end of the available data
                obj.timestep = ts;                              % Update the timestep parameter to reflect the current timestep
                obj.hypothesis{ts} = zeros(size(obj.data{ts},1),1);  % Initilize a blank hypothesis vector
                
                % Add any available labeled data
                if sum(obj.core_support{ts} == 1) >= 1          % If there are labeled instances in the current time step 
                    obj.hypothesis{ts}(obj.core_support{ts} == 1) = obj.labels{ts}(obj.core_support{ts} == 1); % copy the labels into the hypothesis
                end

                % Add info from core supports from previous time step
                if start ~= ts                                  % If its not initialization
                    n_cs = sum(obj.core_support{ts-1} == 2);    % find number of core supports from previous time step
                    obj.data{ts}(end+1:end+n_cs,:) = obj.data{ts-1}(obj.core_support{ts-1}==2,:);             % append the current data with the core supports
                    obj.hypothesis{ts}(end+1:end+n_cs,:) = obj.hypothesis{ts-1}(obj.core_support{ts-1}==2,:); % append the hypothesis to include the classes ofthe core supports
                    obj.labels{ts}(end+1:end+n_cs,:) = obj.labels{ts-1}(obj.core_support{ts-1}==2,:);         % append the labels to include the class of the core support
                    obj.core_support{ts}(end+1:end+n_cs,:) = zeros(n_cs,1);                                   % append the core_supports to accomodate the added coresupports of the previouse time step
                end
                
                obj.step = 1; % labeled and unlabeled data
                
                % Plot unlabeled/labeled (given) data
                if obj.verbose == 2
                    if obj.axes_h == 0
                        obj.axes_h = gca;
                    end
                    
                    scatter(obj.axes_h, obj.data{ts}(:, 1), obj.data{ts}(:, 2), 16, obj.colormap(obj.hypothesis{ts} + 1, :));
                    title(['Labeled/Unlabled data for timestep ' num2str(ts)]);
                    xlim(obj.figure_xlim);
                    ylim(obj.figure_ylim);
                    drawnow;
                    pause(.05)
                end
                
                if isa(obj.plot_callback, 'function_handle')
                    obj.plot_callback(obj.data{ts}, obj.hypothesis{ts}, obj);
                end
                
                unlabeled_ind = obj.CLASSIFY(ts);                              % Complete classification
                obj.step = 2; % predicted labels
                
                % Plot labeled data
                if obj.verbose == 2
                    scatter(obj.axes_h, obj.data{ts}(:, 1), obj.data{ts}(:, 2), 16, obj.colormap(obj.hypothesis{ts} + 1, :));
                    title(['Classified Data for timestep ' num2str(ts)]);
                    xlim(obj.figure_xlim);
                    ylim(obj.figure_ylim);
                    drawnow;
                end
                
                if isa(obj.plot_callback, 'function_handle')
                    obj.plot_callback(obj.data{ts}, obj.hypothesis{ts}, obj);
                end
                
                unlabeled_ind = obj.CSE(unlabeled_ind,ts);                    % Complete core support extraction
                if unlabeled_ind == -1
                    break
                end
                obj.step = 3; % core supports are now in.
                
                obj.PERFORMANCE(unlabeled_ind,ts);                            % Calculate performance
                
                % Plot core supports extracted
                if obj.verbose == 2 && ts ~= start                                       % IF plotting is requested
                    
                    cs_data = obj.data{ts}(obj.core_support{ts}==2,:);
                    cs_labels = obj.hypothesis{ts}(obj.core_support{ts}==2,:);
                    
                    scatter(obj.axes_h, cs_data(:, 1), cs_data(:, 2), 16, obj.colormap(cs_labels + 1, :));
                    title(['Core supports from timestep ' num2str(ts)]);
                    xlim(obj.figure_xlim);
                    ylim(obj.figure_ylim);
                    drawnow;
                    pause(.05);
                end
                
                % User defined plotting
                if isa(obj.plot_callback, 'function_handle')
                    cs_data = obj.data{ts}(obj.core_support{ts}==2,:);
                    cs_labels = obj.hypothesis{ts}(obj.core_support{ts}==2,:);
                    
                    obj.plot_callback(cs_data, cs_labels, obj);
                end
                
                if obj.verbose == 1
                    disp(['Finished iteration ' num2str(ts) ' with performance ' num2str(obj.performance(ts))]);
                end
            end               
            
            filename = strcat('temp_',num2str(cputime+floor(10000*rand(1,1))),'.mat');
            save(filename,'obj')
          
        end
    end
    
    %% Methods - Protected %%
    methods (Access = protected)
    
        %% Classify Current Data %% 
        function unlabeled_ind = CLASSIFY(obj,ts)
            % Sort the data so that labeled data is at the top, unlabeled data to follow
            [obj.hypothesis{ts}, sortID] = sort(obj.hypothesis{ts},'descend');  % sort so unlabeled data is at bottom of list
            obj.data{ts} = obj.data{ts}(sortID,:);                     % sort data to match hypothesis shifts
            obj.labels{ts} = obj.labels{ts}(sortID,:);                 % sort labels to match hypothesis shifts
            obj.core_support{ts} = obj.core_support{ts}(sortID,:);     % sort core_supports to match hypothesis shifts
            unlabeled_ind = obj.hypothesis{ts} == 0;                   % keep track which instances were originally unlabeled so we know which instances to use for performance caclulations

            % Classify Data
            obj.OBJ_learner.set_data(obj.data{ts},obj.hypothesis{ts}); % load the current time step into the learner object
            tic; obj.OBJ_learner.classify;                             % classifiy the data using the learner object
            
            obj.comp_time(ts,1) = toc;                                 % store computation time
            obj.hypothesis{ts} = obj.OBJ_learner.labels;               % copy the classification into the compose object

            if obj.verbose >= 1                                        % if updates are requested
                disp(['Classification took ',num2str(obj.comp_time(ts,1)),' seconds.']); % display classification time
            end
        end    

        %% Extract Core Supports from Current Data %% 
        function unlabeled_ind = CSE(obj, unlabeled_ind,ts)
            % Remove duplicate data instances
            [obj.data{ts}, sortID] = unique(obj.data{ts},'rows');      % check the data by rows and remove duplicate instances keeping track what was removed
            obj.labels{ts} = obj.labels{ts}(sortID,:);                 % remove labels of removed instances
            obj.core_support{ts} = obj.core_support{ts}(sortID,:);     % remove core_support data of removed instances
            obj.hypothesis{ts} = obj.hypothesis{ts}(sortID,:);         % remove hypothesis data of removed instances
            unlabeled_ind = unlabeled_ind(sortID,:);                   % remove unlabeled data indices of removed instances

            % Sort be class
            class = unique(obj.hypothesis{ts});                        % identify number of class
            [obj.hypothesis{ts},sortID] = sort(obj.hypothesis{ts});    % sort by class
            obj.data{ts} = obj.data{ts}(sortID,:);                     % match data with sort
            obj.labels{ts} = obj.labels{ts}(sortID,:);                 % match labels with sort
            obj.core_support{ts} = obj.core_support{ts}(sortID,:);     % match core_supports with sort
            unlabeled_ind = unlabeled_ind(sortID,:);                   % match unlabeled data with sort

            tic;                                                       % start timer
            c_offset = 0;                                              % c_offset is used to keep track how many instances have been analyzed so each class can be returned in the correct spot after cse
            for c = 1:numel(class)                                     % iterate over each class
                class_ind = obj.hypothesis{ts} == class(c);            % find indices of class being analyzed
                obj.OBJ_cse.set_data(obj.data{ts}(class_ind,:));       % set the current class as the data to be analyzed
                inds = obj.OBJ_cse.extract;                            % extract core supports of class being analyzed
if inds == -1
    unlabeled_ind = -1;
    break
end
                obj.core_support{ts}(inds + c_offset) = 2;             % indicate those instances as a core support denoted by a 2
                c_offset = c_offset + obj.OBJ_cse.n_instances;         % update the offset so the classes are stored in correct spot
            end
            obj.comp_time(ts,2) = toc;                                 % store computation time

            if obj.verbose >= 1                                        % If updates are requested
                disp(['Compaction took ',num2str(obj.comp_time(ts,2)),' seconds.']); % display core support extratction time
            end
        end
        
        %% Calculate Performance of Classification %% 
        function PERFORMANCE(obj, unlabeled_ind,ts)
            diff = obj.hypothesis{ts}(unlabeled_ind,:)-obj.labels{ts}(unlabeled_ind,:); % not equal to 0 indicates an instance is different from true label
            obj.performance(ts) = sum(diff==0)/sum(unlabeled_ind)*100; % store the number of correct classifications
        end
    end
end