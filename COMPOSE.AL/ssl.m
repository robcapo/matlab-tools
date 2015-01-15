% This class includes functionality for various semi-supervised learning
% algorithms including: S3VM, Label Propagation, Label Spreading, and
% Cluster and Label. It is intended to be used as a property of an instance
% of the compose class defined in compose.m
%
% Created by Karl Dyer and Rob Capo
classdef ssl < handle
    % SSL is a class of semi-supervised learning classifiers that may be 
    % used for stationary or non-stationary environments. Depending on the
    % classifier choosen a variety of and class balancing techniques are
    % available to reduce SSL problem of assigning all data to one class
    %
    % Created 09.04.2012 by Karl Dyer and Rob Capo
    %                       Rowan University - SPPRL Lab
            
    %% Class Properties - Visible %%
    properties (SetAccess = protected)
        verbose     = 2     % [INTEGER] : Control outputs to screen
                            %   0 : Suppress all output
                            %   1 : Give text updates to command window
                            %   2 : Plot data when dimensionality allows and give text updates to command window
        
        data                % [MATRIX] N instances x D dimensions : Features of data with labeled data grouped at top of matrix
        labels              % [VECTOR] N instances x 1 (integer format) -OR- [MATRIX] N instances x C classes (one-hot format)
                            %   Unlabeled data should be represented by 0 or a vector of 0's
                            %   Note: Output will match the users input format
        
        classifier          % [STRING] : Type of SSL classifier to use
        classifier_opts     % [STRUCT] : Options that correspond with SSL Classifier selected - see individual methods for options
        
        balance             % [STRING] : Type of class balancing to use
        balance_opts        % [STRUCT] : Options that correspond with Balance Function selected - see individual methods for options
    end
    
    %% Class Properties - Hidden %%
    properties (Hidden)
        % All of these properties are calculated from the data when loaded
        % using the set_data method. They are properties that are useful in
        % a variety of classifiers. The automatically calculated
        % class_priors at load may be overwritten by the user if
        % class_priors is passed in classifier_opts or balance_opts
        
        n_features          % [INTEGER] : Number of features in data (i.e. dimensionality of data)
        n_classes           % [INTEGER] : Number of classes different class labels
        n_instances         % [INTEGER] : Number of instances in data
        n_labeled           % [INTEGER] : Number of labeled instances in data
        n_unlabeled         % [INTEGER] : Number of unlabeled instances in data
        
        input_label_format  % [STRING] : Format of labels passed by user - 'integer' OR 'vector'
        input_label_ids     % [VECTOR] : Records the class identifiers of the labels passed by user
        
        label_format        % [STRING] : Current format of labels
    end
    
    %% Class Properties - Hidden & Constant
    properties (Hidden, Constant)
        % The cells below contain text strings that match the SSL 
        % classifiers and class balance methods available in this object
        % If if other classifiers or balancing methods are added to this
        % class these cells must be modified to include those methods
        valid_classifier = {'s3vm'; 'label_prop'; 'label_spread'; 'cluster_n_label'};
        valid_balance = {'none','mass','bid'}; %,'reg'};
    end
    
    %% Initialization and Executing Methods %%
    methods
        %% SSL Constructor %%
        function obj = ssl(verbose)
            if isempty(strfind(version,'7.9'))
                narginchk(0,1);                                                     % Verify correct number of arguments are passed
            end
            if nargin == 1                                                                  % If user specifies verbosity
                validateattributes(verbose,{'numeric'},{'scalar','integer','>=',0,'<=',2})  % Validate verbose and
                obj.verbose = verbose;                                                      % Load verbosity level into object
            end
        end
        
        %% Load Data & Labels %%
        function set_data(obj, data, labels)
            if isempty(strfind(version,'7.9'))
                narginchk(3,3);                                 % Verify correct number of arguments are passed
            end
            
            validateattributes(data,{'numeric'},{'real','nonnan'})      % Validate data
            validateattributes(labels,{'numeric'},{'real','nonnan'})    % Validate labels
            
            if size(data,1) == size(labels,1)       % Validate data and labels have same number of instances
                obj.data   = data;                  % load the data into the object
                obj.labels = labels;                % load the labels into the object
                
                % Propagate Hidden Properties    
                [obj.n_instances, obj.n_features] = size(obj.data); % Obtain size information from data
                obj.n_unlabeled = sum(sum(obj.labels,2)==0);        % Sum across the each row and count how many instances have a zero class meaning unlabled
                obj.n_labeled = obj.n_instances - obj.n_unlabeled;  % The remaining instances must be labeled
                
                % Properties conditional to label format
                if all(obj.labels(:)==0 | obj.labels(:)==1)         % IF labels contains only zeros and ones (i.e. vector format) then
                    obj.input_label_format = 'vector';              % store the input format as 'vector' and
                    obj.input_label_ids = (1:size(obj.labels,2))';  % store the label identifiers as 1 through the size of the matrix
                    obj.label_format = 'vector';                    % set label format flag to indicate type vector
                else                                                % OTHERWISE
                    obj.input_label_format = 'integer';             % store the input format as 'integer' and
                    obj.input_label_ids = unique(obj.labels)';      % store the numbers used to identify each class
                    obj.input_label_ids(obj.input_label_ids==0)=[]; % remove the zero "class" used to indicate no class
                    obj.label_format = 'integer';                   % set label format flag to indicate type integer
                end
                obj.n_classes = numel(obj.input_label_ids);         % store number of classes 
                
                if obj.n_features > 3 && obj.verbose == 2           % If dimensionalty is to large to plot and plotting was requested
                    obj.verbose = 1;                                % change the verbosity level to 1
                    disp('Dimensionality of data to large to plot - verbosity changed to 1') % notify user and
                end
                
            else                                    % If validation fails issue an error
                error('Data must be instances x features and labels must have the same number of instances');
            end
        end
        
        %% Set Classifier Type & Options %%
        function set_classifier(obj, classifier_selection, opts)
            if isempty(strfind(version,'7.9'))  
                narginchk(2,3);                                     % Verify correct number of arguments are passed
            end
            validateattributes(classifier_selection,{'char'},{'nonempty'}); % validate users classifer selection
            if nargin == 3                                                  % if opts parameter is passed
                validateattributes(opts,{'struct','char'},{'nonempty'});    % validate opts 
            else                                                            % if no options are passed
                opts = [];                                                  % make options a blank variable
            end
            
            obj.classifier_opts = [];                                       % ensure options are cleared from any previous classifiers
            
            if sum(strcmp(classifier_selection,obj.valid_classifier)) == 1  % Check user requested classifier against validation list
                obj.classifier = str2func(classifier_selection);            % load classifier selection into object and convert to function handle
                obj.set_default_opts('classifier');                         % load the default options for specified classifier
                if ~isempty(opts)                                           % if the user passed options
                    obj.set_user_opts('classifier',opts);                   % overwrite any defaults with user options that validate
                end
            else
                display([classifier_selection,...
                    ' is an invalid classifier selection - choose from: ']); % if invalid selection notify user
                display(obj.valid_classifier);                               % and provide acceptable selections
            end
        end
        
        %% Set Class Balance Type & Options %%
        function set_balance(obj, balance_selection, opts)
            if isempty(strfind(version,'7.9'))
                narginchk(2,3);                                     % Verify correct number of arguments are passed
            end
            validateattributes(balance_selection,{'char'},{'nonempty'});    % validate users balance selection
            if nargin == 3                                                  % if opts parameter is passed
                validateattributes(opts,{'struct'},{'nonempty'});           % validate opts 
            else                                                            % if no options are passed
                opts = [];                                                  % make options a blank variable
            end
            
            obj.balance_opts = [];                                          % ensure options are cleared from any previous classifiers
            
            if sum(strcmp(balance_selection,obj.valid_balance)) == 1        % Check user requested class balancing method against validation list
                obj.balance = str2func(balance_selection);                  % load balancing method into object and convert to function handle
                obj.set_default_opts('balance');                            % load the default options for specified balance method
                if ~isempty(opts)                                           % if the user passed options
                    obj.set_user_opts('balance',opts);                      % overwrite any defaults with user options that validate
                end
            else
                display([class_balance, ...
                    ' is an invalid class balance selection - choose from: ']); % if invalid selection notify user
                display(obj.valid_balance);                                 % and provide acceptable selections
            end
        end
        
        %% Run Classifier %%
        function classify(obj)
            if isempty(obj.data)                                            % If data hasn't been loaded 
                error('You must load data before classifying');             % issue error to user
            end
            
            if isempty(obj.classifier)                                      % Determine if classifier has been set, if not
                display('Classifier not set - default classifier and options loaded') % issue warning that a default is being used and
                obj.set_classifier('cluster_n_label');                      % load cluster and label as the default classifier
            end;
            
            if isempty(obj.balance)                                         % Determine if classifier has been set, if not
                display('Class balance method was not set - no balancing will occur') % issue warning that a default is being used and
                obj.set_balance('none');                                    % load no balancing as default
            end;
            
            if obj.verbose == 2         % If plotting is requested
                obj.ssl_plot;            % plot the labeled and unlabeled data
            end
            
            obj.classifier(obj);                      % Run the classifier selected
            obj.balance(obj);                         % Run the balance method
            
            if ~strcmp(obj.input_label_format, obj.label_format) % if the current format is different then the input format
                switch obj.label_format                          % use the current format to determine what change to make
                    case 'vector'                                % if it is currently vector
                        obj.vec2int(obj.input_label_ids);         % change it to intger
                    case 'integer'                               % if it is currently integer
                        obj.int2vec;                              % change it to vector
                end
            end
                        
            if obj.verbose == 2         % If plotting is requested
                obj.ssl_plot            % plot the results of the classification
            end
            
        end
        
    end
    
    %% Protected General Methods %%
    methods (Access = protected)
        %% Set Default Options %%
        function set_default_opts(obj,type)
            switch type
                case 'classifier'
                    switch func2str(obj.classifier)
                        case 's3vm'
                            % Defaults loaded by external s3vm application
                        case 'label_prop'
                            obj.classifier_opts.sigma = 1;     % set default sigma to null so algorithm will calculate it
                        case 'label_spread'
                            obj.classifier_opts.alpha = 0.5;    % set default alpha parameter - balance between local and global influences
                            obj.classifier_opts.sigma = [];     % set default sigma to null so algorithm will calculate it
                        case 'cluster_n_label'
                            obj.classifier_opts.kmin = obj.n_classes; % set default minimum number of clusters
                            obj.classifier_opts.k = 2*obj.n_classes;  % set default starting number of clusters
                    end
                case 'balance'
                    switch func2str(obj.balance)
                        case 'none'
                            % No free parameters to set
                        case {'mass', 'bid'}
                            if strcmp(obj.label_format,'vector')                                % IF labels are vector format then
                                obj.balance_opts.class_priors = sum(obj.labels)./obj.n_labeled; % calculate and store class priors from available labeled data 
                            else                                                                % OTHERWISE labels are integer format so
                                obj.balance_opts.class_priors = histc(obj.labels,obj.input_label_ids)./obj.n_labeled; % calculate and store class priors from available labeled data
                            end
                    end
            end
        end
        
        %% Set User Options
        function set_user_opts(obj,type,opts)
            switch type
                case 'classifier'
                    switch func2str(obj.classifier)
                        case 's3vm'                             % s3vm gets handled differently because is is a text string
                            validateattributes(opts,{'char'})   % validate opts
                            obj.classifier_opts = opts;         % load opts into object
                            
                            % TODO - can add parser to validate string options
                            
                        otherwise                               % all classifiers with structure formated option get handled the same
                            user_fields = fieldnames(opts);     % make list of struct fields passed by user
                            for field = 1:numel(user_fields)    % iterate thru each struct field passed
                                if isfield(obj.classifier_opts,user_fields(field))  % If the field passed by user is a valid field
                                    obj.classifier_opts.(user_fields{field}) = opts.(user_fields{field}); % update the option to the users selection
                                else                                                                      % OTHERWISE
                                    disp(['Warning: Option ', user_fields{field},...                      % inform user the field is not a valid parameter
                                        ' is not a valid option for ', func2str(obj.classifier),...
                                        ' classifier'])                   
                                end
                            end
                    end
                case 'balance'
                    switch func2str(obj.balance)
                        case 'none'
                            % No free parameters to set
                        
                        otherwise                               % all classifiers with structure formated option get handled the same
                            user_fields = fieldnames(opts);     % make list of struct fields passed by user
                            for field = 1:numel(user_fields)    % iterate thru each struct field passed
                                if isfield(obj.balance_opts,user_fields(field))  % If the field passed by user is a valid field
                                    obj.balance_opts.(user_fields{field}) = opts.(user_fields{field}); % update the option to the users selection
                                else                                                                      % OTHERWISE
                                    disp(['Warning: Option ', user_fields{field},...                      % inform user the field is not a valid parameter
                                        ' is not a valid option for ', func2str(obj.balance),...
                                        ' class balance method'])                   
                                end
                            end
                    end
            end
        end
        
        %% Convert Labels from Integer Format to Vector Format %%
        function rev_key = int2vec(obj,collapse,nC)
            % This function will create a NxC matrix of class labels from
            % a vector of integers. This can be used to take an integer
            % vector of class labels and develop a set of vectors used for
            % classification that can be used to train a classifier like
            % the MLP or SVM which give posterior probabilites.
            %
            % FUNCTION CALLS:
            % int2vec(obj) - converts obj.labels property from integer 
            %   format to one-hot vector format using default collapse and 
            %   setting nC parameter equal to the number of classes
            %   represented in the current obj.labels property
            %
            % int2vec(obj,collapse) - converts obj.labels property from 
            %   integer format to one-hot vector format but allows control 
            %   over spanning method for unused classes
            %      0  : Maintain a blank column when a class is unused
            %     {1} : Collapse into sequential classes
            %
            % int2vec(obj,...,nC) - converts obj.labels property from 
            %   integer format to one-hot vector format but allows control 
            %   over total number of classes in case a class will be added 
            %   at a later time step. If nC is greater than the current 
            %   number of classes represented additional columns will be 
            %   added to reserve space for future classes
            %
            % Created 09.04.2012 by Karl Dyer
            
            %% Set variable inputs
            if strcmp(obj.label_format,'vector')       % If the labels are already a vectot 
                return;                                % do nothign and return to call
            end
            
            if nargin == 2; nC = numel(unique(obj.labels))-1; end
            if nargin == 1; nC = numel(unique(obj.labels))-1; collapse = 1; end;
                
            %% Construct New Labels
            % Select appropriate spanning method
            if collapse == 0                                        % If the spanning method doesn't want to collapse missing classes
                all_classes = min(obj.labels):max(obj.labels);      % Find the min and max class labels and create a vector of all classes in between
            else                                                    % If collapsing is allowed
                all_classes = unique(obj.labels);                   % Find the uniquw class labels
            end
            all_classes(all_classes == 0) = [];                     % Remove the zero "class label" used to represent unlabeled data
            
            calc_nC = numel(all_classes);                           % Determine the number of classes to be represented
            vect_labels = zeros(size(obj.labels,1),calc_nC);        % initialize a zero vector large enough to accomodate all instances and classes

            % Construct one-hot matrix
            for i = 1:calc_nC
                temp_row = zeros(1,calc_nC); temp_row(i) = 1;
                if sum(obj.labels == all_classes(i)) > 0
                    vect_labels(obj.labels == all_classes(i),:) = repmat(temp_row,sum(obj.labels == all_classes(i)),1);
                end
            end

            rev_key = all_classes; % Export the key that will reverse the transformation
                    
            % Append matrix to accomodate requested extra classes and add
            % another class info to the reversal key
            if calc_nC < nC
                vect_labels = [vect_labels zeros(size(vect_labels,1),nC-calc_nC)];
                rev_key = [rev_key max(rev_key)+1:max(rev_key)+nC-calc_nC];
            end
            
            obj.labels = vect_labels;       % Store the new labels into the object
            obj.label_format = 'vector';    % Change label format flag
        end
        
        %% Convert Labels from Vector Format to Integer Format %%
        function vec2int(obj,rev_key)
            % This function converts one-hot vector labels to an integer
            % representation for more compact storage and ease of use with
            % plotting features. Function may be called:
            %
            % FUNCTION CALLS:
            % vec2int(obj) - converts obj.labels property from one-hot 
            %   vector format to integer format
            % 
            % vec2int(obj,rev_key) - returns obj.label property from 
            %   one-hot vector format to integer format using the same 
            %   class identifiers originally provided by user in rev_key
            %     rev_key : [VECTOR] 1 x C classes - original class labels
            %
            % Created 09.05.2012 by Karl Dyer
            
            if strcmp(obj.label_format,'integer')     % If the format is already integer
                return;                               % do nothign and return to call
            end
            
            if nargin == 2; 
                uncoder = repmat(rev_key,size(obj.labels,1),1);
                %Return to orignal matrix thru point by point
                %multiplication with a replicated matrix of the reversal
                %key - where the one is makes that instance retain that class
                int_labels = sum(obj.labels.*uncoder,2);
            else
                % Find the column index of the the 1 in each row
                [~,int_labels] = max(obj.labels,[],2);
                int_labels(sum(obj.labels,2)==0)=0;
            end
            
            obj.labels = int_labels;        % Store the labels into the object
            obj.label_format = 'integer';   % Change label format flag
        end
        
        %% Plotting Feature %%
        function ssl_plot(obj)
            
            changed = 0;                                        % Initialize label changed flag
            set(0,'DefaultAxesColorOrder',hsv(obj.n_classes));  % create enough colors for each class and set as default colors
            hold all;                                           % make each plot command cycle through the default colors
                       
            if strcmp(obj.label_format,'vector')       % If in vector format
                obj.vec2int                            % convert to integer format for easy plotting
                changed = 1;                           % change flag that format was changed
            end
            
            if sum(sum(obj.labels,2)==0) ~= 0          % If there are any unlabeled instances
                cla;                                   % Clear any current plot data
                switch obj.n_features                  % send to correct plotting calls based on dimensionality
                    case 1
                    case 2
                        legend_entry = cell(1,obj.n_classes+1);                                                 % Initialize legend entry cell
                        scatter(obj.data(obj.labels == 0,1),obj.data(obj.labels == 0,2),'k.');                  % Plot unlabeld data
                        legend_entry{1} = 'Unlabeled Data';                                                     % Make legend entry for unlabeld data
                        for class = 1:obj.n_classes                                                             % Iterate through each class
                            id = obj.input_label_ids(class);                                                    % Get class label
                            scatter(obj.data(obj.labels == id,1),obj.data(obj.labels == id,2),'marker','p');    % Plot labeld data by class
                            legend_entry{class+1} = ['Class ', num2str(id)];                                    % Record class to legend
                        end
                        xlabel('Feature 1')                                                                     % label x axis 
                        ylabel('Feature 2')                                                                     % label y axis 
                        title({['Classifier: ', func2str(obj.classifier)];['Balance Method: ', func2str(obj.balance)]}); % Add a title
                        legend(legend_entry);                                                                   % add a legend
                    case 3
                        legend_entry = cell(1,obj.n_classes+1);                                                 % Initialize legend entry cell
                        scatter3(obj.data(obj.labels == 0,1),obj.data(obj.labels == 0,2),obj.data(obj.labels == 0,3),'k.'); % Plot unlabeld data
                        legend_entry{1} = 'Unlabeled Data';                                                     % Make legend entry for unlabeld data
                        for class = 1:obj.n_classes                                                             % Iterate through each class
                            id = obj.input_label_ids(class);                                                    % Get class label
                            scatter(obj.data(obj.labels == id,1),obj.data(obj.labels == id,2),obj.data(obj.labels == 0,3),'marker','p'); % Plot labeld data by class
                            legend_entry{class+1} = ['Class ', num2str(id)];                                    % Record class to legend
                        end
                        xlabel('Feature 1')                                                                     % label x axis 
                        ylabel('Feature 2')                                                                     % label y axis
                        zlabel('Feature 3')                                                                     % label z axis
                        title({['Classifier: ', func2str(obj.classifier)];['Balance Method: ', func2str(obj.balance)]}); % Add a title
                        legend(legend_entry);                                                                   % add a legend
                end
            else                                       % OTHERWISE all instances are labeled
                switch obj.n_features                  % send to correct plotting calls based on dimensionality
                    case 1
                    case 2
                        for class = 1:obj.n_classes                         % iterate through each class
                            id = obj.input_label_ids(class);                % get class id
                            scatter(obj.data(obj.labels == id,1),obj.data(obj.labels == id,2),'marker','x','sizedata',30,'linewidth',1); % plot data
                        end
                    case 3
                        for class = 1:obj.n_classes                         % iterate through each class
                            id = obj.input_label_ids(class);                % get class id
                            scatter3(obj.data(obj.labels == id,1),obj.data(obj.labels == id,2),obj.data(obj.labels == id,3),'marker','x','sizedata',20,'linewidth',2); % plot data
                        end
                end
            end
            
            if changed == 1                            % If conversion was made for plotting
                obj.int2vec                            % convert back to original format
            end
                
        end
    end
    
    %% S3VM & DEPENDENCIES %%
    % For S3VM obj.classifier_opts use the following call 
    % obj.set_classifier('s3vm','text string of options') - see fold below
    %{ 
     Example : '-v 0 -p [0.5] -t 2 -g 1.5'

     General options:
              -?          - this help
              -v [0..3]   - verbosity level (default 1)
     Learning options:
              -z {c,r,p}  - select between classification (c), regression (r), and 
                            preference ranking (p) (see [Joachims, 2002c])
                            (default classification)          
              -c float    - C: trade-off between training error
                            and margin (default [avg. x*x]^-1)
              -w [0..]    - epsilon width of tube for regression
                            (default 0.1)
              -j float    - Cost: cost-factor, by which training errors on
                            positive examples outweight errors on negative
                            examples (default 1) (see [Morik et al., 1999])
              -b [0,1]    - use biased hyperplane (i.e. x*w+b0) instead
                            of unbiased hyperplane (i.e. x*w0) (default 1)
              -i [0,1]    - remove inconsistent training examples
                            and retrain (default 0)
     Performance estimation options:
              -x [0,1]    - compute leave-one-out estimates (default 0)
                            (see [5])
              -o ]0..2]   - value of rho for XiAlpha-estimator and for pruning
                            leave-one-out computation (default 1.0) 
                            (see [Joachims, 2002a])
              -k [0..100] - search depth for extended XiAlpha-estimator
                            (default 0)
     Transduction options (see [Joachims, 1999c], [Joachims, 2002a]):
              -p [0..1]   - fraction of unlabeled examples to be classified
                            into the positive class (default is the ratio of
                            positive and negative examples in the training data)
     Kernel options:
              -t int      - type of kernel function:
                             0: linear (default)
                             1: polynomial (s a*b+c)^d
                             2: radial basis function exp(-gamma ||a-b||^2)
                             3: sigmoid tanh(s a*b + c)
                             4: user defined kernel from kernel.h
              -d int      - parameter d in polynomial kernel
              -g float    - parameter gamma in rbf kernel
              -s float    - parameter s in sigmoid/poly kernel
              -r float    - parameter c in sigmoid/poly kernel
             -u string   - parameter of user defined kernel
     Optimization options (see [Joachims, 1999a], [Joachims, 2002a]):
              -q [2..]    - maximum size of QP-subproblems (default 10)
              -n [2..q]   - number of new variables entering the working set
                            in each iteration (default n = q). Set n<q to prevent
                            zig-zagging.
              -m [5..]    - size of cache for kernel evaluations in MB (default 40)
                            The larger the faster...
              -e float    - eps: Allow that error for termination criterion
                            [y [w*x+b] - 1] = eps (default 0.001) 
              -h [5..]    - number of iterations a variable needs to be
                            optimal before considered for shrinking (default 100) 
              -f [0,1]    - do final optimality check for variables removed by
                            shrinking. Although this test is usually positive, there
                            is no guarantee that the optimum was found if the test is
                            omitted. (default 1) 
              -y string   -> if option is given, reads alphas from file with given
                             and uses them as starting point. (default 'disabled')
              -# int      -> terminate optimization, if no progress after this
                             number of iterations. (default 100000)
     Output options: 
              -l char     - file to write predicted labels of unlabeled examples 
                            into after transductive learning 
              -a char     - write all alphas to this file after learning (in the 
                            same order as in the training set) 
    %}
    methods (Access = protected)
        %% S3VM Classifier %%
        function s3vm(obj)
            % Runs Joachims svm_light package coded in C to classify data 
            % using semi-supervised methods. The matlab code reformats the 
            % data before passing to executable located int the rood and 
            % then restores the data to a usable matlab format after the 
            % executable runs.
            %
            % FUNCTION CALL:
            % s3vm_light(obj) - uses obj.data and obj.labels to classify
            %   unlabeled data using using Joachims S3VM_light in root
            %
            % TODO
            %   - Add plotting capabilites
            %
            % Created 09.05.2012 by Karl Dyer

            %% Prepare Data for S3VM
            % Ensure labels are in a integer format for proper use with S3VM            
            if strcmp(obj.label_format, 'vector') % If the input format is vector 
                obj.vec2int;                            % convert it to integer
            end
               
            % Ensure there are only two classes and a neutral class for unlabeled data
            if numel(obj.input_label_ids) ~= 2 || sum(obj.labels == 0) < 1                      % If data doesn't have two classes and unlabeld data
                error('S3VM is for two class problems only and must contain unlabeled data')    % throw an error
            end
                                   
            % Convert Labels to [-1 0 +1] class labels
            templabels = zeros(size(obj.labels));                    % Create temporary labels so transfrom is error free
            templabels(obj.labels==obj.input_label_ids(1)) = -1;     % Replace class #1 with a -1 label
            templabels(obj.labels==obj.input_label_ids(2)) =  1;     % Replace class #2 with a +1 label
            obj.labels = templabels; clear('templabels');            % Load templabels into object and delete temp
                        
            %% Format data for use with SVM Light
            namehash = num2str(floor(sum(1e3*rand(1,2))));      % attempt to create a unique filename
            while exist(strcat(namehash,'_train.dat'),'file')       % make sure it does not exist
                namehash = num2str(floor(sum(1e3*rand(1,2))));  % if it does create a new filename
            end
            traindata = obj.s3vm_format(namehash,'train');      % obtain training data
            testdata = obj.s3vm_format(namehash,'test');     % create testing data 
            
            %% Run S3VM Light
            % Determine the operating system matlab is running on and call the correct svm light executable from the root
            if ispc
                eval(['!svm_learn ',obj.classifier_opts,' ',traindata, strcat(' model',namehash)]);
                eval(['!svm_classify ',testdata, strcat(' model',namehash,' predictions', namehash)]);
            elseif ismac
%                 Check executable calls for mac
%                 eval(['!svm_learn ',obj.classifier_opts,' ',traindata, strcat(' model',namehash)]);
%                 eval(['!svm_classify ',testdata, strcat(' model',namehash,' predictions', namehash)]);
            elseif isunix
%                 Check executable calls for unix/linux
%                 eval(['!svm_learn ',obj.classifier_opts,' ',traindata, strcat(' model',namehash)]);
%                 eval(['!svm_classify ',testdata, strcat(' model',namehash,' predictions', namehash)]);
            end
            
            %% Read SVM Light Outputs back into matlab
            Z = zeros(obj.n_unlabeled,1);       % Construct a hypothesis vector same size as unlabeled instances
            i = 1;                              % Initialize line counter to read file in to 1

            fileID = fopen(strcat('predictions',namehash));                             % Attempt to open predictions file
            while fileID == -1                                                          % If the file could not be found
              display('Error reading predictions file from disk - regenerating now');   % display an error
              eval(['!svm_learn ',obj.classifier_opts,' ',traindata, strcat(' model',namehash)]); % then rerun the training
              eval(['!svm_classify ',testdata, strcat(' model',namehash,' predictions', namehash)]);  % and classifying 
              fileID = fopen(strcat('predictions',namehash));                           % then retry opening the file
            end
            one_line = fgetl(fileID);           % read the first line of the predictions file
            while ischar(one_line)              % while there is more information in the file
                Z(i,1) = str2double(one_line);  % convert the current line from text to a double and store in Z vector
                i = i+1;                        % increment for the next line of text
                one_line = fgetl(fileID);       % read a line of the file
            end
            fclose(fileID);                     % close the file so it can be deleted

            % Delete temporary files
            delete(strcat(namehash,'_train.dat'),strcat(namehash,'_test.dat'),...
                strcat('model',namehash),strcat('predictions',namehash));

            %% Transform returned output to users class labels
            templabels = zeros(size(obj.labels));              % Create temporary labels so transfrom is error free
            obj.labels(obj.n_labeled+1:end,:) = Z;             % Recombine the unlabeled instances into the obj.labels
            templabels(obj.labels<0) = obj.input_label_ids(1); % Restore the -1 labels to the class identifier originally presented by user
            templabels(obj.labels>0) = obj.input_label_ids(2); % Restore the +1 labels to the class identifier originally presented by user
            obj.labels = templabels;                           % temp labels is used to ensure the >< thresholds do not pick up classes that have already been changed
                        
        end
        
        %% Matlab Format to S3VM Format Converter %%
        function [fileName] = s3vm_format(obj,namehash,filetype)
            % Converts matlab data to a format for Joachims SVM-Light
            
            switch filetype
                case 'train'
                    fileName = [namehash '_train.dat'];
                    D = obj.data;
                    L = obj.labels;
                    
                case 'test'
                    fileName = [namehash '_test.dat'];
                    D = obj.data(obj.n_labeled+1:end,:);     % isolate unlabeled data for classification
                    L = obj.labels(obj.n_labeled+1:end,:);   % isolate unlabeled data labels 
            end
                        
            if exist(fileName,'file') ~= 0      % If this filename already exists
                delete(fileName);               % delete it
            end
            fileID = fopen(fileName,'a');       % Open a file with this name and allow it to be appended
            
            [instances, features] = size(D);
            
            %Build a single line entry for each instance using the format
            % classlabel 1:feature1 2:feature2 3:feature3 ...
            for i = 1:instances
                templine = num2str(L(i));
                for j = 1:features 
                    if D(i,j) ~= 0
                        templine = strcat(templine,[' ',num2str(j),':',num2str(D(i,j))]);
                    end        
                end
                fprintf(fileID,'%s\n',templine);        %write the line to the file
            end

            fclose('all');                              %close the file
        end
        
    end
    
    %% Label Propagation/Spreading & Dependencies %%
    % For Label Spreading obj.classifier_opts use the following call
    % obj.set_classifier('label_spread',struct('alpha',#))
    methods (Access = protected)
        %% Label Propagation Classifier %%
        function label_prop(obj)
            % Label propagation method based on a 2002 paper by Xiaojin Zhu
            % and Zoubin Ghahramani - "Learning from Labeled and Unlabeled 
            % Data with Label Propagation".
            %
            % FUNCTION CALL:
            % label_prop(obj) - uses obj.data and obj.labels to classify
            %   unlabeled data using label propagation 
            %
            % Created 09.05.2012 by Karl Dyer
            
            %% Prepare Data for Label Propagation
            % Ensure labels are in a vector format for proper use with label_prop            
            if strcmp(obj.label_format, 'integer')      % If the format is integer
                obj.int2vec;                            % convert it to vector
            end
                        
            %% Compute Sigma
            % Sigma is computed by using Kruskal's Algorithm to find the minimum
            % spanning tree of all points. The first tree edge a that connects two
            % components with different labeled points in them is d0. Sigma = d0/3 tp
            % follow the 3-sigma rule of normal distributions
            if isempty(obj.classifier_opts.sigma)       % if the user hasn't specified a sigma value
                sig = obj.kruskal_mst;                % let the classifier determine it
            else
                sig = obj.classifier_opts.sigma;
            end

            %% Compute weight matrix
            % Weight matrix measures euclidean distance between every node, the closer
            % the nodes are the larger the weight assigned to that pair. Weights are
            % controlled by sigma parameter.

            % Gaussian Distance Metric
            W = zeros(obj.n_instances); %Allocate space for Gaussian Weight Matrix

            % Compute upper triangle of matrix
            for i = 1:obj.n_instances 
                for j = i+1:obj.n_instances
                    % Calculates Squared Euclidean Distance over sigma squared
                    W(i,j) = exp(-1*sum((obj.data(i,:)-obj.data(j,:)).^2)/(sig^2));  
                end
            end
            W = W + W'; % Mirror top triangle to lower triangle

            %% Construct Transition Matrix
            % The probabilty to jump from node j to i : T_ij = P(j->i)
            T = W./repmat(sum(W),obj.n_instances,1); %Column Normalize W
            T(isnan(T)) = 0;             % Convert any NaNs to zero elements

            % Normalize the transisiton matrix and truncate small elements (1e-5)
            Tn = pinv(diag(sum(T,2)))*T; % Row Normalize T
            Tn(Tn<1e-5)=0;               % Removing small elements increases speed

            % We refer to differenct segements of the matrices as follows:
            % Tn = [Tll Tlu     Y = [Yl
            %       Tul Tuu]         Yu]
            %
            % Break the normalized transisiton matrix into its components
            Tul = Tn(obj.n_labeled+1:end, 1:obj.n_labeled);     %Grab lower left corner of normalized transition matrix
            Tuu = Tn(obj.n_labeled+1:end, obj.n_labeled+1:end);   %Grab lower right corner of normalized transition matrix

            %% Label Propagation
            % CONVERGENCE THEOREM METHOD:
            % Refernce paper for entire proof
            % YU_pr = (eye(size(Tuu))-Tuu)\Tul*YL; % prob. of belonging to each class
            Z = pinv((eye(size(Tuu))-Tuu))*Tul*obj.labels(1:obj.n_labeled,:);
            obj.labels(obj.n_labeled+1:end,:) = Z;
            
        end
        
        %% Label Spreading Classifier %%
        function label_spread(obj)
            %% Prepare Data for Label Spreading
            % Ensure labels are in a vector format for proper use with label_prop            
            if strcmp(obj.label_format, 'integer')      % If the format is integer 
                obj.int2vec;                            % convert it to vector
            end
            
            %% Compute Sigma
            % Sigma is computed by using Kruskal's Algorithm to find the minimum
            % spanning tree of all points. The first tree edge a that connects two
            % components with different labeled points in them is d0. Sigma = d0/3 tp
            % follow the 3-sigma rule of normal distributions
            if isempty(obj.classifier_opts.sigma)       % if the user hasn't specified a sigma value
                sigma = obj.kruskal_mst;                % let the classifier determine it
            end

            %% Compute weight matrix
            % Weight matrix measures euclidean distance between every node, the closer
            % the nodes are the larger the weight assigned to that pair. Weights are
            % controlled by sigma parameter above.

            % Gaussian Distance Metric
            W = zeros(obj.n_instances); %Allocate space for Gaussian Weight Matrix

            % Compute upper triangle of matrix
            for i = 1:obj.n_instances 
                for j = i+1:obj.n_instances
                    % Calculates Squared Euclidean Distance over sigma squared
                    W(i,j) = exp(-1*(norm(obj.data(i,:)-obj.data(j,:))^2)/(2*sigma^2));  
                end
            end
            W = W + W'; % Mirror top triangle to lower triangle

            %% Compute Laplacian Matrix
            D = diag(sum(W,2));
            S = D^-0.5*W*D^-0.5; % Laplacian Matrix

            %% Label Propagation
            % CONVERGENCE THEOREM METHOD:
            % Refernce paper for entire proof
            obj.labels = (1-obj.classifier_opts.alpha)*(eye(obj.n_instances)-obj.classifier_opts.alpha*S)\obj.labels; % prob. of belonging to each class
            
        end
        
        %% Kruskal Walis Minimum Spanning Tree - Modified %%
        function sigma = kruskal_mst(obj)
            
            % Create temporary integer labels that can be modified without changing the object labels property
            [~,L] = max(obj.labels,[],2);
            L(sum(obj.labels,2)==0)=0;
            
            % Calculate Edge Weights
            E = dist(obj.data,obj.data');       % Calculate distance between each instance and every other instance
            E = E + tril(Inf*ones(size(E,1)));  % Make the lower triangle and diagonal equal to infinity to make search easier

            % Initiate Loop
            c = 1;                              
            tc = 0;                             % Set branch count index to zero
            d0 = Inf;                           % Initialize distance between oposing classes to Infinity
            T = zeros(size(E,1),1);             % Create space in memory for the Tree matrix which tracks which branch the instance belongs to

            while c < size(obj.data,1) && d0 == Inf;   % While there are still instances to be analyzed and two opposing classes haven't been linked
                %% Find Shortest Edge
                [temp,v2] = min(E,[],2);         %Find shortest edge in each row - v1 indicates what column it was found in
                [~,v1] = min(temp,[],1);         %Find shortest edge globally by searching temp - v2 indicates what row is was found in
                v2 = v2(v1);                     %Get the index of the column that the shortest row was found in
                % At this point v1 and v2 indicate which instances or vertices are connected at a minimum span
                
                if (T(v1)~=0 && T(v2)~=0) && T(v1)==T(v2)   % They are part of same tree already so make no connection
                    E(v1,v2)=Inf;                           % Remove that edge from analysis
                    c = c - 1;                              % Undo the effect of the counter increment

                elseif T(v1)==0 && T(v2)==0                 % Neither points are in tree, assign to new branch
                    tc = tc+1;                              % Generate new branch index
                    T([v1 v2],1) = tc;                      % Assign both instances to new branch
                    [L,E,d0] = label_rem(T,L,E,v1,v2);  % Check for possibility to merge branches

                elseif T(v1)~=0 && T(v2)==0                 % If 1st point is in tree, assign 2nd to it
                    T(v2,1) = T(v1,1);                      % Propagate the branch labeled from node 2 to node 1
                    [L,E,d0] = label_rem(T,L,E,v1,v2);  % Check for possibility to merge branches

                elseif T(v1)==0 && T(v2)~=0                 % If 2nd point is in tree, assign 1st to it
                    T(v1,1) = T(v2,1);                      % Propagate the branch labeled from node 2 to node 1
                    [L,E,d0] = label_rem(T,L,E,v1,v2);  % Check for possibility to merge branches

                elseif T(v1)~=0 && T(v2)~=0 && T(v1)~=T(v2) % If node are from opposing classes
                    high = max(T(v1),T(v2));
                    low = min(T(v1),T(v2));
                    T(T==high)=low;
                    [L,E,d0] = label_rem(T,L,E,v1,v2);
                end
                
                c = c + 1;              % Increment the instance counter
            end

            sigma = d0/3;

            % Internal function to Kruskal MST
            function [L,E,d0] = label_rem(T,L,E,v1,v2)

                if L(v1)~= 0 && L(v2)~=0 && L(v1)~=L(v2)
                    d0 = E(v1,v2);
                elseif L(v1)~=0 && L(v2)==0
                    L(v2) = L(v1);
                    d0 = Inf;
                elseif L(v1)==0 && L(v2)~=0
                    L(v1) = L(v2);
                    d0 = Inf;
                else
                    d0 = Inf;
                end

                % Propogate tree label to other nodes in joined tree
                if L(v1)~=0;
                    L(T==T(v1)) = L(v1);
                end;

                % Remove node from analysis
                E(v1,v2)=Inf;
            end
        end    
        
    end
    
    %% Cluster and Label & Dependencies
    % For Label Spreading obj.classifier_opts use the following call
    % obj.set_classifier('cluster_n_label',struct('k',#))
    methods (Access = protected) 
        %% Cluster and Label Classifier %%
        function cluster_n_label(obj)
            
            %% Prepare Data for S3VM
            % Ensure labels are in a integer format for proper use with cluster and label            
            if strcmp(obj.label_format, 'vector')       % If the input format is vector 
                obj.vec2int;                            % convert it to integer
            end
            
            Z = zeros(size(obj.data,1),1);              % Initialize an empty label vector
            missed = 0;                                 % Initialize missed to have any value to enter loop

            %% Cluster all data and make sure each cluster contains as least 1 labeled point
            errCt = 0;
            while ~isempty(missed)                                      % If there are clusters that do not have labeled data 
                [IDX,C] = kmeans(obj.data,obj.classifier_opts.k);            % Generate a set of k clusters
                labeled_clusters = unique(IDX(1:obj.n_labeled));        % Find which clusters have labeled data in them
                missed = setdiff(1:obj.classifier_opts.k,labeled_clusters);  % Determine which clusters are missing labeled data
                
                errCt = errCt + 1;                                                  % Increase error count
                if errCt == 100 && obj.classifier_opts.k > obj.classifier_opts.kmin % If unable to find a suitable model and haven't hit minimum clusters
                    obj.classifier_opts.k = obj.classifier_opts.k-1;                % reduce the number of clusters
                    display(['k reduced to ', num2str(obj.classifier_opts.k)]);     % notify the user
                    errCt = 0;                                                      % reset the error count
                elseif errCt == 100 && obj.classifier_opts.k == obj.classifier_opts.kmin % If the minium clusters has been reached
                    display('Exiting Cluster and Label')                            % notify the user
                    return;                                                         % abort the cluster and label
                end
            end

            YL = obj.labels(1:obj.n_labeled);
            %Use simple majority vote to classify clusters
            for i = 1:obj.classifier_opts.k
                freq = hist(YL(IDX(1:obj.n_labeled)==i),obj.input_label_ids);        % Tally up vote
                [val,c] = max(freq,[],2);                                            % Choose largest
                if sum(freq==val)>1                                                  % If there is a tie
                    d(:,i) = dist(obj.data(IDX(1:obj.n_labeled)==i,:),C(i,:)');      % Measure each instance to the center of the cluster
                    d = sum(d,2);
                    [~,col] = min(d,[],1);                                             % Choose the class that is closest to the center
                    cc = obj.labels(IDX(1:obj.n_labeled)==i,:);
                    c = cc(col,1);
                end
                % handle case when only one class exists in a 2 class prob
                if numel(obj.input_label_ids) == 1
                    Z(IDX == i) = obj.input_label_ids(1);                                % Assign rest of cluster labels
                else
                    Z(IDX == i) = obj.input_label_ids(c);                                % Assign rest of cluster labels
                end
            end
            obj.labels = Z;
        end
        
    end
       
    %% CLASS BALANCE METHODS %%
    methods (Access = protected)
        %% Class Mass Balance %%
        function none(obj)
            obj.vec2int                                     % Get rid of posterior probabilites and assign to a class
        end
                
        function mass(obj)
            class_mass = sum(obj.labels);                                           % total "mass" or probability of each class class 
            balance_factor = obj.balance_opts.class_priors./class_mass;             % scaling factor for unlabeled instances
            obj.labels = obj.labels.*repmat(balance_factor,size(obj.labels,1),1);   % rebalance labels
        end
        
        %% Label Bidding Balance %%
        function bid(obj)
            available_each_class = ceil(obj.balance_opts.class_priors*obj.n_instances); % determine how many labels should be available for each class 
            
            % Could improve by using round to checks which class is over represented
            while sum(available_each_class,2) ~= obj.n_instances           % if available labels does not match number of instances
                random_class = randi(obj.n_classes,1,1);                   % randomly choose one of the classes
                available_each_class(random_class) = available_each_class(random_class) - 1; % remove one label from it
            end

            % Start dealing out labels
            bidded_labels = zeros(obj.n_instances,obj.n_classes);          % initialize matrix for labels after bidding
            
            while sum(available_each_class,2) ~= 0                         % while there are still labels to be bid out
                [temp,r] = max(obj.labels,[],1);                           % find row of the highest bidder (i.e. the highest posterior probability)
                [~,c] = max(temp,[],2);                                    % find column or what class is bidding for
                if available_each_class(1,c) > 0                           % if there are labels still available for that class
                    bidded_labels(r(c),c) = 1;                             % assign the instances to it requested class
                    obj.labels(r(c),:) = -inf;                             % remove that instance from the bidding process
                    available_each_class(1,c) = available_each_class(1,c) - 1; % decrease available labels for that class by one
                else                                                       % OTHERWISE that class has been bidded out so
                    obj.labels(r(c),c) = -Inf;                             % remove the instances probability to belong to that class
                end
            end
            
            obj.labels = bidded_labels;                                    % store the bidded labels into the object
        end
        
        %% Regularization Balance %%
%         function reg(obj)
%         end
%         
    end
end