% This class includes functionality for extracting core supports for the
% COMPOSE algorithm. It is designed to be a property of an instance of the
% compose class defined in compose.m for more information on core support
% extraction, see :
%
% "COMPOSE: A Semi-Supervised Learning Framework for Initially Labeled Non-Stationary Streaming Data" IEEE Transactions on Neural Networks and Learning Systems, Special issue on Learning in Nonstationary and Dynamic Environments
% - and -
% "Core Support Extraction for Learning from Initially Labeled Nonstationary Environments using COMPOSE," World Conference on Computational Intelligence (WCCI 2014)
%
% Written by Karl Dyer and Rob Capo
classdef cse < handle
    %% Class Properties - Visible %%
    properties (SetAccess = protected)
        synthetic_data = 0  % [INTEGER] True/False whether to allow synthetic data
        verbose        = 1  % [INTEGER] Level of feedback displayed during run {default}
                            %    0  : No Information Displayed
                            %   {1} : Command line progress updates
                            %    2  : Plots when possible and Command line progress updates
                            
        data                % [MATRIX] N instances x D dimensions : Features of data 
        
        n_cores        = 1  % [INTEGER] : The number of cores to use for parallel processing segments of code
        
        boundary                % [STRING] : Type of boundary construction to use 
        boundary_opts           % [STRUCT] : Options that correspond to with boundary construction function selected
        boundary_data = struct  % [STRUCT] : Data recorded by the boundary for API usage
        
        user_data = struct;     % [STRUCT] : Extra data can be set during extraction and will be available to user after
        
        parent
        
    end
    
    %% Class Properties - Hidden %%
    properties
        n_instances         % [INTEGER] : Number of instances in data
        n_features          % [INTEGER] : Number of features in data (i.e. dimensionality of data)
    end
    
    %% Class Properties - Hidden & Constant %%
    properties (Hidden, Constant)
        % The cells below contain text strings that match the boundary
        % construction methods available in this object
        % If if other boundary construction methods are added to this class 
        % these cells must be modified to include those methods
        valid_boundary = {'a_shape';'gmm'}
    end
    
    %% Initialization and Executing Methods %%
    methods
        %% CSE Constructor %%
        function  obj = cse(verbose, synthetic_data, parent)
            error(nargchk(1,3,nargin));                            % Verify correct number of arguments are passed
            
            if nargin == 1                                                                  % If user specifies verbosity
                validateattributes(verbose,{'numeric'},{'scalar','integer','>=',0,'<=',2})  % Validate verbose and
                obj.verbose = verbose;                                                      % Load verbosity level into object
            end
            
            if nargin >= 2                                                                         % If user specifies synthetic data use
                validateattributes(synthetic_data,{'numeric'},{'scalar','integer','>=',0,'<=',2})  % Validate synthetic data and
                obj.synthetic_data = synthetic_data;                                               % Load use of synthetic data into object
            end
            
            if nargin == 3
                obj.parent = parent;
            end
        end
        
        %% Set verbosity level
        function obj = set_verbose(verbose)
            if ~isnumeric(verbose), return; end
            if verbose > 2, verbose = 2; end
            if verbose < 0, verbose = 0; end
            
            obj.verbose = floor(verbose);
        end
        
        %% Load Data %%
        function set_data(obj,data)
            error(nargchk(2,2,nargin));                                 % Verify correct number of arguments are passed
            validateattributes(data,{'numeric'},{'real','nonnan'})      % Validate data
            
            obj.data = data;                                            % load the data into the object
            [obj.n_instances, obj.n_features] = size(data);             % load dimensionality of data into object
        end
                
        %% Set Boundary Construction Type & Options %%
        function set_boundary(obj, boundary_selection, opts)
            error(nargchk(2,3,nargin));                                     % Verify correct number of arguments are passed
            validateattributes(boundary_selection,{'char'},{'nonempty'});   % validate users classifer selection
            if nargin == 3                                                  % if opts parameter is passed
                validateattributes(opts,{'struct'},{'nonempty'});           % validate opts 
            else                                                            % if no options are passed
                opts = [];                                                  % make options a blank variable
            end
            
            obj.boundary_opts = [];                                       % ensure options are cleared from any previous classifiers
            
            if sum(strcmp(boundary_selection,obj.valid_boundary)) == 1      % Check user requested classifier against validation list
                obj.boundary = str2func(boundary_selection);              % load classifier selection into object and convert to function handle
                obj.set_default_opts;                                       % load the default options for specified classifier
                if ~isempty(opts)                                           % if the user passed options
                    obj.set_user_opts(opts);                                % overwrite any defaults with user options that validate
                end
            else
                display([boundary_selection,...
                    ' is an invalid boundary construction method - choose from: ']); % if invalid selection notify user
                display(obj.valid_boundary);                               % and provide acceptable selections
            end
        end
       
        %% Extract Core Supports Using Boundary Type Selected %%
        function inds = extract(obj)
            if isempty(obj.data)                              % If data is not present
                error('You must load data before extracting core supports') % throw error
            end
            
            if isempty(obj.boundary)                          % Determine if classifier has been set, if not
                display('Boundary construction type not set - default classifier and options loaded') % issue warning that a default is being used and
                obj.set_boundary('gmm');                      % load cluster and label as the default classifier
            end;
            
            if obj.verbose == 2         % If plotting is requested
                hold on;
                obj.cse_plot([])            % plot the labeled and unlabeled data
                drawnow
            end
            
            inds = obj.boundary(obj);   % Run the boundary construction and extract indices of core supporting instances
            
            
            if obj.verbose == 2 && inds ~= -1         % If plotting is requested
                obj.cse_plot(inds)      % plot the labeled and unlabeled data
                hold off;
            end
        end
    end
    
    %% Protected General Methods %%
    methods (Access = protected)
        %% Set Default Options %%
        function set_default_opts(obj)
            switch func2str(obj.boundary)
                case 'a_shape'
                    obj.boundary_opts.alpha = 2;       % set alpha parameter of alpha shape function 
                    obj.boundary_opts.p = 0.20;        % set the percentage of points to be used for core supports
                case 'gmm'
                    obj.boundary_opts.k = 10;          % set the number of centers to find
                    obj.boundary_opts.p = 0.4;         % set percentage of points to be used for core supports
            end
        end
       
        %% Set User Options %%
        function set_user_opts(obj,opts)
            user_fields = fieldnames(opts);     % make list of struct fields passed by user
            for field = 1:numel(user_fields)    % iterate thru each struct field passed
                if isfield(obj.boundary_opts,user_fields(field))  % If the field passed by user is a valid field
                    obj.boundary_opts.(user_fields{field}) = opts.(user_fields{field}); % update the option to the users selection
                else                                                                      % OTHERWISE
                    disp(['Warning: Option ', user_fields{field},...                      % inform user the field is not a valid parameter
                        ' is not a valid option for ', func2str(obj.boundary),...
                        ' boundary construction method'])                   
                end
            end
        end
        
        %% Special Plotting For CSE Class %%
        function cse_plot(obj, ind)
            if isempty(ind)                     % If there are no indices specified
                ind = 1:size(obj.data,1);       % Set the indices equal to all the data
                color = '.r';                   % set the plot color to red dots 
            else                                % OTHERWISE we are plotting the core supports so 
                color = '+k';                   % set the plot color to black dots
            end
            
            switch obj.n_features                % send to correct plotting calls based on dimensionality
                case 1
                case 2
                    disp(obj.data(ind, 1));
                    plot(obj.data(ind,1),obj.data(ind,2),color);            % Plot unlabeld data
                    xlabel('Feature 1')                                        % label x axis 
                    ylabel('Feature 2')                                        % label y axis 
                    title(['Boundary Constructor: ', func2str(obj.boundary)]); % Add a title
                case 3
                    scatter3(obj.data(ind,1),obj.data(ind,2),obj.data(ind,2),color); % Plot unlabeld data
                    xlabel('Feature 1')                                              % label x axis 
                    ylabel('Feature 2')                                              % label y axis 
                    zlabel('Feature 3')                                              % label z axis
                    title(['Boundary Constructor: ', func2str(obj.boundary)]);       % Add a title
            end
        end
    end
    
    %% Alpha Shape & DEPENDENCIES %%
    methods (Access = protected)
        %% Alpha Shape Core Support Extraction %%
        function support_inds = a_shape(obj)
            
            ashape = obj.alpha_shape;               % Construct Alpha Shape

            if isempty(ashape)                      % If the returned alpha shape is empty
                display('No Alpha Shape could be constructed try different alpha or check data');   % notify user
                support_inds = -1;
                return;                             % return to calling function
            end
            
            if isa(obj.parent, 'compose')
                if isa(obj.parent.plot_callback, 'function_handle')
                    plot_data = obj.data;
                    plot_labels = unique(ashape.simplexes(ashape.include == 1,:), 'rows');

                    obj.parent.step = 2.1;
                    obj.parent.plot_callback(plot_data, plot_labels, obj.parent);
                end
            end
            
            %%%%% ADDED ONION METHOD %%%%%   
            ashape.N_start_instances = size(obj.data, 1);                           % Number of data points to start
            ashape.N_core_supports = ceil(size(obj.data,1)*obj.boundary_opts.p);    % Number of desired core supports based on percentage parameter
            ashape.core_support = ones(ashape.N_start_instances,1);                 % Binary vector indicating instance is core support or not
            ashape.layer = zeros(size(obj.data,1), 1);
            
            too_many_core_supports = true;         % Flag denoting if the target number of coresupports has been obtained
            
            lc = 1;                              % Loop counter for plotting layers of alpha shape compaction            

%             hold on;
%             cmap = colormap(hsv);
%             for cc = 1:size(ashape.simplexes,1)
%                 if ashape.include(cc) == 1;
%                     fill(obj.data(ashape.simplexes(cc,:),1),obj.data(ashape.simplexes(cc,:),2),'c','FaceAlpha',0.3,'FaceColor',cmap(lc,:),'EdgeColor',[0 0 0]);
%                 end
%             end

            %% Start removing layers
            while sum(ashape.core_support) >= ashape.N_core_supports && too_many_core_supports == true
                %% Find edge d-1 simplexes
                Tid = repmat(find(ashape.include==1),size(ashape.simplexes,2),1);           %Identifiy which simplex each d-1 simplex comes from
                
                edges = [];
                nums = 1:size(ashape.simplexes,2);
                for ic = 1:size(ashape.simplexes,2)
                    edges = [edges; ashape.simplexes(ashape.include==1,nums(1:(size(ashape.simplexes,2)-1)))];
                    nums = circshift(nums,[0 1]);    
                end
                
                edges = sort(edges,2);                      %sort the d-1 simplexes so small node is on left in each row
                [edges, Sid] = sortrows(edges);             %sort by rows placing copies of d-1 simplexes in adjacent rows
                Tid = Tid(Sid);                             %sort the simplex identifiers to match
                
                consec_edges = sum(diff(edges),2);          %find which d-1 simplexes are duplicates - a zero in row N indicates row N and N+1 are the same
                consec_edges(find(consec_edges==0)+1) = 0;  %throw a zero mark on the subsequent row (N+1) as well
                
                ashape.include(Tid(consec_edges~=0)) = 0;
                ashape.layer(Tid(consec_edges~= 0)) = lc;
                
                 %% Plot current shape
%                  clf;
%                  hold on;
%                 for cc = 1:size(ashape.simplexes,1)
%                     if ashape.include(cc) == 1;
%                         fill(obj.data(ashape.simplexes(cc,:),1),obj.data(ashape.simplexes(cc,:),2),'c','FaceAlpha',0.3,'FaceColor',[.9 0 0],'EdgeColor',[0 0 0]);
%                     end
%                 end
%                  hold off;
%                  xlim([-6 6]);
%                  ylim([-6 6]);
%                  getframe;
%                  pause

%                 plot(D(:,1),D(:,2),'d','markersize',5,'markerfacecolor',cmap(lc+1,:));
                lc = lc + 1;

                %%
                points_remaining = unique(ashape.simplexes(ashape.include==1));
                if numel(points_remaining) >= ashape.N_core_supports 
                   ashape.core_support(setdiff(1:ashape.N_start_instances,points_remaining),:) = 0;
                else
                   too_many_core_supports = false;
                end
            end
            
            ashape.layer(ashape.include == 1) = lc;
            %%%%% END ONION METHOD %%%%%
            
            if isa(obj.parent, 'compose')
                if isa(obj.parent.plot_callback, 'function_handle')
                    plot_labels = unique(ashape.simplexes(ashape.core_support == 1, :), 'rows');

                    obj.parent.step = 2.2;
                    obj.parent.plot_callback(plot_data, plot_labels, obj.parent);
                end
            end
            
%             D = obj.data(unique(ashape.simplexes(ashape.include == 1)), :);
%             figure;
%             plot(D(:, 1), D(:, 2), '+b');
            
            support_inds = unique(ashape.simplexes(ashape.include == 1));
%             plot(obj.data(support_inds,1),obj.data(support_inds,2),'ko','markerfacecolor','k','markersize',6);
            obj.user_data.ashape = ashape;
        end
               
        %% Alpha Shape Boundary Construction %%
        function [ashape] = alpha_shape(obj)
            validateattributes(obj.boundary_opts.alpha,{'numeric'},{'scalar','>',0})   % Validate alpha paramter is in allowable range
            
            obj.set_data(unique(obj.data,'rows'));    % Remove duplicate instances
            ashape = struct;                          % Allocate memory for structured output

            if obj.n_instances < obj.n_features+1                                      % If a class does not have enought points to construct a tesselation
                fprintf(['Warning::Alpha_Shape::Tesselation_Construction\n',...        % inform the user
                    'Data of dimension %i requires a minimum of %i unique points.',...
                    ' Alpha shape was not constructed for this data\n'],obj.n_features,obj.n_features+1);
                ashape = [];                                                           % set the output to any empty variable
                return;                                                                % and return to calling function
            else                                                                       % OTHERWISE there are enough instance for the number of dimensions
                ashape.simplexes = delaunayn(obj.data,{'Qt','Qbb','Qc','Qz'});         % set the output simplexes to the Delaunay Triangulation
                ashape.include   = zeros(size(ashape.simplexes, 1), 1);                % create a blank vector to indicate whether to include each simplex in the alpha shape
                
                for sID = 1:size(ashape.simplexes,1)                                                % iterate through each simplex     
                    if obj.boundary_opts.alpha > calc_radius(obj.data(ashape.simplexes(sID, :), :)) % If alpha is larger than the radius of the simplex
                        ashape.include(sID) = 1;                                                    % Include that simplex in the alpha shape
                    end
                end
            end
            
            % PLOT OPTIONS
            if obj.verbose == 2
                hold on;
                % Shade d-simplexes inside alpha shape
                switch obj.n_features
                    case 2
                        trimesh(alpha.simplexes(ashape.include == 1),obj.data,'facecolor','k','facealpha',0.05','edgecolor','none','edgealpha',0.5);   
                    case 3
                        tetramesh(alpha.simplexes(ashape.include == 1),obj.data,'facecolor','k','facealpha',0.05','edgecolor','none','edgealpha',0.5);   
                end
                drawnow
            end
            
            % Alpha Shape Boundary Construction Dependency
            function [radius,c] = calc_radius(points)
                %CALC_RADIUS - finds the radius of an D-dimensional sphere from D+1 points
                %           
                % INPUTS
                %   points : [matrix] D+1 points x D dimensions
                %
                % OUPUTS
                %   radius : [double] radius of D-sphere
                %
                % Reference - http://mysite.verizon.net/res148h4j/zenosamples/zs_sphere4pts.html

                %%% Validation and Parser
                nD = size(points,2); % number of dimensions
                validateattributes(nD,{'numeric'},{'>=',2});

                p = inputParser;
                p.addRequired('points',@(x)validateattributes(x,{'numeric'},{'real','nonnan','finite','size',[nD+1,nD]}));

                p.parse(points)
                in = p.Results;         %Struct of inputs ex: in.variable

                %%% Construct Matrix
                M(:,1) = sum(in.points.^2,2);
                M = [M in.points ones(nD+1,1)];

                %%% Calculate minors
                m = zeros(size(M,2),1);
                for mID = 1:size(M,2)
                    temp = M;
                    temp(:,mID) = [];
                    m(mID) = det(temp);
                end

                %%% Calculate center for each dimension
                c = zeros(1,nD);
                for j = 1:nD
                    c(1,j)=((-1)^(j+1))*0.5*m(j+1)/m(1);
                end

                %%% Measure distance from center point to any first given point
                radius = sqrt(sum((c-in.points(1,:)).^2,2));
            end
        end
        
    end
    
    %% GMM & Dependencies %%
    methods (Access = protected)
        %% GMM Mahalbonois Transform Core Support Extraction %%
        function support_inds = gmm(obj)
                       
            [~,M,V,~] = obj.EM_GM;
            M = M';
            
%             inds = zeros(obj.n_instances, 1);
            
            % Create matrix to hold mahal distances of every point to every
            % distribution
            dists = zeros(obj.n_instances, obj.boundary_opts.k);
            
            for i = 1:obj.boundary_opts.k
                MU = M(i, :);
                SIG = V(:, :, i);

                new_data = obj.mahal_transform(MU, SIG);

                for j = 1:size(new_data, 1)
                    dists(j, i) = sqrt(sum(new_data(j, :).^2, 2));

%                     if d < obj.boundary_opts.threshold;
%                         inds(j) = 1;
%                     end
                end
            end
            
            % Find the smallest distance of each point to any of the
            % distributions
            min_dists = min(dists, [], 2);
            
            [b ix] = sort(min_dists);
            
            support_inds = ix(1:round(obj.boundary_opts.p * obj.n_instances));
        end
        
        % GMM Mahal CSE Dependency %
        function z = mahal_transform(obj, mu, s)
            if nargin < 3, s = cov(obj.data); end
            if nargin < 2, mu = mean(obj.data); end

            z = s^(-.5)*(obj.data-repmat(mu, obj.n_instances, 1))';
            z = z';
        end
        
        % GMM Mahal CSE Dependency %
        function D = mahal_dist(obj, Y, mu, SIGMA)
            D = (Y-mu)/SIGMA*(Y-mu)';
        end
        
        %% EM GM %%
        function [W,M,V,L] = EM_GM(obj,ltol,maxiter,Init)
            % [W,M,V,L] = EM_GM(X,k,ltol,maxiter,pflag,Init) 
            % 
            % EM algorithm for k multidimensional Gaussian mixture estimation
            %
            % Inputs:
            %   X(n,d) - input data, n=number of observations, d=dimension of variable
            %   k - maximum number of Gaussian components allowed
            %   ltol - percentage of the log likelihood difference between 2 iterations ([] for none)
            %   maxiter - maximum number of iteration allowed ([] for none)
            %   Init - structure of initial W, M, V: Init.W, Init.M, Init.V ([] for none)
            %
            % Outputs:
            %   W(1,k) - estimated weights of GM
            %   M(d,k) - estimated mean vectors of GM
            %   V(d,d,k) - estimated covariance matrices of GM
            %   L - log likelihood of estimates
            %
            % Written by
            %   Patrick P. C. Tsui,
            %   PAMI research group
            %   Department of Electrical and Computer Engineering
            %   University of Waterloo, 
            %   March, 2006
            %
            % Edited by Rob Capo, August 2012

            %%%% Validate inputs %%%%
            if nargin < 6, Init = []; end
            if nargin < 5 || isempty(maxiter), maxiter = 1000; end
            if nargin < 4 || isempty(ltol), ltol = 0.1; end
            if nargin < 1
                disp('EM_GM needs at least data and a k value input');
                return;
            end
            
            if obj.Verify_X || obj.Verify_K || obj.Verify_ltol(ltol) || obj.Verify_maxiter(maxiter) || obj.Verify_Init(Init)
                return;
            end

            %%%% Initialize W, M, V,L %%%%
            if isempty(Init),  
                [W,M,V] = obj.Init_EM;    
            else
                W = Init.W;
                M = Init.M;
                V = Init.V;
            end
            Ln = obj.Likelihood(W,M,V); % Initialize log likelihood
            Lo = 2*Ln;

            %%%% EM algorithm %%%%
            niter = 0;
            while (abs(100*(Ln-Lo)/Lo)>ltol) && (niter<=maxiter),
                E = obj.Expectation(W,M,V); % E-step    
                [W,M,V] = obj.Maximization(E);  % M-step
                Lo = Ln;
                Ln = obj.Likelihood(W,M,V);
                niter = niter + 1;
            end 
            L = Ln;
        end
        
        % EM GM Dependency %
        function E = Expectation(obj,W,M,V)
            a = (2*pi)^(0.5*obj.n_features);
            S = zeros(1,obj.boundary_opts.k);
            iV = zeros(obj.n_features,obj.n_features,obj.boundary_opts.k);
            for j=1:obj.boundary_opts.k,
                if V(:,:,j)==zeros(obj.n_features,obj.n_features) 
                    V(:,:,j)=ones(obj.n_features,obj.n_features)*eps;
                end
                S(j) = sqrt(det(V(:,:,j)));
                iV(:,:,j) = inv(V(:,:,j));    
            end
            E = zeros(obj.n_instances,obj.boundary_opts.k);
            for i=1:obj.n_instances,    
                for j=1:obj.boundary_opts.k,
                    dXM = obj.data(i,:)'-M(:,j);
                    pl = exp(-0.5*dXM'*iV(:,:,j)*dXM)/(a*S(j));
                    E(i,j) = W(j)*pl;
                end
                E(i,:) = E(i,:)/sum(E(i,:));
            end
        end

        % EM GM Dependency %
        function [W,M,V] = Maximization(obj,E)
            W = zeros(1,obj.n_features, obj.boundary_opts.k); M = zeros(obj.n_features,obj.boundary_opts.k);
            V = zeros(obj.n_features, obj.n_features, obj.boundary_opts.k);
            for i=1:obj.boundary_opts.k,  % Compute weights
                for j=1:obj.n_instances,
                    W(i) = W(i) + E(j,i);
                    M(:,i) = M(:,i) + E(j,i)*obj.data(j,:)';
                end
                M(:,i) = M(:,i)/W(i);
            end
            for i=1:obj.boundary_opts.k,
                for j=1:obj.n_instances,
                    dXM = obj.data(j,:)'-M(:,i);
                    V(:,:,i) = V(:,:,i) + E(j,i)*(dXM*dXM');
                end
                V(:,:,i) = V(:,:,i)/W(i);
            end
            W = W/obj.n_instances;
        end

        % EM GM Dependency %
        function L = Likelihood(obj,W,M,V)
            % Compute L based on K. V. Mardia, "Multivariate Analysis", Academic Press, 1979, PP. 96-97
            % to enchance computational speed
            U = mean(obj.data)';
            S = cov(obj.data);
            L = 0;
            for i=1:obj.boundary_opts.k,               
                L = L + W(i)*(-0.5*obj.n_instances*log(det(2*pi*V(:,:,i))) ...
                    -0.5*(obj.n_instances-1)*(trace(V(:,:,i)\S)+(U-M(:,i))'/V(:,:,i)*(U-M(:,i))));
            end
        end

        % EM GM Dependency %
        function err_X = Verify_X(obj)
            err_X = 1;
            if obj.n_instances<obj.n_features,
                disp('Input data must be n x d!/n');
                return
            end
            err_X = 0;
        end
        
        % EM GM Dependency %
        function err_k = Verify_K(obj)
            err_k = 1;
            if ~isnumeric(obj.boundary_opts.k) || ~isreal(obj.boundary_opts.k) || obj.boundary_opts.k<1,
                disp('k must be a real integer >= 1!/n');
                return
            end
            err_k = 0;
        end

        % EM GM Dependency %
        function err_ltol = Verify_ltol(obj, ltol)
            err_ltol = 1;
            
            if isempty(ltol) || ~isreal(ltol) || ltol<=0,
                disp('ltol must be a positive real number!');
                return;
            end
            err_ltol = 0;
        end

        % EM GM Dependency %
        function err_maxiter = Verify_maxiter(obj, maxiter)
            err_maxiter = 1;
            if ~isreal(maxiter) || maxiter<=0,
                disp('ltol must be a positive real number!');
                return
            end
            err_maxiter = 0;
        end

        % EM GM Dependency %
        function err_Init = Verify_Init(obj, Init)
            err_Init = 1;
            if isempty(Init)
                % Do nothing
            elseif isstruct(Init),
                [~,Wk] = size(Init.W);
                [Md,Mk] = size(Init.M);
                [Vd1,Vd2,Vk] = size(Init.V);
                if Wk~=Mk || Wk~=Vk || Mk~=Vk,
                    disp('k in Init.W(1,k), Init.M(d,k) and Init.V(d,d,k) must equal!/n')
                    return
                end
                if Md~=Vd1 || Md~=Vd2 || Vd1~=Vd2,
                    disp('d in Init.W(1,k), Init.M(d,k) and Init.V(d,d,k) must equal!/n')
                    return
                end
            else
                disp('Init must be a structure: W(1,k), M(d,k), V(d,d,k) or []!');
                return
            end
            err_Init = 0;
        end

        % EM GM Dependency %
        function [W,M,V] = Init_EM(obj)
            [Ci,C] = kmeans(obj.data,obj.boundary_opts.k,'Start','cluster', ...
                'Maxiter',100, ...
                'EmptyAction','drop', ...
                'Display','off'); % Ci(nx1) - cluster indeices; C(k,d) - cluster centroid (i.e. mean)
            while sum(isnan(C))>0,
                [Ci,C] = kmeans(obj.data,obj.boundary_opts.k,'Start','cluster', ...
                    'Maxiter',100, ...
                    'EmptyAction','drop', ...
                    'Display','off');
            end
            M = C';
            Vp = repmat(struct('count',0,'X',zeros(obj.n_instances,obj.n_features)),1,obj.boundary_opts.k);
            for i=1:obj.n_instances, % Separate cluster points
                Vp(Ci(i)).count = Vp(Ci(i)).count + 1;
                Vp(Ci(i)).obj.data(Vp(Ci(i)).count,:) = obj.data(i,:);
            end
            V = zeros(obj.n_features,obj.n_features,obj.boundary_opts.k);
            W = zeros(1,obj.boundary_opts.k);
            for i=1:obj.boundary_opts.k,
                W(i) = Vp(i).count/obj.n_instances;
                V(:,:,i) = cov(Vp(i).obj.data(1:Vp(i).count,:));
            end
        end
    end
    
    %%% Parzen Window & Dependencies %%%
    methods (Access = protected)
        
    end
    
    %%% KNN & Dependencies %%%
    methods (Access = protected)
        
    end
end