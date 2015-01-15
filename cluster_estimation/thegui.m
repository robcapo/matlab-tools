function varargout = thegui(varargin)
% THEGUI MATLAB code for thegui.fig
%      THEGUI, by itself, creates a new THEGUI or raises the existing
%      singleton*.
%
%      H = THEGUI returns the handle to a new THEGUI or the handle to
%      the existing singleton*.
%
%      THEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THEGUI.M with the given input arguments.
%
%      THEGUI('Property','Value',...) creates a new THEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before thegui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to thegui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help thegui

% Last Modified by GUIDE v2.5 17-Dec-2013 22:16:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @thegui_OpeningFcn, ...
                   'gui_OutputFcn',  @thegui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before thegui is made visible.
function thegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to thegui (see VARARGIN)

% Choose default command line output for thegui
handles.output = hObject;


% Global to hold the handle to the status label. This will be set from
% other functions when clustering is complete
global status_h;
status_h = handles.lblStatus;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes thegui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = thegui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txtGmmK_Callback(hObject, eventdata, handles)
% hObject    handle to txtGmmK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGmmK as text
%        str2double(get(hObject,'String')) returns contents of txtGmmK as a double


% --- Executes during object creation, after setting all properties.
function txtGmmK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGmmK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnGenGmm.
function btnGenGmm_Callback(hObject, eventdata, handles)
% hObject    handle to btnGenGmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the user parameters to generate the dataset
mutype = get(handles.popGmmMus, 'String');
covtype = get(handles.popGmmMus, 'String');
priortype = get(handles.popGmmPriors, 'String');

% Generate the data
data = generate_gmm(str2double(get(handles.txtGmmK, 'String')), str2double(get(handles.txtGmmN, 'String')), mutype(get(handles.popGmmMus, 'Value')), covtype(get(handles.popGmmCovs, 'Value')), priortype(get(handles.popGmmPriors, 'Value')));

gd = guidata(hObject);
gd.data = data;

guidata(hObject, gd);

% Plot the data in the GUI
scatter(handles.ax, data(:, 1), data(:, 2));
grid on;


% --- Executes on selection change in popMethod.
function popMethod_Callback(hObject, eventdata, handles)
% hObject    handle to popMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popMethod


% --- Executes during object creation, after setting all properties.
function popMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMinSearch_Callback(hObject, eventdata, handles)
% hObject    handle to txtMinSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMinSearch as text
%        str2double(get(hObject,'String')) returns contents of txtMinSearch as a double


% --- Executes during object creation, after setting all properties.
function txtMinSearch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMinSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMaxSearch_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaxSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxSearch as text
%        str2double(get(hObject,'String')) returns contents of txtMaxSearch as a double


% --- Executes during object creation, after setting all properties.
function txtMaxSearch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaxSearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popGmmMus.
function popGmmMus_Callback(hObject, eventdata, handles)
% hObject    handle to popGmmMus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popGmmMus contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popGmmMus


% --- Executes during object creation, after setting all properties.
function popGmmMus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popGmmMus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popGmmCovs.
function popGmmCovs_Callback(hObject, eventdata, handles)
% hObject    handle to popGmmCovs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popGmmCovs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popGmmCovs


% --- Executes during object creation, after setting all properties.
function popGmmCovs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popGmmCovs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popGmmPriors.
function popGmmPriors_Callback(hObject, eventdata, handles)
% hObject    handle to popGmmPriors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popGmmPriors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popGmmPriors


% --- Executes during object creation, after setting all properties.
function popGmmPriors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popGmmPriors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGmmD_Callback(hObject, eventdata, handles)
% hObject    handle to txtGmmD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGmmD as text
%        str2double(get(hObject,'String')) returns contents of txtGmmD as a double


% --- Executes during object creation, after setting all properties.
function txtGmmD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGmmD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGmmN_Callback(hObject, eventdata, handles)
% hObject    handle to txtGmmN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGmmN as text
%        str2double(get(hObject,'String')) returns contents of txtGmmN as a double


% --- Executes during object creation, after setting all properties.
function txtGmmN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGmmN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnCluster.
function btnCluster_Callback(hObject, eventdata, handles)
% hObject    handle to btnCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Execute clustering

% Determine user selection of method
method = get(handles.popMethod, 'String');
method = method{get(handles.popMethod, 'Value')};

gd = guidata(hObject);

if ~isfield(gd, 'data')
    status('You must generate data first');
    return;
end

data = gd.data;


switch method
    case 'BIC' % GMMs using Bayes information criterion
        status('Starting BIC Method');
        range = [str2double(get(handles.txtMinSearch,'String')), str2double(get(handles.txtMaxSearch,'String'))];
        
        gs = cell(max(range)-min(range)+1, 1);
        
        i = 1;
        
        ks = zeros(max(range)-min(range), 1);
        
        % Test all GMMs in search range and plot at each iteration
        for k = min(range):max(range)
            status(['Fitting with ' num2str(k) ' components..']);
            drawnow;
            g = gmmfit(data, k, get(handles.chkVaryConditions, 'Value'));
            gs{i} = g;
            
            hold on;
            if isfield(gd, 'h');
                if ishandle(gd.h)
                    delete(gd.h);
                end
            end
            h = ezcontour(handles.ax, @(x, y)pdf(g, [x y]), [min(data(:, 1)) max(data(:, 1))], [min(data(:, 2)) max(data(:, 2))]);
            gd.h = h;
            hold off;
            drawnow;
            
            ks(i) = k;
            i = i + 1;
        end
        
        
        
        % Find and plot AIC, BIC, and -log(L) for all k's
        [AICs, ~, BICs, ~, NlogLs] = gmm_stats(gs, size(data, 1));
        
        [~, mAIC] = min(AICs);
        [~, mBIC] = min(BICs);
        [~, mNlogL] = min(NlogLs);
        
        % Plot the values of all statistics
        figure;
        hold on;
        plot(ks, BICs);
        text(ks, BICs, 'B');
        plot(ks, AICs);
        text(ks, AICs, 'A');
        plot(ks, NlogLs);
        text(ks, NlogLs, 'L');
        
        % Plot minima of all stats
        plot(ks(mAIC), AICs(mAIC), 'r*');
        plot(ks(mBIC), BICs(mBIC), 'r*');
        plot(ks(mNlogL), NlogLs(mNlogL), 'r*');
        
        hold off;
        title(['Stats for GMMs tested from ', num2str(min(range)), ' to ', num2str(max(range))]);
        
        ylabel('BICs');
        status('Done!');
        
        % Replot the best GMM
        [~, ind] = min(BICs);
        
        % Plot the GMM with best fit
        g = gs{ind};
        hold on;
        delete(gd.h);
        h = ezcontour(handles.ax, @(x, y)pdf(g, [x y]), [min(data(:, 1)) max(data(:, 1))], [min(data(:, 2)) max(data(:, 2))]);
        gd.h = h;
        hold off;
        
        set(handles.txtOutput, 'String', num2str(ks(ind)));
        
    case 'Gap Statistic' % Gap statistic as detailed in R. Tibshirani, G. Walther, and T. Hastie, ?Estimating the number of clusters in a dataset via the gap statistic,? vol. 63, pp. 411?423, 2000.
        status('Starting Gap Satistic Method');
        
        % Find search range
        range = [str2double(get(handles.txtMinSearch,'String')), str2double(get(handles.txtMaxSearch,'String'))];
        
        gs = cell(max(range)-min(range)+1, 1);
        
        i = 1;
        
        ks = zeros(max(range)-min(range), 1);
        
        % Test all GMMs in search range
        for k = min(range):max(range)
            status(['Fitting with ' num2str(k) ' components..']);
            drawnow;
            g = gmmfit(data, k, get(handles.chkVaryConditions, 'Value'));
            gs{i} = g;
            
            hold on;
            if isfield(gd, 'h');
                if ishandle(gd.h)
                    delete(gd.h); % Delete existing plot if it exists
                end
            end
            
            % Create new plot
            h = ezcontour(handles.ax, @(x, y)pdf(g, [x y]), [min(data(:, 1)) max(data(:, 1))], [min(data(:, 2)) max(data(:, 2))]);
            gd.h = h;
            hold off;
            drawnow;
            
            ks(i) = k;
            i = i + 1;
        end
        
        [gaps, Ws, Wref] = calc_gaps(data, gs);
        
        figure;
        hold on;
        plot(ks, gaps, 'k');
        text(ks, gaps, 'G');
        
        plot(ks, Ws, 'b');
        text(ks, Ws, 'O');
        
        plot(ks, Wref, 'r');
        text(ks, Wref, 'R');
        hold off
        
        title(['Gap Statistic for k = ' num2str(min(range)) ' to ' num2str(max(range))]);
        legend({'Gap statistic', 'Measured Objective Function', 'Reference Objective Function'})
        
        
        status('Done!');
        
        [~, ind] = max(gaps); % Find maxima of gap statistic
        
        % Plot best performing GMM
        g = gs{ind};
        hold on;
        delete(gd.h);
        h = ezcontour(handles.ax, @(x, y)pdf(g, [x y]), [min(data(:, 1)) max(data(:, 1))], [min(data(:, 2)) max(data(:, 2))]);
        gd.h = h;
        hold off;
        
        set(handles.txtOutput, 'String', num2str(ks(ind)));
        
    case 'Hierarchical' % Hierarchical clustering with heuristic moving average filter
        y = linkage(data, 'average');
        dy = y(:, 3) - [0; y(1:end-1, 3)];
        
        diffs = zeros(length(dy) - 1, 1);
        for i = 2:length(dy)
            avgdy = mean(dy(1:i-1));
            
            diffs(i - 1) = (dy(i) - avgdy) / avgdy;
            
        end
        
        figure;
        h1=subplot(311);
        stem(y(:, 3));
        title('Distance of Joins, D(c_i, c_j)');
        
        h2=subplot(312);
        stem(dy);
        title('Differential distance, \Delta');
        
        h3=subplot(313);
        stem([0;diffs]);
        title('Moving Avg of \Delta');
        
        linkaxes([h1 h2 h3], 'x');
        
        [~, maxind] = max(diffs);
        clusters = length(diffs) - maxind + 2;
        
        C = cluster(y, 'cutoff', y(maxind - 6, 3), 'depth', length(y));
%         length(C)
%         axes(handles.ax);
%         hold on;
%         for i = 1:max(C)
%             conv = convhull(data(C == i, 1), data(C == i, 2));
%             plot(data(conv, 1), data(conv, 2), 'k-');
%             drawnow;
%             pause(.5);
%         end
%         hold off;
disp('hierarchical');
        set(handles.txtOutput, 'String', num2str(clusters));
end

guidata(hObject, gd);


% Set the status label on the GUI
function status(sts)
global status_h

set(status_h, 'String', sts);


% --- Executes on button press in chkVaryConditions.
function chkVaryConditions_Callback(hObject, eventdata, handles)
% hObject    handle to chkVaryConditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkVaryConditions


% --- Executes on button press in chkVaryIterative.
function chkVaryIterative_Callback(hObject, eventdata, handles)
% hObject    handle to chkVaryIterative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkVaryIterative



function txtKIterative_Callback(hObject, eventdata, handles)
% hObject    handle to txtKIterative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKIterative as text
%        str2double(get(hObject,'String')) returns contents of txtKIterative as a double


% --- Executes during object creation, after setting all properties.
function txtKIterative_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKIterative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
