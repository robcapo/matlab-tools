function varargout = main_gui(varargin)
% MAIN_GUI MATLAB code for main_gui.fig
%      MAIN_GUI, by itself, creates a new MAIN_GUI or raises the existing
%      singleton*.
%
%      H = MAIN_GUI returns the handle to a new MAIN_GUI or the handle to
%      the existing singleton*.
%
%      MAIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GUI.M with the given input arguments.
%
%      MAIN_GUI('Property','Value',...) creates a new MAIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_gui

% Last Modified by GUIDE v2.5 06-Nov-2013 12:23:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @main_gui_OutputFcn, ...
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


% --- Executes just before main_gui is made visible.
function main_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_gui (see VARARGIN)

% Choose default command line output for main_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slide_cut_Callback(hObject, eventdata, handles)
% hObject    handle to slide_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Plot the cluster connections as the slider moves
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

cla;
hold on;
scatter(handles.x(:, 1), handles.x(:, 2), 'b');
for i = 1:ind
    c = handles.c{i};
    
    for j = 1:size(c, 1)
        plot([handles.x(c(j, 1), 1); handles.x(c(j, 2), 1)], [handles.x(c(j, 1), 2); handles.x(c(j, 2), 2)]);
    end
end
hold off;


% --- Executes during object creation, after setting all properties.
function slide_cut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slide_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btn_gen.
function btn_gen_Callback(hObject, eventdata, handles)
% hObject    handle to btn_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.btn_cluster, 'Enable', 'on');
set(handles.btn_dend, 'Enable', 'on');

n = str2double(get(handles.edit_n, 'String'));
k = str2double(get(handles.edit_k, 'String'));

% Generate data
mus = repmat(10*rand(k, 2), n, 1);
sigmas = repmat([1 0; 0 1], [1 1 n*k]);

x = mvnrnd(mus, sigmas);

% Compute linkage
methods = get(handles.menu_d, 'String');
y = linkage(pdist(x), methods(get(handles.menu_d, 'Value')));

% Generate connection array
c = cell(size(y, 1));

for i = 1:size(y, 1)
    c{i} = get_connections(y, k*n, i);
end



% Set data for other callbacks
handles.x = x;
handles.y = y;
handles.c = c;

set(handles.slide_cut, 'Max', size(y, 1));
set(handles.slide_cut, 'Value', 0);
set(handles.slide_cut, 'SliderStep', [1/size(y, 1) 1/size(y, 1)]);

% Plot data
axes(handles.axh);
cla;
scatter(x(:, 1), x(:, 2));

guidata(hObject, handles);

% Get all connections between observations from linkage matrix
function cons = get_connections(y, n, i)
    instances = unique(get_instances(y, n, i));
    cons = combnk(instances, 2);
    
    
% Get observations belonging to a cluster based on linkage matrix
function instances = get_instances(y, n, i)
    
    if y(i, 1) > n
        y(i, 1)
        instances = get_instances(y, n, mod(y(i, 1), n));
    else
        instances = y(i, 1);
    end
    
    if y(i, 2) > n
        y(i, 2)
        instances = [instances; get_instances(y, n, mod(y(i, 2), n))];
    else
        instances = [instances; y(i, 2)];
    end

% --- Executes on selection change in menu_d.
function menu_d_Callback(hObject, eventdata, handles)
% hObject    handle to menu_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_d contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_d


% --- Executes during object creation, after setting all properties.
function menu_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_k_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k as text
%        str2double(get(hObject,'String')) returns contents of edit_k as a double


% --- Executes during object creation, after setting all properties.
function edit_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_n_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n as text
%        str2double(get(hObject,'String')) returns contents of edit_n as a double


% --- Executes during object creation, after setting all properties.
function edit_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_dend.
function btn_dend_Callback(hObject, eventdata, handles)
% hObject    handle to btn_dend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Generate dendrogram from linkage matrix
figure;
dendrogram(handles.y);


% --- Executes on button press in btn_cluster.
function btn_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to btn_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Compute linkage

x = handles.x;

methods = get(handles.menu_d, 'String');
y = linkage(pdist(x), methods(get(handles.menu_d, 'Value')));

% Generate connection array
c = cell(size(y, 1));

for i = 1:size(y, 1)
    c{i} = get_connections(y, size(x, 1), i);
end



% Set data for other callbacks
handles.y = y;
handles.c = c;

% Set slider range
set(handles.slide_cut, 'Max', size(y, 1));
set(handles.slide_cut, 'Value', 0);
set(handles.slide_cut, 'SliderStep', [1/size(y, 1) 1/size(y, 1)]);

% Plot data
axes(handles.axh);
cla;
scatter(x(:, 1), x(:, 2));

guidata(hObject, handles);
