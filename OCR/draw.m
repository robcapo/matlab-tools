function draw
% DRAW  Creates an Optical Character Recognition classifier that classifies
% user drawn numbers
%
% Created by Rob Capo - Sept 2012
% Rowan University ECE

    global is_pressed big_pix small_pix
    
    % Global to tell whether mouse is clicked or not
    is_pressed = 0;
    
    % big_pix, small_pix:
    % Axes structures with
    %  .box_x      As used in fill() for the drawing space on the axes
    %  .box_y      As used in fill() ..
    %  .blocks     Grid to hold pixel values
    %  .path_tag   String containing path to draw pixels on
    %  .axes_h     Handle for axes object
    
    % Create canvas
    canvas = figure;
    
    
    % Aesthetics
    bgcolor = [.8 .8 .8];
    set(canvas,'numbertitle','off');               
    set(canvas,'name', 'Character Recognition');
    set(canvas,'MenuBar','none');
    set(canvas,'doublebuffer','on');
    set(canvas,'tag','CHARREC');
    set(canvas,'Color', bgcolor);
    set(canvas, 'Units', 'pixels');
    pos = get(canvas, 'Position');
    pos(3) = 390;
    pos(4) = 328;
    set(canvas, 'Position', pos);

    % Define callbacks for mouse move, mouse down, mouse up
    set(canvas,'WindowButtonMotionFcn',@mousemove);
    set(canvas,'WindowButtonDownFcn',@mousedown);
    set(canvas,'WindowButtonUpFcn',@mouseup);
    
    
    big_pix = struct;
    big_pix.box_x    = [.1 .9 .9 .1];
    big_pix.box_y    = [.9 .9 .1 .1];
    big_pix.blocks   = zeros(32, 32);
    big_pix.path_tag = 'big_path';
    big_pix.axes_h = axes('Units', 'pixels',...
                          'position', [10 10 216 288]);
    
    % Create the box for user input
    fill(big_pix.box_x, big_pix.box_y, 'w');
    xlim([min(big_pix.box_x) max(big_pix.box_x)]);
    ylim([min(big_pix.box_y) max(big_pix.box_y)]);
    axis(big_pix.axes_h,'off');
    hold(big_pix.axes_h, 'on');
    
    small_pix = struct;
    small_pix.box_x    = [.1 .9 .9 .1];
    small_pix.box_y    = [.9 .9 .1 .1];
    small_pix.blocks   = zeros(8, 8);
    small_pix.path_tag = 'small_path';
    small_pix.axes_h = axes('Units', 'pixels',...
                          'position', [236 154 144 144]);
                      
    % Create the box for user input
    fill(small_pix.box_x, small_pix.box_y, 'w');
    xlim([min(small_pix.box_x) max(small_pix.box_x)]);
    ylim([min(small_pix.box_y) max(small_pix.box_y)]);
    axis(small_pix.axes_h,'off');
    hold(small_pix.axes_h, 'on');


    % Draw Process button
    uicontrol('Style', 'pushbutton',...
                'String', 'Classify',...
                'units', 'pixels',...
                'Position', [236 109 67 35],...
                'Callback', @process);

    % Draw Clear button
    uicontrol('Style', 'pushbutton',...
                'String', 'Clear',...
                'units', 'pixels',...
                'Position', [313 109 67 35],...
                'Callback', @clear);

    % Label "Draw a number here:"
    uicontrol('Style', 'text', ...
                'String', 'Draw a number here:', ...
                'units', 'pixels', ...
                'Position', [10 298 288 20],...
                'FontSize', 13,...
                'BackgroundColor', bgcolor,...
                'HorizontalAlignment', 'left');
    
    % Label "Classifier Input"
    uicontrol('Style', 'text', ...
                'String', 'Classifier input:', ...
                'units', 'pixels', ...
                'Position', [236 298 144 20],...
                'FontSize', 13,...
                'BackgroundColor', bgcolor,...
                'HorizontalAlignment', 'left');
    
    % Label "Output"
    uicontrol('Style', 'text', ...
                'String', 'Output:', ...
                'units', 'pixels', ...
                'Position', [236 79 144 20],...
                'FontSize', 13,...
                'BackgroundColor', bgcolor,...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'Tag', 'OUTPUTLABEL');
            
    % PC font sizes are larger for some reason.
    % change the font size of the output label depending on OS
    if ispc, outsize = 36; else outsize = 62; end
    
    % Label for number output
    uicontrol('Style', 'text', ...
                'String', '', ...
                'units', 'pixels', ...
                'Position', [236 10 144 59],...
                'FontSize', outsize,...
                'BackgroundColor', bgcolor,...
                'HorizontalAlignment', 'center',...
                'Tag', 'OUTPUT');
    
end

% Called when the mouse is moved over the window
function mousemove(varargin)
    global big_pix is_pressed
    
    % Only draw if the mouse is pressed
    if is_pressed == 1
        % Get current x, y of mouse
        p = get(big_pix.axes_h,'CurrentPoint');
        x = p(1,1);
        y = p(1,2);

        % Detect which of the 32x32 pixels the mouse is in
        [r c] = detect_pixel(x, y, big_pix);

        % Fill the pixels around the current pixel
        [rf cf] = fill_pixel(r, c, big_pix);

        % Set the pixels around the current pixel to 1
        % (currently a 3x3 square brush)
        big_pix.blocks(rf, cf) = 1;
    end
    
end

% Called when the mouse is clicked
function mousedown(varargin)
    global big_pix is_pressed
    
    % Get current x, y of mouse
    p = get(big_pix.axes_h,'CurrentPoint');
    x = p(1,1);
    y = p(1,2);
    
    % Detect which of the 32x32 pixels the mouse is in
    [r c] = detect_pixel(x, y, big_pix);
    
    % Fill the pixels around the current pixel
    [rf cf] = fill_pixel(r, c, big_pix);
    
    % Set the pixels around the current pixel to 1
    % (currently a 3x3 square brush)
    big_pix.blocks(rf, cf) = 1;
    
    % Mouse is down, so set is_pressed to 1
    is_pressed = 1;
end

% Called when mouse is unclicked
function mouseup(varargin)
    global is_pressed
    
    is_pressed = 0;
    
    convert_image
end

% Clears the canvas
function clear(varargin)
    global big_pix small_pix
    
    big_pix.blocks = zeros(size(big_pix.blocks));
    small_pix.blocks = zeros(size(small_pix.blocks));
    
    delete(findobj('Tag', big_pix.path_tag));
    delete(findobj('Tag', small_pix.path_tag));
    
    h_text = findobj('Tag', 'OUTPUT');
    set(h_text, 'String', '');
    
    set(findobj('Tag', 'OUTPUTLABEL'), 'Visible', 'off');
end

% Classifies the user input
function process(varargin)
    global big_pix
    
    % Create 4x4 groups of features to count number of pixels
    features = zeros(8, 8);
    
    % Iterate through bitmap creating ranges
    for row = 1:8
        row_start = (row - 1) * 4 + 1;
        row_end = row_start + 3;
        
        for col = 1:8
            col_start = (col - 1) * 4 + 1;
            col_end = col_start + 3;
            
            % Search 4x4 ranges of bitmap and count number of pixels that
            % are filled in, creating an 8x8 featurespace of 0 - 16
            features(row, col) = sum(sum(big_pix.blocks(row_start:row_end, col_start:col_end)));
        end
        
    end
    
    % Classify the image using a KNN
    load('data');
    
    features = reshape(features', 1, 64);
    
    out = knn(features, train_data, train_labels, 5);
    
    set(findobj('Tag', 'OUTPUTLABEL'), 'Visible', 'on');
    h_text = findobj('Tag', 'OUTPUT');
    set(h_text, 'String', num2str(out));
    
end

%[r c] = detect_pixel(x, y, axs)
% Detects which pixel a cursor is in given a grid size, axes lims, and
% current x, y position
%
% Inputs:
%  x, y             Position of current point
%  axs              Axes structure with:
%   .box_x, .box_y  Axis limits in the form [min max] or as used in fill()
%   .blocks         matrix: zeros(#pixel rows, #pixel cols)
%
% Outputs:
% r, c          Row and column that point is in
function [r c] = detect_pixel(x, y, axs)
    [rs cs] = size(axs.blocks);
    x_lim = [min(axs.box_x) max(axs.box_x)];
    y_lim = [min(axs.box_y) max(axs.box_y)];
    
    % Make sure pixel is within drawing box
    if x >= x_lim(1) && x <= x_lim(2) && y >= y_lim(1) && y <= y_lim(2)
        
        % Get x & y b/w 0 - 1;
        x = (x - x_lim(1)) / (x_lim(2) - x_lim(1));
        y = 1 - ((y - y_lim(1)) / (y_lim(2) - y_lim(1)));
        
        % Get r & c to be a whole number between 0 and max #pixels
        c = ceil(x * cs);
        r = ceil(y * rs);
    else
        r = [];
        c = [];
    end
end

%fill_pixel(r, c, axs)
% Fills the pixel on the axs struct in a given row and column with a 3x3
% brush.
% 
% Inputs:
%  r, c            Row and column of pixel to be filled
%  axs             Axes structure with:
%   .axes_h        Handle to axes
%   .box_x, box_y  Box matrices. Either in the same form as for fill() or
%                   in the xlim/ylim vector form: [x/y_min x/y_max]
%   .blocks        Grid in the form: zeros(#pixel rows, #pixel cols)
%   .path_tag      String containing the tag of the fill path
%  colr            Color to fill pixel
%  siz             Brush size
%
% Outputs:
% rf, cf            Rows and columns filled
function [rf cf] = fill_pixel(r, c, axs, colr, siz)
    if nargin < 4, colr = 'black'; end
    if nargin < 5, siz = 3; end
    if mod(siz, 2) ~= 1, siz = siz + 1; end
    
    if ~isempty(r) && ~isempty(c)
        % Find min and max of box
        x_lim = [min(axs.box_x) max(axs.box_x)];
        y_lim = [min(axs.box_y) max(axs.box_y)];

        [rs cs] = size(axs.blocks);


        % Create intervals corresponding to each corner of a pixel on bmp
        [x_ints y_ints] = meshgrid(linspace(x_lim(1), x_lim(2), cs + 1), linspace(y_lim(2), y_lim(1), rs + 1));

        % Get size of grid
        [gr gc] = size(x_ints);
        
        % Brush span left and brush span right
        bsl = (siz - 1) / 2;
        bsr = bsl + 1;

        % Make sure row and column are padded enough for brush to fit
        r = max(r, bsl + 1);
        c = max(c, bsl + 1);
        r = min(r, gr - bsr);
        c = min(c, gc - bsr);
        
        % Declare outputs of which rows and columns were filled
        rf = r-bsl:r+bsr-1;
        cf = c-bsl:c+bsr-1;

        % Fill in the 3x3 brush space
        fill_x = [x_ints(r - bsl, c - bsl) x_ints(r - bsl, c + bsr) x_ints(r + bsr, c + bsr) x_ints(r + bsr, c - bsl)];
        fill_y = [y_ints(r - bsl, c - bsl) y_ints(r - bsl, c + bsr) y_ints(r + bsr, c + bsr) y_ints(r + bsr, c - bsl)];

        axes(axs.axes_h)
        h = fill(fill_x, fill_y, colr, 'Tag', axs.path_tag);
        set(h, 'EdgeColor', 'none');
        xlim(x_lim);
        ylim(y_lim);
        axis(axs.axes_h,'off');
    else
        rf = [];
        cf = [];
    end
end

%convert_image
% Converts the 32x32 binary grid to a 8x8 grayscale grid
% Executes on a timer
function convert_image(varargin)
    global small_pix big_pix
    
    % Iterate through bitmap creating ranges
    for row = 1:8
        row_start = (row - 1) * 4 + 1;
        row_end = row_start + 3;

        for col = 1:8
            col_start = (col - 1) * 4 + 1;
            col_end = col_start + 3;

            % Search 4x4 ranges of bitmap and count number of pixels that
            % are filled in, creating an 8x8 featurespace of 0 - 16

            section_sum = sum(sum(big_pix.blocks(row_start:row_end, col_start:col_end)));
            if small_pix.blocks(row, col) ~= section_sum
                small_pix.blocks(row, col) = section_sum;
                fill_pixel(row, col, small_pix, ones(1, 3)-sqrt(section_sum/16), 1);
            end
        end
    end
end