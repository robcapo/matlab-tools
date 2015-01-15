function safe_save(varargin)
% SAFE_SAVE('dir', 'var1', 'var2', ..) save files without messups.
%  SAFE_SAVE creates the destination directory if it doesn't exist and
%  appends numbers to the end of your file if the name exists.
%
%  SAFE_SAVE('some/path/')
%   Creates 'some/path/' if it doesn't exist and saves in it:
%    myfile.mat if available, or
%    myfile1.mat if available, or
%    myfile2.mat etc.
%
% By Rob Capo

    fpath = varargin{1};
    slashes = strfind(fpath, '/');
    
    % If we're saving in a directory other than current.
    if ~isempty(slashes)
        folder = fpath(1:slashes(end));
        file = fpath(slashes(end)+1:end);
        
        % Make directory if it doesn't exist
        if ~exist(folder, 'dir'), mkdir(folder); end
    else
        folder = [];
        file = fpath;
    end
    
    % Append default .mat extension if one isn't specified
    if isempty(strfind(file, '.'))
        file = [file '.mat'];
    end
    
    % Append a number to the filename if one exists
    if exist([folder file], 'file')
        i = 1;
        dots = strfind(file, '.');
        filename = file(1:dots(end) - 1);
        ext = file(dots(end):end);
        
        % Increment i until we get to a filename that doesn't already exist.
        while exist([folder filename num2str(i) ext], 'file')
            i = i + 1;
        end
        
        file = [filename num2str(i) ext];
    end
    
    varargin{1} = [folder file];
    
    for i = 1:length(varargin) - 1
        varargin{i} = [' ''' varargin{i} ''','];
    end
    
    varargin{end} = [' ''' varargin{end} ''' '];
    
    func = ['save(' varargin{:} ')'];
    evalin('caller', func);
end
