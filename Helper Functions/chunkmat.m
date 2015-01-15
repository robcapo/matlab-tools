% A = CHUNKMAT(data, chunksize)
% Divides the rows of data matrix into chunks of size chunksize and puts
% them into cell vector A.
%
% Example: A = chunkmat(magic(12), 3) would set A to a cell vector of
% length 4 with each cell containing consecutive 3x12 chunks of the
% magic(12) matrix
function A = chunkmat(data, chunksize)
    n = size(data, 1);
    chunks = ceil(n/chunksize);
    A = cell(chunks, 1);
    
    for i = 1:chunks
        start = (i-1)*chunksize+1;
        fin = start+chunksize-1;
        
        if i == chunks
            fin = n;
        end
        
        A{i} = data(start:fin, :);
    end
end