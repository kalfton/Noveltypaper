function varargout = findseq(A,dim)

% FINDSEQ Find sequences of repeated (adjacent/consecutive) numeric values
%
%   FINDSEQ(A) Find sequences of repeated numeric values in A along the
%              first non-singleton dimension. A should be numeric.
%
%   FINDSEQ(...,DIM) Look for sequences along the dimension specified by the 
%                    positive integer scalar DIM.
%
%   OUT = findseq(...)
%       OUT is a "m by 4" numeric matrix where m is the number of sequences found.
%       
%       Each sequence has 4 columns where:
%           - 1st col.:  the value being repeated
%           - 2nd col.:  the position of the first value of the sequence
%           - 3rd col.:  the position of the last value of the sequence
%           - 4th col.:  the length of the sequence
%       
%   [VALUES, INPOS, FIPOS, LEN] = findseq(...)
%       Get OUT as separate outputs. 
%
%       If no sequences are found no value is returned.
%       To convert positions into subs/coordinates use IND2SUB
%
% 
% Examples:
%
%     % There are sequences of 20s, 1s and NaNs (column-wise)
%     A   =  [  20,  19,   3,   2, NaN, NaN
%               20,  23,   1,   1,   1, NaN
%               20,   7,   7, NaN,   1, NaN]
%
%     OUT = findseq(A)
%     OUT =  
%            20        1          3        3
%             1       14         15        2
%           NaN       16         18        3
%
% Additional features:
% - <a href="matlab: web('http://www.mathworks.com/matlabcentral/fileexchange/28113','-browser')">FEX findseq page</a>
% - <a href="matlab: web('http://www.mathworks.com/matlabcentral/fileexchange/6436','-browser')">FEX rude by us page</a>
%
% See also: DIFF, FIND, SUB2IND, IND2SUB

% Author: Oleg Komarov (oleg.komarov@hotmail.it) 
% Tested on R14SP3 (7.1) and on R2009b. In-between compatibility is assumed.
% 02 jul 2010 - Created
% 05 jul 2010 - Reorganized code and fixed bug when concatenating results
% 12 jul 2010 - Per Xiaohu's suggestion fixed bug in output dimensions when A is row vector
% 26 aug 2010 - Cast double on logical instead of single
% 28 aug 2010 - Per Zachary Danziger's suggestion reorganized check structure to avoid bug when concatenating results

% NINPUTS
error(nargchk(1,2,nargin));

% NOUTPUTS
error(nargoutchk(0,4,nargout));

% IN
if ~isnumeric(A) 
    error('findseq:fmtA', 'A should be numeric')
elseif isempty(A) || isscalar(A)
    varargout{1} = [];
    return
end

% DIM
szA = size(A);
if nargin == 1 || isempty(dim)
    % First non singleton dimension
    dim = find(szA ~= 1,1,'first');
elseif ~(isnumeric(dim) && dim > 0 && rem(dim,1) == 0) || dim > numel(szA)
    error('findseq:fmtDim', 'DIM should be a scalar positive integer <= ndims(A)');
end

% Less than two elements along DIM
if szA(dim) == 1
    varargout{1} = [];
    return
end

% ISVECTOR
if nnz(szA ~= 1) == 1
    A = A(:);
    dim = 1;
    szA = size(A);
end

% Detect 0, NaN, Inf and -Inf
OtherValues    = cell(1,4);
OtherValues{1} = A ==    0;
OtherValues{2} = isnan(A) ;
OtherValues{3} = A ==  Inf;
OtherValues{4} = A == -Inf;
Values         = [0,NaN, Inf,-Inf];

% Remove zeros
A(OtherValues{1}) = NaN;                             

% Make the bread
bread = NaN([szA(1:dim-1),1,szA(dim+1:end)]);

% [1] Get chunks of "normal" values
Out = mainengine(A,bread,dim,szA);

% [2] Get chunks of 0, NaN, Inf and -Inf
for c = 1:4
    if nnz(OtherValues{c}) > 1
        % Logical to double and NaN padding
        OtherValues{c} = double(OtherValues{c});                        
        OtherValues{c}(~OtherValues{c}) = NaN;                          
        % Call mainengine and concatenate results
        tmp = mainengine(OtherValues{c}, bread,dim,szA);
        if ~isempty(tmp)
            Out = [Out; [ones(size(tmp,1),1)*Values(c) tmp(:,2:end)]];  %#ok
        end
    end
end

% Distribute output
if nargout == 1 || nargout == 0  
    varargout{1} = Out;
else
    for k = 1:nargout
        varargout{k} = Out(:,k);
    end
end

end

% MAINENGINE This functions uses run length encoding and retrieve positions 
function Out = mainengine(meat,bread,dim,szMeat)

% Make a sandwich  
sandwich    = cat(dim, bread, meat, bread);

% Find chunks (run length encoding engine)
IDX         = diff(diff(sandwich,[],dim) == 0,[],dim);

% Initial positions and final positions
InPos       = find(IDX  ==  1);
FiPos       = find(IDX  == -1);

% Coordinates of initial/final positions
nd          = ndims(sandwich);                           
[cIn{1:nd}] = ind2sub(szMeat,InPos); 
[cFi{1:nd}] = ind2sub(szMeat,FiPos); 

% Assign output
Out         = [meat(InPos)                ,...    % Values
               InPos                      ,...    % Initial positions 
               FiPos                      ,...    % Final   positions
               cFi{dim} - cIn{dim} + 1   ];...    % Length of the blocks

end
