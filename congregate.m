%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -BRIEF DESCRIPTION:
%       Fetch values from an array of cells (or an array of structures) and create a matrix.
%
%%%% -DETAILED DESCRIPTION:
%       Have you ever tried to create a submatrix from an array of cells? Or extract values of a specific member from an array of strctures? 
%       If so, this function is for you!! Please see the examples below to get more info upon this function.
%
%%%% -EXAMPLES:
%
%       1. Extracting submatrix from an array of cells
%           A = {11,12,13;21,22,23};
%           congregate('A{1:2,2:3}')
%               ans =
%                   12    13
%                   22    23
%
%       2. Converting an array of structures to a 2x3 matrix
%           S(1).a=[1,2,3]';
%           S(2).a=[4,5,6]';
%           congregate('S(:).a')
%                 ans =
%                      1     2     3
%                      4     5     6
%
%       3. Converting an array of structures of structures to a 2x2x3 matrix
%           T(1).a(1).b=[111,112,113]';
%           T(1).a(2).b=[121,122,123]';
%           T(2).a(1).b=[211,212,213]';
%           T(2).a(2).b=[221,222,223]';
%           congregate('T(:).a(:).b')
%                 ans(:,:,1) =
%                    111   121
%                    211   221
%                 ans(:,:,2) =
%                    112   122
%                    212   222
%                 ans(:,:,3) =
%                    113   123
%                    213   223
%
%
%%%% -OPTIONS:
%       1. 'isFillWithNan' = true (Default) means, absent fields would be regarded as NaN 
%       2. 'squeezeLevel': determines the level of squeezing of the input. Squeezing is defined as the removal of all singleton dimensions in a matrix/cell.
%               'squeezeLevel' = 0 is the default, which implies no squeezing
%               'squeezeLevel' = 1 implies moderate squeezing, where the content inside cells are squeezed.
%               'squeezeLevel' = 2 implies high level squeezing, where the content inside both cells and matrices are squeezed.
%
%%%% -DETAILED NOTES:
%       1. Here number of dimensions = number of brackets appearing in 'inString' + dimension of last field  
%       2. Single ':' in structExprStr implies output is a column vector of cells  
%       3. s{[1,2;3,4],[1,2]}.a will result in 1D output
%
%%%% -PROGRAMMING NOTES:
%       1. There are some design errors in MATLAB software, as follows
%           a). MATLAB ignores, tralining singleton dimensions. 
%           b). 'a(2:5)' and 'a{:}' are outputted as row-vectors, instead of column vecotors
%           c). ndims of column vector is 2 (instead of 1)!!. 
%       Solution:
%           The problem (a) is rectified using 'nDimPerScope'.  To overcome (c), I implemented ndims_().
%
%%%% -TROUBLESHOOTS:
%       1. One way to check if this function works or not, is to execute your input by replacing all variable indices with '1's 
%           i.e.  s{1}.s1(1).s2{1}.a
%
%%%% -TODO:
%       1. Include accepting 'end' keyword.
%       2. Nonworking Input: congregate_('in{:,2,3}{1}.F')% 
%
%%%% -VERSIONS:
%       Version 1 release: 2019 January 31
%
%%%% -AUTHORS:
%       Anver Hisham <anverhisham@gmail.com>
%
%%%% -License
%       This work is under a Creative Commons Share Alike (CC-SA) license
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Main function
function out = congregate(inString,isFillWithNan,squeezeLevel)
   if(isempty_('isFillWithNan')) isFillWithNan = true; end
   if(isempty_('squeezeLevel')) squeezeLevel = 0; end
   %% -USER calling.
   %    i.e., Just bring the USER's indexing variables to this scope, rename it (also rename inString), and pass it to next scope  
       inString = strtrim(inString);
       %% Fetch variables from USER scope
       [variableNames,remainder] = regexp(inString,'(?<=^|\(|\[|\{|,|:)[a-zA-Z]\w*','match','split'); %% variableNames{1} is base object
       remainder(1)=[];
       for i=1:numel(variableNames)
          userVariables{i} = evalin('caller',variableNames{i}); 
          userVariablesNewname{i} = ['userVariables{',num2str(i),'}'];
       end
       inStringNew = strjoin_(userVariablesNewname,remainder);
       [out,nDimPerScope] = congregate_Recursive_(inStringNew,isFillWithNan);
       out = cellOfCell2mat_(out,squeezeLevel,nDimPerScope);
end



%% NOTE: Calling of this function is restricted to only congregate_()
%   'nDimsCurrent' : Local variable indicating whether current scope deserve how many dimensions
%            NaN   => Discard this value in cellOfCell2mat_   
%            0     => 0 dimension in cellOfCell2mat_   
function [out,nDimPerScope] = congregate_Recursive_(inString,isFillWithNan,recursiveCallingScope) 
       if(nargin<3) recursiveCallingScope = 1; end                      %% Avoided using isempty_() for speed
       userVariables = evalin('caller','userVariables');
       % Get 's(1).s1(p)' from 's(1).s1(p).s2(q).a'
       [from, to] = findComplicatedParentheses_(inString,1);
       %  TODO_Verify(2016/11/3):   Removed && to~=numel(inString) from below expression
       if(~isempty(to) && ~isempty_(inString(1:to)))  %% TODO_Verify(2016/11/3): Changed from isempty_(['[',inString(1:to),']']) to isempty_(['inString(1:to)'])
           inStringnew = inString;
           inStringnew(from) = '(';  inStringnew(to) = ')';
           tmpIn = eval(inStringnew(1:to));                             %% Changing the base object
           [userVariables{1},nDimsCurrent] = congregate_reshape_(tmpIn); %% Changing vector 'tempIn' to intended matrix 'userVariables{1}'
           out = cell(size(userVariables{1}));
           for i=1:numel(out)
               if(isa(userVariables{1},'cell'))
                inStringNew = ['userVariables{1}{',num2str(i),'}',inString(to+1:end)];
               else
                inStringNew = ['userVariables{1}(',num2str(i),')',inString(to+1:end)];
               end
               [out{i},nDimsSub{i}] = congregate_Recursive_(inStringNew,isFillWithNan,recursiveCallingScope+1);
           end
           nDimsSub = nanmax(cellOfCell2mat_(nDimsSub)',[],2);
           nDimPerScope = [nDimsCurrent;nDimsSub];      %% Since sometimes 'inStringNew' would be empty cells
       else
          if(recursiveCallingScope==1 && ~isempty(to) && isempty_(['[',inString(1:to),']']))  %% TODO_Verify(2016/11/3): Changed from isempty_(['[',inString(1:to),']']) to isempty_(inString(1:to))
             error(['Error: Non-existent parent object  ']);
          end
          if(eval(['~isempty_(''',inString,''')']))
             out = eval(inString);
          else                  %% Missing or Empty field
             if(isFillWithNan) 
                 out = NaN;
             else
                 error(['Non-existant field: ',inString]);
             end
          end
          nDimPerScope = ndims_(out);
       end
    %% This function is written in order to make sure that 'a{:,:,:}' is considered as a 3D matrix instead of doing vectorization! 
    %   'nDimsCurrent': 
    %       0   => tmpIn is a single object without any ':'
    %       NaN => Ignore this value by any clients
    %           for each ':', and '1:1' is taken as extra dimension
    % PROGRAMMING NOTE:
    %   1. MATLAB create 'a(:)' and 'a{:}' as a row-vector, but ideally it should be column vectors. So here we are always sticking with the ideal.  
    % TODO_Future: 
    %   Program to handle every complicated expressions too 
    function [tmpIn,nDimsCurrent] = congregate_reshape_(tmpIn)
        indicesInputStr = regexp(strtrim(inString(from+1:to-1)),',','split');   %% Get entries inside last bracket 
        nDimsCurrent = 0;
%                 if(numel(indicesInputStr)>1)      %% TODO: Clarify: Why did I placed this constraint in the first place? 
            clear finalSizePerDim;
            for i=1:numel(indicesInputStr)
                if(regexp(indicesInputStr{i},'^\s*:\s*$'))
                  if(numel(indicesInputStr)==1)     %% If 'a(:)'
                    finalSizePerDim(i) = numel(eval(inString(1:from-1)));
                  else
                    finalSizePerDim(i) = size(eval(inString(1:from-1)),i);
                  end
                  nDimsCurrent = nDimsCurrent + 1;
                else
                  finalSizePerDim(i) = numel(eval(indicesInputStr{i}));
                  if(finalSizePerDim(i)>1 || ~isempty(regexp(indicesInputStr{i},':')))    
                    nDimsCurrent = nDimsCurrent + 1;
                  end
                end
            end
            tmp1 = reshape(tmpIn,[finalSizePerDim,1]);     %% '1' is appended to avoid the error when numel(finalSizePerDim)==1   
            tmpIn = tmp1;
%                 end
%                 nDimsCurrent = numel(finalSizePerDim);
    end            
end


%% Function to find '([1,2,3])' in 's(1).s1([1,2,3])' (Note from & to corresponds to bracket location)
%       '(1,2)', '('abcd')' are not complicated parentheses.
%   'count' determines how many outputs needed. If count>1, then an array of cells is outputted   
function [from,to] = findComplicatedParentheses_(str,count)
    if(nargin<2) count = Inf; end                                               %% Avoided using isempty_() for speed
    if(count>1)
        % Count the number of brackets in str
        i = 1; to_last = 0;
        while i<=count && to_last < numel(str)
           [from{i},to{i}] = findComplicatedParentheses_(str(to_last+1:end),1);
           if(isempty(from{i})) 
               break;           
           else
               from{i} = from{i} + to_last;
               to{i} = to{i} + to_last;
               to_last = to{i};
               i = i+1;
           end
        end
    else
        indices = regexp(str,'{|(');
        for i=1:numel(indices)
           from = indices(i);
           to = findMatchingParentheses_(str,from);
           if(isempty(regexp(str(from+1:to-1),'^((\s*\d+\s*)|(''.*''))$')))    %% If number or a string
               return;
           end
        end
        from = []; to = [];
    end
end






%% SAME AS MATLAB 'isempty', but input can be anything (including structure fields and non-existing variables) %%%%%%%
%   NOTE: input must be name of a variable (string format)!
function out = isempty_(inputString)
    if(~isempty(strfind(inputString,':')) && ~isempty(regexp(inputString,':\d*}(?=\s*$)','once')))      %% For 'a{:}', we need to ensure that, this function is not tricked with a{1} alone.  NOTE: 'contains' not supported in old(Glenn) matlab versions :-(
        to = regexp(inputString,'}(?=\s*$)','once');
        from = findMatchingParentheses_(inputString,to);
        if(findComplicatedParentheses_(inputString(from:to),1))
            inputString = ['{',inputString,'}'];
        end
    end
    try
        out = evalin('caller',['isempty(',inputString,')']);
    catch
       out = true;                
   end
end





%% {{[1,2],[3,4,5]},{[11,12]'}}  ->  []  
%   Brief Description:
%       Each {} introduces an extra dimension (eventhough it's singleton dimension). To avoid this feature, input 'nDimPerScope' carefully. 
%   -'squeezeLevel'
%       0: No Squeezing. Expect redundant dimensions for row-vectors   
%       1: 'in' (which is a matrix of cells) and all cells within is (in all scope) are squeezed 
%       2: Both 'in' and all the numerical matrix content inside any cells are squeezed
%   - 'nDimPerScope'
%       This would make sure that, ith scope has got nDimPerScope(i) dimensions, by inserting singleton dimensions. (Highly useful in congregate_Recursive)_ 
%       nDimPerScope=nan  =>  Ignore this parameter
function out = cellOfCell2mat_(in,squeezeLevel,nDimPerScope)
    if(nargin < 2 || isempty_('squeezeLevel')) squeezeLevel = 2; end
    if(nargin < 3 || isempty_('nDimPerScope')) nDimPerScope = [NaN]; end
    if(numel(nDimPerScope)==1) nDimPerScope = [nDimPerScope;NaN]; end
    %
    if(isa(in,'cell'))
        outTemp = cell(size(in));
        if(squeezeLevel>0)
            outTemp = squeezeAndTransposeToRowVector_(outTemp);
        end
        offsetDim = ndims_(outTemp);
        for iCell = 1:numel(in)
            outTemp{iCell} = cellOfCell2mat_(in{iCell},squeezeLevel,nDimPerScope(2:end));
            outTemp{iCell} = shiftdim(outTemp{iCell},-offsetDim);
        end
        out = cell2mat_(outTemp);         %% TODO_Next: Write another program called cell2matRaw_() which wouldn't loose any dimensions. Also clearly define 'squeezeLevel'   
        if(~isnan(nDimPerScope(1)) && ndims_(outTemp)<ndims_(out))   %% Insert singleton dimensions based on 'nDimPerScope' 
            if(nDimPerScope(1)>0 && ndims_(outTemp)>nDimPerScope(1))
                error('Error(): Inputted nDimPerScope is not matching with');
            end
            if(ndims_(outTemp)<nDimPerScope(1))   %% Adding extra dimensions to the current scope
                newDimIndices = [1:ndims_(outTemp), ndims_(out)+1:ndims_(out)+1+(nDimPerScope(1)-ndims_(outTemp)-1), ...
                     ndims_(outTemp)+1:ndims_(out)];
                out = permute(out,newDimIndices);
            elseif(nDimPerScope(1)==0)              %% Deleting current dimension, by appending it the subscope dimensions 
                outSize = size(out);
                outSize(ndims_(outTemp)+1) = prod(outSize(1:ndims_(outTemp)+1));
                outSize(1:ndims_(outTemp)) = [];
                out = reshape(out,[outSize,1]);
            end
        end
    else
       if(squeezeLevel>1)
        out = squeezeAndTransposeToRowVector_(in);
       else
        out = in; 
       end
    end
end


%% Same as MATLAB strjoin, but this one concatenate equal sized in1 & in2
function out = strjoin_(in1,in2)
    if(numel(in1)==numel(in2))
        in1{end+1} = '';
    end
    out = strjoin(in1,in2);
end


%% Find matching parentheses
%   findMatchingParentheses_('a(bcd)',2)  ->  6  
%   findMatchingParentheses_('a(bcd)',6)  ->  2  
% TODO: Publish it
% https://www.mathworks.com/matlabcentral/cody/problems/111-find-matching-parenthesis
function index = findMatchingParentheses_(str,inIndex)
    parentheses = str(inIndex);
    switch(parentheses)
        case '{', antiParentheses = '}';
        case '(', antiParentheses = ')';
        case '[', antiParentheses = ']';
        case '}', antiParentheses = '{';
        case ')', antiParentheses = '(';
        case ']', antiParentheses = '[';
        otherwise, error('Index in the inputted str doesn''t corresponds to a bracket!!');
    end
    isParentheses = zeros(numel(str),1);
    ind1 = regexp(str,parentheses);
    ind2 = regexp(str,antiParentheses);
    isParentheses(ind1) = 1;
    isParentheses(ind2) = -1;
    if(regexp(parentheses,'{|(|['))
        index = find(cumsum(isParentheses(inIndex:end))==0,1,'first');
        index = inIndex-1 + index;
    else
        index = find(cumsum(isParentheses(1:inIndex),1,'reverse')==0,1,'last');
    end
end


%% Same as MATLAB inbuilt 'squeeze' except that this function also squeezes row vector to column vector 
function out =  squeezeAndTransposeToRowVector_(in)
    out = squeeze(in);
    if(numel(size(in))==2 && size(in,1)==1)
        out = in.';
    end
end


%% -Returns true, if input is a column vector
function out = isColumnVector_(in)
    if(numel(size(in))==2 && size(in,2)==1)   %% If in is a column vector
        out = true;
    else
        out = false;
    end
end

%% Same as MATLAB inbuilt 'cell2mat' except that this function compensate for incompatible sizes of c{i} by placing 'Nan'  
%   'in' is a matrix of cells, with each cell containing a matrix along high dimensions.
% -Algorithm:
%   1. Find the lowest dimension 'currentDim' with length >1
%   2. Split the 'in' along the dimension 'currentDim', and call cell2mat recursively for each sub 'in' 
%   3. Concatenate all the results got in the above step. 
function out = cell2mat_(in)
    nDim = numel(size(in));
    currentDim = find(size(in)>1,1);
    if(~isempty(currentDim))
        exString = repmat(':,',1,nDim);     %% Create string as ':,:,12,:,:'
        exString(end) = [];         %% Removing trailing ','
        outSub = cell(size(in,currentDim),1);
        for i = 1:size(in,currentDim)
            exStringCurrent = [exString(1:2*(currentDim-1)),num2str(i),exString(2*currentDim:end)];
            outSub{i} = cell2mat_(eval(['in(',exStringCurrent,')']));
        end
        out = cat_(currentDim,outSub{:});
    else
        out = in{1};
    end
end




%% Same as MATLAB inbuilt 'cell2mat' except that this function concatenate incompatible sized by appending 'Nan's 
%   cat_(1,A,B,'-appendValue',123) will append '123' for missing values 
function out = cat_(dim,varargin)
    optionIndex = contains_(varargin,'-appendValue');
    if(~isempty(optionIndex))
        optionIndex = optionIndex{1}(1);
        appendValue = varargin{optionIndex+1};
        varargin(optionIndex:optionIndex+1) = [];
    else
        appendValue = nan;
    end
    sizePerInput = [];
    for i = 1:numel(varargin)
        if(~isempty(sizePerInput) && size(sizePerInput,1) ~= numel(size(varargin{i})))
            sizePerInput = cat_(2,sizePerInput,size(varargin{i})'); %% Used cat_ instead of cat, since varargin{n} can have higher dimensions than varargin{1},varargin{2},... varargin{n-1}
        else
            sizePerInput = cat(2,sizePerInput,size(varargin{i})');    %% To avoid inifinite recursion
        end
    end
    maxSizePerDim = nanmax(sizePerInput,[],2);  
    in = cell(numel(varargin),1);
    for i = 1:numel(varargin)
        desiredSizePerDim = maxSizePerDim;
        desiredSizePerDim(dim) = size(varargin{i},dim);
        if(isnan(appendValue))
            in{i} = augmentMatrix_(varargin{i},desiredSizePerDim);                %% For speeding up
        else
            in{i} = augmentMatrix_(varargin{i},desiredSizePerDim,appendValue);
        end
    end
    out = cat(dim,in{:});
end


%% Function to compare different types of objects
%   'indexPerin2' is an array of cells, with indexPerin2{i}=indices => in1{indices}==in2{i}  
function indexPerin2 = contains_(in1,in2,isCaseInsensitive,isRegex)
    if(nargin < 3 || isempty_('isCaseInsensitive')) isCaseInsensitive = true; end
    if(nargin < 4) isRegex = false; end
    if(~isa(in1,'cell')) in1 = {in1}; end
    if(~isa(in2,'cell')) in2 = {in2}; end
    indexPerin2 = cell(numel(in2),1);
    for i=1:numel(in2)
        if(isa(in2{i},'char'))
            if(isRegex)
                fh = @regexp;
            elseif(isCaseInsensitive) 
                fh=@strcmpi; 
            else
                fh=@strcmp; 
            end
        else
            fh=@isequal;
        end
        for j=1:numel(in1)
            if(strcmp(class(in1{j}),class(in2{i})))
                if(fh(in1{j},in2{i}))
                    indexPerin2{i} = [indexPerin2{i};j];
                end
            end
        end
    end
    if(numel(indexPerin2)==1 && isempty(indexPerin2{1}))    %% Make it fully empty output.
        indexPerin2 = cell(numel(in2),0);
    end
end


%% Enlarge a numeric matrix by adding 'NaN's 
%  NOTES:
%       1. NaN values 'outSize' indicate, don't touch those dimensions.
%       1. 'isIgnoreSmallerOutSize' = true  =>  Program ignore ith dimension for which outSize(i) < size(in,i)
%   Inputs: 
%       1. 'in' must be a matrix, shouldn't contain any cells 
%       2. 'outSize' must be a column vector, and in 'Matlab size' format (Not 'Standard Size', ie. [3,1] instead of [3])
%               NaN in 'outSize' indicates "do not touch" that dimension
%       3. 'outIndicesPerDim' must be an array of cells. outIndicesPerDim{i} indicates the indices in output to consider while filling along ith dimension.
%               A cell containing ':' indicates, no reordering along corresponding dimension.
%   Warnings:
%       1. Changing default inputs affect cat_()
function out = augmentMatrix_(in,outSize,appendValue,isIgnoreSmallerOutSize,outIndicesPerDim)
    if(nargin<3 || isempty_('appendValue')) appendValue = nan; end
    if(nargin<4 || isempty_('isIgnoreSmallerOutSize')) isIgnoreSmallerOutSize = false; end
    % -Rectify inputted sizes
    outSize = convertRowVectorsToColumnVectors_(outSize);
    inSize = size(in)';
    if(numel(inSize)<numel(outSize))
        inSize = augmentMatrix_(inSize,size(outSize),1);
    elseif(numel(inSize)>numel(outSize))
        outSize = augmentMatrix_(outSize,size(inSize),1);
    end
    outSize(isnan(outSize)) = inSize(isnan(outSize));
    if(isIgnoreSmallerOutSize)
        outSize(outSize<inSize) = inSize(outSize<inSize);
    end
    if(nargin<5 || isempty_('outIndicesPerDim')) outIndicesPerDim = num2cell(repmat(':',numel(inSize),1)); end
    for iDim = 1:numel(inSize)
        if(strcmp(outIndicesPerDim{iDim},':'))
            outIndicesPerDim{iDim} = [1:inSize(iDim)];
        end
    end
    % -Verification
    assert(all(inSize==numel_(outIndicesPerDim{:})),'Error: outSize must match to the corresponding number of element in ''outIndicesPerDim'' !!')
    assert(numel(outSize)>=2,'Error: ''outSize'' must be in MATLAB Format!!');
    assert(all(inSize<=outSize),'Error: Desired out size is less than ''in'' size!!');
    % -Assigning to output
    out = appendValue*ones(outSize');
    out(outIndicesPerDim{:}) = in;
end

%% -Function to convert all row vectors to column vectors..
% -NOTES:
%      1. #Outputs = #Inputs
%      2. If input is not a vector, then don't touch it.
%      3. Inputting multi-dimensional matrix would be resulted in throwing error
% -LOG:
%      1. (2017/02/27): Splitted the existing function convertRowVectorsToColumnVectors_() into 2 function:
%           convertRowVectorsToColumnVectorsRecursively_() and convertRowVectorsToColumnVectors_()
function varargout = convertRowVectorsToColumnVectors_(varargin)
    %% -For multiple inputs, call the function recursively for each input
    if(nargin>1)
        for iInput=1:nargin
            varargout{iInput} = convertRowVectorsToColumnVectors_(varargin{iInput});
        end
    else
        varargout{1} = varargin{1};
        if(isa(varargout{1},'cell') || isa(varargout{1},'double') || isa(varargout{1},'logical'))
            if(size(varargout{1},1)==1 && size(varargout{1},2)>1) 
                varargout{1} = varargout{1}';
            end
        end
    end
end                

%% Same as MATLAB ndims, but here it returns '1' for column vector 
function nDim = ndims_(in)
    nDim = ndims(in);
    if(nDim == 2 && size(in,2)==1)
       nDim = 1; 
    end
end        


%% -Same as MATLAB built in numel function except that this one supports multiple inputs
function out = numel_(varargin)
    if(nargin>1)
       out = nan(nargin,1);
       for i=1:nargin
        out(i) = numel_(varargin{i});
       end
    else
       out = numel(varargin{1});
    end
end







