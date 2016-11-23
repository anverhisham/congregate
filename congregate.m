%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -BRIEF DESCRIPTION:
%   
%
%%%% -DETAILED DESCRIPTION:
%       
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1.  
%           s(1).a=[1,2,3];
%           s(2).a=[11,12];
%           congregate('s(:).a')
%
%%%% -NOTES:
% 
%
%%%% -TODO:
%       1. Handle single element input
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%
% PROGRAMMING NOTES:
%   1. One way to check if this function works or not, is to execute your input by replacing all variable indices with '1's 
%       i.e.  s{1}.s1(1).s2{1}.a
%   2. 'isFillWithNan' = true (Default) means, absent fields would be regarded as NaN 
%   3. Here number of dimensions = number of brackets appearing in 'inString' + dimension of last field  
%   4. Single ':' in structExprStr implies output is a column vector of cells  
%   5. s{[1,2;3,4],[1,2]}.a will result in 1D output
function out = congregate(inString,isFillWithNan,recursiveCallingScope)
   if(isempty_('isFillWithNan')) isFillWithNan = true; end
   if(isempty_('recursiveCallingScope')) recursiveCallingScope = 0; end
   %% -USER calling
   if(~recursiveCallingScope)    
       inString = strtrim(inString);
       %% Fetch variables from USER scope
       [variableNames,remainder] = regexp(inString,'(?<=^|\(|\[|\{|,|:)[a-zA-Z]\w*','match','split'); %% variableNames{1} is base object
       remainder(1)=[];
       for i=1:numel(variableNames)
          userVariables{i} = evalin('caller',variableNames{i}); 
          userVariablesNewname{i} = ['userVariables{',num2str(i),'}'];
       end
       inStringNew = strjoin_(userVariablesNewname,remainder);
       out = congregate(inStringNew,isFillWithNan,recursiveCallingScope+1);
       out = cellOfCell2mat_(out,2);
   %% -Recursive calling    
   else
       userVariables = evalin('caller','userVariables');
       % Get 's(1).s1(p)' from 's(1).s1(p).s2(q).a'
       [from, to] = findComplicatedParentheses_(inString,1);
       %  TODO_Verify(2016/11/3):   Removed && to~=numel(inString) from below expression
       if(~isempty(to) && ~isempty_(inString(1:to)))  %% TODO_Verify(2016/11/3): Changed from isempty_(['[',inString(1:to),']']) to isempty_(['inString(1:to)'])
           tmp = eval(['{',inString(1:to),'}']);   %% Changing the base object
           if(numel(tmp)>1)
               userVariables{1} = tmp;                              %% To ensure that we don't forget to place NaN for any missing field when isFillWithNan=false
           else
               userVariables{1} = eval(['[',inString(1:to),']']);   %% Changing the base object                       
           end
           out = cell(size(userVariables{1}));
           for i=1:numel(out)
               if(isa(userVariables{1},'cell'))
                inStringNew = ['userVariables{1}{',num2str(i),'}',inString(to+1:end)];
               else
                inStringNew = ['userVariables{1}(',num2str(i),')',inString(to+1:end)];
               end
               out{i} = congregate(inStringNew,isFillWithNan,recursiveCallingScope+1);
           end
       else
          if(recursiveCallingScope==1 && ~isempty(to) && isempty_(inString(1:to)))  %% TODO_Verify(2016/11/3): Changed from isempty_(['[',inString(1:to),']']) to isempty_(['inString(1:to)'])
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
       end
   end            
end   


        %% Function to find '([1,2,3])' in 's(1).s1([1,2,3])' (Note from & to corresponds to bracket location)
        %       '(1,2)', '('abcd')' are not complicated parentheses.
        %   'count' determines how many outputs needed. If count>1, then an array of cells is outputted   
        function [from,to] = findComplicatedParentheses_(str,count)
            if(isempty_('count')) count = Inf; end
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
           if(regexp(inputString,':\d*}(?=\s*$)'))      %% For 'a{:}', we need to ensure that, this function is not tricked with a{1} alone.
            to = regexp(inputString,'}(?=\s*$)');
            from = findMatchingParentheses_(inputString,to);
            if(findComplicatedParentheses_(inputString(from:to),1))
                inputString = ['{',inputString,'}'];
            end
           end
           try
               in = evalin('caller',inputString);
               out = isempty(in);
           catch
               out = true;
           end
        end
        
        
        
                
        %% {{[1,2],[3,4,5]},{[11,12]'}}  ->  []  
        %   -'squeezeLevel'
        %       0: No Squeezing. Expect redundant dimensions for row-vectors   
        %       1: 'in' (which is a matrix of cells) and all cells within is (in all scope) are squeezed 
        %       2: Both 'in' and all the numerical matrix content inside any cells are squeezed 
        function out = cellOfCell2mat_(in,squeezeLevel)
            if(isempty_('squeezeLevel'))
                squeezeLevel = 2;
            end
            if(isa(in,'cell'))
                outTemp = cell(size(in));
                if(squeezeLevel>0)
                    outTemp = squeezeAndTransposeToRowVector_(outTemp);
                end
                if(isColumnVector_(outTemp))
                    offsetDim = 1;
                else
                    offsetDim = numel(size(outTemp));
                end
                for iCell = 1:numel(in)
                    outTemp{iCell} = cellOfCell2mat_(in{iCell},squeezeLevel);
                    outTemp{iCell} = shiftdim(outTemp{iCell},-offsetDim);
                end
                out = cell2mat_(outTemp);
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
        %   findMatchingParentheses('a(bcd)',2)  ->  6  
        %   findMatchingParentheses('a(bcd)',6)  ->  1  
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
            optionIndex = contains_(varargin,'-appendValue',true);
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
                in{i} = enlarge_(varargin{i},desiredSizePerDim',appendValue);
            end
            out = cat(dim,in{:});
        end
        
        
        %% Enlarge a numeric matrix by adding 'NaN's 
        %   Warning: 'in' must be a matrix, shouldn't contain any cells 
        function out = enlarge_(in,outSize,appendValue)
            if(isempty_('appendValue')) appendValue = nan; end
            out = appendValue*ones(outSize);
            exString = '';
            for iDim=1:numel(size(in))
                exString = [exString,'1:',num2str(size(in,iDim)),','];
            end
            for iDim=numel(size(in))+1:numel(outSize)
                exString = [exString,'1,'];            
            end
            exString(end) = [];                 %% Removing trailing ','
            eval(['out(',exString,')=in;']);    %% This is the limiting part, when we try to enlarge matrix of cells 
        end
 
        
        %% Function to compare different types of objects
        %   'indexPerin2' is an array of cells, with indexPerin2{i}=indices => in1{indices}==in2{i}  
        function indexPerin2 = contains_(in1,in2,isCaseInsensitive,isRegex)
            if(isempty_('isCaseInsensitive')) isCaseInsensitive = false; end
            if(isempty_('isRegex')) isRegex = true; end
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
        








