function UTIL_showStruct(structure, name, valuesOnly, terminator)
% showStruct prints out the values in a single structure
% using the format <name> '.' <fieled> ' = ' <value><terminator>
% e.g. 
% showStruct(params,'params')         % standard usage
%showStruct(params,'params', 0, ';')  % adds final comma
%showStruct(params,'params', 1)      % omits structure and field names

if nargin<3, valuesOnly=0; end
if nargin<4, terminator=[]; end

fields=fieldnames(eval('structure'));
for i=1:length(fields)
    % y is the contents of this field
    y=eval([ 'structure.' fields{i}]);
    if ischar(y),
        % strings
        if valuesOnly
            fprintf('\n%s', y)
        else
            fprintf('\n%s.%s=\t''%s%s''', name, fields{i},y,terminator)
        end
    elseif isnumeric(y)
        % arrays
        if length(y)>1
            % matrices and vectors
            [r c]=size(y);
            if r>c, y=y'; end   % make row vector from column vector
            [r c]=size(y);
            
            if r>10
                % large matrix
                fprintf('\n%s.%s=\t%g x %g matrix',name, fields{i}, r, c)
                
            elseif r>1
                % small matrix
                for row=1:r
                    fprintf('\n%s.%s(%1.0f)=\t%s;', name, fields{i},row, num2str(y(row,:)))
                end
                
            elseif c>20
                % long row vector
                fprintf('\n%s.%s=\t %g...   [%g element array]%s',name, fields{i}, y(1),c, terminator)
                
            else
                fprintf('\n%s.%s=\t[%s]%s', name, fields{i},num2str(y),terminator)
            end
        else
            % single valued arrays
            if valuesOnly
                fprintf('\n%s%s', num2str(y), terminator)
            else
                fprintf('\n%s.%s=\t%s%s', name, fields{i},num2str(y), terminator)
            end

        end % length (y)
    elseif iscell(y)
        fprintf('\n%s.%s=\t cell array', name, fields{i})
            
    elseif isstruct(y)
        fprintf('\n%s.%s=\t structure', name, fields{i})
    end     % isstr/ numeric

end         % field
fprintf('\n')