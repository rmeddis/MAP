function UTIL_showStructureSummary(structure, name, maxNoArrayValues)
% showStructureSummary prints out the values of a single structure
% The header is the structure name and each row is a field
% e.g. showStructureSummary(params,'params')
% This not the same as 'UTIL_showstruct'


if nargin<3
	maxNoArrayValues=20;
end

fprintf('\n%s:', name)

fields=fieldnames(eval('structure'));
% for each field in the structure
for i=1:length(fields)
	y=eval([ 'structure.' fields{i}]);
	if isstr(y),
		% strings
		fprintf('\n%s=\t''%s''',  fields{i},y)
	elseif isnumeric(y)
		% arrays
		if length(y)>1
			% vectors
			[r c]=size(y);
			if r>c, y=y'; end

			[r c]=size(y);
			if r>1
				%                 fprintf('\n%s.%s=\t%g x %g matrix',name, fields{i}, r, c)
				fprintf('\n%s=\t%g x %g matrix',fields{i}, r, c)

			elseif c<maxNoArrayValues
				%                     fprintf('\n%s=\t[%s]',  fields{i},num2str(y))
				fprintf('\n%s=',  fields{i})
				fprintf('\t%g',y)

			else
				fprintf('\n%s=\t %g...   [%g element array]', fields{i}, y(1),c)
			end
		else
			% single valued arrays
			%             fprintf('\n%s.%s=\t%s;', name, fields{i},num2str(y))
			fprintf('\n%s=\t%s', fields{i},num2str(y))
		end
	elseif iscell(y)
		fprintf('\n%s=\t cell array', fields{i})

	elseif isstruct(y)
		fprintf('\n%s=\t structure', fields{i})
	end,

end,
fprintf('\n')

