function UTIL_printTabTable(M, headers, format)
% printTabTable prints a matrix as a table with tabs
%headers are optional
%headers=strvcat('firstname', 'secondname')
%  UTIL_printTabTable([1 2; 3 4],strvcat('a1','a2'));

if nargin<3
    format='%g';
end

if nargin>1
    [r c]=size(headers);
    for no=1:r
        fprintf('%s\t',headers(no,:))
    end
    fprintf('\n')
end

[r c]=size(M);

for row=1:r
    for col=1:c
        fprintf('%s',num2str(M(row,col),format))
        if col<c
            fprintf('\t')
        end
    end
    fprintf('\n')
end

