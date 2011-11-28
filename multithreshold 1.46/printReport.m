function printReport(fileName, printTracks)

% End of run report (no args)
% *or*
% reprint previous report from file

global experiment stimulusParameters betweenRuns withinRuns statsModel
global LevittControl expGUIhandles
global paramChanges

% NB all globals are used even though this is not obvious
global inputStimulusParams OMEParams DRNLParams
global IHC_VResp_VivoParams IHCpreSynapseParams  AN_IHCsynapseParams
global MacGregorParams MacGregorMultiParams  filteredSACFParams
global experiment % used by calls from multiThreshold only
global IHC_cilia_RPParams

printReportGuide.structures=1;
printReportGuide.showPsychometric=0;
printReportGuide.HorizontalTracks=1;


if nargin==0
    % print new report
    printReportGuide.showTracks=experiment.printTracks;
    printReportGuide.fileName=[];
    % save this data (just in case)
    saveFileName=['savedData/mostRecentResults'];
    experiment.minElapsed=etime(clock, betweenRuns.timeNow)/60;
    save(saveFileName, 'experiment', 'stimulusParameters',...
        'betweenRuns', 'withinRuns', 'statsModel', 'expGUIhandles')
    disp(['data saved as: ' saveFileName]);

else
    % reprint request (i.e print out old data)
    printReportGuide.fileName=fileName;
    load(printReportGuide.fileName)
    saveFileName=printReportGuide.fileName;
    if nargin>1
        printReportGuide.showTracks=printTracks;
    else
        printReportGuide.showTracks=experiment.printTracks;
    end
end

fprintf('******** multithreshold ')

x=pwd; disp(['version ' x(end-3:end)])
fprintf('\nName:\t%s ', experiment.name)
fprintf('\nparadigm:\t%s ', experiment.paradigm)
fprintf('\nEar:\t%s ', experiment.ear)

method=experiment.threshEstMethod;
if stimulusParameters.includeCue && ~strcmp(method(1:6),'2I2AFC')
    method=[method '/ withCue'];
end
fprintf('\nmethod:\t%s ', method)
fprintf('\ndate:\t%s ', experiment.date)
fprintf('\n\n')

if isempty(betweenRuns.thresholds)
    disp('no thresholds found')
end

% prepare results as matrices ready to print
%  first sort the actual sequence into a more readable sequence
[idx1, idx2, var1values, var2values]=...
    sortVariables(betweenRuns.variableList1, betweenRuns.variableList2,...
    betweenRuns.var1Sequence, betweenRuns.var2Sequence);

header1=betweenRuns.variableName1;
header2=betweenRuns.variableName2;
header1 = strrep(header1, 'none', ' '); % none is not a useful header
header2 = strrep(header2, 'none', ' '); % none is not a useful header
headers=strvcat([header1 '/'], header2);

disp('thresholds')
global resultsTable
resultsTable= sortTablesForPrinting(idx1,idx2,...
    var1values,var2values, betweenRuns.thresholds);
msg=printTabTable(resultsTable,  headers);
addToMsg(msg,0)
if ~isempty(paramChanges)
    fprintf('\n')
    fprintf('%s\n', char(paramChanges))
end

% sort tracks into the same order
betweenRuns.levelTracks=betweenRuns.levelTracks(idx1);
betweenRuns.responseTracks=betweenRuns.responseTracks(idx1);
betweenRuns.bestThresholdTracks=betweenRuns.bestThresholdTracks(idx1);
betweenRuns.levelTracks=betweenRuns.levelTracks(idx2);
betweenRuns.responseTracks=betweenRuns.responseTracks(idx2);
betweenRuns.bestThresholdTracks=betweenRuns.bestThresholdTracks(idx2);

if printReportGuide.structures
    maxNoArrayValues=30;
    showStructureSummary(stimulusParameters, 'stimulusParameters', ...
        maxNoArrayValues)
    showStructureSummary(experiment, 'experiment',maxNoArrayValues)
    showStructureSummary(betweenRuns, 'betweenRuns',maxNoArrayValues)
    showStructureSummary(withinRuns, 'withinRuns')

    switch experiment.threshEstMethod
        case {'2I2AFC++', '2I2AFC+++'}
            showStructureSummary(LevittControl, 'LevittControl', ...
                maxNoArrayValues)
    end

    switch experiment.ear
        case {'statsModelLogistic','statsModelRareEvent'}
            showStructureSummary(statsModel, 'statsModel', ...
                maxNoArrayValues)
    end
end

if printReportGuide.showTracks
    % NB this procedure can only be used if all the tracks are present and
    % of equal length
    bigTable=[]; header=[];
    disp(' '); 
    disp('Leveltracks starting from 1 response before the first reversal')
    for i=1:length(betweenRuns.levelTracks)
        if printReportGuide.HorizontalTracks
            printTabTable(betweenRuns.levelTracks{i});
        end
        header=strvcat(header, 'level');
    end

    disp(' '); 
    disp('Response tracks  from 1 response before the first reversal')
    for i=1:length(betweenRuns.responseTracks)
        if printReportGuide.HorizontalTracks
            printTabTable(betweenRuns.responseTracks{i});
        end
        header=strvcat(header, 'resp');
    end

    disp(' '); disp('threshold tracks starting from the first reversal')
    for i=1:length(betweenRuns.bestThresholdTracks)
        if printReportGuide.HorizontalTracks
        end
        printTabTable(betweenRuns.bestThresholdTracks{i});
        header=strvcat(header, 'mean');
    end

end

switch experiment.ear
    case  {'MAPmodelMultiCh', 'MAPmodelSingleCh'}
        % show all parameters but do not compute the model
        nm=UTIL_paramsList(whos);
        for i=1:length(nm)
            try
                eval(['UTIL_showStruct(' nm{i} ', ''' nm{i} ''')'])
            catch
            end
        end
end

if experiment.saveData
    fprintf('\n')
    disp('To reprint this report with tracks use:')
    disp([ 'printReport(''' saveFileName ''',1)'])
end

% print final summary (repeat of above)
fprintf('\n')
fprintf('\n')
disp('thresholds')
msg=printTabTable(sortTablesForPrinting...
    (idx1,idx2,var1values,var2values, betweenRuns.thresholds),  headers);
addToMsg(msg,0)
fprintf('\n')

if length(var1values)==1 && length(var2values)==1 ...
        && experiment.maxTrials>49
    [psy, levelsBinVector, binFrequencies]= ...
        psychometricFunction(withinRuns.levelsPhaseTwo,...
        withinRuns.responsesPhaseTwo, experiment.psyBinWidth);
    disp('Psychometric function')
    fprintf(' level  \tfreq\tprob\n')
    fprintf('%6.0f\t%6.2f\t%6.0f\n',[levelsBinVector; binFrequencies; psy])

    switch experiment.threshEstMethod
        % only one value required for level change
        case {'MaxLikelihood',  'oneIntervalUpDown'};
%             fprintf('\n')
%             fprintf('k \t %6.2f\n',logistic.bestK)
%             fprintf('g  \t%7.5f\n',rareEvent.bestGain)
    end
    fprintf('\n')

end

fprintf('\nparadigm:\t%s\n ', experiment.paradigm)
if ~isempty(paramChanges)
    fprintf('%s\n', char(paramChanges))
end

% ------------------------------------------------- sortTablesForPrinting
function table= sortTablesForPrinting(idx1,idx2, var1values,var2values, x)
% table converts a vector to a table
% after sorting according to idx1 and idx2
% the table is completed by adding row and column values
x=x(idx1);
x=x(idx2);
xMatrix=reshape(x,length(var1values),length(var2values));

table=[[-1000 var2values]; [var1values' xMatrix]];

% --------------------------------------------------- showStructureSummary
function showStructureSummary(structure, name, maxNoArrayValues)
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
                % fprintf('\n%s.%s=\t%g x %g matrix',name, fields{i}, r, c)
                fprintf('\n%s=\t%g x %g matrix',fields{i}, r, c)

            elseif c<maxNoArrayValues
                %  fprintf('\n%s=\t[%s]',  fields{i},num2str(y))
                fprintf('\n%s=',  fields{i})
                fprintf('\t%g',y)

            else
                fprintf('\n%s=\t %g...   [%g element array]', ...
                    fields{i}, y(1),c)
            end
        else
            % single valued arrays
            % fprintf('\n%s.%s=\t%s;', name, fields{i},num2str(y))
            fprintf('\n%s=\t%s', fields{i},num2str(y))
        end
    elseif iscell(y)
        fprintf('\n%s=\t cell array', fields{i})

    elseif isstruct(y)
        fprintf('\n%s=\t structure', fields{i})
    end,

end,
fprintf('\n')


% ------------------------------------------------------- printTabTable
function strings= printTabTable(M, headers)
% printTabTable prints a matrix as a table with tabs
%headers are optional
%headers=strvcat('firstname', 'secondname')
%  printTabTable([1 2; 3 4],strvcat('a1','a2'));
stringCount=1; strings{stringCount}=[];

if nargin>1
    [r c]=size(headers);
    for no=1:r
        % print all headers in a row
        fprintf('%s\t',headers(no,:))
        strings{stringCount}=...
            sprintf('%s\t',headers(no,:)); stringCount=stringCount+1;
    end
    fprintf('\n')
end

[r c]=size(M);

for row=1:r
    string=[];
    for col=1:c
        if row==1 && col==1 && M(1,1)==-1000
            %   Print nothing (tab follows below)
        else
            fprintf('%s',num2str(M(row,col)))
            string=[string ' ' sprintf('%s',num2str(M(row,col)))];
        end
        if col<c
            fprintf('\t')
            %strings{stringCount}=sprintf('\t'); stringCount=stringCount+1;
        end
    end % col
    strings{stringCount}=string; stringCount=stringCount+1;
    fprintf('\n')
end % row

% ------------------------------------------------------- sortVariables
function [idx1, idx2, var1values, var2values]= ...
    sortVariables(var1values, var2values, var1Sequence, var2Sequence)

[x idx1]= sort(var1Sequence);
var1Sequence= x;
var2Sequence= var2Sequence(idx1);
depVarName= 'th';

[x idx2]=sort(var2Sequence);
var2Sequence=x;
var1Sequence=var1Sequence(idx2);

var1values=sort(var1values);
var2values=sort(var2values);
