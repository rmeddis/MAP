function paramsFound=UTIL_paramsList(myWhos)
% UTIL_paramsList looks for structures with names ending in Params.
% nm=UTIL_paramsList(whos);
% for i=1:length(nm), eval(['showStruct(' nm{i} ', ''' nm{i} ''')']),end

% find structures ending with params
paramsUsedCount=0;
for i=1:length(myWhos)
    var=myWhos(i).name;
    if length(var)>5
        tag=var(end-5:end);
        if strcmp(tag,'Params')
            paramsUsedCount=paramsUsedCount+1;
            paramsUsed{paramsUsedCount}=var;
        end
    end
end



orderedList={'controlParams', 'globalStimParams', 'inputStimulusParams',...
     'OMEParams', 'DRNLParams', ...
    'IHC_cilia_RPParams', 'IHCpreSynapseParams', 'AN_IHCsynapseParams', ...
    'MacGregorMultiParams', 'MacGregorParams'};


% check that they belong to the approved list
paramsFoundcount=0;
for i=1:length(orderedList)
    for j=1:length(paramsUsed)
        usedName=paramsUsed{j};
        if strcmp(orderedList{i},paramsUsed{j})
            paramsFoundcount=paramsFoundcount+1;
            paramsFound{paramsFoundcount}=usedName;
        end
    end
end

% return the list of names

