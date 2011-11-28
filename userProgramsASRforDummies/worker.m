%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a function the makes sure all of the jobs get completed in an
% orderly fashion. It is responsible for keeping unruly nodes in line.
% NC 2010/1011

function worker(workFolder)

%unit testing
% close all; clear all; clc
% workFolder = 'D:\exps\_bar';


%main script
isloaded = 0;
while ~isloaded
    if numel(dir(fullfile(workFolder,'jobLock.txt'))) %Check to see if lock already in place
        pTime = randi(30);
        disp(['Worker function locked out -> waiting for ' num2str(pTime) ' seconds until next retry']);
        pause(pTime);
    else
        load(fullfile(workFolder,'jobObject.mat'))
        isloaded = 1;
    end
end



x=obj;
clear obj;
x.initMAP; %Need to alert it to the path

personalWork = 0;
while(any(x.todoStatus==0))        
    x.lockJobList;
    x = x.loadSelf; %Reload incase changed
    rJobs = 8+randi(8);%Grab 1st 9-16 open jobs
    todoNow = find(~x.todoStatus,rJobs,'first');     
    x.todoStatus(todoNow) = 1; %Flag it (them) as pending
    x.storeSelf; %store pending flag as quickly as possible to minimise race condition impact
    x.unlockJobList;
    disp(['Grabbed ' num2str(numel(todoNow)) ' sound files for current job.'])
        
    % ---  DO WORK  ---
    for tt=1:numel(todoNow)
        if ~numel(dir(fullfile(x.opFolder,strrep(x.wavList(todoNow(tt)).name, '.wav','.map'))));
            x.genFeat(x.wavList(todoNow(tt)).name);
        else
            disp(['File ' x.wavList(todoNow(tt)).name ' already processed'])
        end
        
    end
    % --- END OF WORK ---
    
    
    x.lockJobList;
    x = x.loadSelf; %Reload incase changed while processing (probably has)
    x.todoStatus(todoNow) = 2; %Flag as complete
    x.storeSelf; %Update as done immediately    
    x.unlockJobList;
        
    clc
    personalWork = personalWork+1;
    disp( ['This process has completed ' num2str(personalWork) ' lists'] )
    x.checkStatus        
end

disp('-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-')
disp( ' > COMPLETED CURRENT JOB' )
disp( ['  In the folder ' workFolder '  .....'] )
disp( ['  This process completed ' num2str(personalWork) ' lists'] )
disp('-*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-')


