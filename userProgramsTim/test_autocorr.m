addpath('..\utilities');
addpath('..\MAP')

%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=1;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=1;               % shows SACF (probability only)
showMapOptions.showEfferent=1;          % tracks of AR and MOC
showMapOptions.surfProbability=1;       % 2D plot of HSR response 
showMapOptions.surfSpikes=1;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk

UTIL_showMAP(showMapOptions, {})
