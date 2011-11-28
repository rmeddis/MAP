function plotProfile(fgName, bgName)

addpath (['..' filesep 'profiles'])

%% plot profile
if nargin<1
    fgName = 'profile_JSAN_R';
    bgName = '';
end

cmd=['foreground = ' fgName ';'];
eval(cmd)

if nargin==2
    cmd=['background = ' bgName ';'];
    eval(cmd)
else
    bgName='';
end

% absolute thresholds
figure(90), clf
set(gcf, 'name', 'Profile')
subplot(2,1,2)
semilogx(foreground.BFs,foreground.LongTone,'ko-','lineWidth',2); hold on
semilogx(foreground.BFs,foreground.ShortTone,'bo-','lineWidth',2); hold on
if ~isempty(bgName)
    semilogx(background.BFs,background.LongTone,'ko:'); hold on
    semilogx(background.BFs,background.ShortTone,'bo:'); hold on
end
ylim([0 100])

% TMC
for BFno=1:length(foreground.TMCFreq)
    subplot(2,6,BFno)
    % SL
% plot(foreground.Gaps,foreground.TMC(BFno,:)-foreground.LongTone(BFno),'r','lineWidth',3), hold on
    plot(foreground.Gaps,foreground.TMC(BFno,:),'b','lineWidth',3), hold on
    ylim([-10 110])
    xlim([0.01 0.1])
%     grid on
    if BFno==1
        ylabel('masker dB SL')
        xlabel('    gap (s)')
    end
    title([num2str(foreground.TMCFreq(BFno)) ' Hz'])
    set(gca,'XTick',[ 0.02:0.02:0.1],'xTickLabel', { '', '0.04', '', '',  '0.1'})
end

if ~isempty(bgName)
    for BFno=1:length(background.TMCFreq)
        BF = background.TMCFreq(BFno);
        idx = find(BF == foreground.TMCFreq);
        if ~isempty(idx);
            
            subplot(2,6,idx)
            % SL
% plot(background.Gaps,background.TMC(BFno,:)-background.LongTone(BFno),'k:')
            plot(background.Gaps,background.TMC(BFno,:),'k:')
            ylim([-10 110])
            xlim([0.01 0.1])
        end
    end
end

% IFMCs
for BFno=1:length(foreground.IFMCFreq)
    freq=foreground.MaskerRatio'*foreground.IFMCFreq(BFno);
    subplot(2,1,2)
    semilogx(freq,foreground.IFMCs(BFno,:),'r','lineWidth',3), hold on
    ylim([-20 100])
    xlim([100 12000])
%     grid on
end
xlabel('frequency (Hz)')
ylabel('masker dB / probe dB')
set(gca,'XTick',foreground.IFMCFreq)
set(gca,'Ytick', [-20 0 50 100])

if ~isempty(bgName)
    for BFno=1:length(background.IFMCFreq)
        freq=background.MaskerRatio'*background.IFMCFreq(BFno);
        subplot(2,1,2)
        semilogx(freq,background.IFMCs(BFno,:),'k:')
        ylim([0 100])
        xlim([100 12000])
    end
end
set(get(gca,'title'),'interpreter','None')

title([fgName ' / ' bgName])
% mydate=datestr(now); idx=findstr(':',mydate); mydate(idx)='_';

% fileName= ['savedData/' mydate ];
% 
% save (fileName)
% set(gcf,'name', mydate)
% disp(fileName)
