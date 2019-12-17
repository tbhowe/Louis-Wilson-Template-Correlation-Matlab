% ############ Script to run moving template correlation ################
%(Louie & Wilson, 2001)

% data required are two cell arrays, one containing spike rasters for the
% RUN period, the other containing spike rasters for the REM period. The
% unit identities and order must be preserved across the two arrays, so
% that RUN{i} is the same unit as REM{i}.

%========= Data loading =================================================


target='r3test'; % name of experiment folder

% determine whether script is running on native machine (PC) or HPC (Linux)
% and assign working directories accordingly

if ispc
    
    pat = 'C:\CODE\LouieV3test\data';
    cd 'C:\PROJECTS\LouieV3\code';
else
    
    path(path,'/panfs/panasas01/phph/th17624/Louie_v3/functions') % check this!
    home = getenv('HOME');
    cd ([home '/Louie_v3/functions'])
    pat = [home '/Louie_v3/data'];
    
end
if ~exist([pat filesep target],'dir')
    mkdir([pat filesep target])
end

load([ pat filesep target filesep 'runspikes.mat']);
load([ pat filesep target filesep 'remspikes.mat']);
%========= mode selection ================================================

mode=2; % 1= debug, 2= normal run

% DEBUG (1) runs only a few bootstrap repetitions, and has low resolution
% in the time-scale dimension

% WORKER (2) runs a full complement of bootstraps, with a high resolution
% in the time-scale dimension

%========= Correlation parameters =======================================
SF_list=0.5:0.5:2.5; % vector of scaling factors for the RUN template

%========  run correlation ==============================================
[scaled_CT, zedmap, shuffles] = LW_template_corr(remspikes,runspikes, SF_list, mode);

 save([ pat filesep target filesep target '_output.mat'],'zedmap','scaled_CT','shuffles');
% ========  BASIC PLOTTING ===============================================
figure
imagesc(scaled_CT)
colormap jet
caxis([0 1])
xlabel('time')
ylabel('SF')
yticklabs=linspace(min(SF_list),max(SF_list),5);
ylabels={};
for i=1:numel(yticklabs)
    ylabels{i}=num2str(yticklabs(i));
end
yticko=linspace(1, size(scaled_CT,1),5);
yticks(yticko)
yticklabels(ylabels);
set(gca,'TickDir','out'); 
set(gca,'box','off')
title('z-scored C(t,SF)')
  set(gcf, 'Position',  [100, 100, 1600, 400])
