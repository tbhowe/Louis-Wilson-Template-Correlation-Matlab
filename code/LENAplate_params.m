function[Rem_RMS,Run_RMS,Xbar,Ybar,stdx,stdy]=LENAplate_params(Remwin,Runwin)

% if nargin <2
%     disp('needs rem and run templates!')
%     return
% else
% end

clear stdx stdy Ncells RemNBins RunNbins Binsum binterm
%% Define RMS of each cell rate
Ncells=numel(Runwin(:,1));
RunNBins=numel(Runwin(1,:));
RemNBins=numel(Remwin(1,:));

% for iCell=1:Ncells;
% Run_RMS(iCell,1)=rms( (Runwin(iCell,:))  );
% end

Run_RMS=rms(Runwin');

% for iCell=1:Ncells;
% Rem_RMS(iCell,1)=rms(  (Remwin(iCell,:))  );
% end

Rem_RMS=rms(Remwin');


%% Define Xbar and Ybar (see Louie & Wilson 2001 methods!)
clear Xc Binsum;

% for iCell=1:Ncells
%     
%     Xc=Rem_RMS(iCell);
%     
%     Binsum(iCell)=sum(Remwin(iCell,:)/Xc);
% end

Binsum=sum(Remwin')./Rem_RMS;
Binsum(isnan(Binsum))=0;
k=1/(Ncells*RemNBins);
Xbar = k*sum(Binsum);

clear Binsum ;

% for iCell=1:Ncells
%     
%     Yc=Run_RMS(iCell);
%     
%     Binsum(iCell)=sum(Runwin(iCell,:)/Yc);
% end
Binsum=sum(Runwin')./Run_RMS;
Binsum(isnan(Binsum))=0;
Ybar= k*sum(Binsum);

        

%% Define stdx and stdy
%check what the value of k should be in each case... it's not that obvious
%from the paper..
clear Xc binterm;

for iCell=1:Ncells
    Xc=Rem_RMS(iCell);
    binterm(iCell)=sum(   ((Remwin(iCell,:)/Xc) - Xbar).^2);
end

binterm(isnan(binterm))=0;
stdx=sqrt( (1/(Ncells*RemNBins)) * sum(binterm));

clear Xc Yc binterm;

for iCell=1:Ncells
    Yc=Run_RMS(iCell);
    binterm(iCell)=sum(   ((Runwin(iCell,:)/Yc) - Ybar).^2);
end
binterm(isnan(binterm))=0;

stdy=sqrt( (1/(Ncells*RemNBins)) * sum(binterm));
