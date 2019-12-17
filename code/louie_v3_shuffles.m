function  scaled_CT = louie_v3_shuffles(remcounts,runspikes, SF_list)

% ######### temporary variables (remove once functionalised) ############


% ######### correct time offset from zero in both windows ###############

% correct start point of RUN window to zero;
Ncells=size(remcounts,1);
run_tmin=min(cellfun(@min, runspikes'));
for icell=1:numel(runspikes)
    runspikes_correg{icell}=runspikes{icell}-run_tmin;
    if min(min(diff(runspikes{icell})))
    mindiff(icell)=min(min(diff(runspikes{icell})));
    else
        mindiff(icell)=NaN;
    end
end
 mindiff=min(mindiff(~isnan(mindiff))); % minimum inter-spike interval
 est_Fs=10^(floor(log10(mindiff)));
run_tmax=ceil(max(cellfun(@max,runspikes_correg')));


% create smoothing kernel
tb_sec= 0:run_tmax; % timebase in seconds
sigma=1.5;
bw           = mean(diff(tb_sec));
edges       = -3*sigma:bw:3*sigma; 
mu_    = 0; 
kernel = exp(-0.5 * ((edges - mu_)./(sigma/bw)).^2) ./ (sqrt(2*pi) .* (sigma/bw));
kernel = kernel*bw; % Scale kernel integral to 1 

% IFR of REM window:
for iCell=1:numel(remcounts(:,1))
    currentcell=remcounts(iCell,:);
    PDF    = conv(currentcell,kernel);
    kernal_center = ceil(length(edges)/2);
 remSMOO(iCell,:)=PDF(kernal_center:length(PDF)-kernal_center); % smoothed REM window
end

% ########## template correlation #####################################

% some variables for C(t)
RemNBins=numel(remSMOO(1,:)); % number of bins in REM window
k=1/(Ncells*RemNBins); % constant for Ct algorithm

RemDim=size(remSMOO);
tb_run=0:0.25:run_tmax;
run_Fs=mean(diff(tb_run));



%cycles timestep, starting at half length of Remwin/SF



% initialise scaled CT variable
scaled_CT=zeros(numel(SF_list),numel(tb_run));
% ======= LOOP THROUGH SCALING FACTORS =======================
for iscale=1:numel(SF_list)
     SF=SF_list(iscale);
    halfrange=0.5*(RemDim(2)/SF);
ts_init=FindClosestIndex(tb_run,halfrange);
    

    
Ct=zeros(size(tb_run)); % initialise C(t)
    % ======= LOOP THROUGH TIMEBASE ==============================
tic
for count=1:numel(tb_run)

if count < ts_init
    Ct(count)=0; % positions before this point are invalid!
    runcounts=zeros(size(remcounts));
else
    rwin_mintime=tb_run(count)-halfrange;
    rwin_maxtime=tb_run(count)+halfrange;
    run_edges=linspace(rwin_mintime,rwin_maxtime,RemDim(2)+1);
    
    runcounts=zeros(size(remcounts));
    for icell=1:numel(runspikes)
%         thiscell=runspikes_correg{icell};
%         thiscell=thiscell(thiscell >= rwin_mintime & thiscell < rwin_maxtime);
%         runcounts(icell,:)=histcounts(thiscell, run_edges);
runcounts(icell,:)=histcounts(runspikes_correg{icell}(runspikes_correg{icell} >= rwin_mintime & runspikes_correg{icell} < rwin_maxtime), run_edges);
    end
end

% make iFR of run window

runSMOO=zeros(size(remSMOO));
    for iCell=1:numel(runcounts(:,1))
    currentcell=runcounts(iCell,:);
    PDF    = conv(currentcell,kernel);
    kernal_center = ceil(length(edges)/2);
 runSMOO(iCell,:)=PDF(kernal_center:length(PDF)-kernal_center);
    end
 
 % get template correlation parameters
[Rem_RMS,Run_RMS,Xbar,Ybar,stdx,stdy]=LENAplate_params(remSMOO,runSMOO);

SumBin=zeros(Ncells,1);

% -------------------- execute C(t) for this window pair ---------------
for iCell=1:Ncells
%     clear evaluated
    Xc=Rem_RMS(iCell);
    Yc=Run_RMS(iCell);
    evaluated=zeros(1,RemNBins);
%     for iBin=1:RemNBins
%         
%     evaluated(iBin)=   ((Remwin(iCell,iBin)/Xc) - Xbar)*((Runwin(iCell,iBin)/Yc) -Ybar);
%     end

%vectorisation of the iBin forloop. No fucking difference!
    evaluated= ( (remSMOO(iCell,:)./Xc)-Xbar ).* ( ( runSMOO(iCell,:)./Yc) - Ybar);
    evaluated(isnan(evaluated))=0;
    SumBin(iCell)=sum(evaluated);
end
SumBin(isnan(SumBin))=0;
SumCell=sum(SumBin);
% clear Sumbin
Ct(count)= (k*SumCell)/(stdx*stdy);

%============= end timestamp loop =============================
end
%============= end scaling factor loop ========================
scaled_CT(iscale,:)=Ct;
end


end
