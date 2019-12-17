function [scaled_CT, zedmap, shuffles] = LW_template_corr(remspikes,runspikes, SF_list, mode)

% ############ Script to run moving template correlation ################
%(Louie & Wilson, 2001)

% data required are two cell arrays, one containing spike rasters for the
% RUN period, the other containing spike rasters for the REM period. The
% unit identities and order must be preserved across the two arrays, so
% that RUN{i} is the same unit as REM{i}.

%========= Correlation parameters =======================================

if mode == 1
% SF_list=0.5:0.5:2.5; % vector of scaling factors for the RUN template
nshuffs=2;  % number of bootstraps for each shuffle type

elseif mode == 2
%     SF_list=0.5:0.5:2.5; 
% SF_list=0.3:0.1:3; % vector of scaling factors for the RUN template
nshuffs=50;  % number of bootstraps for each shuffle type
else
 error('please select mode')
end
%========= Create REM template ==========================================

% correct start point of REM to zero, get window length
REM_tmin=min(cell2mat(cellfun(@min, remspikes', 'Un',0)));
for icell=1:numel(runspikes)
    remspikes_correg{icell}=remspikes{icell}-REM_tmin;
end
% REM template dimensions
Rem_length=ceil(max(cell2mat(cellfun(@max, remspikes_correg', 'Un',0)))); % time
Ncells=numel(remspikes); % number of cells

% create firing rate histogram with 1s bin width

remcounts=zeros(Ncells,Rem_length);
rem_edges=0:Rem_length; % 1s time bins for histogram
for icell=1:numel(remspikes)
    if ~isempty(remspikes_correg{icell})
    remcounts(icell,:)=histcounts(remspikes_correg{icell},rem_edges); % convert this to "histcounts" function if not limited to v2013b!
    else
    end
end

% ======== Correlation algorithm ========================================

scaled_CT = louie_v3(remspikes,runspikes, SF_list); % run the "real" correlation pass
[BINshuffs,SWAPshuffs,SHIFTshuffs,COLshuffs]=shuffle_generator(remcounts,nshuffs);

% ======= Bootstrapping ==================================================

% code runs each shuffle group sequentially, and within each shuffle group
% executes inside a parfor loop for use on an HPC.

% binwise shuffle

parfor fake=1:numel(BINshuffs)
shuff_scaledBIN{fake}=louie_v3_shuffles(BINshuffs{fake},runspikes,SF_list);
end

% cellwise shuffle

parfor fake=1:numel(SWAPshuffs)
shuff_scaledSWAP{fake}=louie_v3_shuffles(SWAPshuffs{fake},runspikes,SF_list);
end

% column shuffle
parfor fake=1:numel(COLshuffs)
shuff_scaledCOL{fake}=louie_v3_shuffles(COLshuffs{fake},runspikes,SF_list);
end

% "shift" shuffle
parfor fake=1:numel(SHIFTshuffs)
shuff_scaledSHIFT{fake}=louie_v3_shuffles(SHIFTshuffs{fake},runspikes,SF_list);
end
delete gcp;
% ========= Z-scoring of real result against shuffled distributions ======
zed_extract = @(C, k) cellfun(@(c)c(k), C) ;


%initialise variables
    zedded_BIN=zeros(size(scaled_CT));
    zedded_SWAP=zeros(size(scaled_CT));
    zedded_COL=zeros(size(scaled_CT));
    zedded_SHIFT=zeros(size(scaled_CT));
    zedmap=zeros(size(scaled_CT));
    
for icell=1:size(scaled_CT,1)
for ibin=1:size(scaled_CT,2) %loop through bins of the corr result
    
    idx=sub2ind(size(scaled_CT),icell,ibin);
    
    % zscore for BIN shuffle
    this_BINdist=zed_extract(shuff_scaledBIN, idx);
    this_BINdist(isnan(this_BINdist))=0;
    mu=mean(this_BINdist);
    sigma=std(this_BINdist);
    zedded_BIN(idx)=(scaled_CT(idx)-mu)/sigma;
    
    % zscore for SWAP shuffle
    this_SWAPdist=zed_extract(shuff_scaledSWAP, idx);
    this_SWAPdist(isnan(this_SWAPdist))=0;
    mu=mean(this_SWAPdist);
    sigma=std(this_SWAPdist);
    zedded_SWAP(idx)=(scaled_CT(idx)-mu)/sigma;
    
    % zscore for COL shuffle
    this_COLdist=zed_extract(shuff_scaledCOL, idx);
    this_COLdist(isnan(this_COLdist))=0;
    mu=mean(this_COLdist);
    sigma=std(this_COLdist);
    zedded_COL(idx)=(scaled_CT(idx)-mu)/sigma;
    
    % zscore for SHIFT shuffle
    this_SHIFTdist=zed_extract(shuff_scaledSHIFT, idx);
    this_SHIFTdist(isnan(this_SHIFTdist))=0;
    mu=mean(this_SHIFTdist);
    sigma=std(this_SHIFTdist);
    zedded_SHIFT(idx)=(scaled_CT(idx)-mu)/sigma;
    
    zedmap(idx)=min([zedded_BIN(idx),zedded_SWAP(idx),zedded_COL(idx),zedded_SHIFT(idx)]);
end
end
zedmap=reshape(zedmap,size(scaled_CT));
shuffles.BIN=shuff_scaledBIN;
shuffles.COL=shuff_scaledCOL;
shuffles.SHIFT=shuff_scaledSHIFT;
shuffles.SWAP=shuff_scaledSWAP;
% save([ pat filesep target filesep target '_output.mat'],'zedmap','scaled_CT');
% sprintf('zscored correlation script complete')
end

