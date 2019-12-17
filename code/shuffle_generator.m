
function [BINshuffs,SWAPshuffs,SHIFTshuffs,COLshuffs]=shuffle_generator(remcounts,nshuffs)

% rem_edges=0:Rem_length; % 1s time bins for histogram
% for icell=1:numel(remspikes)
%     remcounts(icell,:)=histcounts(remspikes_correg{icell},rem_edges);
% end

% shuffles the binned data n times

BINshuffs=cell(nshuffs,1);
SWAPshuffs=cell(nshuffs,1);
SHIFTshuffs=cell(nshuffs,1);
COLshuffs=cell(nshuffs,1);

for i = 1:nshuffs
    BINshuffs{i}=binshufflev1(remcounts);
    SWAPshuffs{i}=rowshuffler1(remcounts);
    SHIFTshuffs{i}=circshuffv1(remcounts);
    COLshuffs{i}=columnshufflev1(remcounts);
end
end
