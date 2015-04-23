% sta over time
% sandbox for direction-speed STA
load notevars.mat
% clearvars -except stimtrain_dir spktri stimtrain_spd nstimdir nstimspd seg_dur nTypes
pairs=[{'dirL_nospd'},{[1 6]},{{[2,3],[1]}};{'dirH_nospd'},{[1 6]},{{[1],[2,3]}};...
        {'spdL_nodir'},{[2 7]},{{[2,3],[1]}};{'spdH_nodir'},{[2 7]},{{[1],[2,3]}};...
        {'dirL_spdL'},{[3 9 10 11 12]},{{[3],[2],[2],[3],[2]}};{'dirH_spdH'},{[3 9]},{{[2],[3]}};... %  {'dirL_spdL'},{[3 5 9 10 11 12]},{{[3],[3],[2],[2],[3],[2]}}
        {'dirL_spdH'},{[4 8 11 12]},{{[3],[2],[2],[3]}};{'dirH_spdL'},{[4 8 10]},{{[2],[3],[3]}}]; %{'dirH_spdL'},{[4 5 8 10]},{{[2],[2],[3],[3]}}
% for varpair=1:size(pairs,1)/2
numtriTypes=size(pairs,1);

% 
% stimspk=cell(size(pairs,1),3); % dir stim, spd stim, spikes trials x time.
% for i=1:numtriTypes
%     spkseg_tmp=[];
%     for trI=1:size(pairs{i,2},2)            %trialtype index
%         tr=pairs{i,2}(trI);                 %trialtype
%         for segI=1:size(pairs{i,3}{trI},2)  %segment index
%             seg=pairs{i,3}{trI}(segI);      %segment
%             stim{i,1}=[stimspk{i,1};stimtrain_dir{1,tr}(:,(seg-1)*seg_dur+1:seg*seg_dur)];
%             stim{i,2}=[stimspk{i,2};stimtrain_spd{1,tr}(:,(seg-1)*seg_dur+1:seg*seg_dur)];
%             spkseg_tmp=[spkseg_tmp,spktri{1,tr}(seg,:)];
%         end
%     end
%     ntr=size(spkseg_tmp,2);
%     for k=1:ntr
%         stimspk.spk{i}(k,1:length(spkseg_tmp{k}))=round(spkseg_tmp{k});
%     end
% end
%
% generate dir x spd x time array for STAs
stimarray_spd=[];
stimarray_dir=[];

nTypes=length(nstimspd);
for i=1:nTypes
    stimarray_spd=[stimarray_spd,reshape(nstimspd{i},1,[])];
    stimarray_dir=[stimarray_dir,reshape(rad2deg(unwrap(deg2rad(nstimdir{i}))),1,[])];
end
%
%
for i=1:length(stimarray_dir)
    if stimarray_dir(i)<-180
        stimarray_dir(i)=stimarray_dir(i)+360;
    end
    if stimarray_dir(i)>180
        stimarray_dir(i)=stimarray_dir(i)-360;
    end
    if stimarray_dir(i)<-150
        stimarray_dir(i)=stimarray_dir(i)+180;
    end
end
nBins=50;
[h, binc_spd]=hist(stimarray_spd,nBins);
binsize_spd=binc_spd(2)-binc_spd(1);
%
[h,binc_dir]=hist(stimarray_dir,nBins);
binsize_dir=binc_dir(2)-binc_dir(1);%%
        
tic
rep=1;
tL=400;
nStimDims=2;
sta=zeros(nBins,nBins,tL,rep,numtriTypes);
STA=zeros(nBins,nBins,tL,rep);
STC=zeros(nBins,nBins,tL,rep);

%
for i=1:numtriTypes
    i
    stim=cell(1,2);
    for trI=1:size(pairs{i,2},2)            %trialtype index
        tr=pairs{i,2}(trI);                 %trialtype
        for segI=1:size(pairs{i,3}{trI},2)  %segment index
            seg=pairs{i,3}{trI}(segI);      %segment
            stim{1,1}=[stim{1,1};stimtrain_dir{1,tr}(:,(seg-1)*seg_dur+1:seg*seg_dur)];
            stim{1,2}=[stim{1,2};stimtrain_spd{1,tr}(:,(seg-1)*seg_dur+1:seg*seg_dur)];
        end
    end
    stim2d=zeros(nBins,nBins,size(stim{1,2},1),size(stim{1,2},2));
    for tr=1:size(stim{1,2},1)
        for hs=1:nBins
            inds1=find(stim{1,2}(tr,:)<binc_spd(hs)+binsize_spd/2);
            inds2=find(stim{1,2}(tr,:)>binc_spd(hs)-binsize_spd/2);
            inds_spd=intersect(inds1,inds2);
            for hd=1:nBins
                inds1=find(stim{1,1}(tr,:)<binc_dir(hd)+binsize_dir/2);
                inds2=find(stim{1,1}(tr,:)>binc_dir(hd)-binsize_dir/2);
                inds_dir=intersect(inds1,inds2);
                inds_ds=intersect(inds_dir,inds_spd);
    %             [inds_r,inds_c]=ind2sub(size(stimspk{i,1}),inds_ds);
                stim2d(hd,hs,tr,inds_ds)=1;
            end     
        end
    end
    
    spkseg_tmp=[];
    for trI=1:size(pairs{i,2},2)            %trialtype index
        tr=pairs{i,2}(trI);                 %trialtype
        for segI=1:size(pairs{i,3}{trI},2)  %segment index
            seg=pairs{i,3}{trI}(segI);      %segment
            spkseg_tmp=[spkseg_tmp,spktri{1,tr}(seg,:)];
        end
    end
    
    ntr=size(spkseg_tmp,2);
    for k=1:ntr
        spk{1,i}(k,1:length(spkseg_tmp{k}))=round(spkseg_tmp{k});
    end
    
    for rp=1:rep
        ntri=size(spk{1,i},1);
        index=randperm(ntri);
        ind_use=index(1:round(0.8*length(index)));       
        [sta1,STA1,STC1] = get_2sta(spk{1,i}(ind_use,:), stim2d(:,:,ind_use,:));
        sta(:,:,:,rp,i)=sta1;
    end
end
toc
clear stim
% save('stim2d.mat','stim2d','-v7.3');

%
% tlags=1:150;
% 
% for i=1:numtriTypes
%     for l=tlags
%         for t=1:size(stimspk{i,3},4)
%             validspks=find(stimspk{i,4});
%         end
%     end
% end

% rep=1;
% tL=400;
% nStimDims=2;
% sta=zeros(nBins,nBins,tL,rep,ntriTypes);
% STA=zeros(nBins,nBins,tL,rep);
% STC=zeros(nBins,nBins,tL,rep);
% 
% for j=2:7
%      spkseg_tmp=[];
%     for trI=1:size(pairs{j,2},2)            %trialtype index
%         tr=pairs{j,2}(trI);                 %trialtype
%         for segI=1:size(pairs{j,3}{trI},2)  %segment index
%             seg=pairs{j,3}{trI}(segI);      %segment
%             spkseg_tmp=[spkseg_tmp,spktri{1,tr}(seg,:)];
%         end
%     end
%     ntr=size(spkseg_tmp,2);
%     for k=1:ntr
%         spk{1,j}(k,1:length(spkseg_tmp{k}))=round(spkseg_tmp{k});
%     end
%     for rp=1:rep
%         ntri=size(spk{1,j},1);
%         index=randperm(ntri);
%         ind_use=index(1:round(0.8*length(index)));       
%         [sta1,STA1,STC1] = get_2sta(spk{1,j}(ind_use,:), stim2d{j,1}(:,:,ind_use,:));
%         sta(:,:,:,rp,j)=sta1;
%     end
%     j
% end
clearvars -except sta 
save('ga031015_STA.mat','sta','-v7.3')

%%
lags=[-199:200];
numtriTypes=8;

maxsta=squeeze(max(max(max(mean(sta,4),[],3),[],2)));
figure;
for t=7
    normsta(:,:,:,t)=mean(sta(:,:,:,:,t),4)./maxsta(t);

% for t=1:numtriTypes
    for i=101:1:201
        imagesc(-50:(150/50):50,-100:(200/50):100,normsta(:,:,i,t))
        title([num2str(lags(i)),', type',num2str(t)])
        colorbar
        pause(0.05)
    end
end