% sta over time
% sandbox for direction-speed STA
% load notevars.mat
clearvars -except experiment pairs stimtrain_dir spktri stimtrain_spd nstimdir nstimspd seg_dur nTypes numtriTypes 
% for varpair=1:size(pairs,1)/2

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
    stimarray_dir=[stimarray_dir,reshape(nstimdir{i},1,[])];
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
% nBins=60;
% [h, binc_spd]=hist(stimarray_spd,nBins);
% binsize_spd=binc_spd(2)-binc_spd(1);
% %
% [h,binc_dir]=hist(stimarray_dir,nBins);
% binsize_dir=binc_dir(2)-binc_dir(1);
%%
        
tic
rep=1;
tL=400;
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
    nBins=60;
    [hspd(i,:), binc_spd(i,:)]=hist(stim{1,2},nBins);
    binsize_spd=binc_spd(i,2)-binc_spd(i,1);
    %
    [hdir(i,:),binc_dir(i,:)]=hist(stim{1,1},nBins);
    binsize_dir=binc_dir(i,2)-binc_dir(i,1);
    
    stim2d=zeros(nBins,nBins,size(stim{1,2},1),size(stim{1,2},2));
    for tr=1:size(stim{1,2},1)
        for hs=1:nBins
            inds1=find(stim{1,2}(tr,:)<binc_spd(i,hs)+binsize_spd/2);
            inds2=find(stim{1,2}(tr,:)>binc_spd(i,hs)-binsize_spd/2);
            inds_spd=intersect(inds1,inds2);
            for hd=1:nBins
                inds1=find(stim{1,1}(tr,:)<binc_dir(i,hd)+binsize_dir/2);
                inds2=find(stim{1,1}(tr,:)>binc_dir(i,hd)-binsize_dir/2);
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
        sta(:,:,:,rp,i)=-sta1;
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
clearvars -except pairs sta binc_dir binc_spd experiment
save([experiment,'_STA.mat'],'stdlev','sta','binc_dir','binc_spd','hdir','hspd','pairs','-v7.3')



%% if just loading sta
% load([experiment,'_linfilt_spd.mat'],'stdlev')
lags=[-199:200];
numtriTypes=size(sta,5);
numBins=size(sta,1);
maxsta=squeeze(max(max(max(mean(sta,4),[],3),[],2)));
minsta=squeeze(min(min(min(mean(sta,4),[],3),[],2)));
tsd=repmat(binc_dir',1,numBins);
tss=repmat(binc_spd,numBins,1);
colors=distinguishable_colors(8);
h1=figure;
h2=figure;
h3=figure;
for i=1:numtriTypes
    for t=1:size(sta,3)
        sta_spd_exp(:,t,i)=sum(sta(:,:,t,1,i).*tss,1);
        sta_dir_exp(:,t,i)=sum(sta(:,:,t,1,i).*tsd,2);
        sta_dir(t,i)=sum(sta_dir_exp(:,t,i));
        sta_spd(t,i)=sum(sta_spd_exp(:,t,i));
    end
%     if abs(min(min(sta_dir)))>abs(max(max(sta_dir)))...
%             ||abs(min(min(sta_spd)))>abs(max(max(sta_spd)))
%         sta=-sta;
%         for t=1:size(sta,3)
%             sta_spd_exp(:,t,i)=sum(sta(:,:,t,1,i).*tss,1);
%             sta_dir_exp(:,t,i)=sum(sta(:,:,t,1,i).*tsd,2);
%             sta_dir(t,i)=sum(sta_dir_exp(:,t,i));
%             sta_spd(t,i)=sum(sta_spd_exp(:,t,i));
%         end
%     end
    figure(h1);
    subplot(numtriTypes/2,2,i)
    imagesc(lags,binc_dir,sta_dir_exp(:,:,i))
    xlabel('time');ylabel('direction')
    colorbar
    title(pairs{i,1},'Interpreter','none')
    figure(h2);
    subplot(numtriTypes/2,2,i)
    imagesc(lags,binc_spd,sta_spd_exp(:,:,i))
    colorbar
    title(pairs{i,1},'Interpreter','none')
    xlabel('time');ylabel('speed')
    title(pairs{i,1})
    figure(h3);subplot 211
    plot(lags,sta_dir(:,i),'Color',colors(i,:),'LineWidth',2);hold all
    figure(h3);subplot 212
    plot(lags,sta_spd(:,i),'Color',colors(i,:),'LineWidth',2);hold all
    
end


figure(h1)
suptitle('direction STA over time')
figure(h2)
suptitle('speed STA over time')
figure(h3)
subplot 211
legend(pairs{:,1})
set(legend,'Interpreter','none')
xlabel('time');ylabel('direction')
figure(h3)
subplot 212
legend(pairs{:,1})
set(legend,'Interpreter','none')
xlabel('time');ylabel('speed')


    %
%% for t=4
    figure;
%
for tt=1:8
%     t=[-105 -65];
%     t=t+200;
    if tt==1||tt==2
        tstr=['SD',sprintf('%2.0f',stdlev(tt,1)),' deg / ','no speed variance'];
    elseif tt==3||tt==4
        tstr=['no direction variance / SD ',sprintf('%2.0f',stdlev(tt,2)),' dps'];
    else
        tstr=['SD',sprintf('%2.0f',stdlev(tt,1)),' deg / ',sprintf('%2.0f',stdlev(tt,2)),' dps'];
    end
%     title(tstr)
    writerObj=VideoWriter(['ga031015strf_',pairs{tt,1},'.avi']);
    set(writerObj,'FrameRate',30,'Quality',100)
    open(writerObj)
    for i=50:1:200  
        ct=i-50;
        normsta(:,:,:,tt)=mean(sta(:,:,:,:,tt),4)./max(max(max(sta(:,:,:,:,tt),[],3),[],2));
%         diffsta=normsta(:,:,i,tt)'-normsta(:,:,200,tt)';

%         imagesc(binc_spd,binc_dir,diffsta./max(max(diffsta)))
%         subplot(2,4,tt)
%         diffsta=normsta(:,:,t(1),tt)'-normsta(:,:,t(2),tt)';
%         imagesc(binc_dir,binc_spd,diffsta)
        normsta(60,60,i,tt)=1;

        imagesc(binc_dir(1:47),binc_spd(8:50),(fliplr(normsta(1:47,8:50,i,tt)')))
        caxis([0 0.85])
%         title([num2str(abs(lags(t(i)))),'ms pre-spike'])
        colorbar
%         title(tstr)
        title(sprintf('%2.0f',lags(i)),'FontSize',16)
        ylabel('Speed','FontSize',16)
        xlabel('Direction','FontSize',16)
        F=getframe(gcf);
        writeVideo(writerObj,F);
        pause(0.01)

    end
    close(writerObj)
end
%%   suptitle(tstr)
figure;
for tt=1:numtriTypes
    t=[-105 -65];
    t=t+200;
%     for i=1:length(t)
        normsta(:,:,:,tt)=mean(sta(:,:,:,:,tt),4)./maxsta(tt);
%         diffsta=normsta(:,:,i,tt)'-normsta(:,:,200,tt)';

%         imagesc(binc_spd,binc_dir,diffsta./max(max(diffsta)))
        subplot(2,4,tt)
        
        diffsta=normsta(:,:,t(2),tt)'-normsta(:,:,t(1),tt)';
        diffsta(end,end)
        imagesc(fliplr(binc_dir(1:47)),binc_spd(8:50),(diffsta(8:50,1:47)))
        if tt==1||tt==5
            ylabel('Speed (dps)','FontSize',16)
        end
        if tt==5||tt==6||tt==7||tt==8
            xlabel('Direction (deg)','FontSize',16)
        end

            
%         imagesc(binc_dir,binc_spd,normsta(:,:,t(i),tt)')

        
%         title([num2str(abs(lags(t(i)))),'ms pre-spike'])
caxis([-0.6 0.5])
        colorbar
%         title(lags(t(i)))
%         pause(0.05)
%     end
    if tt==1||tt==2
        tstr=['SD ',sprintf('%2.0f',stdlev(tt,1)),' deg / ','0 dps'];
    elseif tt==3||tt==4
        tstr=['SD 0 deg / ',sprintf('%2.0f',stdlev(tt,2)),' dps'];
    else
        tstr=['SD ',sprintf('%2.0f',stdlev(tt,1)),' deg / ',sprintf('%2.0f',stdlev(tt,2)),' dps'];
    end
    title(tstr,'FontSize',16)
end

