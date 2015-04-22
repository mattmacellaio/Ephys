
clear all;

if ispc
    rt='Z:/';
else
    rt='/Volumes/Hub/';
end
path(path,[rt,'MT/MATLAB/bing_ana']);
path(path,[rt,'RMVideo_backup']);
trialsdir=[rt,'MT/Data/ga031015/maestro'];
savedir=[rt,'MT/Data/ga031015/mat/'];
experiment='ga031015a1_dirspd_simul45';
first=1;
last=1763;


neuron_idx=1;

% %
try
    load([savedir,'data.mat'])
catch
    [data] = load_maestrospks(trialsdir,experiment,first,last);
    save([savedir,'data.mat'],'data')
end
%
seg_dur=500;  % important if we change the length of the stimulus.
% tdur=2*seg_dur;
tdur=3*seg_dur; %-250 to remove transient at beginning

nTypes = size(data,1);
nTags = 1;            % tags ??
nSegs = 0;
nstimdir=[];
stdpdir = []; 
spikearray=[]; 
spiketimes = cell(nTypes*nTags,1);
nSegs=0;
for i=1:nTypes
    temp = data(i,1);
    temp= [temp.targ.patdir];
    nReps = size(temp,2);
    nSegs = nSegs + nTags*nReps;
end 
celltimes = cell(nSegs,1);
%
segdata.nSpks = zeros(nTypes,nTags,nReps);
segdata.target = zeros(tdur,nTypes,nTags,nReps);
segdata.bspks = segdata.target;
segdata.spiketimes = zeros(200,nTypes,nTags,nReps); % 200 is bigger than the max spikes num
segdata.nfiles = zeros(1,nTypes);
tagval=[]; %tagvals identifies the trial condition for all segments within trials (so nTypes*nSegs*nReps long);
count=0;
%
for i=[1:4, 6:nTypes] %skip dir045_HTL_LTL for now: too many issues with putting short pert seg in matrices
        clear ptemp stemp;
        temp = data(i,1);
        tempdir = [temp.targ.patdir];
        temph=[temp.targ.hpatvel];
        tempv=[temp.targ.vpatvel];
        nReps = size(tempdir,2);
        tempspd=[data(i,1).targ.patspeed];
        stemp = data(i,1);
        stemp = stemp.spks; 
        for k= 1:nTags
            
%this was recorded using old version of pert_adapt_dir_spd45: motion and perts start together at
%401 for dir045_HTL_LTL (and similar trials) but perts start at 651 for
%dir045_HLL_spd (and similar trials)
%             minT=901;
            if strcmp(data(i).trname{1}(9),'T')
                minT=401; %401+250 
            else
                minT=651;
            end
%             if strcmp(data(i).trname{1},'dir045_HTL_LTL')
%                 maxT=tdur-200+minT-1; %2nd pert 300 instead of 500
%             else
                maxT=tdur+minT-1;
%             end

            %prefDir = [data(1).targ.predir]; % direction that has been added to target motion
            %prefDir = prefDir(1);
            ptempdir = tempdir(minT:maxT,1:nReps)';
            pdirarray = zeros(nReps,size(ptempdir,2));
            npdirarray = pdirarray;
            ptempspd = tempspd(minT:maxT,1:nReps)';
            pspdarray = zeros(nReps,size(ptempspd,2));
            npspdarray = pspdarray;
            
            bspks = pdirarray;
            %bspks = zeros(nReps,length(minT:maxT));
            stimes=[];
            for j=1:nReps
                pdirarray(j,:)=ptempdir(j,:);
%                 npdirarray(j,:) = pdirarray(j,:);
                npdirarray(j,:) = pdirarray(j,:) - repmat(mean(pdirarray(j,:)),1,size(ptempdir,2)); 
                pspdarray(j,:)=ptempspd(j,:);
%                 npspdarray(j,:) = pspdarray(j,:); %without removing mean speed
                npspdarray(j,:) = pspdarray(j,:) - repmat(mean(pspdarray(j,:)),1,size(ptempspd,2)); 
                % get spikes
                
                stmp = round(stemp(j,:));
                stmp = stmp(stmp>minT);
                stmp = stmp(stmp<maxT);
                nspks = length(stmp);
                stimes(j,1:nspks) = stmp;
                binspks = zeros(1,length(minT:maxT));
                %binspks([stmp-minT])==1;
                binspks([stmp-round(double(minT))])=1;
                bspks(j,1:length(minT:maxT)) = binspks;
                count=count+1;
                celltimes{count} = stmp-round(double(minT)); %store times relative to segment onset
                segdata.nSpks(i,k,j) = nspks;
                tagval(count) = i;
            end
            %stdpdir(i,k,1:nReps) = sqrt(var(reshape(npdirarray,1,nReps*size(ptemp,2))); 
            % make an array the same length as stim and spike arrays with the std levels of the perturbations
            stdpdir = [stdpdir repmat(sqrt(var(reshape(rad2deg(unwrap(deg2rad(npdirarray))),1,nReps*size(ptempdir,2)))),1,nReps)];
            nstimdir{i} = npdirarray;
            nstimspd{i}=npspdarray;
            spikearray = [spikearray ; bspks];
            %spiketimes{(i-1)*nTags+k} = stimes; 
            segdata.target(1:length(minT:maxT),i,k,1:nReps) = reshape(pdirarray',length(minT:maxT),1,1,nReps); 
            segdata.bspks(1:length(minT:maxT),i,k,1:nReps) = reshape(bspks',length(minT:maxT),1,1,nReps); 
            segdata.spiketimes(1:size(stimes,2),i,k,1:nReps) = reshape(stimes',size(stimes,2),1,1,nReps);
            segdata.nfiles(i) = nReps;

        end
end 
%
path(path,[rt,'/MattM/Code/population decoding'])
%%
% for i=[1:4, 6:nTypes]
%     figure(1)
%     subplot(ceil(nTypes/2),2,i)
%     plot(nstimdir{i}')
%     ylim([mean(mean(nstimdir{i}'))-90 mean(mean(nstimdir{i}))+90])
%     title([data(i,1).trname])
% 
% %
%     figure(2)
%     subplot(ceil(nTypes/2),2,i)
%     plot(nstimspd{i}')
%     ylim([0 40])
%     title([data(i,1).trname])
% 
% 
% figure(3)
% 
%     subplot(ceil(nTypes/2),2,i)
%     plot(rad2deg(circ_std(deg2rad(nstimdir{i})).^2))
%     title(data(i,1).trname{1})
%     ylim([0 50])
%     if i==3
%         ylabel('Variance of Stimulus Direction')
%     end
%     figure(4)
%     subplot(ceil(nTypes/2),2,i)
%     plot(var(nstimspd{i}))
%     ylim([0 100])
%     title(data(i,1).trname{1})
%     if i==3
%         ylabel('Variance of Stimulus Speed')
%     end
% 
% end
%%
%
try
    load([savedir,'nexfile.mat'])
catch
  [nexFile] = readNexFile()
  save([savedir,'nexfile.mat'],'nexFile')
end
ts=nexFile.markers{1}.timestamps;
%events
ev=nexFile.markers{1}.values{1}.strings;
n_ev=length(ev);
 
%
%trial start
st=zeros(1,n_ev);
for i=1:n_ev
    st(i)=str2num(ev{i})==2;
end
st1=find(st==1);

%trial end
st=zeros(1,n_ev);
for i=1:n_ev
    st(i)=str2num(ev{i})==3;
end
ed1=find(st==1);

%trial saved
st=zeros(1,n_ev);
for i=1:n_ev
    st(i)=str2num(ev{i})==6;
end
savd=find(st==1);


trialnames=cell(1,length(st1));
maestronames=cell(1,length(st1));
ind=0;
for i=1:length(st1)
    while (str2num(ev{st1(i)+ind})~=0)
        tmp=native2unicode(str2num(ev{st1(i)+ind}));
        ind=ind+1;
        trialnames{i}=[trialnames{i} tmp];
    end
    ind1=1;
    while (str2num(ev{st1(i)+ind+ind1})~=0)
        tmp=native2unicode(str2num(ev{st1(i)+ind+ind1}));
        ind1=ind1+1;
        maestronames{i}=[maestronames{i} tmp];
    end
    ind=0;
end

%
stsn=zeros(1,1);
idx=1;
for i=1:length(ed1)
    tmp=ed1(ed1>st1(i));
    tped=tmp(1);
    if str2num(ev{tped-1})==6
        stsn(idx)=i;
        idx=idx+1;
    end
end

% stsn=stsn(3:end);       % we have to change this for different day's data
% we get 116 trials for the results


% pay attention
 
 %
 
  stt=zeros(1,length(stsn));
edt=zeros(1,length(stsn));

for i=1:length(stsn)
    stt(i)=ts(st1(stsn(i)));
    tmp=ed1(ed1>st1(stsn(i)));
    edt(i)=ts(tmp(1));
end
%

startmarkers=zeros(2,length(stsn));
markdou2=zeros(2,length(stsn));
markdou3=zeros(5,length(stsn));

mark2=nexFile.events{1}.timestamps;
dou2=nexFile.events{2}.timestamps;
dou3=nexFile.events{3}.timestamps;

for i=1:length(stt)
    tmp=mark2(mark2>stt(i));
    tmp=tmp(tmp<edt(i));
    startmarkers(:,i)=tmp;
    
    tmp=dou2(dou2>stt(i));
    tmp=tmp(tmp<edt(i));
    markdou2(:,i)=tmp;
    
    tmp=dou3(dou3>stt(i));
    tmp=tmp(tmp<edt(i));
    if size(tmp,1)==4
        tmp=[tmp(1);0;tmp(2:end)]; %for dir045_HTL_HTL trials missing motion/no-pert seg
    end
    markdou3(:,i)=tmp;
end

%
spikes=nexFile.neurons{neuron_idx}.timestamps;   % pay attention here, choose the unit you want

%
spk_tr=cell(3,length(stsn));
spk_fortrans=cell(1,length(stsn));
for i=1:length(stsn)
    %1st pert
    tmp=spikes;
    tmp=tmp(tmp>markdou3(3,i)); 
    tmp=tmp(tmp<markdou3(4,i));
    spk_tr{1,i}=1000.*(tmp-markdou3(3,i));
    
    %2nd pert
    tmp=spikes;
    tmp=tmp(tmp>markdou3(4,i));
    tmp=tmp(tmp<markdou3(5,i));
    spk_tr{2,i}=1000*(tmp-markdou3(4,i));
    
    %3rd pert
    tmp=spikes;
    tmp=tmp(tmp>markdou3(5,i));
    tmp=tmp(tmp<markdou2(2,i));
    spk_tr{3,i}=1000.*(tmp-markdou3(5,i));
    
%     all pert segs
    tmp=spikes;
%     tmp=tmp(tmp>markdou3(3,i));
    tmp=tmp(tmp>markdou3(3,i)); 
    tmp=tmp(tmp<markdou2(2,i));
    spk_fortrans{i}=1000.*(tmp-markdou3(3,i));
end
% have a look of the spk_tr, if the firing rate seems too low, maybe we will
% not get good results.

  %

spk_use=spk_tr;      % here get the trials we want to use
tri_use=trialnames(stsn); % same aim as spk_use

tri_num=zeros(1,length(stsn));
tname={};
for i=1:length(stsn)
    tmp=tri_use{i};
    tname = unique([tname tmp]);
end

for i=1:length(stsn)
    tri_num(i)=find(strcmp(tname,tri_use{i}));
end

%  the firing rate of the neuron= spk_tr, trialnames=tri_num. 
% perts of the stimulus = pert.
 
%


%%

numlevels=2; %HTL, LTH
numconditions=2; %dir, spd

ntype=unique(tri_num);
tris=length(tri_num);
n_eachtri=zeros(1,length(ntype));
spktri=cell(1,length(ntype));
transpk=cell(1,length(ntype));
ntype_short=[1:4, 6:nTypes];
%trial order: HLL_dir, HLL_spd, HTL_HTL, HTL_LTH, HTL_LTL (skipped), 
%LHH_dir, LHH_spd, LTH_HTL, LTH_LTH, LTH_LTL, LTL_HTL, LTL_LTH
for i=1:tris
%     n_eachtri(tri_num(i))=n_eachtri(tri_num(i))+1;
%     spktri{tri_num(i),(4*n_eachtri(tri_num(i))-3):(4*n_eachtri(tri_num(i)))}=spk_use(:,i);
      spktri{tri_num(i)}=[spktri{tri_num(i)} spk_use(:,i)];
      transpk{tri_num(i)}=[transpk{tri_num(i)} spk_fortrans(:,i)];
end

for i=1:tris
    for j=1:length(ntype)
       if tri_num(i)==ntype(j)
           n_eachtri(j)=n_eachtri(j)+1;
       end
    end
end

%
%
cumn=1;
stimtrain_dir=cell(1,length(ntype));
%
for i=ntype_short
    tp=n_eachtri(i);
    for j=1:tp
        stimtrain_dir{i}=[stimtrain_dir{i}; nstimdir{i}(j,:)];
%         cumn=cumn+1;
    end
end
%
stimtrain_spd=cell(1,length(ntype));

cumn=1;
for i=ntype_short
    tp=n_eachtri(i);
    for j=1:tp
        stimtrain_spd{i}=[stimtrain_spd{i}; nstimspd{i}(j,:)];
%         cumn=cumn+1;
    end
end
%% bin spikes to calc gain
tnoise=cell(1,2);
spkcounts=cell(1,2);
bintransspk=cell(1,12);
bin=20;
nbin=1:(seg_dur-bin+1);
i=1;
spkseg_trans=cell(1,12); 
%alltrials 1pert-alltrials 2pert-alltrials 3pert: should be
%allperts-eachtrial
for j=1:12
    spkseg_trans{j}=[spktri{j}(1,:) spktri{j}(2,:) spktri{j}(3,:)];
end
skseg_trans=cell(1,12);
i=1
for j=1:12
    tmp=spkseg_trans{i,j};
    ntr=length(tmp);
    for k=1:ntr
        %only use spikes after end of initial transient
%             inds_posttrans=find(spkseg_dir{i,j}{k}>transient_endtime);
        skseg_trans{i,j}(k,1:length(spkseg_trans{i,j}{k}))=round(spkseg_trans{i,j}{k});
%             skseg_dir{i,j}(k,1:length(spkseg_trans))=round(spkseg_dir{i,j}{k});
    end

    tmp=skseg_trans{i,j}';
    spkcounts{i,j}=zeros(size(skseg_trans{i,j}));
    for k=1:size(tmp,2)
        counts=tmp(:,k);
        for tt=1:max(counts);
            spkcounts{i,j}(tt,k) = length(find(counts==tt));    
        end	
        for b=1:length(nbin)
            tp=counts(counts>b);
            tp=tp(tp<b+bin-1);
            bintransspk{i,j}(b,k)=length(tp);
        end
    end
end


%%
% spklag=55; %-10 to get 20ms bin surrounding peak at 65
% 
% for j=1:12
%     bintransspk{i,j}=bintransspk{i,j}.*1000./bin;
%     for t=1:size(stimtrain_dir{1,j},2)
%         bin_i=t+spklag;
%         p=polyfit(stimtrain_dir{1,j}(:,t),bintransspk{1,j}(:,bin_i));
%         gain(j,t)=p(1);
%     end
% end
%% sandbox for direction-speed STA
pairs=[{'dirL_nospd'},{[1 6]},{{[2,3],[1]}};{'dirH_nospd'},{[1 6]},{{[1],[2,3]}};...
        {'spdL_nodir'},{[2 7]},{{[2,3],[1]}};{'spdH_nodir'},{[2 7]},{{[1],[2,3]}};...
        {'dirL_spdL'},{[3 9 10 11 12]},{{[3],[2],[2],[3],[2]}};{'dirH_spdH'},{[3 9]},{{[2],[3]}};... %  {'dirL_spdL'},{[3 5 9 10 11 12]},{{[3],[3],[2],[2],[3],[2]}}
        {'dirL_spdH'},{[4 8 11 12]},{{[3],[2],[2],[3]}};{'dirH_spdL'},{[4 8 10]},{{[2],[3],[3]}}]; %{'dirH_spdL'},{[4 5 8 10]},{{[2],[2],[3],[3]}}
% for varpair=1:size(pairs,1)/2
    
stimspk=cell(size(pairs,1),3); % dir stim, spd stim, spikes trials x time.
numtriTypes=size(pairs,1);
%
for i=1:numtriTypes
    spkseg_tmp=[];
    for trI=1:size(pairs{i,2},2)            %trialtype index
        tr=pairs{i,2}(trI);                 %trialtype
        for segI=1:size(pairs{i,3}{trI},2)  %segment index
            seg=pairs{i,3}{trI}(segI);      %segment
            stimspk{i,1}=[stimspk{i,1};stimtrain_dir{1,tr}(:,(seg-1)*seg_dur+1:seg*seg_dur)];
            stimspk{i,2}=[stimspk{i,2};stimtrain_spd{1,tr}(:,(seg-1)*seg_dur+1:seg*seg_dur)];
            spkseg_tmp=[spkseg_tmp,spktri{1,tr}(seg,:)];
        end
    end
    ntr=size(spkseg_tmp,2);
    for k=1:ntr
        stimspk{i,3}(k,1:length(spkseg_tmp{k}))=round(spkseg_tmp{k});
    end
end

% blorp
% dirspd_sta


%
%% direction linfilt
skseg(:,1)=stimspk(:,3);
stimseg(:,1)=stimspk(:,1);

tnoise=cell(numtriTypes,1);
spkcounts=cell(numtriTypes,1);
sldspk=cell(numtriTypes,1);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:numtriTypes
   j=1;
        tnoise{i,j}=stimseg{i,j}';
        
        tmp=skseg{i,j}';
        spkcounts{i,j}=zeros(size(tnoise{i,j}));
        for k=1:size(tmp,2)
            counts=tmp(:,k);
            for tt=1:max(counts);
	            spkcounts{i,j}(tt,k) = length(find(counts==tt));    
            end	
            for b=1:length(nbin)
                tp=counts(counts>b);
                tp=tp(tp<b+bin-1);
                sldspk{i,j}(b,k)=length(tp);
            end
        end
        
end
%
for i=1:numtriTypes
    j=1;
        sldspk{i,j}=sldspk{i,j}.*1000./bin;
end
%


%%
shift =50;
shifts=[0 shift];
cutoff=25;
ft=250;
maxK=20;
pshift=1; 
lag=[1:ft];

for i=1:numtriTypes
    j=1;
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
end

mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error

amp_filter=zeros(numtriTypes,1); 
std_filter=zeros(numtriTypes,1);
m_stdsti=zeros(numtriTypes,1);
figure;
colors=distinguishable_colors(numtriTypes);
for i=[1,2,5:8]
    hold on;
%     subplot(1,2,1);
     j=1;
        stim=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [filtdir, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(stim,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if i==numtriTypes
                break;
            end
            continue;
        end
        myfilter=[filtdir(pshift,myidx).filter_allt]*1000;
        [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stimseg{i,j}'));
           errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',colors(i,:)); 
           hold on; 
       
        linfilt_dir(i,j).results=filtdir;
        linfilt_dir(i,j).eyeest=eye_est;
        linfilt_dir(i,j).index2=index2_list;
        linfilt_dir(i,j).residuals=residuals;
        linfilt_dir(i,j).ccoef=ccoef;
        linfilt_dir(i,j).error=error;
        linfilt_dir(i,j).idx=myidx;
%         filter(i,j).coef=mean(ccoef.allt(pshift,myidx));
end
    set(gca,'FontSize',16,'Box','Off','TickDir','Out');
    xlim([-50 250]);
%     legend(['20 deg ccoef: ',num2str(filter(1).coef)],['60 deg ccoef: ',num2str(filter(2).coef)]);
    legend(pairs{:,1})
    xlabel('Time lag (ms)');
    ylabel('Filter Amplitude');
    axis square;
%     title('bu080312 filter');
    box off;
saveas(gcf,[experiment,'linfilt_dir.fig'])
save([experiment,'_linfilt_dir.mat'],'linfilt_dir','pairs')
%

% amp_filter=abs(amp_filter);
% subplot(1,2,2);
% 
% id1=[1];
% color=[0 0 1; 0 1 0];
% for i=1:length(id1)
%     errorbar(m_stdsti(:,id1(i)),amp_filter(:,id1(i)),std_filter(:,id1(i)),'Color',color(i,:));
%     hold on;
% end
% legend('dir45','dir-45')
% set(gca,'FontSize',16,'Box','Off','TickDir','Out');
% xlabel('std of stimulus (deg)');
% ylabel('Filter Amplitude');
% axis square;
% box off;

% speed linfilt
skseg(:,1)=stimspk(:,3);
stimseg(:,1)=stimspk(:,2);

tnoise=cell(numtriTypes,1);
spkcounts=cell(numtriTypes,1);
sldspk=cell(numtriTypes,1);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:numtriTypes
   j=1;
        tnoise{i,j}=stimseg{i,j}';
        
        tmp=skseg{i,j}';
        spkcounts{i,j}=zeros(size(tnoise{i,j}));
        for k=1:size(tmp,2)
            counts=tmp(:,k);
            for tt=1:max(counts);
	            spkcounts{i,j}(tt,k) = length(find(counts==tt));    
            end	
            for b=1:length(nbin)
                tp=counts(counts>b);
                tp=tp(tp<b+bin-1);
                sldspk{i,j}(b,k)=length(tp);
            end
        end
        
end
%
for i=1:numtriTypes
    j=1;
        sldspk{i,j}=sldspk{i,j}.*1000./bin;
end
%


%
% shift =50;
% shifts=[0 shift];
% cutoff=30;
% ft=250;
% maxK=20;
% pshift=1; 
% lag=[1:ft];

for i=1:numtriTypes
    j=1;
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
end

mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error

amp_filter=zeros(numtriTypes,1); 
std_filter=zeros(numtriTypes,1);
m_stdsti=zeros(numtriTypes,1);
figure;
colors=distinguishable_colors(numtriTypes);
for i=[3:8] 
    hold on;
%     subplot(1,2,1);
     j=1;
        stim=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [filtspd, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(stim,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if i==numtriTypes
                break;
            end
            continue;
        end
        myfilter=[filtspd(pshift,myidx).filter_allt]*1000;
        [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stimseg{i,j}'));
           errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',colors(i,:)); 
           hold on; 
       
        linfilt_spd(i,j).results=filtspd;
        linfilt_spd(i,j).eyeest=eye_est;
        linfilt_spd(i,j).index2=index2_list;
        linfilt_spd(i,j).residuals=residuals;
        linfilt_spd(i,j).ccoef=ccoef;
        linfilt_spd(i,j).error=error;
        linfilt_spd(i,j).idx=myidx;
%         filter(i,j).coef=mean(ccoef.allt(pshift,myidx));
end
    set(gca,'FontSize',16,'Box','Off','TickDir','Out');
    xlim([-50 250]);
%     legend(['20 deg ccoef: ',num2str(filter(1).coef)],['60 deg ccoef: ',num2str(filter(2).coef)]);
    legend(pairs{:,1})
    xlabel('Time lag (ms)');
    ylabel('Filter Amplitude');
    axis square;
%     title('bu080312 filter');
    box off;
saveas(gcf,[experiment,'linfilt_spd.fig'])
save([experiment,'_linfilt_spd.mat'],'linfilt_spd','pairs')

%% speed STA and linfilt

skseg=skseg_spd;
stimseg=stimseg_spd;
rep=10;
tL=400;
sta=zeros(tL,rep,1,2);
STA=zeros(tL,rep,1,2);
STC=zeros(tL,rep,1,2);


for i=1:1
    for j=1:2
        for rp=1:rep
            ntri=size(skseg{i,j},1);
            index=randperm(ntri);
            ind_use=index(1:round(0.8*length(index)));
            [sta1,STA1,STC1] = get_sta(skseg{i,j}(ind_use,:), stimseg{i,j}(ind_use,:));
            sta(:,rp,i,j)=sta1;
        end
    end
end

%%%%
% stseg{2,1}=stseg{2,1}(20:end,:);
% skseg{2,1}=skseg{2,1}(20:end,:);
%%%%

lags=-199:200;
flip=-1;

%
figure;
% subplot(1,2,1);
errorbar(lags,flip*mean(sta(:,:,1,1),2),std(sta(:,:,1,1)'));
hold on;
errorbar(lags,flip*mean(sta(:,:,1,2),2),std(sta(:,:,1,2)'),'r');
set(gca,'FontSize',16,'Box','Off','TickDir','Out');
legend('4 dps','8 dps');
xlabel('Time lag (ms)');
ylabel('Spike Triggered Average (dps)');
xlim([-200 200]);
axis square;
% title('bu080312 predir sta');
box off;

% hold on;
% subplot(1,2,2);
% errorbar([-199:200],mean(sta(:,:,2,1)'),std(sta(:,:,1,1)'));
% hold on;
% errorbar([-199:200],mean(sta(:,:,2,2)'),std(sta(:,:,1,2)'),'r');
% set(gca,'FontSize',16,'Box','Off','TickDir','Out');
% legend('20 deg','60 deg');
% xlabel('Time lag (ms)');
% ylabel('Spike Triggered Average (deg)');
% xlim([-200 200]);
% axis square;
% % title('bu080312 45dir sta');
% box off;

for j=1:2
    [maxval, maxind]=max(flip*mean(sta(:,:,1,j),2));
    latency=lags(maxind)
end

%

% we will get the filters of the spikes to the stimulus
%  here we can use the stseg as the targets and the skseg as the spikes.
%
tnoise=cell(1,2);
spkcounts=cell(1,2);
sldspk=cell(1,2);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:1
    for j=1:2
        tnoise{i,j}=stimseg{i,j}';
        
        tmp=skseg{i,j}';
        spkcounts{i,j}=zeros(size(tnoise{i,j}));
        for k=1:size(tmp,2)
            counts=tmp(:,k);
            for tt=1:max(counts);
	            spkcounts{i,j}(tt,k) = length(find(counts==tt));    
            end	
            for b=1:length(nbin)
                tp=counts(counts>b);
                tp=tp(tp<b+bin-1);
                sldspk{i,j}(b,k)=length(tp);
            end
        end
        
    end
end

for i=1:1
    for j=1:2
        sldspk{i,j}=sldspk{i,j}.*1000./bin;
    end
end



%
shift =50;
shifts=[0 shift];
cutoff=30;
ft=250;
maxK=3;
pshift=1; 
lag=[1:ft];

for i=1:1
    for j=1:2
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
    end
end
%
mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error
%
amp_filter=zeros(1,2); 
std_filter=zeros(1,2);
m_stdsti=zeros(1,2);
figure;
i=1;
    hold on;
    subplot(1,2,1);
    for j=1:2
        stim=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [filtspd, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(stim,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if j==2
                break;
            end
            continue;
        end
            myfilter=[filtspd(pshift,myidx).filter_allt]*1000;
            [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stimseg{i,j}'));
           errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',mycolor(-30+(31*j),:)); 
           hold on; 
       
        linfilt_spd(i,j).results=filtspd;
        linfilt_spd(i,j).eyeest=eye_est;
        linfilt_spd(i,j).index2=index2_list;
        linfilt_spd(i,j).residuals=residuals;
        linfilt_spd(i,j).ccoef=ccoef;
        linfilt_spd(i,j).error=error;
        linfilt_spd(i,j).idx=myidx;
%         filter(i,j).coef=mean(ccoef.allt(pshift,myidx));
    end
    set(gca,'FontSize',16,'Box','Off','TickDir','Out');
    xlim([-50 250]);
%     legend(['20 deg ccoef: ',num2str(filter(1).coef)],['60 deg ccoef: ',num2str(filter(2).coef)]);
    legend('var 4','var 8')
    xlabel('Time lag (ms)');
    ylabel('Filter Amplitude');
    axis square;
%     title('bu080312 filter');
    box off;
    


amp_filter=abs(amp_filter);
subplot(1,2,2);

id1=[1];
color=[0 0 1; 0 1 0];
for i=1:length(id1)
    errorbar(m_stdsti(:,id1(i)),amp_filter(:,id1(i)),std_filter(:,id1(i)),'Color',color(i,:));
    hold on;
end
legend('dir45','dir-45')
set(gca,'FontSize',16,'Box','Off','TickDir','Out');
xlabel('std of stimulus (dps)');
ylabel('Filter Amplitude');
axis square;
box off;


%% to calc filters of non-changing speed

skseg=skseg_dir;
stimseg=stimseg_spd_off;

tnoise=cell(1,2);
spkcounts=cell(1,2);
sldspk=cell(1,2);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:1
    for j=1:2
        tnoise{i,j}=stimseg{i,j}';
        
        tmp=skseg{i,j}';
        spkcounts{i,j}=zeros(size(tnoise{i,j}));
        for k=1:size(tmp,2)
            counts=tmp(:,k);
            for tt=1:max(counts);
	            spkcounts{i,j}(tt,k) = length(find(counts==tt));    
            end	
            for b=1:length(nbin)
                tp=counts(counts>b);
                tp=tp(tp<b+bin-1);
                sldspk{i,j}(b,k)=length(tp);
            end
        end
        
    end
end

for i=1:1
    for j=1:2
        sldspk{i,j}=sldspk{i,j}.*1000./bin;
    end
end



%
shift =50;
shifts=[0 shift];
cutoff=30;
ft=250;
maxK=3;
pshift=1; 
lag=[1:ft];

for i=1:1
    for j=1:2
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
    end
end
%
mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error
%
amp_filter=zeros(1,2); 
std_filter=zeros(1,2);
m_stdsti=zeros(1,2);
figure;
i=1;
    hold on;
    subplot(1,2,1);
    for j=1:2
        stim=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [filt1spd, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(stim,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if j==2
                break;
            end
            continue;
        end
            myfilter=[filt1spd(pshift,myidx).filter_allt]*1000;
            [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stimseg{i,j}'));
           errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',mycolor(-30+(31*j),:)); 
           hold on; 
       
        linfilt_1spd(i,j).results=filt1spd;
        linfilt_1spd(i,j).eyeest=eye_est;
        linfilt_1spd(i,j).index2=index2_list;
        linfilt_1spd(i,j).residuals=residuals;
        linfilt_1spd(i,j).ccoef=ccoef;
        linfilt_1spd(i,j).error=error;
        linfilt_1spd(i,j).idx=myidx;
%         filter(i,j).coef=mean(ccoef.allt(pshift,myidx));
    end
    set(gca,'FontSize',16,'Box','Off','TickDir','Out');
    xlim([-50 250]);
%     legend(['20 deg ccoef: ',num2str(filter(1).coef)],['60 deg ccoef: ',num2str(filter(2).coef)]);
    legend('var 4','var 8')
    xlabel('Time lag (ms)');
    ylabel('Filter Amplitude');
    axis square;
%     title('bu080312 filter');
    box off;
    


amp_filter=abs(amp_filter);
subplot(1,2,2);

id1=[1];
color=[0 0 1; 0 1 0];
for i=1:length(id1)
    errorbar(m_stdsti(:,id1(i)),amp_filter(:,id1(i)),std_filter(:,id1(i)),'Color',color(i,:));
    hold on;
end
legend('dir45','dir-45')
set(gca,'FontSize',16,'Box','Off','TickDir','Out');
xlabel('std of stimulus (dps)');
ylabel('Filter Amplitude');
axis square;
box off;
%% to calc filters of non-changing dir

skseg=skseg_spd;
stimseg=stimseg_dir_off;

tnoise=cell(1,2);
spkcounts=cell(1,2);
sldspk=cell(1,2);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:1
    for j=1:2
        tnoise{i,j}=stimseg{i,j}';
        
        tmp=skseg{i,j}';
        spkcounts{i,j}=zeros(size(tnoise{i,j}));
        for k=1:size(tmp,2)
            counts=tmp(:,k);
            for tt=1:max(counts);
	            spkcounts{i,j}(tt,k) = length(find(counts==tt));    
            end	
            for b=1:length(nbin)
                tp=counts(counts>b);
                tp=tp(tp<b+bin-1);
                sldspk{i,j}(b,k)=length(tp);
            end
        end
        
    end
end

for i=1:1
    for j=1:2
        sldspk{i,j}=sldspk{i,j}.*1000./bin;
    end
end



%
shift =50;
shifts=[0 shift];
cutoff=30;
ft=250;
maxK=3;
pshift=1; 
lag=[1:ft];

for i=1:1
    for j=1:2
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
    end
end
%
mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error
%
amp_filter=zeros(1,2); 
std_filter=zeros(1,2);
m_stdsti=zeros(1,2);
figure;
i=1;
    hold on;
    subplot(1,2,1);
    for j=1:2
        stim=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [filt1dir, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(stim,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if j==2
                break;
            end
            continue;
        end
            myfilter=[filt1dir(pshift,myidx).filter_allt]*1000;
            [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stimseg{i,j}'));
           errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',mycolor(-30+(31*j),:)); 
           hold on; 
       
        linfilt_1dir(i,j).results=filt1dir;
        linfilt_1dir(i,j).eyeest=eye_est;
        linfilt_1dir(i,j).index2=index2_list;
        linfilt_1dir(i,j).residuals=residuals;
        linfilt_1dir(i,j).ccoef=ccoef;
        linfilt_1dir(i,j).error=error;
        linfilt_1dir(i,j).idx=myidx;
%         filter(i,j).coef=mean(ccoef.allt(pshift,myidx));
    end
    set(gca,'FontSize',16,'Box','Off','TickDir','Out');
    xlim([-50 250]);
%     legend(['20 deg ccoef: ',num2str(filter(1).coef)],['60 deg ccoef: ',num2str(filter(2).coef)]);
    legend('var 20','var 60')
    xlabel('Time lag (ms)');
    ylabel('Filter Amplitude');
    axis square;
%     title('bu080312 filter');
    box off;
    


amp_filter=abs(amp_filter);
subplot(1,2,2);

id1=[1];
color=[0 0 1; 0 1 0];
for i=1:length(id1)
    errorbar(m_stdsti(:,id1(i)),amp_filter(:,id1(i)),std_filter(:,id1(i)),'Color',color(i,:));
    hold on;
end
legend('dir45','dir-45')
set(gca,'FontSize',16,'Box','Off','TickDir','Out');
xlabel('std of stimulus (dps)');
ylabel('Filter Amplitude');
axis square;
box off;
%%
% bin spikes

bin=20;
spktrain=cell(size(transpk));
spkcount=cell(size(transpk));

for i=1:length(spktrain)
    temp=transpk{i};
    trnb=length(temp);
    spkcount{i}=zeros(tdur,trnb);
    for j=1:trnb
        tmp=temp{j};
        if ~isempty(tmp)
        if round(tmp(1))==0
            tmp(1)=1;
        end
        spkcount{i}(round(tmp),j)=1;
        for kk=1:tdur-bin+1
            tp=tmp(tmp>=kk);
            tp=tp(tp<=kk+bin-1);
            spktrain{i}(kk,j)=length(tp);
        end
        end
    end
end

for i=1:length(spktrain)
    spktrain{i}=spktrain{i}.*1000./bin;
end
%% plot psth 
figure
for i=1:nTypes
    subplot(ceil(nTypes/2),2,i)
    errorbar(mean(spktrain{i},2),std(spktrain{i},[],2)./sqrt(size(spktrain{i},2)));hold all
    line([seg_dur seg_dur],[0 max(mean(spktrain{i},2))],'Color','r')
    line([2*seg_dur 2*seg_dur],[0 max(mean(spktrain{i},2))],'Color','r')
    ylim([0 150])
    title(tname{i},'Interpreter','none')
end
    
%% calculate gain: delta r/delta s
for i=1:nTypes
    %get stim for each trial:nstimdir=trials x 1500
    
    %get response lagged by latency:spktrain = 1500-bin+1 x trials
    %scatter r vs s: bin s or plot r for each s?
    %for each timepoint across all trials
    %linear fit of r vs s
    %gain (timepoint)=slope(linear fit)
    
    
end

%% generate nonlinearity and spike predictions
%what to do about high and low var filters vs. HTL/LTH stimulus?
% speed
trialnums=[2,4];
for trialtype=1:length(trialnums)
    %using mean of high and low variance filters as proxy for basic filter
    meanfilt_spd=mean([mean([linfilt_spd(1,1).results(pshift,:).filter_allt],2), mean([linfilt_spd(1,2).results(pshift,:).filter_allt],2)],2);
    meanfilt_dir=mean([mean([linfilt_dir(1,1).results(pshift,:).filter_allt],2), mean([linfilt_dir(1,2).results(pshift,:).filter_allt],2)],2);

    for trial=1:size(nstimspd{trialnums(trialtype)},1)
        spk_est_spd{trialtype}(trial,:)=conv(meanfilt_spd',nstimspd{trialnums(trialtype)}(trial,:));
        spk_est_1dir{trialtype}(trial,:)=conv(meanfilt_dir',nstimdir{trialnums(trialtype)}(trial,:));
    end

    linest_spd{trialtype}=reshape(spk_est_spd{trialtype}(:,1:size(spktrain{trialnums(trialtype)},1)),1,[]);
    spkrate_spd{trialtype}=reshape(spktrain{trialnums(trialtype)}',1,[]);
    
    linest_1dir{trialtype}=reshape(spk_est_1dir{trialtype}(:,1:size(spktrain{trialnums(trialtype)},1)),1,[]);
    linest_combspd{trialtype}=linest_spd{trialtype}+linest_1dir{trialtype};
    
    [R,P]=corrcoef(linest_spd{trialtype},spkrate_spd{trialtype});
    R(1,2)
    [R,P]=corrcoef(linest_combspd{trialtype},spkrate_spd{trialtype});
    R(1,2)
%
%     inds0=find(abs(linest_spd{trialtype})<1);
    clear meanest_spd meanspkrate_spd;
    nbin=1000;
    [odest,myidx]=sort(linest_combspd{trialtype});
    num_sordata=floor(length(linest_combspd{trialtype})/nbin);
    for i=1:nbin
        meanest_spd{trialtype}(i)=mean(odest(((i-1)*num_sordata+1):i*num_sordata));
        meanspkrate_spd{trialtype}(i)=mean(spkrate_spd{trialtype}(myidx(((i-1)*num_sordata+1):i*num_sordata)));
    end

    figure;scatter(meanest_spd{trialtype},meanspkrate_spd{trialtype})

    for i=1:length(spkrate_spd{trialtype})
        N=hist(linest_combspd{trialtype}(i),meanest_spd{trialtype});
        spk_non_est_spd{trialtype}(i)=meanspkrate_spd{trialtype}(N==1);
    end

    [R,P]=corrcoef(spk_non_est_spd{trialtype},spkrate_spd{trialtype});
    R(1,2)


    non_est_spd{trialtype}=reshape(spk_non_est_spd{trialtype},size(spktrain{trialnums(trialtype)}'));
    figure;plot(mean(spktrain{trialnums(trialtype)}',1))
    hold on;plot(mean(non_est_spd{trialtype},1),'r')
    title(tname(trialnums(trialtype)))
end

%% direction
clear spk_est_dir linest_dir spkrate_dir
trialnums=[1,3]; %for nstimdir or low, high var if using stiseg_dir
for trialtype=1:length(trialnums)
    meanfilt_spd=mean([mean([linfilt_spd(1,1).results(pshift,:).filter_allt],2), mean([linfilt_spd(1,2).results(pshift,:).filter_allt],2)],2);
    meanfilt_dir=mean([mean([linfilt_dir(1,1).results(pshift,:).filter_allt],2), mean([linfilt_dir(1,2).results(pshift,:).filter_allt],2)],2);
    
    %should I use nstimdir or stiseg here?
    for trial=1:size(nstimdir{trialnums(trialtype)},1)
        spk_est_dir{trialtype}(trial,:)=conv(meanfilt_dir,nstimdir{trialnums(trialtype)}(trial,:));
        spk_est_1spd{trialtype}(trial,:)=conv(meanfilt_spd,nstimspd{trialnums(trialtype)}(trial,:));

    end

    linest_dir{trialtype}=reshape(spk_est_dir{trialtype}(:,1:size(spktrain{trialnums(trialtype)},1)),1,[]);
    spkrate_dir{trialtype}=reshape(spktrain{trialnums(trialtype)}',1,[]);
    
    linest_1spd{trialtype}=reshape(spk_est_1spd{trialtype}(:,1:size(spktrain{trialnums(trialtype)},1)),1,[]);
    linest_combdir{trialtype}=linest_dir{trialtype}+linest_1spd{trialtype};
    
    [R,P]=corrcoef(linest_dir{trialtype},spkrate_dir{trialtype});
    R(1,2)
    [R,P]=corrcoef(linest_combdir{trialtype},spkrate_dir{trialtype});
    R(1,2)
% end
%
clear meanest meanorieye;
nbin=1000;
% for trialtype=1:4
    [odest,myidx]=sort(linest_combdir{trialtype});
    num_sordata=floor(length(linest_combdir{trialtype})/nbin);
    for i=1:nbin
        meanest_dir{trialtype}(i)=mean(odest(((i-1)*num_sordata+1):i*num_sordata));
        meanspkrate_dir{trialtype}(i)=mean(spkrate_dir{trialtype}(myidx(((i-1)*num_sordata+1):i*num_sordata)));
    end

    figure;scatter(meanest_dir{trialtype},meanspkrate_dir{trialtype})

    for i=1:length(spkrate_dir{trialtype})
        N=hist(linest_combdir{trialtype}(i),meanest_dir{trialtype});
        spk_non_est_dir{trialtype}(i)=meanspkrate_dir{trialtype}(N==1);
    end
    [R,P]=corrcoef(spk_non_est_dir{trialtype},spkrate_dir{trialtype});
    R(1,2)


    non_est_dir{trialtype}=reshape(spk_non_est_dir{trialtype},size(spktrain{trialnums(trialtype)}'));

    figure;plot(mean(spktrain{trialnums(trialtype)}',1))
    hold on;plot(mean(non_est_dir{trialtype},1),'r')
end

%
%% info

seg_dur =690;
non_est_dir=eye(:,HTL(index2_list));
stim=tnoise.v(:,HTL(index2_list));



all_win_dur = 60; % duration of the sliding time window over which we'll measure spike count, direction

shift_max = 120;
shift_min = 70;
shift_step = 10;
 
all_nBins = 20;  % 8:4:20;  %number of adaptive bins

shift_vec = shift_min:shift_step:shift_max;
nshifts = length(shift_vec);



tnoise_info = [stim(:,1:size(non_est_dir,2))]';
pnoise_info = non_est_dir';
tnoise_shuffle = tnoise_info(randperm(size(tnoise_info,1)),randperm(size(tnoise_info,2)));
pnoise_shuffle = pnoise_info(randperm(size(pnoise_info,1)),randperm(size(pnoise_info,2)));

I=cell(length(all_win_dur),length(all_nBins)); Ierr=cell(length(all_win_dur),length(all_nBins)); Ifrac=cell(length(all_win_dur),length(all_nBins)); IR=cell(length(all_win_dur),length(all_nBins)); IRerr=cell(length(all_win_dur),length(all_nBins)); IRfrac=cell(length(all_win_dur),length(all_nBins));
I_shuffle=cell(length(all_win_dur),length(all_nBins)); Ierr_shuffle=cell(length(all_win_dur),length(all_nBins)); Ifrac_shuffle=cell(length(all_win_dur),length(all_nBins)); IR_shuffle=cell(length(all_win_dur),length(all_nBins)); IRerr_shuffle=cell(length(all_win_dur),length(all_nBins)); IRfrac_shuffle=cell(length(all_win_dur),length(all_nBins));
Pjoint=cell(length(all_win_dur),length(all_nBins),nshifts); %,ntbins); 
Pjoint_shuffle=Pjoint; 
PjointR = Pjoint; 
PjointR_shuffle=Pjoint; 


for iwindur=1:length(all_win_dur)
    win_dur=all_win_dur(iwindur);

    for ibins=1:length(all_nBins)

        nBins=all_nBins(ibins);


        for ishift = 1:nshifts

            fprintf('ishift = %d\n',ishift);
            shift_val = shift_vec(ishift);

            nWins = size(non_est_dir,2)-win_dur-shift_val+1;
            step = 8; 
            ntbins = nWins-(step-1);         

            vShifts = ones(1,nWins)*shift_val;
            vShifts(1,1:seg_dur) = shift_val;

            twins=zeros(nWins,size(tnoise_info,2));
            pwins=twins;

            twins_shuffle=zeros(nWins,size(tnoise_info,2));
            pwins_shuffle=twins_shuffle;


            for n=1:nWins

                twins(n,:) = mean(tnoise_info(n:n+win_dur-1,:),1);
                pwins(n,:) = mean(pnoise_info(n+vShifts(n):n+vShifts(n)+win_dur-1,:),1);

%                     pwins(n,:)=smooth(pwins(n,:),20);

                twins_shuffle(n,:) = mean(tnoise_shuffle(n:n+win_dur-1,:),1);
                pwins_shuffle(n,:) = sum(pnoise_shuffle(n+vShifts(n):n+vShifts(n)+win_dur-1,:),1);

%                     pwins_shuffle(n,:)=smooth(pwins_shuffle(n,:),20);
            end

%                 for n=1:size(twins,2)
%                     pwins(:,n)=smooth(pwins(:,n),20);
%                     pwins_shuffle(:,n)=smooth(pwins_shuffle(:,n),20);
%                 end

            % Adaptively bin values
            tcuts = []; tcuts_shuffle=[];
            pcuts = tcuts; pcuts_shuffle=[];

            tbins=[]; tocc=[]; pbins=[]; pocc=[];

            tbins_shuffle=[]; tocc_shuffle=[]; pbins_shuffle=[]; pocc_shuffle=[];

            for j=1:nWins
                [tcuts(j,:),tbins(j,:),tocc(j,:)] = adaptive_bins(twins(j,:),nBins);

                 [pcuts(j,:),pbins(j,:),pocc(j,:)] = adaptive_bins(pwins(j,:),nBins);
          %      pbins(j,:)=pwins(j,:)+1;

                [tcuts_shuffle(j,:),tbins_shuffle(j,:),tocc_shuffle(j,:)] = adaptive_bins(twins_shuffle(j,:),nBins);

                 [pcuts_shuffle(j,:),pbins_shuffle(j,:),pocc_shuffle(j,:)] = adaptive_bins(pwins_shuffle(j,:),nBins);
          %      pbins_shuffle(j,:)=pwins_shuffle(j,:)+1;


            end

            tbins_shuffle=tbins_shuffle(randperm(size(tbins_shuffle,1)),randperm(size(tbins_shuffle,2)));
            pbins_shuffle=pbins_shuffle(randperm(size(pbins_shuffle,1)),randperm(size(pbins_shuffle,2)));

            %steps = 1:step:nWins-(step);
            steps = 1:nWins-(step-1);

            for j = 1:ntbins
                istep = steps(j);
               % keyboard    
                 nyBins=max(max(pbins));
                [I{iwindur,ibins}(ishift,j),Ierr{iwindur,ibins}(ishift,j),Ifrac{iwindur,ibins}(ishift,j),IR{iwindur,ibins}(ishift,j),IRerr{iwindur,ibins}(ishift,j),IRfrac{iwindur,ibins}(ishift,j),Pjoint{iwindur,ibins,ishift}{j}(:,:), PjointR{iwindur,ibins,ishift}{j}(:,:)] =  ...
                        calc_info_P_joint(reshape(tbins(istep:istep+step-1,:),1,numel(tbins(istep:istep+step-1,:))), ...
                        reshape(pbins(istep:istep+step-1,:),1,numel(pbins(istep:istep+step-1,:))),nBins,nyBins,[1 0.95 0.9 0.85 0.8 0.5],30);

                nyBins=max(max(pbins_shuffle));
                [I_shuffle{iwindur,ibins}(ishift,j),Ierr_shuffle{iwindur,ibins}(ishift,j),Ifrac_shuffle{iwindur,ibins}(ishift,j),IR_shuffle{iwindur,ibins}(ishift,j),IRerr_shuffle{iwindur,ibins}(ishift,j),IRfrac_shuffle{iwindur,ibins}(ishift,j),Pjoint_shuffle{iwindur,ibins,ishift}{j}(:,:), PjointR_shuffle{iwindur,ibins,ishift}{j}(:,:)] =  ...
                        calc_info_P_joint(reshape(tbins_shuffle(istep:istep+step-1,:),1,numel(tbins(istep:istep+step-1,:))), ...
                        reshape(pbins_shuffle(istep:istep+step-1,:),1,numel(pbins(istep:istep+step-1,:))),nBins,nyBins,[1 0.95 0.9 0.85 0.8 0.5],30);

            end

        end

    end



end

%

figure;imagesc(I{1})
tlag=2;

figure;h=errorbar([1:length(I_shuffle{1})]+(all_win_dur(1))/2,I{1}(tlag,:),Ierr{1}(tlag,:),'r');
hold on
h2=errorbar([1:length(I_shuffle{1})]+(all_win_dur(1))/2,I_shuffle{1}(tlag,:),Ierr_shuffle{1}(tlag,:),'k');
xlabel('Time wrt beginning of perturbation (ms)');
ylabel('bits');







%%



R=cell(1);
S=cell(1);
for i=1:size(mysldspk_sm,2)
    [R{i},S{i}]=corrcoef(non_est_dir(:,i),mysldspk_sm(:,i));
end 



myr=zeros(1,size(mysldspk_sm,2));

for i=1:size(mysldspk_sm,2)
    myr(i)=R{i}(1,2);
end

figure;plot(myr)

%%
save experiment
save([experiment,'_infodata'],'stitrain','spkcount','tname')

      

      
      



