clear all;

if ispc
    rt='Z:/';
else
    rt='/Volumes/Hub/';
end
path(path,[rt,'MT/MATLAB/bing_ana']);
path(path,[rt,'RMVideo_backup']);
trialsdir=[rt,'MT/Data/ga032715/maestro'];
savedir=[rt,'MT/Data/ga032715/mat/'];
experiment='ga032714a_dir15_5spds'; %yes I screwed this date up in recording

first=1;
last=1944;



% %
try
    load([savedir,'data.mat'])
catch
    [data] = load_maestrospks(trialsdir,experiment,first,last);
    save([savedir,'data.mat'],'data')
end

%
seg_dur=400;  % important if we change the length of the stimulus.
tdur=6*seg_dur; %-250 to remove transient at beginning

nTypes = size(data,1);
nTags = 1;            % we should have 4 tags ??
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
% segdata.spiketimes = zeros(200,nTypes,nTags,nReps); % 200 is bigger than the max spikes num
segdata.nfiles = zeros(1,nTypes);
tagval=[]; %tagvals identifies the trial condition for all segments within trials (so nTypes*nSegs*nReps long);
count=0;

for i=[1:nTypes] %skip dir045_HTL_LTL for now: too many issues with putting short pert seg in matrices
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
            

            minT=201;  
            maxT=tdur+minT+200-1;

            %prefDir = [data(1).targ.predir]; % direction that has been added to target motion
            %prefDir = prefDir(1);
            nstimdir{i} = tempdir(minT:maxT,1:nReps)';
%             pdirarray = zeros(nReps,size(ptempdir,2));
%             npdirarray = pdirarray;
            nstimspd{i} = tempspd(minT:maxT,1:nReps)';
%             pspdarray = zeros(nReps,size(ptempspd,2));
%             npspdarray = pspdarray;
            
%             bspks = pdirarray;
            %bspks = zeros(nReps,length(minT:maxT));
%             stimes=[];
%             for j=1:nReps
%                 pdirarray(j,:)=ptempdir(j,:);
% %                 npdirarray(j,:) = pdirarray(j,:);
%                 %centered around mean direction 135
%                 npdirarray(j,:) = pdirarray(j,:) - repmat(mean(pdirarray(j,:)),1,size(ptempdir,2)); 
%                 pspdarray(j,:)=ptempspd(j,:);
%                 %speeds here are 4,8,16,32,64 so mean speed analysis???
% %                 npspdarray(j,:) = pspdarray(j,:); %without removing mean speed
%                 npspdarray(j,:) = pspdarray(j,:) - repmat(mean(pspdarray(j,:)),1,size(ptempspd,2)); 
%                 % get spikes
%                 
%                 stmp = round(stemp(j,:));
%                 stmp = stmp(stmp>minT);
%                 stmp = stmp(stmp<maxT);
%                 nspks = length(stmp);
%                 stimes(j,1:nspks) = stmp;
%                 binspks = zeros(1,length(minT:maxT));
%                 %binspks([stmp-minT])==1;
%                 binspks([stmp-round(double(minT))])=1;
%                 bspks(j,1:length(minT:maxT)) = binspks;
%                 count=count+1;
%                 celltimes{count} = stmp-round(double(minT)); %store times relative to segment onset
%                 segdata.nSpks(i,k,j) = nspks;
%                 tagval(count) = i;
%             end
%             %stdpdir(i,k,1:nReps) = sqrt(var(reshape(npdirarray,1,nReps*size(ptemp,2))); 
%             % make an array the same length as stim and spike arrays with the std levels of the perturbations
%             stdpdir = [stdpdir repmat(sqrt(var(reshape(rad2deg(unwrap(deg2rad(npdirarray))),1,nReps*size(ptempdir,2)))),1,nReps)];
%             nstimdir{i} = npdirarray;
%             nstimspd{i}=npspdarray;
%             spikearray = [spikearray ; bspks];
            %spiketimes{(i-1)*nTags+k} = stimes; 
%             segdata.target(1:length(minT:maxT),i,k,1:nReps) = reshape(pdirarray',length(minT:maxT),1,1,nReps); 
%             segdata.bspks(1:length(minT:maxT),i,k,1:nReps) = reshape(bspks',length(minT:maxT),1,1,nReps); 
%             segdata.spiketimes(1:size(stimes,2),i,k,1:nReps) = reshape(stimes',size(stimes,2),1,1,nReps);
%             segdata.nfiles(i) = nReps;

        end
end 
%
path(path,[rt,'/MattM/Code/population decoding'])

%
%
try
    load([savedir,'nexfile.mat'])
catch
  [nexFile] = readNexFile()
  save([savedir,'nexfile.mat'],'nexFile')
end
%
%primary cell: unit 2
neuron_idx=2;
%secondary cell: unit 1 spikes safely iso'd and extracted from trial ~500-700 to end
% neuron_idx=1;

ts=nexFile.markers{1}.timestamps;
ev=nexFile.markers{1}.values{1}.strings;
n_ev=length(ev);
 
st=zeros(1,n_ev);
for i=1:n_ev
    st(i)=str2num(ev{i})==2;
end
st1=find(st==1);


st=zeros(1,n_ev);
for i=1:n_ev
    st(i)=str2num(ev{i})==3;
end
ed1=find(st==1);


st=zeros(1,n_ev);
for i=1:n_ev
    st(i)=str2num(ev{i})==6;
end
savd=find(st==1);


trialnames=cell(1,length(st1));
maestronames=cell(1,length(st1));
ind=0;
for i=1:length(ed1)
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
%
for i=1:length(stsn)
    triname{i}=trialnames{stsn(i)};
end
%
prefSpd=data(1,1).targ(1,1).prefSpd; %all speeds multiplied by Target Velocity Scale
prefDir=data(1,1).targ(1,1).prefDir; %all directions rotated by Target Velocity Rot

tri_num=zeros(1,length(triname));
tname=unique(triname);

for i=1:length(triname)    
    tri_num(i)=find(strcmp(tname,triname{i}));
    ind=strfind(triname{i},'tria1'); %because this trialname was screwed up
    if ind>0
        triname{i}=[' trial1',triname{i}(7:end)];
%         tname{tri_num(i)}=[' trial1',tname{tri_num(i)}(7:end)];
    end
    
    undind=8;
    dind=strfind(triname{i},'d');
    trispds(i)=str2num(triname{i}(undind+1:dind-1))*prefSpd;
end

% fix this after spds are 02/04/08/16/32
% spds=[16 2 32 4 8];
spds=str2num(num2str(unique(trispds)));


for i=1:length(tri_num)
    trispdind(i)=find(spds==trispds(i));
end

stt=zeros(1,length(stsn));
endt=zeros(1,length(stsn));

for i=1:length(stsn)
    stt(i)=ts(st1(stsn(i)));
    tmp=ed1(ed1>st1(stsn(i)));
    endt(i)=ts(tmp(1));
end

%
startmarkers=zeros(2,length(stsn));
markdou2=zeros(2,length(stsn));
markdou3=zeros(7,length(stsn));
%
mark2=nexFile.events{1}.timestamps;
dou2=nexFile.events{2}.timestamps;
dou3=nexFile.events{3}.timestamps;

for i=1:length(stt)
    tmp=mark2(mark2>stt(i));
    tmp=tmp(tmp<endt(i));
    startmarkers(:,i)=tmp;
    
    tmp=dou2(dou2>stt(i));
    tmp=tmp(tmp<endt(i));
    markdou2(:,i)=tmp;
    
    tmp=dou3(dou3>stt(i));
    tmp=tmp(tmp<endt(i));
    markdou3(:,i)=tmp;
end

%
trialdirs=[225 270 315 330 345 0 15 30 45 90 135 180];
%trial direction indices by segment
segdirs={[315,345,225],[90,135,270],[330,0,15],[45,30,180]};
for i=1:4
    for j=1:3
        segdirinds{i}(j)=find(trialdirs==segdirs{i}(j));
    end
end
%for each trialtype
%for each seg
%get mean nstimdir and mean nstimspd to check

spikes=nexFile.neurons{neuron_idx}.timestamps;
%
spk_tr=cell(3,length(stsn));
for i=1:length(stsn)
    %1st seg
    tmp=spikes;
    tmp=tmp(tmp>markdou3(2,i)); 
    tmp=tmp(tmp<markdou3(3,i)); %+0.1 for each 100 ms post-stim end
    spk_tr{1,i}=1000.*(tmp-markdou3(2,i));
    
    %2nd seg
    tmp=spikes;
    tmp=tmp(tmp>markdou3(4,i));
    tmp=tmp(tmp<markdou3(5,i));
    spk_tr{2,i}=1000*(tmp-markdou3(4,i));
    
    %3rd seg
    tmp=spikes;
    tmp=tmp(tmp>markdou3(6,i));
    tmp=tmp(tmp<markdou3(7,i));
    spk_tr{3,i}=1000.*(tmp-markdou3(6,i));    
end
% have a look of the spk_tr, if the firing rate seems too low. maybe we will
% not get good results.


%{segdirinds{tri_type(i)}(seg),tri_spd} to index

        %
spk_nf=cell(length(trialdirs),length(spds));
tri_type=ceil(tri_num/5);

for i=1:length(tri_num)
    for j=1:size(spk_tr,1); %for each seg
        spdind=trispdind(i);
        triind=tri_type(i);
        spk_nf{segdirinds{triind}(j),spdind}=[spk_nf{segdirinds{triind}(j),spdind}; {spk_tr{j,i}}];
    end
end

%why is this not trials x 400?
tt=400;
spk_tp=cell(size(trialdirs,2),size(spds,2));
%
for d=1:length(trialdirs)
    for s=1:length(spds)
        ntri=length(spk_nf{d,s});
        spk_tp{d,s}=zeros(ntri,tt);
        for j=1:ntri
            spk_tp{d,s}(j,1:tt)=0;
            clear tp;
            tp=round(spk_nf{d,s}{j});
            if numel(tp)~=0
              if tp(1)==0
                   tp(1)=1;
                end
                spk_tp{d,s}(j,tp)=1;
            end
        end
    end
end
%
% spk_pt=cell(size(trialdirs,2),size(spds,2));
% for d=1:length(trialdirs)
%     for s=1:length(spds)
%         ntri=length(spk_tp{d,s});
%         for j=1:ntri
%             spk_pt{d,s}=[spk_pt{d,s}; spk_tp{d,s}{j}];
%         end
%     end
% end

%

% ind1=1;
% ind2=1;
% 
% for i=1:length(tname)
% %     ind=find(i==spd2_inds);
% %     if ~isempty(ind)
% %         spk_nf_2spd{ind}=spk_nf{i};
% %         spk_pt_2spd{ind}=spk_pt{i};
% %         spk_tp_2spd{ind}=spk_tp{i};
% %         ind2=ind2+1;
% %     else
%         spk_nf{ind1}=spk_nf{i};
%         spk_pt{ind1}=spk_pt{i};
%         spk_tp{ind1}=spk_tp{i};
%         ind1=ind1+1;
% %     end
% end

% rasters
% 
numdirs=size(spk_nf,1);
trialdirs_rot(1:5)=trialdirs(1:5)-360;
trialdirs_rot(6:12)=trialdirs(6:12);

mycolor=colormap(hsv);
% for spd=1:size(spk_nf,2);   
%     figure;
%     hold all
%     for dir=1:numdirs
%         numtri=length(spk_nf{dir,spd});            
%         for trial=1:length(spk_nf{dir,spd})
%             for spk=1:length(spk_nf{dir,spd}{trial})
%                 line([spk_nf{dir,spd}{trial}(spk),spk_nf{dir,spd}{trial}(spk)],[trialdirs_rot(dir)+15*trial/numtri,trialdirs_rot(dir)+15*trial/numtri+0.7*15/numtri],'Color',mycolor(floor(dir*64/numdirs),:),'LineWidth',1)
%             end
%         end     
% 
%     end
%     xlim([0 tt+300])
%     ylim([-180 180])
%     set(gca,'YTick',-180:45:180,'YTickLabel',-180:45:180)
% end
%
% % tuning curves
numdirs=size(spk_tp,1);
numspds=size(spk_tp,2);
% % trialdirs=linspace(0,360,numdirs);
% 
% for dir=1:numdirs
%     for spd=1:numspds
%         for trial=1:size(spk_nf{dir,spd})
%             spkct{dir,spd}(trial)=length(spk_nf{dir,spd}{trial});
%         end
%         spk_ct_mean(dir,spd)=mean(sum(spk_tp{dir,spd},2)); %cell of total spks per trial per dir/spd needed?
%         spk_cumct{dir,spd}=cumsum(spk_tp{dir,spd});
%     end
% end
% 
% 
% figure;hold all
% for i=1:5
%     plot(trialdirs_rot,spk_ct_mean(:,i)./max(max(spk_ct_mean)))
% end
% legend(cellstr(num2str(spds')));

%% bin
binsize=20;
% figure;
for dir=1:numdirs
    for spd=1:numspds
        for t=1:size(spk_tp{dir,spd},2)-binsize
%             bint=(t-1)*binsize+1;
            for trial=1:size(spk_tp{dir,spd},1)
                fr_bin{dir,spd}(trial,t)=sum(spk_tp{dir,spd}(trial,t:t+binsize))*2; %/200ms*1000ms=spks/s
            end
        end
%     plot(mean(fr_bin{dir,spd},1));
%     hold all
    end
end

%%
data_x=cell(1,numdirs);
data_y=cell(1,numdirs);
fracs=[1 0.9 0.8 0.5];
n_x=12;
n_y=20;
nreps=20;

for dir=1:numdirs
    for spd=1:numspds
        data_x{dir}=[data_x{dir};ones(size(fr_bin{dir,spd})).*spds(spd)];
        data_y{dir}=[data_y{dir};fr_bin{dir,spd}];
    end
    xdata=data_x{dir}';
    ydata=data_y{dir}';
    %%wtf why this reshaping problem?
    info_forarup
%     [I_spd{dir},I_spd_err_std{dir},I_spd_err_frac,I_spdR,I_spdR_err_std,I_spdR_err_frac,Pjoint, PjointR] = calc_info_P_joint(data_x{dir},data_y{dir},n_x,n_y,fracs,nreps);
end
%% stimulus-specific info from butts/goldman 2006 and mutual info of dir and spkct
timebin=20;
maxK=50;
for t=1:tt-timebin
    t
    for K=1:maxK
        clear pr pthetar pvr Hvr prv pv Hthetar prtheta ptheta tmp maxr rcount allresps
        i=1;
        for dir=1:numdirs
            for spd=1:numspds 
                ind= sub2ind(size(spkct),dir,spd);
                trialinds=(randperm(length(spkct{dir,spd}),floor(0.7*length(spkct{dir,spd}))));
                spkct_bs{ind}=[];
                for trial=trialinds;
                    spkct_bs{ind}=[spkct_bs{ind},length(find(spk_nf{dir,spd}{trial}<=(t+timebin)))];
                end
                tmp(i)=max(spkct_bs{ind});
                i=i+1;
            end
        end

        maxr=max(tmp);
        
        for dir=1:numdirs
            trialstim=cell(2,numdirs*numspds);
            allresps=cell(1,numdirs*numspds);
            for spd=1:numspds
                ind=sub2ind(size(spkct),dir,spd);
                for trial=1:length(spkct_bs{ind})
                    trialstim{1,ind}(end+1)=trialdirs_rot(dir);
                    trialstim{2,ind}(end+1)=spds(spd);
                end
                allresps{ind}=[allresps{ind},spkct_bs{ind}];    
                pstim(ind)=length(spkct_bs{ind});
            end
            pstim=pstim/sum(pstim);
            Hv=-sum(pstim.*log2(pstim));

            for r=1:maxr+1
                r_inds=find(allresps{ind}==r-1);
                rcount(r)=length(r_inds);
%                 stims_byspkct=trialstim{ind}(r_inds);
%                 pvr(dir,r,:)=hist(spds_byspkct,spds)/(eps+length(spds_byspkct));
%                 Hvr(dir,r)=-nansum(pvr(dir,r,:).*log2(pvr(dir,r,:)));
%                 isp(dir,r)=Hv-Hvr(spd,r);
            end
            pr(ind,:)=rcount./sum(rcount);
        end
        for dir=1:numdirs
            ptn=0;
            for spd=1:numspds
            end
            for spd=1:numspds
                ind= sub2ind(size(spkct),dir,spd);
                spdinds=find(trialstim{2,ind}==spds(spd));
                dirinds=find(trialstim{1,ind}==trialdirs_rot(dir));
                stiminds=intersect(spdinds,dirinds);
                r_bystim=allresps{ind}(stiminds);
                prstim(ind,:)=hist(r_bystim,0:maxr)/(eps+length(r_bystim));
    %             issi(spd,dir)=nansum(prtheta(spd,:,dir).*isp(spd,:));
                ptn=ptn+pstim(ind).*prstim(ind,:);
            end
            for spd=1:numspds    
                ind= sub2ind(size(spkct),dir,spd);                
                info(ind)=pstim(ind)*nansum(prstim(ind,:).*log2(prstim(ind,:)./(eps+ptn)));
    %             eff(spd,dir)=(-nansum(pr(spd,:).*log2(pr(spd,:)))-(-nansum(prtheta(spd,:,dir).*log2(prtheta(spd,:,dir)))))...
    %                 /(-nansum(pr(spd,:).*log2(pr(spd,:))));
            end
        end
        mutinfo_stim(t,K)=sum(info,1);
        eff_stim(t,K)=squeeze(mutinfo_spd(t,K))./(-nansum(pr.*log2(pr),2));
        
        clear pr pthetar pvr Hvr prv pv Hthetar prtheta ptheta tmp maxr isp rcount allresps
        i=1;
        for dir=1:numdirs
            for spd=1:numspds 
                trialinds=(randperm(length(spkct{dir,spd}),floor(0.7*length(spkct{dir,spd}))));
                spkct_bs{dir,spd}=[];
                for trial=trialinds;
                    spkct_bs{dir,spd}=[spkct_bs{dir,spd},length(find(spk_nf{dir,spd}{trial}<=(t+timebin)))];
                end
                tmp(i)=max(spkct_bs{dir,spd});
                i=i+1;
            end
        end

        maxr=max(tmp);

        for spd=1:numspds
            trialdir{spd}=[];
            allresps{spd}=[];
            for dir=1:numdirs
                for trial=1:length(spkct_bs{dir,spd})
                    trialdir{spd}(end+1)=trialdirs_rot(dir);
                end
                allresps{spd}=[allresps{spd},spkct_bs{dir,spd}];    
                ptheta(dir)=length(spkct_bs{dir,spd});
            end
            ptheta=ptheta/sum(ptheta);
            Htheta=-sum(ptheta.*log2(ptheta));

            for r=1:maxr+1
                r_inds=find(allresps{spd}==r-1);
                rcount(r)=length(r_inds);
                dirs_byspkct=trialdir{spd}(r_inds);
                pthetar(spd,r,:)=hist(dirs_byspkct,trialdirs_rot)/(eps+length(dirs_byspkct));
                Hthetar(spd,r)=-nansum(pthetar(spd,r,:).*log2(pthetar(spd,r,:)));
%                 isp(spd,r)=Htheta-Hthetar(spd,r);
            end
            pr(spd,:)=rcount./sum(rcount);
        end
        for spd=1:numspds
            ptn=0;
            for dir=1:numdirs
                dirinds=find(trialdir{spd}==trialdirs_rot(dir));
                r_bydir=allresps{spd}(dirinds);
                prtheta(spd,:,dir)=hist(r_bydir,0:maxr)/(eps+length(r_bydir));
                ptn=ptn+ptheta(dir).*prtheta(spd,:,dir);
    %             issi(spd,dir)=nansum(prtheta(spd,:,dir).*isp(spd,:));
            end
            for dir=1:numdirs
                info(spd,dir)=ptheta(dir)*nansum(prtheta(spd,:,dir).*log2(prtheta(spd,:,dir)./(eps+ptn)));
    %             eff(spd,dir)=(-nansum(pr(spd,:).*log2(pr(spd,:)))-(-nansum(prtheta(spd,:,dir).*log2(prtheta(spd,:,dir)))))...
    %                 /(-nansum(pr(spd,:).*log2(pr(spd,:))));
            end
        end
        mutinfo_dir(t,K,:)=sum(info,2);
        eff_dir(t,K,:)=squeeze(mutinfo_dir(t,K,:))./(-nansum(pr.*log2(pr),2));

        clear pr pthetar pvr Hvr prv pv Hthetar prtheta ptheta tmp maxr rcount allresps
        i=1;
        for dir=1:numdirs
            for spd=1:numspds 
                trialinds=(randperm(length(spkct{dir,spd}),floor(0.7*length(spkct{dir,spd}))));
                spkct_bs{dir,spd}=[];
                for trial=trialinds;
                    spkct_bs{dir,spd}=[spkct_bs{dir,spd},length(find(spk_nf{dir,spd}{trial}<=(t+timebin)))];
                end
                tmp(i)=max(spkct_bs{dir,spd});
                i=i+1;
            end
        end

        maxr=max(tmp);
        
        for dir=1:numdirs
            trialspd{dir}=[];
            allresps{dir}=[];
            for spd=1:numspds
                for trial=1:length(spkct_bs{dir,spd})
                    trialspd{dir}(end+1)=spds(spd);
                end
                allresps{dir}=[allresps{dir},spkct_bs{dir,spd}];    
                pv(spd)=length(spkct_bs{dir,spd});
            end
            pv=pv/sum(pv);
            Hv=-sum(pv.*log2(pv));

            for r=1:maxr+1
                r_inds=find(allresps{dir}==r-1);
                rcount(r)=length(r_inds);
                spds_byspkct=trialspd{dir}(r_inds);
                pvr(dir,r,:)=hist(spds_byspkct,spds)/(eps+length(spds_byspkct));
                Hvr(dir,r)=-nansum(pvr(dir,r,:).*log2(pvr(dir,r,:)));
%                 isp(dir,r)=Hv-Hvr(spd,r);
            end
            pr(dir,:)=rcount./sum(rcount);
        end
        for dir=1:numdirs
            ptn=0;
            for spd=1:numspds
            end
            for spd=1:numspds
                spdinds=find(trialspd{dir}==spds(spd));
                r_byspd=allresps{dir}(spdinds);
                prv(dir,:,spd)=hist(r_byspd,0:maxr)/(eps+length(r_byspd));
    %             issi(spd,dir)=nansum(prtheta(spd,:,dir).*isp(spd,:));
                ptn=ptn+pv(spd).*prv(dir,:,spd);
            end
            for spd=1:numspds            
                info(spd,dir)=pv(spd)*nansum(prv(dir,:,spd).*log2(prv(dir,:,spd)./(eps+ptn)));
    %             eff(spd,dir)=(-nansum(pr(spd,:).*log2(pr(spd,:)))-(-nansum(prtheta(spd,:,dir).*log2(prtheta(spd,:,dir)))))...
    %                 /(-nansum(pr(spd,:).*log2(pr(spd,:))));
            end
        end
        mutinfo_spd(t,K,:)=sum(info,1);
        eff_spd(t,K,:)=squeeze(mutinfo_spd(t,K,:))./(-nansum(pr.*log2(pr),2));
    end
end

%
% figure;subplot 211;hold all
% for spd=1:3
%     plot(trialdirs_rot,issi(spd,:))
% end
% legend(num2str(noiselevels'))
% title([experiment(1:9),' stimulus-specific info'])
% figure;hold all
% for spd=1:3
%     plot(0:maxr,isp(spd,:))
% end
% legend(num2str(noiselevels'))
% title('specific information of response')
%%

figure;subplot 211;hold all
for spd=1:numspds
    errorbar(1:size(mutinfo_dir,1),mean(mutinfo_dir(:,:,spd),2),std(mutinfo_dir(:,:,spd),[],2))
end
ylabel('bits')
xlabel('time (ms)')
title('mutual information between direction and spike count')
legend(cellstr(num2str(spds')),'Location','EastOutside');

colors=distinguishable_colors(numdirs);
subplot 212;

figure;
hold all
for dir=1:numdirs
    errorbar(1:size(mutinfo_spd,1),mean(mutinfo_spd(1:size(mutinfo_spd,1),:,dir),2),std(mutinfo_spd(1:size(mutinfo_spd,1),:,dir),[],2),'Color',colors(dir,:))
end
ylabel('bits')
xlabel('time (ms)')
title('mutual information between speed and spike count')
legend(cellstr(num2str(trialdirs_rot')),'Location','EastOutside');