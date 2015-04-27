


clear all;

if ispc
    rt='Z:/';
else
    rt='/Volumes/Hub/';
end
path(path,[rt,'MT/MATLAB/bing_ana']);
path(path,[rt,'RMVideo_backup']);
path(path,[rt,'MattM/Code/population decoding']);

trialsdir=[rt,'MT/Data/ga022015/maestro'];
savedir=[rt,'MT/Data/ga022015/mat/'];
experiment='ga022015d_pert_adapt_dirspd45';
first=1;
last=1000;

seg_dur=500;  % important if we change the length of the stimulus.
tdur=3*seg_dur;

neuron_idx=1;

%

try
    load([savedir,'data.mat'])
catch
    [data] = load_maestrospks(trialsdir,experiment,first,last);
    save([savedir,'data.mat'],'data')
end

% save([savedir,'data.mat'],'data')
seg_dur=500;  % important if we change the length of the stimulus.
tdur=3*seg_dur;
nTypes = size(data,1);
nTags = 1;            % we should have 4 tags
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

segdata.nSpks = zeros(nTypes,nTags,nReps);
segdata.target = zeros(tdur,nTypes,nTags,nReps);
segdata.bspks = segdata.target;
segdata.spiketimes = zeros(200,nTypes,nTags,nReps); % 200 is bigger than the max spikes num
segdata.nfiles = zeros(1,nTypes);
tagval=[]; %tagvals identifies the trial condition for all segments within trials (so nTypes*nSegs*nReps long);
count=0;
for i=1:nTypes
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
%                 minT = cell2mat(data(i).tagStart(k));
%                 maxT = cell2mat(data(i).tagEnd(k));
%                 minT=data(i).tagStart(k);
%                 minT=minT{1};
%                 maxT=data(i).tagEnd(k);
%                 maxT=maxT{1};
            minT=401;
            maxT=tdur+minT-1;

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
%                 if i<=12
%                     pdirarray(j,:)=ptemp(j,:);
%                 else
%                     pdirarray(j,:) = mod(unwrap(ptemp(j,:)*pi/180)*180/pi,360);
%                 end 
%                pdirarray(j,:) = mod(unwrap(ptemp(j,:)*pi/180)*180/pi,360);
                pdirarray(j,:)=ptempdir(j,:);
%                npdirarray(j,:) = pdirarray(j,k,:) - prefDir; 
                npdirarray(j,:) = pdirarray(j,:) - repmat(mean(pdirarray(j,:)),1,size(ptempdir,2)); 
                pspdarray(j,:)=ptempspd(j,:);
%                npdirarray(j,:) = pdirarray(j,k,:) - prefDir; 
                npspdarray(j,:) = pspdarray(j,:) - repmat(mean(pdirarray(j,:)),1,size(ptempspd,2)); 
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
            stdpdir = [stdpdir repmat(sqrt(var(reshape(unwrap(npdirarray),1,nReps*size(ptempdir,2)))),1,nReps)];
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
%
% for i=1:nTypes
%     figure(1)
%     subplot(4,2,i)
%     plot(nstimdir{i}')
%     ylim([mean(mean(nstimdir{i}'))-90 mean(mean(nstimdir{i}))+90])
%     title([data(i,1).trname])
%     figure(2)
%     subplot(4,2,i)
%     plot(nstimspd{i}')
%     ylim([0 40])
% 
% title([data(i,1).trname])
% figure(3)
% 
%     subplot(4,2,i)
%     plot(rad2deg(circ_std(deg2rad(nstimdir{i})).^2))
%     title(data(i,1).trname{1})
%     ylim([0 50])
%     if i==3
%         ylabel('Variance of Stimulus Direction')
%     end
%     figure(4)
%     subplot(4,2,i)
%     plot(var(nstimspd{i}))
%     ylim([0 100])
%     title(data(i,1).trname{1})
%     if i==3
%         ylabel('Variance of Stimulus Speed')
%     end
% 
% end


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
markdou3=zeros(4,length(stsn));

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
    markdou3(:,i)=tmp;
end

%
spikes=nexFile.neurons{neuron_idx}.timestamps;   % pay attention here, choice the unit you want

%
spk_tr=cell(3,length(stsn));
spk_fortrans=cell(1,length(stsn));
for i=1:length(stsn)
    %1st pert
    tmp=spikes;
    tmp=tmp(tmp>markdou3(2,i));
    tmp=tmp(tmp<markdou3(3,i));
    spk_tr{1,i}=1000.*(tmp-markdou3(2,i));
    
    %2nd pert
    tmp=spikes;
    tmp=tmp(tmp>markdou3(3,i));
    tmp=tmp(tmp<markdou3(4,i));
    spk_tr{2,i}=1000*(tmp-markdou3(3,i));
    
    %3rd pert
    tmp=spikes;
    tmp=tmp(tmp>markdou3(4,i));
    tmp=tmp(tmp<markdou2(2,i));
    spk_tr{3,i}=1000.*(tmp-markdou3(4,i));
    
    %all pert segs
    tmp=spikes;
    tmp=tmp(tmp>markdou3(2,i));
    tmp=tmp(tmp<markdou2(2,i));
    spk_fortrans{i}=1000.*(tmp-markdou3(2,i));
end
% have a look of the spk_tr, the firing rate seems too low. maybe we will
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

% now we got three important things, the first is the firing rate of the
% neuron, spk_tr,   and the second is the trialnames, tri_num. And the
% perts of the stimulus, pert.
 
%


%

numlevels=2; %HTL, LTH
numconditions=2; %dir, spd

ntype=unique(tri_num);
tris=length(tri_num);
n_eachtri=zeros(1,length(ntype));
spktri=cell(1,length(ntype));
transpk=cell(1,length(ntype));
stitrain_dir=cell(1,length(ntype));

%trial order: HLL_dir, HLL_spd, LHH_dir, LHH_spd
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

cumn=1;
stimtrain_dir=cell(1,length(ntype));
%
for i=ntype
    tp=n_eachtri(i);
    for j=1:tp
        stimtrain_dir{i}=[stimtrain_dir{i}; nstimdir{i}(j,:)];
%         cumn=cumn+1;
    end
end
%
stimtrain_spd=cell(1,length(ntype));

cumn=1;
for i=ntype
    tp=n_eachtri(i);
    for j=1:tp
        stimtrain_spd{i}=[stimtrain_spd{i}; nstimspd{i}(j,:)];
%         cumn=cumn+1;
    end
end
%
%% sandbox for direction-speed STA
pairs=[{'dirL_nospd'},{[1 3]},{{[2,3],[1]}};{'dirH_nospd'},{[1 3]},{{[1],[2,3]}};...
        {'spdL_nodir'},{[2 4]},{{[2,3],[1]}};{'spdH_nodir'},{[2 4]},{{[1],[2,3]}}];
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
    stdlev(i,1)=rad2deg(circ_std(deg2rad(reshape(stimspk{i,1},1,[])),[],[],2));
    stdlev(i,2)=std(reshape(stimspk{i,2},1,[]),[],2);

    ntr=size(spkseg_tmp,2);
    for k=1:ntr
        stimspk{i,3}(k,1:length(spkseg_tmp{k}))=round(spkseg_tmp{k});
    end
end
save([experiment,'_stimspk.mat'],'stimspk','stdlev')
% blorp
% dirspd_sta


%
% direction linfilt
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
%
colors=distinguishable_colors(numtriTypes);
trinums=[1,2];
%


%
for i=trinums
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
%         myfilter=[filtdir(pshift,myidx).filter_allt]*1000;
%         [amp_filter(j,i) tempt]=min(mean(myfilter,2));
%        
%         tmpp=std(myfilter');
%         if length(myidx)>1
%            std_filter(j,i)=tmpp(tempt);
%         end
%            m_stdsti(j,i)=mean(std(stimseg{i,j}'));
%            errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',colors(i,:)); 
%            hold on; 
       
        linfilt_dir(i,j).results=filtdir;
        linfilt_dir(i,j).eyeest=eye_est;
        linfilt_dir(i,j).index2=index2_list;
        linfilt_dir(i,j).residuals=residuals;
        linfilt_dir(i,j).ccoef=ccoef;
        linfilt_dir(i,j).error=error;
        linfilt_dir(i,j).idx=myidx;
%         filter(i,j).coef=mean(ccoef.allt(pshift,myidx));
end

% set(gca,'FontSize',12,'Box','Off','TickDir','Out');
% xlim([-50 250]);
% %     legend(['20 deg ccoef: ',num2str(filter(1).coef)],['60 deg ccoef: ',num2str(filter(2).coef)]);
% %
% [legh,~,~,~]=legend(legcell{trinums});   
% set(legh,'Interpreter','none','Location','EastOutside','FontSize',14)
% lc=get(legh,'Children');
% for i=[1 3 5 7 9 11]
%     ts=get(get(lc(i),'Children'),'Children');
%     set(ts(2),'LineWidth',2)
% end
% %
% xlabel('Time lag (ms)');
% ylabel('Filter Amplitude');
% axis square;
% title('Direction');
% box off;
% saveas(gcf,[experiment,'linfilt_dir.fig'])
save([experiment,'_linfilt_dir.mat'],'linfilt_dir','pairs','experiment')



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
%
trinums=3:4;
ct=1;

%
for i=trinums
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
%         myfilter=[filtspd(pshift,myidx).filter_allt]*1000;
%         [amp_filter(j,i) tempt]=min(mean(myfilter,2));
%        
%         tmpp=std(myfilter');
%         if length(myidx)>1
%            std_filter(j,i)=tmpp(tempt);
%         end
%            m_stdsti(j,i)=mean(std(stimseg{i,j}'));
%            errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',colors(i,:)); 
%            hold on; 
       
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
    %
% [legh,~,~,~]=legend(legcell); 
% legh=findobj(gcf,'Type','axes','Tag','legend');
% % set(legh,'Interpreter','none','Location','EastOutside','FontSize',14)
% lc=get(legh,'Children');
% for i=[1 3 5 7 9 11]
%     ts=get(get(lc(i),'Children'),'Children');
%     set(ts(2),'LineWidth',3)
% end
% %
%     xlabel('Time lag (ms)');
%     ylabel('Filter Amplitude (1/s)');
%     axis square;
%     title('Speed');
%     box off;
%     %
% saveas(gcf,[experiment,'linfilt_spd.fig'])
save([experiment,'_linfilt_spd.mat'],'linfilt_spd','pairs','stdlev','experiment')

% to replot separately

h1=figure;
subplot 211
hold all
numtriTypes=size(pairs,1);
pshift=1; 
ft=250;
lag=[1:ft];
colors=distinguishable_colors(numtriTypes);
inds=1:2;
j=1;
for i=1:numtriTypes
    if i==1||i==2
        legcell{i}=[sprintf('%2.0f',stdlev(i,1)),' deg / 0 dps'];
    elseif i==3||i==4
        legcell{i}=['0 deg / ',sprintf('%2.0f',stdlev(i,2)),' dps'];
%     else
%         legcell{i}=[sprintf('%2.0f',stdlev(i,1)),' deg / ',sprintf('%2.0f',stdlev(i,2)),' dps'];
    end
end
%
for i=inds
    myfilter=[linfilt_dir(i).results(pshift,:).filter_allt]*1000;
    [amp_filter(j,i) tempt]=min(mean(myfilter,2));
    tmpp=std(myfilter');
%     if length(myidx)>1
%         std_filter(j,i)=tmpp(tempt);
%     end
%     m_stdsti(j,i)=mean(std(stimseg{i,j}'));
    errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',colors(i,:));
end
legh=legend(legcell{inds});
lc=get(legh,'Children');
set(legh,'Location','SouthEast','FontSize',14)
for i=[1 3]
    ts=get(get(lc(i),'Children'),'Children');
    set(ts(2),'LineWidth',3)
end
title('Direction Filters','FontSize',18)

ylabel('Filter Amplitude (spikes/s ^2deg)','FontSize',18)

xlabel('Time Lag (ms)','FontSize',16)

figure(h1);
subplot 212
hold all
colors=distinguishable_colors(numtriTypes);
inds=[3,4];
for i=inds
    myfilter=[linfilt_spd(i).results(pshift,:).filter_allt]*1000;
    [amp_filter(j,i) tempt]=min(mean(myfilter,2));
    tmpp=std(myfilter');
%     if length(myidx)>1
%         std_filter(j,i)=tmpp(tempt);
%     end
%     m_stdsti(j,i)=mean(std(stimseg{i,j}'));
    errorbar(lag,mean(myfilter(lag,:),2),std(myfilter(lag,:),[],2),'Color',colors(i,:));
end
legh=legend(legcell{inds});
lc=get(legh,'Children');
set(legh,'Location','SouthEast','FontSize',14)
for i=[1 3]
    ts=get(get(lc(i),'Children'),'Children');
    set(ts(2),'LineWidth',3)
end
title('Speed Filters','FontSize',18)

ylabel('Filter Amplitude (spikes/s-deg)','FontSize',18)

xlabel('Time Lag (ms)','FontSize',16)



% save experiment
% save([experiment,'_infodata'],'stitrain','spkcount','tname')

      
      
      



