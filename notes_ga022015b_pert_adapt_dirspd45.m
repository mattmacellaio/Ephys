


clear all;

if ispc
    rt='Z:/';
else
    rt='/Volumes/Hub/';
end
path(path,[rt,'MT/MATLAB/bing_ana']);
path(path,[rt,'RMVideo_backup']);

trialsdir=[rt,'MT/Data/ga022015/maestro'];
savedir=[rt,'MT/Data/ga022015/mat/'];
experiment='ga022015d_pert_adapt_dirspd45';
first=1;
last=1000;

seg_dur=500;  % important if we change the length of the stimulus.
tdur=3*seg_dur;

neuron_idx=1;

%

[data] = load_maestrospks(trialsdir,experiment,first,last);

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
                npdirarray(j,:) = pdirarray(j,:);  % - repmat(mean(pdirarray(j,:)),1,size(ptemp,2)); 
                pspdarray(j,:)=ptempspd(j,:);
%                npdirarray(j,:) = pdirarray(j,k,:) - prefDir; 
                npspdarray(j,:) = pspdarray(j,:);  % - repmat(mean(pdirarray(j,:)),1,size(ptemp,2)); 
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
%%
for i=1:nTypes
    figure(1)
    subplot(4,2,i)
    plot(nstimdir{i}')
    ylim([mean(mean(nstimdir{i}'))-90 mean(mean(nstimdir{i}))+90])
    title([data(i,1).trname])
    figure(2)
    subplot(4,2,i)
    plot(nstimspd{i}')
    ylim([0 40])

title([data(i,1).trname])
figure(3)

    subplot(4,2,i)
    plot(rad2deg(circ_std(deg2rad(nstimdir{i})).^2))
    title(data(i,1).trname{1})
    ylim([0 50])
    if i==3
        ylabel('Variance of Stimulus Direction')
    end
    figure(4)
    subplot(4,2,i)
    plot(var(nstimspd{i}))
    ylim([0 100])
    title(data(i,1).trname{1})
    if i==3
        ylabel('Variance of Stimulus Speed')
    end

end


%%

  [nexFile] = readNexFile()
  save([savedir,'nexfile.mat'],'nexFile')
ts=nexFile.markers{1}.timestamps;
%events
ev=nexFile.markers{1}.values{1}.strings;
n_ev=length(ev);
 
%%
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
%%

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

%%
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


%%

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

%
%% direction stimulus
cumn=1;
for i=1:length(ntype)
    tp=n_eachtri(i);
    for j=1:tp
        stitrain_dir{i}=[stitrain_dir{i}; nstimdir{i}(j,:)];
%         cumn=cumn+1;
    end
end

%

spkseg_dir=cell(1,numlevels); % 2 is the 2 perts.
stiseg_dir=cell(1,numlevels);   % [ L  ,  H ]


%low variance        
spkseg_dir{1}=[spkseg_dir{1} spktri{1}(2,:)];
spkseg_dir{1}=[spkseg_dir{1} spktri{1}(3,:)];
spkseg_dir{1}=[spkseg_dir{1} spktri{3}(1,:)];
%high variance
spkseg_dir{2}=[spkseg_dir{2} spktri{1}(1,:)];
spkseg_dir{2}=[spkseg_dir{2} spktri{3}(2,:)];
spkseg_dir{2}=[spkseg_dir{2} spktri{3}(3,:)];


stiseg_dir{1}=[stiseg_dir{1}; stitrain_dir{1}(:,seg_dur+1:2*seg_dur)];
stiseg_dir{1}=[stiseg_dir{1}; stitrain_dir{1}(:,2*seg_dur+1:3*seg_dur)];
stiseg_dir{1}=[stiseg_dir{1}; stitrain_dir{3}(:,1:seg_dur)];

stiseg_dir{2}=[stiseg_dir{2}; stitrain_dir{1}(:,1:seg_dur)];
stiseg_dir{2}=[stiseg_dir{2}; stitrain_dir{3}(:,seg_dur+1:2*seg_dur)];
stiseg_dir{2}=[stiseg_dir{2}; stitrain_dir{3}(:,2*seg_dur+1:3*seg_dur)];

stseg_dir=stiseg_dir;
skseg_dir=cell(1,numlevels);
for i=1:1
    for j=1:numlevels
        tmp=spkseg_dir{i,j};
        ntr=length(tmp);
        for k=1:ntr
            skseg_dir{i,j}(k,1:length(spkseg_dir{i,j}{k}))=round(spkseg_dir{i,j}{k});
        end
    end
end

%%%%
% stseg{2,1}=stseg{2,1}(20:end,:);
% skseg{2,1}=skseg{2,1}(20:end,:);
%%%%


%% speed stimulus
stitrain_spd=cell(1,length(ntype));

cumn=1;
for i=1:length(ntype)
    tp=n_eachtri(i);
    for j=1:tp
        stitrain_spd{i}=[stitrain_spd{i}; nstimspd{i}(j,:)];
%         cumn=cumn+1;
    end
end

spkseg_spd=cell(1,numlevels); % 2 is the 2 perts.
stiseg_spd=cell(1,numlevels);   % [ L  ,  H ]


        
spkseg_spd{1}=[spkseg_spd{1} spktri{2}(2,:)];
spkseg_spd{1}=[spkseg_spd{1} spktri{2}(3,:)];
spkseg_spd{1}=[spkseg_spd{1} spktri{4}(1,:)];

spkseg_spd{2}=[spkseg_spd{2} spktri{2}(1,:)];
spkseg_spd{2}=[spkseg_spd{2} spktri{4}(2,:)];
spkseg_spd{2}=[spkseg_spd{2} spktri{4}(3,:)];


stiseg_spd{1}=[stiseg_spd{1}; stitrain_spd{2}(:,seg_dur+1:2*seg_dur)];
stiseg_spd{1}=[stiseg_spd{1}; stitrain_spd{2}(:,2*seg_dur+1:3*seg_dur)];
stiseg_spd{1}=[stiseg_spd{1}; stitrain_spd{4}(:,1:seg_dur)];

stiseg_spd{2}=[stiseg_spd{2}; stitrain_spd{2}(:,1:seg_dur)];
stiseg_spd{2}=[stiseg_spd{2}; stitrain_spd{4}(:,seg_dur+1:2*seg_dur)];
stiseg_spd{2}=[stiseg_spd{2}; stitrain_spd{4}(:,2*seg_dur+1:3*seg_dur)];

stseg_spd=stiseg_spd;
skseg_spd=cell(1,numlevels);
for i=1:1
    for j=1:numlevels
        tmp=spkseg_spd{i,j};
        ntr=length(tmp);
        for k=1:ntr
            skseg_spd{i,j}(k,1:length(spkseg_spd{i,j}{k}))=round(spkseg_spd{i,j}{k});
        end
    end
end
%%
% get the raster figure;
% using the skseg;

skseg=skseg_dir;
stseg=stseg_dir;

lft=[0.05 0.375 0.70];
b_btm=[0.70 0.375 0.05];  
wth=0.25;
figure;
for i=1:1
    for j=1:2
        tempd=skseg_dir{i,j};
        hgt=0.25/size(tempd,1);
        hold on;
        k=1;
            btm=b_btm(j)+(k-1)*hgt;
             subplot('position',[lft(i) btm wth hgt]);
%            subplot(size(tempd,1),1,1);
            tpd=tempd(k,tempd(k,:)~=0);
            plot([tpd;tpd],[2*ones(size(tpd));zeros(size(tpd))],'r','LineWidth',2);
            axis([0 seg_dur 0 2])
            set(gca,'TickDir','out')
            set(gca,'YTick', [])
%             set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
            set(gca,'Color',get(gcf,'Color')) % match figure background
            set(gca,'YColor',get(gcf,'Color'))
            box off
            xlabel('time (ms)') 
            
        for k=2:size(tempd,1)
            btm=b_btm(j)+(k-1)*hgt;
%              subplot('position',[lft(i) btm wth hgt]);
           subplot(size(tempd,1),1,k);
            tpd=tempd(k,tempd(k,:)~=0);
            plot([tpd;tpd],[2*ones(size(tpd));zeros(size(tpd))],'r','LineWidth',2);
            hold on;
            axis([0 seg_dur 0 2])
            set(gca,'XTick',[])
            set(gca,'YTick', []) % don't draw y-axis ticks
%             set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
            set(gca,'Color',get(gcf,'Color')) % match figure background
            set(gca,'YColor',get(gcf,'Color')) % hide the y axis
            set(gca,'XColor',get(gcf,'Color'))
            box off
            if k==size(tempd,1)
                title(['raster of dir',num2str(i),'pert',num2str(j*20)]);
            end
        end
    end
end
% suptitle('raster of bu112912c pert')



%%

skseg=skseg_dir;
stseg=stseg_dir;
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
            [sta1,STA1,STC1] = get_sta(skseg{i,j}(ind_use,:), stseg{i,j}(ind_use,:));
            sta(:,rp,i,j)=sta1;
        end
    end
end

%%%%
% stseg{2,1}=stseg{2,1}(20:end,:);
% skseg{2,1}=skseg{2,1}(20:end,:);
%%%%

figure;
% subplot(1,2,1);
errorbar([-199:200],-mean(sta(:,:,1,1)'),std(sta(:,:,1,1)'));
hold on;
errorbar([-199:200],-mean(sta(:,:,1,2)'),std(sta(:,:,1,2)'),'r');
set(gca,'FontSize',16,'Box','Off','TickDir','Out');
legend('20 deg','60 deg');
xlabel('Time lag (ms)');
ylabel('Spike Triggered Average (deg)');
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



%

% we will get the filters of the spikes to the stimulus
%  here we can use the stseg as the targets and the skseg as the spikes.

tnoise=cell(1,2);
spkcounts=cell(1,2);
sldspk=cell(1,2);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:1
    for j=1:2
        tnoise{i,j}=stseg{i,j}';
        
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
maxK=2;
pshift=1; 
rag=[1:ft];

for i=1:1
    for j=1:2
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
    end
end

mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error

amp_filter=zeros(1,2); 
std_filter=zeros(1,2);
m_stdsti=zeros(1,2);
figure;
i=1;
    hold on;
    subplot(1,2,1);
    for j=1:2
        sti=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [linfilt_results, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(sti,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if j==2
                break;
            end
            continue;
        end
            myfilter=[linfilt_results(pshift,myidx).filter_allt]*1000;
            [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stseg_dir{i,j}'));
           errorbar(rag,-mean(myfilter(rag,:),2),std(myfilter(rag,:),[],2),'Color',mycolor(-30+(31*j),:)); 
           hold on; 
       
        filter(i,j).results=linfilt_results;
        filter(i,j).eyeest=eye_est;
        filter(i,j).index2=index2_list;
        filter(i,j).residuals=residuals;
        filter(i,j).ccoef=ccoef;
        filter(i,j).error=error;
        filter(i,j).idx=myidx;
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
xlabel('std of stimulus (deg)');
ylabel('Filter Amplitude');
axis square;
box off;

%%

skseg=skseg_spd;
stseg=stseg_spd;
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
            [sta1,STA1,STC1] = get_sta(skseg{i,j}(ind_use,:), stseg{i,j}(ind_use,:));
            sta(:,rp,i,j)=sta1;
        end
    end
end

%%%%
% stseg{2,1}=stseg{2,1}(20:end,:);
% skseg{2,1}=skseg{2,1}(20:end,:);
%%%%

figure;
% subplot(1,2,1);
errorbar([-199:200],-mean(sta(:,:,1,1)'),std(sta(:,:,1,1)'));
hold on;
errorbar([-199:200],-mean(sta(:,:,1,2)'),std(sta(:,:,1,2)'),'r');
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



%

% we will get the filters of the spikes to the stimulus
%  here we can use the stseg as the targets and the skseg as the spikes.

tnoise=cell(1,2);
spkcounts=cell(1,2);
sldspk=cell(1,2);
bin=20;
nbin=1:(seg_dur-bin+1);

for i=1:1
    for j=1:2
        tnoise{i,j}=stseg{i,j}';
        
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
maxK=2;
pshift=1; 
rag=[1:ft];

for i=1:1
    for j=1:2
        temp=sldspk{i,j};
        sldspk{i,j}=sldspk{i,j}-repmat(mean(sldspk{i,j},1),size(sldspk{1,1},1),1);
        tnoise{i,j}=tnoise{i,j}-repmat(mean(tnoise{i,j},1),size(tnoise{1,1},1),1);
    end
end

mycolor=colormap;
clear linfilt_results eye_est index2_list residuals cceof error

amp_filter=zeros(1,2); 
std_filter=zeros(1,2);
m_stdsti=zeros(1,2);
figure;
i=1;
    hold on;
    subplot(1,2,1);
    for j=1:2
        sti=tnoise{i,j}(bin:end,:); %(1:end-bin+1,:);
        spk=sldspk{i,j};
        [linfilt_results, eye_est, index2_list, residuals, ccoef, error]=fget_linfilt(sti,spk,shift,ft,cutoff,maxK);
        myidx=1:2;  %=find(ccoef.allt(pshift,:)>=0.01);
        if isempty(myidx)
            if j==2
                break;
            end
            continue;
        end
            myfilter=[linfilt_results(pshift,myidx).filter_allt]*1000;
            [amp_filter(j,i) tempt]=min(mean(myfilter,2));
       
        tmpp=std(myfilter');
        if length(myidx)>1
           std_filter(j,i)=tmpp(tempt);
        end
           m_stdsti(j,i)=mean(std(stseg_dir{i,j}'));
           errorbar(rag,-mean(myfilter(rag,:),2),std(myfilter(rag,:),[],2),'Color',mycolor(-30+(31*j),:)); 
           hold on; 
       
        filter(i,j).results=linfilt_results;
        filter(i,j).eyeest=eye_est;
        filter(i,j).index2=index2_list;
        filter(i,j).residuals=residuals;
        filter(i,j).ccoef=ccoef;
        filter(i,j).error=error;
        filter(i,j).idx=myidx;
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


%%

load('C:\Data\online_data\bu051613a_dir32.mat')
psh_glob.unitIDs
opti_ori=225;


uniidnum=1;
rate=psh_glob.meanrate(:,1,uniidnum);
dir   = [0 45 90 135 180 225 270 315]/180*pi; 
stidir=[opti_ori-45 opti_ori opti_ori+45]/180*pi;
maxrate=max(rate);

mydir=[6,5,7];

hold on;
subplot(2,3,4);
polar([dir dir(1)], [rate(1:8)' rate(1)]);
hold on;polar(stidir(2)*ones(size(1:maxrate)),1:maxrate,'r--');
hold on;polar(dir(mydir),rate(mydir)','ro');


%%
% get the info data for leslie.

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

% stitrain is the sti for this data.
%%



save experiment
save([experiment,'_infodata'],'stitrain','spkcount','tname')

      
      
      



