close all
try
    load([savedir,'data.mat'])
catch
    [data] = load_maestrospks(trialsdir,experiment,first,last);
    save([savedir,'data.mat'],'data')
end

%
seg_dur=400;  % important if we change the length of the stimulus.
tdur=6*seg_dur; 

nTypes = size(data,1);
nTags = 1;           
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
  [nexFile] = readNexFile(nexdir)
  save([savedir,'nexfile.mat'],'nexFile')
end
%

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
numspds=length(spds);
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
spk_nf=cell(length(trialdirs),numspds);
tri_type=ceil(tri_num/8);

for i=1:length(tri_num)
    for j=1:size(spk_tr,1); %for each seg
        spdind=trispdind(i);
        triind=tri_type(i);
        spk_nf{segdirinds{triind}(j),spdind}=[spk_nf{segdirinds{triind}(j),spdind}; {spk_tr{j,i}}];
    end
end
% wavelet_info
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

% %to shift responses to make cell non-separable
% for i=1:size(spk_tp,2)
%     spk_tp(:,i)=circshift(spk_tp(:,i),i);
% end

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

for dir=1:numdirs
    for spd=1:numspds
%         for trial=1:size(spk_nf{dir,spd})
%             spkct{dir,spd}(trial)=length(spk_nf{dir,spd}{trial});
%         end
        spk_ct_mean(dir,spd)=mean(sum(spk_tp{dir,spd},2)); %cell of total spks per trial per dir/spd needed?
    end
end
figure;imagesc(trialdirs_rot,spds,spk_ct_mean');colormap('jet');colorbar
set(gca,'YTick',[(2:13.5:13.5*8+2).*str2num(num2str(prefSpd))],'YTickLabel',spds)
xlabel('Direction');ylabel('Speed')
%%
savedir=[rt,'MT/MATLAB/matt_ana/Info/',experiment(1:8)];

mkdir(savedir);
save([savedir,'/',experiment,'2dtunecurve.mat'],'spk_ct_mean');
%
saveas(gcf,[savedir,'/',experiment,'_unit ',num2str(neuron_idx),'_2dtunecurve.fig'])

figure;hold all
for i=1:numspds
    plot(trialdirs_rot,spk_ct_mean(:,i))%./max(max(spk_ct_mean)))
end
legend(cellstr(num2str(spds')));
title([experiment,' unit',num2str(neuron_idx),'Direction Tuning'])
saveas(gcf,[savedir,'/',experiment,'_unit ',num2str(neuron_idx),'_dirtune.fig'])
%
colors=distinguishable_colors(numdirs);
figure;
for i=1:numdirs
    semilogx(spds,squeeze(spk_ct_mean(i,:)),'Color',colors(i,:));hold all %./max(max(spk_ct_mean))
end
set(gca,'XTick',[1 2 4 8 16 32 64 96],'XTickLabel',[1 2 4 8 16 32 64 96])
legend(cellstr(num2str(trialdirs_rot')),'Location','EastOutside')
saveas(gcf,[savedir,'/',experiment,'_unit ',num2str(neuron_idx),'_spdtune.fig'])
%
% bin
binsize=20;
colorsdir=distinguishable_colors(numdirs);
colorsspd=distinguishable_colors(numspds);
isi=cell(numdirs,numspds);
ct_bin=cell(numdirs,numspds);
cumct_bin=cell(numdirs,numspds);
figure(100);
clf
figure(101);
clf
for dir=1:numdirs
    for spd=1:numspds
        for trial=1:size(spk_tp{dir,spd},1)
            for t=1:size(spk_tp{dir,spd},2)-binsize
%             bint=(t-1)*binsize+1;
                ct_bin{dir,spd}(trial,t+binsize)=sum(spk_tp{dir,spd}(trial,t:t+binsize)); %/200ms*1000ms=spks/s
                tunecurve_2d(dir,spd,t+binsize)=mean(ct_bin{dir,spd}(:,t+binsize));   
            end
            isi{dir,spd}=[isi{dir,spd};spk_nf{dir,spd}{trial}(2:end)-spk_nf{dir,spd}{trial}(1:end-1)];
        end
        cumct_bin{dir,spd}=cumsum(spk_tp{dir,spd},2);  
        figure(100)
        subplot(2,4,spd)
        plot(dir+mean(ct_bin{dir,spd},1),'Color',colorsdir(dir,:));
        hold all
        title([num2str(spds(spd)),' dps'])
        if sum(spd==[1 5])
            ylabel('Spikes/bin');
        end
        figure(101)
        subplot(3,4,dir)
        plot(spd+mean(ct_bin{dir,spd},1),'Color',colorsspd(spd,:));
        hold all
        title([num2str(trialdirs_rot(dir)),' degrees'])
        if sum(dir==[1 5 9])
            ylabel('Spikes/bin');
        end
        
    end
end
maxfig(100,1)
saveas(gcf,[savedir,'/',experiment,'_unit ',num2str(neuron_idx),'_psth_spd.fig'])
maxfig(101,1)
saveas(gcf,[savedir,'/',experiment,'_unit ',num2str(neuron_idx),'_psth_dir.fig'])
%%
close all
figure
writerObj=VideoWriter([savedir,'/',experiment,'2dtunecurve.avi']);
set(writerObj,'FrameRate',10,'Quality',30)
open(writerObj)
for t=1:400
    imagesc(trialdirs_rot,spds,tunecurve_2d(:,:,t)',[0 max(max(max(tunecurve_2d)))]);colormap('jet');colorbar
    set(gca,'YTick',[2:13.5:13.5*8+2],'YTickLabel',spds)
    xlabel('Direction','FontSize',16);ylabel('Speed','FontSize',16);title([num2str(t),' ms'],'FontSize',16)
    set(gcf, 'Position', [100, 100, 1000, 800]);    
    F=getframe(gcf);
    writeVideo(writerObj,F);
end

close(writerObj)
bork
    
%% %%
for kind=1:3
    if kind==1
        clear I_dir I_dir_comb I_dir_comb_std I_dir_mean I_dir_xspd I_dirspd_joint...
        I_dirspd_joint_1d I_dirspd_joint_1d_shuffle I_dirspd_joint_shuffle I_spd I_spd_comb I_spd_comb_std I_spd_mean I_spd_xdir Iinf Iinf_shuffle...
        Ixy Ixy_shuffle Ixyz Ixyz_shuffle Ntrialsfrac OUTH OUTM Pjoint Pjoint_shuffle Ptemp Ptemp2 Ptempy PtempyGx Ptempz PtempzGxy Px Pxy Sinf Sinf_shuffle...
        Snoise_inf Snoise_inf_shuffle Sy Sy_given_x Sy_given_x_shuffle Sy_shuffle Sz Sz_given_xy Sz_given_xy_shuffle Sz_given_y Sz_shuffle Xdata Xdata_r...
        Ydata Ydata_r Zdata Zdata_r
        response=cumct_bin;
        tag='from cumulative spike count';    
        info_jointdirspd
    elseif kind==2
        clear I_dir I_dir_comb I_dir_comb_std I_dir_mean I_dir_xspd I_dirspd_joint...
        I_dirspd_joint_1d I_dirspd_joint_1d_shuffle I_dirspd_joint_shuffle I_spd I_spd_comb I_spd_comb_std I_spd_mean I_spd_xdir Iinf Iinf_shuffle...
        Ixy Ixy_shuffle Ixyz Ixyz_shuffle Ntrialsfrac OUTH OUTM Pjoint Pjoint_shuffle Ptemp Ptemp2 Ptempy PtempyGx Ptempz PtempzGxy Px Pxy Sinf Sinf_shuffle...
        Snoise_inf Snoise_inf_shuffle Sy Sy_given_x Sy_given_x_shuffle Sy_shuffle Sz Sz_given_xy Sz_given_xy_shuffle Sz_given_y Sz_shuffle Xdata Xdata_r...
        Ydata Ydata_r Zdata Zdata_r
        response=ct_bin;
        tag=' from binned spike count';
        info_jointdirspd
    elseif kind==3
        clear I_dir I_dir_comb I_dir_comb_std I_dir_mean I_dir_xspd I_dirspd_joint...
        I_dirspd_joint_1d I_dirspd_joint_1d_shuffle I_dirspd_joint_shuffle I_spd I_spd_comb I_spd_comb_std I_spd_mean I_spd_xdir Iinf Iinf_shuffle...
        Ixy Ixy_shuffle Ixyz Ixyz_shuffle Ntrialsfrac OUTH OUTM Pjoint Pjoint_shuffle Ptemp Ptemp2 Ptempy PtempyGx Ptempz PtempzGxy Px Pxy Sinf Sinf_shuffle...
        Snoise_inf Snoise_inf_shuffle Sy Sy_given_x Sy_given_x_shuffle Sy_shuffle Sz Sz_given_xy Sz_given_xy_shuffle Sz_given_y Sz_shuffle Xdata Xdata_r...
        Ydata Ydata_r Zdata Zdata_r
        response=isi;
        tag=' from ISI';
        info_jointdirspd
    end
end