%
if kind==1
    savedir=[rt,'MT/MATLAB/matt_ana/Info/',experiment(1:8),'/cumct/'];
elseif kind==2
    savedir=[rt,'MT/MATLAB/matt_ana/Info/',experiment(1:8),'/binct/'];
elseif kind==3
    savedir=[rt,'MT/MATLAB/matt_ana/Info/',experiment(1:8),'/isi/'];
end
numfracreps=20;

% %spd info
% data_x=cell(1,numdirs);
% data_y=cell(1,numdirs);
% fracs=[1 0.9 0.8 0.5];
% colors=distinguishable_colors(numdirs);
% 
% for dir=1:numdirs
%     for spd=1:numspds
%         clear tmp
%         data_x{dir}=[data_x{dir};ones(size(response{dir,spd})).*spds(spd)];
%         data_y{dir}=[data_y{dir};response{dir,spd}];
%         for i=1:size(response{dir,spd},2)
%             tmp(:,i) =response{dir,spd}(randperm(size(response{dir,spd},1)),i);
%         end
%         data_y_shuffle=[data_y_shuffle;tmp];
%     end
%     end
% end
% h1=figure;
% h2=figure;
% nBins_x=numspds;
% %
% nBins_y=30;
% 
% for ind=1:numdirs
%     xdata=data_x{ind}';
%     ydata=data_y{ind}';
%     stimval=trialdirs_rot(ind);
%     
%     info_forarup
%     %alt:calc_info_P_joint but so many problems with data_x and
%     %data_y: no 0s allowed in response? max(data) must be less than n_
%     %(number of bins)? wtf
% %     for t=1:size(data_x{dir},2)
% %         n_x=5;
% %         n_y=max(data_y{dir}(:,t));
% %         [I_spd{dir},I_spd_err_std{dir},I_spd_err_frac,I_spdR,I_spdR_err_std,I_spdR_err_frac,Pjoint, PjointR] = calc_info_P_joint(data_x{dir}(:,t),data_y{dir}(:,t),n_x,n_y,fracs,numfracreps);
% %     end
%     I_spd{ind}=Iinf;
%     I_spd_poiss_shuffle{ind}=Iinf_1shuffle;
% 
% 
% end
% 
% %
% 
% figure(h1)
% legend(cellstr(num2str(trialdirs_rot')))
% %
% figure(h2)
% title(['I(v,r)',tag])
% I_spd_comb=[];
% I_spd_comb_std=[];
% 
% for i=1:numspds
%     I_spd_comb=[I_spd_comb,I_spd{i}(:,2)];
%     I_spd_comb_std=[I_spd_comb_std,I_spd{i}(:,3)];
% 
% end
% I_spd_mean(:,1)=mean(I_spd_comb,2);
% I_spd_mean(:,2)=mean(I_spd_comb_std,2);
% 
% plot(I_spd_mean(:,1),'k')
% legend([cellstr(num2str(trialdirs_rot'));{'mean'}])
% try
%     saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_spd.fig'])
% catch
%     mkdir(savedir);
%     saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_spd.fig'])
% end
% 
% 
% close all
% %
% 
% 
% %% dir info
% data_x=cell(1,numspds);
% data_y=cell(1,numspds);
% fracs=[1 0.9 0.8 0.5];
% 
% for dir=1:numdirs
%     for spd=1:numspds
%         data_x{spd}=[data_x{spd};ones(size(response{dir,spd})).*trialdirs_rot(dir)];
%         data_y{spd}=[data_y{spd};response{dir,spd}];
%     end
% end
% h1=figure;
% h2=figure;
% nBins_x=numdirs;
% nBins_y=30;
% for ind=1:numspds
%     xdata=data_x{ind}';
%     ydata=data_y{ind}';
%     stimval=spds(ind);
%     info_forarup
%     %alt:calc_info_P_joint but so many problems with data_x and
%     %data_y: no 0s allowed in response? max(data) must be less than n_
%     %(number of bins)? wtf
% %     for t=1:size(data_x{dir},2)
% %         n_x=5;
% %         n_y=max(data_y{dir}(:,t));
% %         [I_spd{dir},I_spd_err_std{dir},I_spd_err_frac,I_spdR,I_spdR_err_std,I_spdR_err_frac,Pjoint, PjointR] = calc_info_P_joint(data_x{dir}(:,t),data_y{dir}(:,t),n_x,n_y,fracs,numfracreps);
% %     end
%     I_dir{ind}=Iinf;
%     I_dir_poiss_shuffle{ind}=Iinf_1shuffle;
% end
% I_dir_comb=[];
% I_dir_comb_std=[];
% 
% for i=1:numspds
%     I_dir_comb=[I_dir_comb,I_dir{i}(:,2)];
%     I_dir_comb_std=[I_dir_comb_std,I_dir{i}(:,3)];
% 
% end
% I_dir_mean(:,1)=mean(I_dir_comb,2);
% I_dir_mean(:,2)=mean(I_dir_comb_std,2);
%     
% 
% figure(h1)
% legend(cellstr(num2str(spds')))
% title(['I(theta,r) at data fracs'])
% 
% figure(h2)
% plot(I_dir_mean(:,1),'k')
% legend([cellstr(num2str(spds'));{'mean'}])
% title(['I(theta,r)',tag])
% saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_dir.fig'])
% close all
% %
% save([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'.mat'],'I_dir','I_spd','I_dir_mean','I_spd_mean','response','I_dir_poiss_shuffle','I_spd_poiss_shuffle')

%% info for each variable combined across other variable

% spd info
data_x=[];
data_y=[];
data_y_shuffle=[];
fracs=[1 0.9 0.8 0.5];
nBins_x=numspds;
nBins_y=30;

for dir=1:numdirs
    for spd=1:numspds
        clear tmp
        data_x=[data_x;ones(size(response{dir,spd})).*spds(spd)];
        data_y=[data_y;response{dir,spd}];
        for i=1:size(response{dir,spd},2)
            tmp(:,i) =response{dir,spd}(randperm(size(response{dir,spd},1)),i);
        end
        data_y_shuffle=[data_y_shuffle;tmp];
    end
end

ind=1;

% 
h1=figure;
h2=figure;
colors=distinguishable_colors(numdirs);
xdata=data_x';
ydata=data_y';
ydata_1shuffle=data_y_shuffle';
info_forarup
I_spd_xdir=Iinf;
I_spd_xdir_poiss_shuffle=Iinf_1shuffle;

figure(h2)
title(['I(v,r), all directions',tag])
hold all;
% load([experiment,'_mutinfo_Unit',num2str(neuron_idx),'.mat'],'I_spd_mean')
% plot(I_spd_mean,'k')
legend('all directions','average')
try
    saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_spd(alldirs).fig'])
catch
    mkdir(savedir);
    saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_spd(alldirs).fig'])
end
close all

%
%dir info
data_x=[];
data_y=[];
data_y_shuffle=[];

fracs=[1 0.9 0.8 0.5];
nBins_x=numdirs;
nBins_y=30;

for dir=1:numdirs
    for spd=1:numspds
        clear tmp
        data_x=[data_x;ones(size(response{dir,spd})).*trialdirs_rot(dir)];
        data_y=[data_y;response{dir,spd}];
        for i=1:size(response{dir,spd},2)
            tmp(:,i) =response{dir,spd}(randperm(size(response{dir,spd},1)),i);
        end
        data_y_shuffle=[data_y_shuffle;tmp];
    end
end
h1=figure;
h2=figure;
colors=distinguishable_colors(numspds);
ind=1;
xdata=data_x';
ydata=data_y';
ydata_1shuffle=data_y_shuffle';
info_forarup
I_dir_xspd=Iinf;
I_dir_xspd_poiss_shuffle=Iinf_1shuffle;
%
save([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'_combined.mat'],'I_dir_xspd','I_spd_xdir','I_spd_xdir_poiss_shuffle','I_dir_xspd_poiss_shuffle')

figure(h2)
title(['I(theta,r), all speeds',tag])
hold all;
% load([experiment,'_mutinfo_Unit',num2str(neuron_idx),'.mat'],'I_dir_mean')
% plot(I_dir_mean,'k')
legend('all speeds','average')


saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_dir(allspds).fig'])

% info for 1-dimensional joint dir-spd distribution

data_x=[];
data_y=[];
data_y_shuffle=[];

fracs=[1 0.9 0.8 0.5];

for dir=1:numdirs
    for spd=1:numspds
        clear tmp
        triind=sub2ind([numdirs,numspds],dir,spd);       
        data_x=[data_x;ones(size(response{dir,spd})).*triind];
        data_y=[data_y;response{dir,spd}];
        for i=1:size(response{dir,spd},2)
            tmp(:,i) =response{dir,spd}(randperm(size(response{dir,spd},1)),i);
        end
        data_y_shuffle=[data_y_shuffle;tmp];
        
    end
end

% h1=figure;
% h2=figure;
% ind=1;
% colors=distinguishable_colors(numdirs);
xdata=data_x';
ydata=data_y';

nBins_x=numdirs*numspds;
nBins_y=30;

%
h1=figure;
h2=figure;
ydata_1shuffle=data_y_shuffle';
info_forarup
I_dirspd_joint_1d=Iinf;
I_dirspd_joint_1d_shuffle=Iinf_shuffle;
I_dirspd_joint_1d_poiss_shuffle=Iinf_1shuffle;



errorbar([1:size(I_dirspd_joint_1d_shuffle,1)]+tShift,I_dirspd_joint_1d_shuffle(:,2),I_dirspd_joint_1d_shuffle(:,3),'Color',[0.5 0.5 0.5]);
title(['Mutual info of 1-dimensional stimulus and response',tag])
saveas(h2,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_dirspd_joint1d.fig'])
save([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'_dirspd_joint1d.mat'],'I_dirspd_joint_1d','I_dirspd_joint_1d_shuffle','I_dirspd_joint_1d_poiss_shuffle')
%% info for joint dir-spd distribution

data_x=[];
data_y=[];
data_z=[];
data_z_shuffle=[];

fracs=[1 0.9 0.8 0.5];

for dir=1:numdirs
    for spd=1:numspds
        clear tmp
%         triind=sub2ind([numdirs,numspds],dir,spd); 
        data_x=[data_x;ones(size(response{dir,spd})).*trialdirs_rot(dir)];
        data_y=[data_y;ones(size(response{dir,spd})).*spds(spd)];
        data_z=[data_z;response{dir,spd}];
        for i=1:size(response{dir,spd},2)
            tmp(:,i) =response{dir,spd}(randperm(size(response{dir,spd},1)),i);
        end
        data_z_shuffle=[data_z_shuffle;tmp];
        
    end
end
h1=figure;
h2=figure;
ind=1;
% colors=distinguishable_colors(numdirs);
xdata=data_x';
ydata=data_y';
zdata=data_z';

nBins_x=12;
nBins_y=5;
nBins_z=30;
zdata_1shuffle=data_z_shuffle';
%
info_forarup2d

%
I_dirspd_joint=Iinf;
I_dirspd_joint_shuffle=Iinf_shuffle;
I_dirspd_joint_poiss_shuffle=Iinf_1shuffle;

%

save([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'_dirspd_joint.mat'],'I_dirspd_joint','I_dirspd_joint_shuffle','I_dirspd_joint_poiss_shuffle')

h=figure;
%% load([experiment,'_mutinfo_Unit',num2str(neuron_idx),'.mat'])
errorbar([1:size(I_dirspd_joint,1)]+tShift,I_dirspd_joint(:,2),I_dirspd_joint(:,3),'k');
hold on
errorbar([1:size(I_dirspd_joint_shuffle,1)]+tShift,I_dirspd_joint_shuffle(:,2),I_dirspd_joint_shuffle(:,3),'Color',[0.5 0.5 0.5]);

load([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'_combined.mat'],'I_dir_xspd','I_spd_xdir')
errorbar([1:size(I_dir_xspd,1)]+tShift,I_dir_xspd(:,2),I_dir_xspd(:,3),'r');
errorbar([1:size(I_spd_xdir,1)]+tShift,I_spd_xdir(:,2),I_spd_xdir(:,3),'b');
errorbar([1:size(I_dir_xspd,1)]+tShift,I_dir_xspd(:,2)+I_spd_xdir(:,2),mean([I_dir_xspd(:,3),I_spd_xdir(:,3)],2),'Color',[1, 0.1034,0.7241]); %purple
[legh,objh,~,~]=legend('I((\theta,v),r)','I(\theta,v),r) shuffled','I(\theta,r), all v','I(v,r), all \theta','\Sigma I(\theta,r) I(v,r)');
% set(objh,'LineWidth',2)
legh=findobj(gcf,'Type','axes','Tag','legend');
% set(legh,'Location','SouthOutside','FontSize',18)
% lc=get(legh,'Children');
% for i=[1 3 5 7 9]
%     ts=get(get(lc(i),'Children'),'Children');
%     set(ts(2),'LineWidth',3)
% end
% set(gca,'FontSize',18)
% xlim([-50 450])
% 
if kind==1
    ylabel('Information from cumulative spike count (bits)')
elseif kind==2
    ylabel('Information from binned spike count (bits)')
elseif kind==3
    ylabel('Information from ISI (bits)')
end
% ylabel('Information from binned spike count (bits)','FontSize',12)
% yh=findobj(gcf,'Type','axes','Tag','ylabel');
% set(yh,'FontSize',18)
% tag='';
% title([' Unit ', num2str(neuron_idx),tag],'Interpreter','none','FontSize',14)
saveas(h,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_dirspd_joint together.fig'])
%
% h=figure;
% 
% errorbar([1:size(I_dirspd_joint,1)]+tShift,I_dirspd_joint(:,2),I_dirspd_joint(:,3),'k');
% hold on
% errorbar([1:size(I_dirspd_joint_shuffle,1)]+tShift,I_dirspd_joint_shuffle(:,2),I_dirspd_joint_shuffle(:,3),'Color',[0.5 0.5 0.5]);
% 
% load([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'.mat'],'I_dir','I_spd','I_dir_mean','I_spd_mean')
% 
% errorbar([1:size(I_dir_mean,1)]+tShift,I_dir_mean(:,1),I_dir_mean(:,2),'y');
% errorbar([1:size(I_spd_mean,1)]+tShift,I_spd_mean(:,1),I_spd_mean(:,2),'g');
% errorbar([1:size(I_dir_mean,1)]+tShift,I_dir_mean(:,1)+I_spd_mean(:,1),mean([I_dir_mean(:,2),I_spd_mean(:,2)],2),'Color',[0,.7,0.7]); %teal
% 
% [legh,objh,OUTH,OUTM]=legend('I({\theta,v},r)','I({\theta,v},r) shuffled','Mean of v-separated I(\theta,r)','Mean of \theta-separated I(v,r)','Sum of mean separated I(\theta,r) and I(v,r)');
% set(legh,'Interpreter','none','Location','SouthOutside')
% % set(OBJH,'LineWidth',2)
% title(['Information, ' experiment, ' Unit ', num2str(neuron_idx),tag],'Interpreter','none')
% 
% saveas(h,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_dirspd_joint separated.fig'])
%%
h=figure;
% load([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'.mat'],'I_dir_poiss_shuffle','I_spd_poiss_shuffle')
load([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'_combined.mat'],'I_spd_xdir_poiss_shuffle','I_dir_xspd_poiss_shuffle')
load([savedir,experiment,'_mutinfo_Unit',num2str(neuron_idx),'_dirspd_joint.mat'],'I_dirspd_joint_shuffle','I_dirspd_joint_poiss_shuffle')

errorbar([1:size(I_dirspd_joint_poiss_shuffle,1)]+tShift,I_dirspd_joint_poiss_shuffle(:,2),I_dirspd_joint_poiss_shuffle(:,3),'k');
hold on
errorbar([1:size(I_dirspd_joint_shuffle,1)]+tShift,I_dirspd_joint_shuffle(:,2),I_dirspd_joint_shuffle(:,3),'Color',[0.5 0.5 0.5]);

errorbar([1:size(I_dir_xspd_poiss_shuffle,1)]+tShift,I_dir_xspd_poiss_shuffle(:,2),I_dir_xspd_poiss_shuffle(:,3),'r');
errorbar([1:size(I_spd_xdir_poiss_shuffle,1)]+tShift,I_spd_xdir_poiss_shuffle(:,2),I_spd_xdir_poiss_shuffle(:,3),'b');
errorbar([1:size(I_dir_xspd_poiss_shuffle,1)]+tShift,I_dir_xspd_poiss_shuffle(:,2)+I_spd_xdir_poiss_shuffle(:,2),mean([I_dir_xspd_poiss_shuffle(:,3),I_spd_xdir_poiss_shuffle(:,3)],2),'Color',[1, 0.1034,0.7241]); %purple
[legh,objh,~,~]=legend('I((\theta,v),r)','I(\theta,v),r) shuffled','I(\theta,r), all v','I(v,r), all \theta','\Sigma I(\theta,r) I(v,r)');
% set(objh,'LineWidth',2)
legh=findobj(gcf,'Type','axes','Tag','legend');
% set(legh,'Location','SouthOutside','FontSize',18)
% lc=get(legh,'Children');
% for i=[1 3 5 7 9]
%     ts=get(get(lc(i),'Children'),'Children');
%     set(ts(2),'LineWidth',3)
% end
% set(gca,'FontSize',18)
% xlim([-50 450])
% 
if kind==1
    ylabel('Information from cumulative spike count (bits)')
elseif kind==2
    ylabel('Information from binned spike count (bits)')
elseif kind==3
    ylabel('Information from ISI (bits)')
end
% ylabel('Information from binned spike count (bits)','FontSize',12)
% yh=findobj(gcf,'Type','axes','Tag','ylabel');
% set(yh,'FontSize',18)
% tag='';
% title([' Unit ', num2str(neuron_idx),tag],'Interpreter','none','FontSize',14)
saveas(h,[savedir,experiment,'_Unit ',num2str(neuron_idx),'_I_dirspd_joint_poiss together.fig'])
