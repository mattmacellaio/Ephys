
% load xdata, ydata (stims), (zdata=response) the arrays for which you wish to compute the mutual information
% arrays are assumed to by time x trials
% appropriately windowed so that xdata and ydata have length nWins
% create shuffled versions
nWins = size(xdata,1);
xdata_shuffle = xdata(randperm(size(xdata,1)),randperm(size(xdata,2)));
ydata_shuffle = ydata(randperm(size(ydata,1)),randperm(size(ydata,2)));
zdata_shuffle = zdata(randperm(size(zdata,1)),randperm(size(zdata,2)));
% for i=1:size(xdata,1)
%     xdata_1shuffle(i,:) = xdata(i,randperm(size(xdata,2)));
%     ydata_1shuffle(i,:) = ydata(i,randperm(size(ydata,2)));
%     zdata_1shuffle(i,:) = zdata(i,randperm(size(zdata,2)));
% end

% Bin data
% test different values, the number you can use depends on the data sample

[x_binned,x_cuts,x_occ] = adaptbin(reshape(xdata,1,numel(xdata)),nBins_x);
x_binned = reshape(x_binned,nWins,size(xdata,2));
[x_binned_shuffle,x_cuts_shuffle,x_occ_shuffle] = adaptbin(reshape(xdata_shuffle,1,[]),nBins_x);
x_binned_shuffle = reshape(x_binned_shuffle,nWins,size(xdata,2));

[y_binned,y_cuts,y_occ] = adaptbin(reshape(ydata,1,numel(ydata)),nBins_y);
y_binned = reshape(y_binned,nWins,size(ydata,2));
[y_binned_shuffle,y_cuts_shuffle,y_occ_shuffle] = adaptbin(reshape(ydata_shuffle,1,[]),nBins_y);
y_binned_shuffle = reshape(y_binned_shuffle,nWins,size(ydata,2));

[z_binned,z_cuts,z_occ] = adaptbin(reshape(zdata,1,numel(zdata)),nBins_z);
z_binned = reshape(z_binned,nWins,size(zdata,2));
[z_binned_shuffle,z_cuts_shuffle,z_occ_shuffle] = adaptbin(reshape(zdata_shuffle,1,[]),nBins_z);
z_binned_shuffle = reshape(z_binned_shuffle,nWins,size(zdata,2));
[z_binned_1shuffle,z_cuts_1shuffle,z_occ_1shuffle] = adaptbin(reshape(zdata_1shuffle,1,[]),nBins_z);
z_binned_1shuffle = reshape(z_binned_1shuffle,nWins,size(zdata,2));

    
datafrac=[1 .95 .90 0.85 .80 .50];

ntrials = size(xdata,2);
Sz = zeros(length(datafrac),numfracreps,nWins);%
Sz_given_xy = Sz;
Sz_given_y=Sz;
Ixyz = Sz;
Sz_shuffle = Sz;
Sz_given_xy_shuffle = Sz;
Ixyz_shuffle = Sz;
Sz_1shuffle = Sz;
Sz_given_xy_1shuffle = Sz;
Ixyz_1shuffle = Sz;

for datafracind=1:length(datafrac)
    fprintf('data fraction index:%d\n',datafracind);
 %                 tic;
    if datafracind==1
       nrepstop=1;
    else
       nrepstop=numfracreps;
    end

    Ntrialsfrac=floor(datafrac(datafracind).*ntrials);


    for nrep=1:nrepstop
        fprintf('rep %d of %d\n',nrep,nrepstop);
        trialvec=randperm(ntrials);
        trialvec=trialvec(1:Ntrialsfrac);
        Xdata = x_binned(:,trialvec);
        Xdata_r = x_binned_shuffle(:,trialvec);
        Ydata = y_binned(:,trialvec);
        Ydata_r = y_binned_shuffle(:,trialvec);
        Zdata = z_binned(:,trialvec);
        Zdata_r = z_binned_shuffle(:,trialvec);
        Zdata_1r = z_binned_1shuffle(:,trialvec);
        
%keep but probably not working anymore with changes
%         Pjoint = zeros(nWins,nBins_x,nBins_x);
%         Pjoint_shuffle = Pjoint;
%         for m=1:nWins
%             for n=1:size(Xdata,2); % number of examples       
%                 x = Xdata(m,n);
%                 y = Ydata(m,n);
%                 Pjoint(m,x,y) = Pjoint(m,x,y) + 1/size(Xdata,2);
%                 x = Xdata_r(m,n);
%                 y = Ydata_r(m,n);
%                 Pjoint_shuffle(m,x,y) = Pjoint_shuffle(m,x,y) + 1/size(Xdata_r,2);
%             end
%             Ptemp = reshape(Pjoint(m,:,:),nBins_x,nBins_x);
%             Ptempy = sum(Ptemp,1);
%             PtempzGx = Ptemp./(sum(Ptemp,2)*ones(1,nBins_x) +eps);
%             Px = sum(Ptemp,2);
%             Px = Px./sum(Px);
%             Sy(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
%             Sy_given_x(datafracind,nrep,m) = Px'*-sum(PtempzGx.*log2(PtempzGx+eps),2);
%             Ixy(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
%           
% 
%             Ptemp = reshape(Pjoint_shuffle(m,:,:),nBins_x,nBins_x);
%             Ptempy = sum(Ptemp,1);
%             PtempzGx = Ptemp./(sum(Ptemp,2)*ones(1,nBins_x) +eps);
%             Px = sum(Ptemp,2);
%             Px = Px./sum(Px);
%             Sy_shuffle(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
%             Sy_given_x_shuffle(datafracind,nrep,m) = Px'*-sum(PtempzGx.*log2(PtempzGx+eps),2);
%             Ixy_shuffle(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
%         end
%keep to here
%2dstim 
        Pjoint = zeros(nWins,nBins_x,nBins_y,nBins_z);
        Pjoint_shuffle = Pjoint;
        Pjoint_1shuffle = Pjoint;

        for m=1:nWins
            for n=1:size(Xdata,2); % number of examples       
                x = Xdata(m,n);
                y = Ydata(m,n);
                z = Zdata(m,n);
                Pjoint(m,x,y,z) = Pjoint(m,x,y,z) + 1/size(Xdata,2);%p of response and stim across time
                x = Xdata_r(m,n);
                y = Ydata_r(m,n);
                z = Zdata_r(m,n);
                Pjoint_shuffle(m,x,y,z) = Pjoint_shuffle(m,x,y,z) + 1/size(Xdata_r,2);
                x = Xdata(m,n);
                y = Ydata(m,n);
                z = Zdata_1r(m,n);
                Pjoint_1shuffle(m,x,y,z) = Pjoint_1shuffle(m,x,y,z) + 1/size(Xdata,2);
                
            end
            %data
            Ptemp2 = squeeze(Pjoint(m,:,:,:)); %p of response and stim
            Ptempz = reshape(sum(sum(Ptemp2,1),2),1,[]);%p of response, %sum=1
            PtempzGxy = Ptemp2./(repmat(sum(Ptemp2,3),[1,1,nBins_z]) +eps);%normalize, sum=1*nBins
            Pxy = sum(Ptemp2,3);
            Pxy = Pxy./sum(sum(Pxy)); %sum=1
            Sz(datafracind,nrep,m) = -sum(Ptempz.*log2(Ptempz+eps));
            Sz_given_xy(datafracind,nrep,m) = sum(sum(Pxy.*sum(PtempzGxy.*log2(PtempzGxy+eps),3)));
            temp=[];
            for zind=1:nBins_z
                temp(:,:,zind)=Pxy.*Ptempz(zind);
            end
            Ixyz(datafracind,nrep,m) = sum(sum(sum(Ptemp2.*log2(Ptemp2./(temp+eps) +eps))));
            
            %all shuffled
            Ptemp2 = squeeze(Pjoint_shuffle(m,:,:,:)); %p of response and stim
            Ptempz = reshape(sum(sum(Ptemp2,1),2),1,[]);%p of response, %sum=1
            PtempzGxy = Ptemp2./(repmat(sum(Ptemp2,3),[1,1,nBins_z]) +eps);%normalize, sum=1*nBins
            Pxy = sum(Ptemp2,3);
            Pxy = Pxy./sum(sum(Pxy)); %sum=1
            Sz_shuffle(datafracind,nrep,m) = -sum(Ptempz.*log2(Ptempz+eps));
            Sz_given_xy_shuffle(datafracind,nrep,m) = sum(sum(Pxy.*sum(PtempzGxy.*log2(PtempzGxy+eps),3)));
            temp=[];
            for zind=1:nBins_z
                temp(:,:,zind)=Pxy.*Ptempz(zind);
            end
            Ixyz_shuffle(datafracind,nrep,m) = sum(sum(sum(Ptemp2.*log2(Ptemp2./(temp+eps) +eps))));
            
            %just trial-shuffled
            Ptemp2 = squeeze(Pjoint_1shuffle(m,:,:,:)); %p of response and stim
            Ptempz = reshape(sum(sum(Ptemp2,1),2),1,[]);%p of response, %sum=1
            PtempzGxy = Ptemp2./(repmat(sum(Ptemp2,3),[1,1,nBins_z]) +eps);%normalize, sum=1*nBins
            Pxy = sum(Ptemp2,3);
            Pxy = Pxy./sum(sum(Pxy)); %sum=1
            Sz_1shuffle(datafracind,nrep,m) = -sum(Ptempz.*log2(Ptempz+eps));
            Sz_given_xy_1shuffle(datafracind,nrep,m) = sum(sum(Pxy.*sum(PtempzGxy.*log2(PtempzGxy+eps),3)));
            temp=[];
            for zind=1:nBins_z
                temp(:,:,zind)=Pxy.*Ptempz(zind);
            end
            Ixyz_1shuffle(datafracind,nrep,m) = sum(sum(sum(Ptemp2.*log2(Ptemp2./(temp+eps) +eps))));
%
%            
        end; %nWins

     end %nrep	

     if datafracind==1 & numfracreps>1 % need to replicate answers for data fraction of 1
        Sz(datafracind,2:numfracreps,1:nWins) =repmat(Sz(datafracind,1,1:nWins),1,numfracreps-1);
        Sz_given_xy(datafracind,2:numfracreps,1:nWins)=repmat(Sz_given_xy(datafracind,1,1:nWins),1,numfracreps-1);
        Ixyz(datafracind,2:numfracreps,1:nWins)=repmat(Ixyz(datafracind,1,1:nWins),1,numfracreps-1);
        
        Sz_shuffle(datafracind,2:numfracreps,1:nWins) =repmat(Sz_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Sz_given_xy_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sz_given_xy_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Ixyz_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Ixyz_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        
        Sz_1shuffle(datafracind,2:numfracreps,1:nWins) =repmat(Sz_1shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Sz_given_xy_1shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sz_given_xy_1shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Ixyz_1shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Ixyz_1shuffle(datafracind,1,1:nWins),1,numfracreps-1);

     end
end %datafrac

nstop=length(datafrac)-1;
%Extrapolate to inf data size
% Use the intercept (2nd value)
for tt = 1:nWins

    Sinf(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sz(1:nstop,:,tt),2)',1);  
    Sinf(tt,3) = std(Sz(length(datafrac),:,tt),0,2);
    Snoise_inf(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sz_given_xy(1:nstop,:,tt),2)',1);
    Snoise_inf(tt,3) = std(Sz_given_xy(length(datafrac),:,tt),0,2);
    Iinf(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixyz(1:nstop,:,tt),2)',1);    
    Iinf(tt,3) = std(Ixyz(length(datafrac),:,tt),0,2);
    
    Sinf_shuffle(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sz_shuffle(1:nstop,:,tt),2)',1);  
    Sinf_shuffle(tt,3) = std(Sz_shuffle(length(datafrac),:,tt),0,2);
    Snoise_inf_shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sz_given_xy_shuffle(1:nstop,:,tt),2)',1);
    Snoise_inf_shuffle(tt,3) = std(Sz_given_xy_shuffle(length(datafrac),:,tt),0,2);
    Iinf_shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixyz_shuffle(1:nstop,:,tt),2)',1);    
    Iinf_shuffle(tt,3) = std(Ixyz_shuffle(length(datafrac),:,tt),0,2);
    
    Sinf_1shuffle(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sz_1shuffle(1:nstop,:,tt),2)',1);  
    Sinf_1shuffle(tt,3) = std(Sz_1shuffle(length(datafrac),:,tt),0,2);
    Snoise_inf_1shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sz_given_xy_1shuffle(1:nstop,:,tt),2)',1);
    Snoise_inf_1shuffle(tt,3) = std(Sz_given_xy_1shuffle(length(datafrac),:,tt),0,2);
    Iinf_1shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixyz_1shuffle(1:nstop,:,tt),2)',1);    
    Iinf_1shuffle(tt,3) = std(Ixyz_1shuffle(length(datafrac),:,tt),0,2);


%     Sinf_shuffle(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sy_shuffle(1:nstop,:,tt),2)',1);
%     Sinf_shuffle(tt,3) = std(Sy_shuffle(length(datafrac),:,tt),0,2);
%     Snoise_inf_shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sy_given_x_shuffle(1:nstop,:,tt),2)',1);
%     Snoise_inf_shuffle(tt,3) = std(Sy_shuffle(length(datafrac),:,tt),0,2);
%     Iinf_shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixy_shuffle(1:nstop,:,tt),2)',1);
%     Iinf_shuffle(tt,3) = std(Ixy_shuffle(length(datafrac),:,tt),0,2);
    
end
tShift=0;
% check finite size correction
figure(h1)
if nWins>100
    samptime=100;
elseif nWins==1
    samptime=1;
else
    samptime=40;
end

plot(1./datafrac,mean(Ixyz(:,:,samptime),2),'x');
hold on
idata = polyval(Iinf(samptime,1:2),[0 1./datafrac]);
plot([0 1./datafrac],idata,'r');
if kind==1
    saveas(gcf,[experiment,'inffit_joint2d.fig'])
end
figure(h2)
if nWins==1
    errorbar([1:size(Iinf,1)]+tShift,Iinf(:,2),Iinf(:,3),'k');
    hold on
    errorbar([1:size(Iinf_shuffle,1)]+tShift,Iinf_shuffle(:,2),Iinf_shuffle(:,3),'Color',[0.5 0.5 0.5]);
else
    errorbar([1:size(Iinf,1)]+tShift,Iinf(:,2),Iinf(:,3),'k');
    hold on
    errorbar([1:size(Iinf_shuffle,1)]+tShift,Iinf_shuffle(:,2),Iinf_shuffle(:,3),'Color',[0.5 0.5 0.5]);
end

% h=plot([1:length(Iinf)]+tShift,Iinf(:,2),'g');
xlabel('Time (ms)');
ylabel('bits');

