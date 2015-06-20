% load xdata, ydata the arrays for which you wish to compute the mutual information
% arrays are assumed to be time x trials
% appropriately windowed so that xdata and ydata have length nWins
% create shuffled versions
nWins = size(xdata,1);
xdata_shuffle = xdata(randperm(size(xdata,1)),randperm(size(xdata,2)));
ydata_shuffle = ydata(randperm(size(ydata,1)),randperm(size(ydata,2)));

% Bin data
% actually this is x_binned, bin_centers, bin_edges (missing bin_occupancy)
[x_binned,x_cuts,x_occ] = adaptbin(reshape(xdata,1,numel(xdata)),nBins_x);
x_binned = reshape(x_binned,nWins,size(xdata,2));
[y_binned,y_cuts,y_occ] = adaptbin(reshape(ydata,1,numel(ydata)),nBins_y);
y_binned = reshape(y_binned,nWins,size(ydata,2));

[x_binned_shuffle,x_cuts_shuffle,x_occ_shuffle] = adaptbin(reshape(xdata_shuffle,1,numel(xdata)),nBins_x);
x_binned_shuffle = reshape(x_binned_shuffle,nWins,size(xdata,2));
[y_binned_shuffle,y_cuts_shuffle,y_occ_shuffle] = adaptbin(reshape(ydata_shuffle,1,numel(ydata)),nBins_y);
y_binned_shuffle = reshape(y_binned_shuffle,nWins,size(ydata,2));

[y_binned_1shuffle,y_cuts_1shuffle,y_occ_1shuffle] = adaptbin(reshape(ydata_1shuffle,1,numel(ydata)),nBins_y);
y_binned_1shuffle = reshape(y_binned_1shuffle,nWins,size(ydata,2));

    
datafrac=[1 .95 .90 0.85 .80 .70 .60 .50 .40 .30 .20 .10];

ntrials = size(xdata,2);
Sy = zeros(length(datafrac),numfracreps,nWins);
Sy_given_x = Sy;
Ixy = Sy;
Sy_shuffle = Sy;
Sy_given_x_shuffle = Sy;
Ixy_shuffle = Sy;
Sy_1shuffle = Sy;
Sy_given_x_1shuffle = Sy;
Ixy_1shuffle = Sy;

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
        Ydata = y_binned(:,trialvec);
        Xdata_r = x_binned_shuffle(:,trialvec);
        Ydata_r = y_binned_shuffle(:,trialvec);
        Ydata_1r = y_binned_1shuffle(:,trialvec);

        Pjoint = zeros(nWins,nBins_x,nBins_y);
        Pjoint_shuffle = Pjoint;
        Pjoint_1shuffle = Pjoint;

        for m=1:nWins
            for n=1:size(Xdata,2); % number of examples       
                x = Xdata(m,n);
                y = Ydata(m,n);
                Pjoint(m,x,y) = Pjoint(m,x,y) + 1/size(Xdata,2);
                x = Xdata_r(m,n);
                y = Ydata_r(m,n);
                Pjoint_shuffle(m,x,y) = Pjoint_shuffle(m,x,y) + 1/size(Xdata_r,2);
                x = Xdata(m,n);
                y = Ydata_1r(m,n);
                Pjoint_1shuffle(m,x,y) = Pjoint_1shuffle(m,x,y) + 1/size(Xdata,2);
            end
            Ptemp = reshape(Pjoint(m,:,:),nBins_x,nBins_y); %nbins x nbins %sum=1
            Ptempy = sum(Ptemp,1); %sum=1   1 x nbins
            PtempyGx = Ptemp./(repmat(sum(Ptemp,2),[1,nBins_y]) +eps); %sum=1*nBins
            Px = sum(Ptemp,2);
            Px = Px./sum(Px); %sum=1
            Sy(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
            Sy_given_x(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
            Ixy(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(Px*Ptempy+eps) +eps)));
          

            Ptemp = reshape(Pjoint_shuffle(m,:,:),nBins_x,nBins_y);
            Ptempy = sum(Ptemp,1);
            PtempyGx = Ptemp./(repmat(sum(Ptemp,2),[1,nBins_y]) +eps); %sum=1*nBins
            Px = sum(Ptemp,2);
            Px = Px./sum(Px);
            Sy_shuffle(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
            Sy_given_x_shuffle(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
            Ixy_shuffle(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
            
            Ptemp = reshape(Pjoint_1shuffle(m,:,:),nBins_x,nBins_y);
            Ptempy = sum(Ptemp,1);
            PtempyGx = Ptemp./(repmat(sum(Ptemp,2),[1,nBins_y]) +eps); %sum=1*nBins
            Px = sum(Ptemp,2);
            Px = Px./sum(Px);
            Sy_1shuffle(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
            Sy_given_x_1shuffle(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
            Ixy_1shuffle(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));


        end; %nWins

     end %nrep	

     if datafracind==1 & numfracreps>1 % need to replicate answers for data fraction of 1
        Sy(datafracind,2:numfracreps,1:nWins) =repmat(Sy(datafracind,1,1:nWins),1,numfracreps-1);
        Sy_given_x(datafracind,2:numfracreps,1:nWins)=repmat(Sy_given_x(datafracind,1,1:nWins),1,numfracreps-1);
        Ixy(datafracind,2:numfracreps,1:nWins)=repmat(Ixy(datafracind,1,1:nWins),1,numfracreps-1);
        
        Sy_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sy_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Sy_given_x_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sy_given_x_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Ixy_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Ixy_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
     
        Sy_1shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sy_1shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Sy_given_x_1shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sy_given_x_1shuffle(datafracind,1,1:nWins),1,numfracreps-1);
        Ixy_1shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Ixy_1shuffle(datafracind,1,1:nWins),1,numfracreps-1);

     end
end %datafrac

nstop=8;
%Extrapolate to inf data size
% Use the intercept (2nd value)
for tt = 1:nWins

    Sinf(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sy(1:nstop,:,tt),2)',1);  
    Sinf(tt,3) = std(Sy(length(datafrac),:,tt),0,2);
    Snoise_inf(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sy_given_x(1:nstop,:,tt),2)',1);
    Snoise_inf(tt,3) = std(Sy_given_x(length(datafrac),:,tt),0,2);
    Iinf(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixy(1:nstop,:,tt),2)',1);    
    Iinf(tt,3) = std(Ixy(length(datafrac),:,tt),0,2);
    
    Sinf_shuffle(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sy_shuffle(1:nstop,:,tt),2)',1);
    Sinf_shuffle(tt,3) = std(Sy_shuffle(length(datafrac),:,tt),0,2);
    Snoise_inf_shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sy_given_x_shuffle(1:nstop,:,tt),2)',1);
    Snoise_inf_shuffle(tt,3) = std(Sy_shuffle(length(datafrac),:,tt),0,2);
    Iinf_shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixy_shuffle(1:nstop,:,tt),2)',1);
    Iinf_shuffle(tt,3) = std(Ixy_shuffle(length(datafrac),:,tt),0,2);
    
    Sinf_1shuffle(tt,1:2) =  polyfit(1./datafrac(1:nstop),mean(Sy_1shuffle(1:nstop,:,tt),2)',1);
    Sinf_1shuffle(tt,3) = std(Sy_1shuffle(length(datafrac),:,tt),0,2);
    Snoise_inf_1shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Sy_given_x_1shuffle(1:nstop,:,tt),2)',1);
    Snoise_inf_1shuffle(tt,3) = std(Sy_1shuffle(length(datafrac),:,tt),0,2);
    Iinf_1shuffle(tt,1:2) = polyfit(1./datafrac(1:nstop),mean(Ixy_1shuffle(1:nstop,:,tt),2)',1);
    Iinf_1shuffle(tt,3) = std(Ixy_1shuffle(length(datafrac),:,tt),0,2);

end

% check finite size correction
figure(h1)
if nWins>100
    samptime=100;
elseif nWins==1
    samptime=1;
else
    samptime=40;
end

plot(1./datafrac,mean(Ixy(:,:,samptime),2),'x');
hold on
idata = polyval(Iinf(samptime,1:2),[0 1./datafrac]);
plot([0 1./datafrac],idata,'Color',colors(ind,:));
if kind==2
    saveas(gcf,[experiment,'inffit_',figtag,'.fig'])
end
%
if kind==2
    tShift=20;
else
    tShift=0;
end

figure(h2)
if nWins==1
   h=errorbar([1:size(Iinf,1)],Iinf(:,2),Iinf(:,3),'Color',colors(ind,:));
    hold on
    errorbar([1:size(Iinf,1)],Iinf_shuffle(:,2),Iinf_shuffle(:,3),'Color',colors(ind,:));
else
    h=errorbar([1:size(Iinf,1)]+tShift,Iinf(:,2),Iinf(:,3),'Color',colors(ind,:));
    hold on
    errorbar([1:size(Iinf,1)]+tShift,Iinf_shuffle(:,2),Iinf_shuffle(:,3),'Color',colors(ind,:));
end
% h=plot([1:length(Iinf)]+tShift,Iinf(:,2),'g');
xlabel('Time (ms)');
ylabel('bits');

