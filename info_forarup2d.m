
% load xdata, ydata (stims), (zdata=response) the arrays for which you wish to compute the mutual information
% arrays are assumed to by time x trials
% appropriately windowed so that xdata and ydata have length nWins
% create shuffled versions
nWins = size(xdata,1);
xdata_shuffle = xdata(randperm(size(xdata,1)),randperm(size(xdata,2)));
ydata_shuffle = ydata(randperm(size(ydata,1)),randperm(size(ydata,2)));
zdata_shuffle = zdata(randperm(size(zdata,1)),randperm(size(zdata,2)));
%add zs

% Bin data
% test different values, the number you can use depends on the data sample
nBins_x=10;
nBins_y=10;
nBins_z=10;

[x_binned,x_cuts,x_occ] = adaptbin(reshape(xdata,1,numel(xdata)),nBins_x);
x_binned = reshape(x_binned,nWins,size(x_binned,2));
[x_binned_shuffle,x_cuts_shuffle,x_occ_shuffle] = adaptbin(reshape(xdata_shuffle,1,numel(xdata)),nBins_x);
x_binned_shuffle = reshape(x_binned_shuffle,nWins,size(x_binned,2));

[y_binned,y_cuts,y_occ] = adaptbin(reshape(ydata,1,numel(ydata)),nBins_y);
y_binned = reshape(y_binned,nWins,size(y_binned,2));
[y_binned_shuffle,y_cuts_shuffle,y_occ_shuffle] = adaptbin(reshape(ydata_shuffle,1,numel(ydata)),nBins_y);
y_binned_shuffle = reshape(y_binned_shuffle,nWins,size(y_binned,2));

[z_binned,z_cuts,z_occ] = adaptbin(reshape(zdata,1,numel(zdata)),nBins_z);
z_binned = reshape(z_binned,nWins,size(z_binned,2));
[z_binned_shuffle,z_cuts_shuffle,z_occ_shuffle] = adaptbin(reshape(zdata_shuffle,1,numel(zdata)),nBins_z);
z_binned_shuffle = reshape(z_binned_shuffle,nWins,size(z_binned,2));

    
datafrac=[1 .95 .90 0.85 .80 .50];
numfracreps=30;

ntrials = size(xdata,2);
Sz = zeros(length(datafrac),numfracreps,nWins);%
Sz_given_x = Sy;
Ixy = Sy;
Sy_shuffle = Sy;
Sy_given_x_shuffle = Sy;
Ixy_shuffle = Sy;
%Ixyz
%Ixyz_shuffle

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
        %Zdata
        %Zdata_r
	Ydata_r = y_binned_shuffle(:,trialvec);
%keep
        Pjoint = zeros(nWins,nBins_x,nBins_x);
        Pjoint_shuffle = Pjoint;
        for m=1:nWins
            for n=1:size(Xdata,2); % number of examples       
                x = Xdata(m,n);
	        y = Ydata(m,n);
	        Pjoint(m,x,y) = Pjoint(m,x,y) + 1/size(Xdata,2);
	        x = Xdata_r(m,n);
		y = Ydata_r(m,n);
	        Pjoint_shuffle(m,x,y) = Pjoint_shuffle(m,x,y) + 1/size(Xdata_r,2);
            end
	    Ptemp = reshape(Pjoint(m,:,:),nBins_x,nBins_x);
	    Ptempy = sum(Ptemp,1);
	    PtempyGx = Ptemp./(sum(Ptemp,2)*ones(1,nBins_x) +eps);
            Px = sum(Ptemp,2);
            Px = Px./sum(Px);
	    Sy(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
	    Sy_given_x(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
	    Ixy(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
          

            Ptemp = reshape(Pjoint_shuffle(m,:,:),nBins_x,nBins_x);
	    Ptempy = sum(Ptemp,1);
            PtempyGx = Ptemp./(sum(Ptemp,2)*ones(1,nBins_x) +eps);
            Px = sum(Ptemp,2);
            Px = Px./sum(Px);
            Sy_shuffle(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
            Sy_given_x_shuffle(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
	    Ixy_shuffle(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
%keep to here
%2dstim 
        Pjoint = zeros(nWins,nBinsx,nBinsy,nBinsz);
        Pjoint_shuffle = Pjoint;
        for m=1:nWins
            for n=1:size(Xdata,2); % number of examples       
                x = Xdata(m,n);
	        y = Ydata(m,n);
            Z=Zdata(m,n);
            
	        Pjoint(m,x,y,z) = Pjoint(m,x,y,z) + 1/size(Xdata,2);%p of response and stim across time
	        x = Xdata_r(m,n);
		y = Ydata_r(m,n);
        z=Zdata_r(m,n);
	        Pjoint_shuffle(m,x,y,z) = Pjoint_shuffle(m,x,y,z) + 1/size(Xdata_r,2);
            end
	    Ptemp = reshape(Pjoint(m,:,:,:),nBinsx,nBinsy,nBinsz);%p of response and stim
	    Ptempz = sum(sum(Ptemp,1),1);%p of response
	    PtempyGx = Ptemp./(sum(Ptemp,2)*ones(1,nBins_x) +eps);%normalize
            Px = sum(Ptemp,2);
            Px = Px./sum(Px);
	    Sy(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
	    Sy_given_x(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
	    Ixy(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
          %kind of ok to here

            Ptemp = reshape(Pjoint_shuffle(m,:,:),nBins_x,nBins_x);
	    Ptempy = sum(Ptemp,1);
            PtempyGx = Ptemp./(sum(Ptemp,2)*ones(1,nBins_x) +eps);
            Px = sum(Ptemp,2);
            Px = Px./sum(Px);
            Sy_shuffle(datafracind,nrep,m) = -sum(Ptempy.*log2(Ptempy+eps));
            Sy_given_x_shuffle(datafracind,nrep,m) = Px'*-sum(PtempyGx.*log2(PtempyGx+eps),2);
	    Ixy_shuffle(datafracind,nrep,m) = sum(sum(Ptemp.*log2(Ptemp./(sum(Ptemp,2)*sum(Ptemp,1)+eps) +eps)));
        end; %nWins

     end %nrep	

     if datafracind==1 & numfracreps>1 % need to replicate answers for data fraction of 1
        Sy(datafracind,2:numfracreps,1:nWins) =repmat(Sy(datafracind,1,1:nWins),1,numfracreps-1);
        Sy_given_x(datafracind,2:numfracreps,1:nWins)=repmat(Sy_given_x(datafracind,1,1:nWins),1,numfracreps-1);
        Ixy(datafracind,2:numfracreps,1:nWins)=repmat(Ixy(datafracind,1,1:nWins),1,numfracreps-1);
        
	Sy_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sy_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
	Sy_given_x_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Sy_given_x_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
	Ixy_shuffle(datafracind,2:numfracreps,1:nWins)=repmat(Ixy_shuffle(datafracind,1,1:nWins),1,numfracreps-1);
     end
end %datafrac

nstop=length(datafrac)-1;
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
    
end

% check finite size correction
figure
plot(1./datafrac,mean(Ixy(:,:,15),2),'x');
hold on
idata = polyval(Iinf(15,1:2),[0 1./datafrac]);
plot([0 1./datafrac],idata,'r');

figure
h=errorbar([1:length(Iinf)]+tShift,Iinf(:,2),Iinf(:,3),'b');
hold on
h=plot([1:length(Iinf)]+tShift,Iinf(:,2),'g');
xlabel('Time (ms)');
ylabel('bits');

