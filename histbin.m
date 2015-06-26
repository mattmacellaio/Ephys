
function [binned_data,bin_cents,bin_cuts,bin_occupancy]= histbin(data,Nbins)
%based on adaptbin

data_size=size(data);
data=data(:);
[bin_occupancy,bin_cents]=hist(data,Nbins);
bin_cuts = ceil([0:1/Nbins:1]*max(data));
binned_data = ones(length(data),1);
% everything starts in bin 1

for k=2:Nbins
    test = (data > bin_cuts(k-1)) & (data <= bin_cuts(k)); %binary
    binned_data(test==1) = k;
    NN(k) = sum(test==1);
end
NN(1) = sum(binned_data==1);
bin_occupancy=NN;

for k=1:Nbins
    index=binned_data==k;
    bin_cents(k)=mean(data(index));
end

binned_data=reshape(binned_data,data_size);

