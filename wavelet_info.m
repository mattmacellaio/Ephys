%for using WIgui: implement after spikes have already been extracted and
%organized from nexfile/maestro files

%class_id: 1 x trials vector of doubles across all 
%spiketimes: 1 x trialtypes vector of cells with spiketimes for each
%trialtype in ms (1 x repeats)
tri_type=ceil(tri_num/5);

%spd info
class_id=[];
clear spiketimes handles
for i=1:length(tri_num)
    for j=1:size(spk_tr,1); %for each seg
        spdind=trispdind(i);
        triind=tri_type(i);
        if ~isempty(spk_tr{j,i}')
            class_id(end+1)=segdirinds{triind}(j);
            spiketimes{length(class_id)}=spk_tr{j,i}';
        end
    end
end

%  class_id=class_id(1:500);
% spiketimes=spiketimes(1:500);

save([experiment,'WIdata.mat'],'class_id','spiketimes')
% WIgui %raster throws errors if any cells in spiketimes are empty arrays: no spikes that trial
WIscript_caller 
%

figure;imagesc(handles.decode.SPKCNTconfusionmatrix);colorbar
title(['confusion matrix from spike count for speed: performance ',num2str(handles.decode.SPKCNTperf)])
figure;imagesc(handles.decode.WVconfusionmatrix);colorbar
title(['confusion matrix from wavelet info for speed: performance ',num2str(handles.decode.WVperf)])

%dir info
class_id=[];
clear spiketimes handles
for i=1:length(tri_num)
    for j=1:size(spk_tr,1); %for each seg
        spdind=trispdind(i);
        triind=tri_type(i);
        if ~isempty(spk_tr{j,i}')
            class_id(end+1)=spdind;
            spiketimes{length(class_id)}=spk_tr{j,i}';
        end
    end
end

%  class_id=class_id(1:500);
% spiketimes=spiketimes(1:500);

save([experiment,'WIdata.mat'],'class_id','spiketimes')
% WIgui %raster throws errors if any cells in spiketimes are empty arrays: no spikes that trial
WIscript_caller 
%

figure;imagesc(handles.decode.SPKCNTconfusionmatrix);colorbar
title(['confusion matrix from spike count for direction: performance ',num2str(handles.decode.SPKCNTperf)])
figure;imagesc(handles.decode.WVconfusionmatrix);colorbar
title(['confusion matrix from wavelet info for direction: performance ',num2str(handles.decode.WVperf)])

%joint info
class_id=[];
clear spiketimes handles
for i=1:length(tri_num)
    for j=1:size(spk_tr,1); %for each seg
        spdind=trispdind(i);
        triind=tri_type(i);
        if ~isempty(spk_tr{j,i}')
            class_id(end+1)=sub2ind([length(trialdirs),length(spds)],segdirinds{triind}(j),spdind);
            spiketimes{length(class_id)}=spk_tr{j,i}';
        end
    end
end

%  class_id=class_id(1:500);
% spiketimes=spiketimes(1:500);

save([experiment,'WIdata.mat'],'class_id','spiketimes')
% WIgui %raster throws errors if any cells in spiketimes are empty arrays: no spikes that trial
WIscript_caller 
%

figure;imagesc(handles.decode.SPKCNTconfusionmatrix);colorbar
title(['confusion matrix from spike count for speed and direction: performance ',num2str(handles.decode.SPKCNTperf)])
figure;imagesc(handles.decode.WVconfusionmatrix);colorbar
title(['confusion matrix from wavelet info for speed and direction: performance ',num2str(handles.decode.WVperf)])
