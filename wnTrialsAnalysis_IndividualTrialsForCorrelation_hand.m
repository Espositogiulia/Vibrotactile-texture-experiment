clear;
clc;
%% Load the merged datasets
[lwdataHand(1).header,lwdataHand(1).data]=CLW_load('merged_epochs wn1 hand');
[lwdataHand(2).header,lwdataHand(2).data]=CLW_load('merged_epochs wn2 hand');
[lwdataHand(3).header,lwdataHand(3).data]=CLW_load('merged_epochs wn3 hand');
[lwdataHand(4).header,lwdataHand(4).data]=CLW_load('merged_epochs wn4 hand');
[lwdataHand(5).header,lwdataHand(5).data]=CLW_load('merged_epochs wn5 hand');
[lwdataHand(6).header,lwdataHand(6).data]=CLW_load('merged_epochs wn6 hand');
[lwdataHand(7).header,lwdataHand(7).data]=CLW_load('merged_epochs wn7 hand');
[lwdataHand(8).header,lwdataHand(8).data]=CLW_load('merged_epochs wn8 hand');

for i=1:8;
    lwdataHand(i).data=squeeze(lwdataHand(i).data);
end;

%%
%assign conditions

%stim_type : wn1-8
lwdataHand(1).stim_type=1;
lwdataHand(2).stim_type=2;
lwdataHand(3).stim_type=3;
lwdataHand(4).stim_type=4;
lwdataHand(5).stim_type=5;
lwdataHand(6).stim_type=6;
lwdataHand(7).stim_type=7;
lwdataHand(8).stim_type=8;


%% Frequency axis bins>frequency
tpxHand=1:size(lwdataHand(1).data,3);
tpxHand=(tpxHand-1)*lwdataHand(1).header.xstep;

%% Baseline correction and Z-score computation
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:8;
    tp=lwdataHand(i).data;
    lwdataHand(i).bl_data=tp;
  %  lwdataHand(i).Z_data=tp;
    for epochpos=1:size(tp,1);
        for channelpos=1:size(tp,2);
            for dx=1:size(tp,3);
                dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
                dxbins(dxbins<1)=[];
                dxbins(dxbins>size(tp,3))=[];
                tpy=tp(epochpos,channelpos,dxbins);
                mean_tpy=mean(tpy);
                SD_tpy=std(tpy);
                lwdataHand(i).bl_data(epochpos,channelpos,dx)=lwdataHand(i).data(epochpos,channelpos,dx)-mean_tpy;
               % lwdataHand(i).Z_data(epochpos,channelpos,dx)=(lwdataHand(i).data(epochpos,channelpos,dx)-mean_tpy)/SD_tpy;
            end;
        end;
    end;
end;

%% group-level average amplitude + average across scalp channels
% 
for i=1:8
    lwdataHand(i).avg=squeeze(mean(lwdataHand(i).data,1));
    lwdataHand(i).avg_bl=squeeze(mean(lwdataHand(i).bl_data,1));
   % lwdataHand(i).avg_Z=squeeze(mean(lwdataHand(i).Z_data,1));
    lwdataHand(i).chanavg_avg=squeeze(mean(lwdataHand(i).avg,1));
    lwdataHand(i).chanavg_avg_bl=squeeze(mean(lwdataHand(i).avg_bl,1));
    lwdataHand(i).chanavg_bl=squeeze(mean(lwdataHand(i).bl_data,2));
   % lwdataHand(i).chanavg_Z=squeeze(mean(lwdataHand(i).Z_data,2));   
    lwdataHand(i).chanavg=squeeze(mean(lwdataHand(i).data,2));
end;

%% zscore of average amplitude 
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:8;
    num_channels=size(lwdataHand(i).data,2);
    num_bins=size(lwdataHand(i).data,3);
    %avg
    for channelpos=1:num_channels;
        for dx=1:num_bins;
            dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
            dxbins(dxbins<1)=[];
            dxbins(dxbins>num_bins)=[];
            tpy=lwdataHand(i).avg(channelpos,dxbins);
            mean_tpy=mean(tpy);
            SD_tpy=std(tpy);
            lwdataHand(i).bl_avg(channelpos,dx)=lwdataHand(i).avg(channelpos,dx)-mean_tpy;
          %  lwdataHand(i).Z_avg(channelpos,dx)=(lwdataHand(i).avg(channelpos,dx)-mean_tpy)/SD_tpy;
        end;
    end;
   % lwdataHand(i).chanvg_Z_avg=mean(lwdataHand(i).Z_avg,1);
    %chanavg_avg
    for dx=1:num_bins;
        dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
        dxbins(dxbins<1)=[];
        dxbins(dxbins>num_bins)=[];
        tpy=lwdataHand(i).chanavg_avg(dxbins);
        mean_tpy=mean(tpy);
        SD_tpy=std(tpy);
        lwdataHand(i).bl_chanavg_avg(dx)=lwdataHand(i).chanavg_avg(dx)-mean_tpy;
       % lwdataHand(i).Z_chanavg_avg(dx)=(lwdataHand(i).chanavg_avg(dx)-mean_tpy)/SD_tpy;
    end;
end;
%% Base and oddball frequencies and corresponding bins
max_frequency=40;

base_freq=1:floor(250/8);
base_freq=base_freq*8;
oddball_freq=1:floor(250/1.6);
oddball_freq=oddball_freq*1.6;

%remove oddball harmonics that are also base harmonics
idx=[];
for i=1:length(base_freq);
    idx=[idx find(oddball_freq==base_freq(i))];
end;
oddball_freq(idx)=[];

%remove frequencies above max_frequency
base_freq=base_freq(find(base_freq<=max_frequency));
oddball_freq=oddball_freq(find(oddball_freq<=max_frequency));


%convert to bin positions
%base_freq=[];
%oddball_freq=[];
for i=1:length(base_freq);
    [a,b]=min(abs(tpxHand-base_freq(i)));
    base_freq_dx(i)=b;
end;
for i=1:length(oddball_freq);
    [a,b]=min(abs(tpxHand-oddball_freq(i)));
    oddball_freq_dx(i)=b;
end;

results.base_freq=base_freq;
results.oddball_freq=oddball_freq;
results.base_freq_dx=base_freq_dx;
results.oddball_freq_dx=oddball_freq_dx;
%% ttests to determine whether or not there is a periodic response

for i=1:8;
    %baseline_corrected amplitude at all base frequency harmonics, averaged across all channels
    %base frequencies
    tp_base=lwdataHand(i).chanavg_bl(:,base_freq_dx);
    mean_tp_base=mean(tp_base,2);
    [a,b]=ttest(mean_tp_base,0,'Tail','right');
    lwdataHand(i).ttest.base_chan_avg_bl=tp_base;
    lwdataHand(i).ttest.mean_base_chan_avg_bl=mean_tp_base;
    lwdataHand(i).ttest.base_significant=a;
    lwdataHand(i).ttest.base_pvalue=b;
    %oddball frequencies
    tp_oddball=lwdataHand(i).chanavg_bl(:,oddball_freq_dx);
    mean_tp_oddball=mean(tp_oddball,2);
    [a,b]=ttest(mean_tp_oddball,0,'Tail','right');
    lwdataHand(i).ttest.oddball_chan_avg_bl=tp_oddball;
    lwdataHand(i).ttest.mean_oddball_chan_avg_bl=mean_tp_oddball;
    lwdataHand(i).ttest.oddball_significant=a;
    lwdataHand(i).ttest.oddball_pvalue=b;
 
end;
%% ttests at each frequency

for i=1:8;
    lwdataHand(i).ttest_base_freq=[];
    lwdataHand(i).ttest_oddball_freq=[];
    tp=lwdataHand(i).chanavg_bl;
    %loop through base_frequencies
    for j=1:length(results.base_freq_dx);
        tp2=tp(:,results.base_freq_dx(j));
        lwdataHand(i).ttest_base_freq(j)=ttest(tp2,0,'Tail','right');
    end;
    %loop through oddball_frequencies
    for j=1:length(results.oddball_freq_dx);
        tp2=tp(:,results.oddball_freq_dx(j));
        lwdataHand(i).ttest_oddball_freq(j)=ttest(tp2,0,'Tail','right');
    end;
end;
%% plot spectra (averaged across subjects and channels with results of ttest
for i = 1:8
maxYamp(i) = ceil(max(lwdataHand(i).chanavg_avg_bl(1,12:1921))*100)/100; % to get ylim for each subplot
end

positions = zeros(8,4); % to position each subplot 
positions(1,2) = 0.895;
for i = 1:8
    positions(i,1) = 0.1;
    positions(i,3) = 0.4;
    positions(i,4) = 0.1;
    if i>1
    positions(i,2) = positions(i-1,2)-0.118;
    end
end

maxYampAll = max(maxYamp);
figure;
for i=1:8;
    %subplot(8,1,i);
    axes('Position',positions(i,:));
    h(i) = plot(tpxHand,lwdataHand(i).chanavg_avg_bl,'k-');
    
    box off
    xlim([0 40]);
    hold on
    %highlight significant base frequencies
   
    sig_base_freqs=results.base_freq(find(lwdataHand(i).ttest_base_freq));
    sig_base_amps=lwdataHand(i).chanavg_avg_bl(results.base_freq_dx(find(lwdataHand(i).ttest_base_freq)));
    scatter(sig_base_freqs,sig_base_amps,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
    
    stem(tpxHand(results.base_freq_dx(find(lwdataHand(i).ttest_base_freq))),...
        lwdataHand(i).chanavg_avg_bl(results.base_freq_dx(find(lwdataHand(i).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
    
    %highlight significant oddball frequencies
    sig_oddball_freqs=results.oddball_freq(find(lwdataHand(i).ttest_oddball_freq));
    sig_oddball_amps=lwdataHand(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdataHand(i).ttest_oddball_freq)));
    scatter(sig_oddball_freqs,sig_oddball_amps,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
    
    stem(tpxHand(results.oddball_freq_dx(find(lwdataHand(i).ttest_oddball_freq))),...
        lwdataHand(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdataHand(i).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    
    set(h(i),'LineWidth',1);
    set(gca,'xlim',[0.25 40]);
    set(gca,'ylim',[0 maxYampAll+0.01]);
    set(gca,'TickLength',[0 0]);
    a = get(gca,'XTickLabel');
    set(gca,'XTick',[10 20 30 40]);
    set(gca,'XTickLabel',a,'fontsize',16)
    set(gca,'YTick',[0 maxYampAll]);
    if i <8
        set(gca, 'XTickLabel', [])
    end
    %hold off;
end;

axes('visible','off')
xlabel('Frequency (Hz)','FontSize',24','FontWeight','bold','Position',[0.5,-0.09,1])
ylabel('Amplitude (\muV)','FontSize',24,'FontWeight','bold','Position',[-0.1,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')

%hold off;
%% find frequencies having Z score > threshold
Z_threshold=1.64;
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins


%compute average amplitude spectrum across all conditions
tp=[];
for i=1:8;
    tp(:,i)=lwdataHand(i).chanavg_avg;
end;
avg_tp=mean(tp,2);

%compute Z score using neighbouring bins
num_bins=length(avg_tp);
for dx=1:num_bins;
    dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
    dxbins(dxbins<1)=[];
    dxbins(dxbins>num_bins)=[];
    tpy=avg_tp(dxbins);
    mean_tpy=mean(tpy);
    SD_tpy=std(tpy);
    bl_avg_tp(dx)=avg_tp(dx)-mean_tpy;
   % Z_avg_tp(dx)=(avg_tp(dx)-mean_tpy)/SD_tpy;
end;
results.avgALL_chanavg_avg=avg_tp;
results.bl_avgALL_chanavg_avg=bl_avg_tp;
%results.Z_avgALL_chanavg_avg=Z_avg_tp;

%find "relevant" base frequencies (Z>threshold)
% sig_base_freq=zeros(size(results.base_freq_dx));
% for i=1:length(results.base_freq_dx);
%     if results.Z_avgALL_chanavg_avg(results.base_freq_dx(i))>=Z_threshold
%         sig_base_freq(i)=1;
%     end;
% end;
    
%find "relevant" oddball frequencies (Z>threshold)
% sig_oddball_freq=zeros(size(results.oddball_freq_dx));
% for i=1:length(results.oddball_freq_dx);
%     if results.Z_avgALL_chanavg_avg(results.oddball_freq_dx(i))>=Z_threshold
%         sig_oddball_freq(i)=1;
%     end;
% end;

%% Topographical plot (baseline-corrected data)

%find harmonics that were significant (ttest) in at least one condition
idx_base=[];
idx_oddball=[];

% load idx_base and idx_oddball from folder (they are the same as those
% selected in the main analysis)
results.idx_base=idx_base;
results.idx_oddball=idx_oddball;

%average (baseline-corrected) signal across the relevant harmonics
for i=1:8;
    lwdataHand(i).base_signal=lwdataHand(i).bl_avg(:,idx_base);
    lwdataHand(i).oddball_signal=lwdataHand(i).bl_avg(:,idx_oddball);
    lwdataHand(i).mean_base_signal=mean(lwdataHand(i).base_signal,2);
    lwdataHand(i).mean_oddball_signal=mean(lwdataHand(i).oddball_signal,2);
end;

%average (baseline-corrected) signal across the relevant harmonics for each
%participant
for i=1:8;
    lwdataHand(i).base_signal_individual=lwdataHand(i).chanavg_bl(:,idx_base);
    lwdataHand(i).oddball_signal_individual=lwdataHand(i).chanavg_bl(:,idx_oddball);
    lwdataHand(i).mean_base_signal_individual=mean(lwdataHand(i).base_signal_individual,2);
    lwdataHand(i).mean_oddball_signal_individual=mean(lwdataHand(i).oddball_signal_individual,2);
    
    for i = 1:8
mean_oddball_hand(i).mean_oddball_hand =  lwdataHand(i).mean_oddball_signal_individual;
end

end;

save('mean_oddball_hand_2410.mat','mean_oddball_hand');

%% average (baseline-corrected) signal across the relevant harmonics for each
%participant and without averaging across channels

for i=1:8;
    lwdataHand(i).base_signal_individual_chan=lwdataHand(i).bl_data(:,:,idx_base);
    lwdataHand(i).oddball_signal_individual_chan=lwdataHand(i).bl_data(:,:,idx_oddball);
    lwdataHand(i).mean_base_signal_individual_chan=mean(lwdataHand(i).base_signal_individual_chan,3);
    lwdataHand(i).mean_oddball_signal_individual_chan=mean(lwdataHand(i).oddball_signal_individual_chan,3);
end;

%% weighted average of oddball responses
for i = 1:8
for j=1:17
    lwdataHand(i).hand_noise_oddball_hand_template(j)=sum(lwdataHand(i).mean_oddball_signal_individual_chan(j,:)'.*hand_template);
end
end

for i = 1:8
    hand_noise_oddball_hand_template(i,:) = lwdataHand(i).hand_noise_oddball_hand_template;
end


% %% average (baseline-corrected) signal across the clusters of harmonics
% %across participants and without averaging across channels
% for i=1:8;
%     lwdataHand(i).oddball_signal_individual_chan_mean_cluster1=mean(mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([1:4])),3),1);
%     lwdataHand(i).oddball_signal_individual_chan_mean_cluster2=mean(mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([5:8])),3),1);
%     lwdataHand(i).oddball_signal_individual_chan_mean_cluster3=mean(mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([9:12])),3),1);
%     lwdataHand(i).oddball_signal_individual_chan_mean_cluster4=mean(mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([13:16])),3),1);
%     lwdataHand(i).oddball_signal_individual_chan_mean_cluster5=mean(mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([17:20])),3),1);
% end;
% 
% % average clusters across sequences
% 
% for i = 1:8
% all_seq_cluster1(i,:,:) = lwdataHand(i).oddball_signal_individual_chan_mean_cluster1(:,:);
% end
% 
% all_seq_cluster1 = squeeze(all_seq_cluster1);
% mean_cluster1 = mean(all_seq_cluster1);
% 
% for i = 1:8
% all_seq_cluster2(i,:,:) = lwdataHand(i).oddball_signal_individual_chan_mean_cluster2(:,:);
% end
% 
% all_seq_cluster2 = squeeze(all_seq_cluster2);
% mean_cluster2 = mean(all_seq_cluster2);
% 
% for i = 1:8
% all_seq_cluster3(i,:,:) = lwdataHand(i).oddball_signal_individual_chan_mean_cluster3(:,:);
% end
% 
% all_seq_cluster3 = squeeze(all_seq_cluster3);
% mean_cluster3 = mean(all_seq_cluster3);
% 
% for i = 1:8
% all_seq_cluster4(i,:,:) = lwdataHand(i).oddball_signal_individual_chan_mean_cluster4(:,:);
% end
% 
% all_seq_cluster4 = squeeze(all_seq_cluster4);
% mean_cluster4 = mean(all_seq_cluster4);
% 
% for i = 1:8
% all_seq_cluster5(i,:,:) = lwdataHand(i).oddball_signal_individual_chan_mean_cluster5(:,:);
% end
% 
% all_seq_cluster5 = squeeze(all_seq_cluster5);
% mean_cluster5 = mean(all_seq_cluster5);
% 
%         figure;
%  CLW_topoplot_vector(lwdataHand(i).header,mean_cluster1);
%  
%     figure;
%  CLW_topoplot_vector(lwdataHand(i).header,mean_cluster2);
%  
%     figure;
%  CLW_topoplot_vector(lwdataHand(i).header,mean_cluster3);
%  
%     figure;
%  CLW_topoplot_vector(lwdataHand(i).header,mean_cluster4);
%  
%     figure;
%  CLW_topoplot_vector(lwdataHand(i).header,mean_cluster5);
% 
% 
% 
% %average (baseline-corrected) signal across the clusters of harmonics
% %for each participant and without averaging across channels
% for i=1:8;
%     lwdataHand(i).oddball_signal_individual_chan_cluster1=mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([1:4])),3);
%     lwdataHand(i).oddball_signal_individual_chan_cluster2=mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([5:8])),3);
%     lwdataHand(i).oddball_signal_individual_chan_cluster3=mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([9:12])),3);
%     lwdataHand(i).oddball_signal_individual_chan_cluster4=mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([13:16])),3);
%     lwdataHand(i).oddball_signal_individual_chan_cluster5=mean(lwdataHand(i).bl_data(:,:,oddball_freq_dx([17:20])),3);
% end;

%% NEXT, without WN2
% mean_oddball_wn2excl_hand(1).mean_oddball_hand = ...
%     mean_oddball_hand(1).mean_oddball_hand;
% 
% mean_oddball_wn2excl_hand(2).mean_oddball_hand = ...
%     mean_oddball_hand(3).mean_oddball_hand;
% 
% mean_oddball_wn2excl_hand(3).mean_oddball_hand = ...
%     mean_oddball_hand(4).mean_oddball_hand;
% 
% mean_oddball_wn2excl_hand(4).mean_oddball_hand = ...
%     mean_oddball_hand(5).mean_oddball_hand;
% 
% mean_oddball_wn2excl_hand(5).mean_oddball_hand = ...
%     mean_oddball_hand(6).mean_oddball_hand;
% 
% mean_oddball_wn2excl_hand(6).mean_oddball_hand = ...
%     mean_oddball_hand(7).mean_oddball_hand;
% 
% mean_oddball_wn2excl_hand(7).mean_oddball_hand = ...
%     mean_oddball_hand(8).mean_oddball_hand;


%%
scatter([1:8],Dsum(1,:));  %,'FaceColor',[24, 116, 205]/255);
ylabel('Envelope dissimuilarity (AU)','FontSize',18);
xlabel('Spectrotemporal sequence','FontSize',18);
l =get(gca,'XTickLabel');
set(gca,'XTickLabel',l,'FontSize',16);
% hold on
% er = errorbar([1:8],grand_avg_oddball_hand,std_oddball_hand,...
%     'LineStyle','none','Color',[0 0 0]);

%% load the 8 individual A and B Hilbert-transformed stimuli
for n=1:8;
    [Ah Ad]=CLW_load(['hilTranA' num2str(n)]);
    [Bh Bd]=CLW_load(['hilTranB' num2str(n)]);
    
    A=squeeze(Ad);
    B=squeeze(Bd);
    
    D(n,:)=abs(B-A);
    Dsum(n)=sum(D(n,:));
    Dsum_onset(n)=sum(D(n,1:2210));
end;
%

wn1hand = mean_oddball_hand(1).mean_oddball_hand;
wn2hand = mean_oddball_hand(2).mean_oddball_hand;
wn3hand = mean_oddball_hand(3).mean_oddball_hand;
wn4hand = mean_oddball_hand(4).mean_oddball_hand;
wn5hand = mean_oddball_hand(5).mean_oddball_hand;
wn6hand = mean_oddball_hand(6).mean_oddball_hand;
wn7hand = mean_oddball_hand(7).mean_oddball_hand;
wn8hand = mean_oddball_hand(8).mean_oddball_hand;
oddball_allHand = [wn1hand wn2hand wn3hand wn4hand wn5hand wn6hand wn7hand wn8hand];
%oddball_hand_excl_Wn2 = [wn1hand wn3hand wn4hand wn5hand wn6hand wn7hand wn8hand];


%%
filename = 'analysis_individual_wn_trials_hand3010.mat';
save (filename, '-v7.3');
