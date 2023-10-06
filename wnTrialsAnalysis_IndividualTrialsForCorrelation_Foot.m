clear;
clc;
%% Load the merged datasets
[lwdataFoot(1).header,lwdataFoot(1).data]=CLW_load('merged_epochs wn1 foot');
[lwdataFoot(2).header,lwdataFoot(2).data]=CLW_load('merged_epochs wn2 foot');
[lwdataFoot(3).header,lwdataFoot(3).data]=CLW_load('merged_epochs wn3 foot');
[lwdataFoot(4).header,lwdataFoot(4).data]=CLW_load('merged_epochs wn4 foot');
[lwdataFoot(5).header,lwdataFoot(5).data]=CLW_load('merged_epochs wn5 foot');
[lwdataFoot(6).header,lwdataFoot(6).data]=CLW_load('merged_epochs wn6 foot');
[lwdataFoot(7).header,lwdataFoot(7).data]=CLW_load('merged_epochs wn7 foot');
[lwdataFoot(8).header,lwdataFoot(8).data]=CLW_load('merged_epochs wn8 foot');

for i=1:8;
    lwdataFoot(i).data=squeeze(lwdataFoot(i).data);
end;

%%
%assign conditions

%stim_type : wn1-8
lwdataFoot(1).stim_type=1;
lwdataFoot(2).stim_type=2;
lwdataFoot(3).stim_type=3;
lwdataFoot(4).stim_type=4;
lwdataFoot(5).stim_type=5;
lwdataFoot(6).stim_type=6;
lwdataFoot(7).stim_type=7;
lwdataFoot(8).stim_type=8;


%% Frequency axis bins>frequency
tpxFoot=1:size(lwdataFoot(1).data,3);
tpxFoot=(tpxFoot-1)*lwdataFoot(1).header.xstep;

%% Baseline correction and Z-score computation
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:8;
    tp=lwdataFoot(i).data;
    lwdataFoot(i).bl_data=tp;
   % lwdataFoot(i).Z_data=tp;
    for epochpos=1:size(tp,1);
        for channelpos=1:size(tp,2);
            for dx=1:size(tp,3);
                dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
                dxbins(dxbins<1)=[];
                dxbins(dxbins>size(tp,3))=[];
                tpy=tp(epochpos,channelpos,dxbins);
                mean_tpy=mean(tpy);
                SD_tpy=std(tpy);
                lwdataFoot(i).bl_data(epochpos,channelpos,dx)=lwdataFoot(i).data(epochpos,channelpos,dx)-mean_tpy;
               % lwdataFoot(i).Z_data(epochpos,channelpos,dx)=(lwdataFoot(i).data(epochpos,channelpos,dx)-mean_tpy)/SD_tpy;
            end;
        end;
    end;
end;

%% group-level average amplitude + average across scalp channels
% 
for i=1:8
    lwdataFoot(i).avg=squeeze(mean(lwdataFoot(i).data,1));
    lwdataFoot(i).avg_bl=squeeze(mean(lwdataFoot(i).bl_data,1));
  %  lwdataFoot(i).avg_Z=squeeze(mean(lwdataFoot(i).Z_data,1));
    lwdataFoot(i).chanavg_avg=squeeze(mean(lwdataFoot(i).avg,1));
    lwdataFoot(i).chanavg_avg_bl=squeeze(mean(lwdataFoot(i).avg_bl,1));
    lwdataFoot(i).chanavg_bl=squeeze(mean(lwdataFoot(i).bl_data,2));
  %  lwdataFoot(i).chanavg_Z=squeeze(mean(lwdataFoot(i).Z_data,2));   
    lwdataFoot(i).chanavg=squeeze(mean(lwdataFoot(i).data,2));
end;

%% zscore of average amplitude 
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:8;
    num_channels=size(lwdataFoot(i).data,2);
    num_bins=size(lwdataFoot(i).data,3);
    %avg
    for channelpos=1:num_channels;
        for dx=1:num_bins;
            dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
            dxbins(dxbins<1)=[];
            dxbins(dxbins>num_bins)=[];
            tpy=lwdataFoot(i).avg(channelpos,dxbins);
            mean_tpy=mean(tpy);
            SD_tpy=std(tpy);
            lwdataFoot(i).bl_avg(channelpos,dx)=lwdataFoot(i).avg(channelpos,dx)-mean_tpy;
           % lwdataFoot(i).Z_avg(channelpos,dx)=(lwdataFoot(i).avg(channelpos,dx)-mean_tpy)/SD_tpy;
        end;
    end;
  %  lwdataFoot(i).chanvg_Z_avg=mean(lwdataFoot(i).Z_avg,1);
    %chanavg_avg
    for dx=1:num_bins;
        dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
        dxbins(dxbins<1)=[];
        dxbins(dxbins>num_bins)=[];
        tpy=lwdataFoot(i).chanavg_avg(dxbins);
        mean_tpy=mean(tpy);
        SD_tpy=std(tpy);
        lwdataFoot(i).bl_chanavg_avg(dx)=lwdataFoot(i).chanavg_avg(dx)-mean_tpy;
     %   lwdataFoot(i).Z_chanavg_avg(dx)=(lwdataFoot(i).chanavg_avg(dx)-mean_tpy)/SD_tpy;
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
    [a,b]=min(abs(tpxFoot-base_freq(i)));
    base_freq_dx(i)=b;
end;
for i=1:length(oddball_freq);
    [a,b]=min(abs(tpxFoot-oddball_freq(i)));
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
    tp_base=lwdataFoot(i).chanavg_bl(:,base_freq_dx);
    mean_tp_base=mean(tp_base,2);
    [a,b]=ttest(mean_tp_base,0,'Tail','right');
    lwdataFoot(i).ttest.base_chan_avg_bl=tp_base;
    lwdataFoot(i).ttest.mean_base_chan_avg_bl=mean_tp_base;
    lwdataFoot(i).ttest.base_significant=a;
    lwdataFoot(i).ttest.base_pvalue=b;
    %oddball frequencies
    tp_oddball=lwdataFoot(i).chanavg_bl(:,oddball_freq_dx);
    mean_tp_oddball=mean(tp_oddball,2);
    [a,b]=ttest(mean_tp_oddball,0,'Tail','right');
    lwdataFoot(i).ttest.oddball_chan_avg_bl=tp_oddball;
    lwdataFoot(i).ttest.mean_oddball_chan_avg_bl=mean_tp_oddball;
    lwdataFoot(i).ttest.oddball_significant=a;
    lwdataFoot(i).ttest.oddball_pvalue=b;
 
end;
%% ttests at each frequency

for i=1:8;
    lwdataFoot(i).ttest_base_freq=[];
    lwdataFoot(i).ttest_oddball_freq=[];
    tp=lwdataFoot(i).chanavg_bl;
    %loop through base_frequencies
    for j=1:length(results.base_freq_dx);
        tp2=tp(:,results.base_freq_dx(j));
        lwdataFoot(i).ttest_base_freq(j)=ttest(tp2,0,'Tail','right');
    end;
    %loop through oddball_frequencies
    for j=1:length(results.oddball_freq_dx);
        tp2=tp(:,results.oddball_freq_dx(j));
        lwdataFoot(i).ttest_oddball_freq(j)=ttest(tp2,0,'Tail','right');
    end;
end;
%% plot spectra (averaged across subjects and channels with results of ttest

% figure;
% for i=1:8;
%     subplot(1,8,i);
%     plot(tpxFoot,lwdataFoot(i).chanavg_avg_bl);
%     xlim([0 40]);
%     hold on
%     %highlight significant base frequencies
%     sig_base_freqs=results.base_freq(find(lwdataFoot(i).ttest_base_freq));
%     sig_base_amps=lwdataFoot(i).chanavg_avg_bl(results.base_freq_dx(find(lwdataFoot(i).ttest_base_freq)));
%     scatter(sig_base_freqs,sig_base_amps,'k','filled');
%     %highlight significant oddball frequencies
%     sig_oddball_freqs=results.oddball_freq(find(lwdataFoot(i).ttest_oddball_freq));
%     sig_oddball_amps=lwdataFoot(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdataFoot(i).ttest_oddball_freq)));
%     scatter(sig_oddball_freqs,sig_oddball_amps,'r','filled');
%     hold off;
%     
%     
% end;

for i = 1:8
maxYamp(i) = ceil(max(lwdataFoot(i).chanavg_avg_bl(1,12:1921))*100)/100; % to get ylim for each subplot
end

positions = zeros(8,4); % to position each subplot 
positions(1,2) = 0.895;
for i = 1:8
    positions(i,1) = 0.58;
    positions(i,3) = 0.4;
    positions(i,4) = 0.1;
    if i>1
    positions(i,2) = positions(i-1,2)-0.118;
    end
end

% titlePositions = zeros(8,3);
% titlePositions(:,1) = 40.513;
% titlePositions(:,3) = 17.321;
% titlePositions(:,2) = 0.03;

maxYampAll = max(maxYamp);
%figure;
for i=1:8;
    %subplot(8,1,i);
    axes('Position',positions(i,:));
    h(i) = plot(tpxFoot,lwdataFoot(i).chanavg_avg_bl,'k-');
%     title('prova','Position',titlePositions(i,:));
    box off
    xlim([0 40]);
    hold on
    %highlight significant base frequencies
   
    sig_base_freqs=results.base_freq(find(lwdataFoot(i).ttest_base_freq));
    sig_base_amps=lwdataFoot(i).chanavg_avg_bl(results.base_freq_dx(find(lwdataFoot(i).ttest_base_freq)));
    scatter(sig_base_freqs,sig_base_amps,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
    
    stem(tpxFoot(results.base_freq_dx(find(lwdataFoot(i).ttest_base_freq))),...
        lwdataFoot(i).chanavg_avg_bl(results.base_freq_dx(find(lwdataFoot(i).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
    
    %highlight significant oddball frequencies
    sig_oddball_freqs=results.oddball_freq(find(lwdataFoot(i).ttest_oddball_freq));
    sig_oddball_amps=lwdataFoot(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdataFoot(i).ttest_oddball_freq)));
    scatter(sig_oddball_freqs,sig_oddball_amps,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
    
    stem(tpxFoot(results.oddball_freq_dx(find(lwdataFoot(i).ttest_oddball_freq))),...
        lwdataFoot(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdataFoot(i).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    
    set(h(i),'LineWidth',1);
    set(gca,'xlim',[0.25 40]);
    %set(gca,'ylim',[0 maxYamp(i)+0.01]);
    set(gca,'ylim',[0 maxYampAll+0.01]);
    set(gca,'TickLength',[0 0]);
    set(gca,'XTick',[10 20 30 40]);
    a = get(gca,'XTickLabel');
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

hold off;

%% find frequencies having Z score > threshold
Z_threshold=1.64;
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins


%compute average amplitude spectrum across all conditions
tp=[];
for i=1:8;
    tp(:,i)=lwdataFoot(i).chanavg_avg;
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
%     
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
for i=1:8;
    idx_base=[idx_base results.base_freq_dx(find(lwdataFoot(i).ttest_base_freq))];
    idx_oddball=[idx_oddball results.oddball_freq_dx(find(lwdataFoot(i).ttest_oddball_freq))];
end;
idx_base=sort(unique(idx_base));
idx_oddball=sort(unique(idx_oddball));

results.idx_base=idx_base;
results.idx_oddball=idx_oddball;

%average (baseline-corrected) signal across the relevant harmonics
for i=1:8;
    lwdataFoot(i).base_signal=lwdataFoot(i).bl_avg(:,idx_base);
    lwdataFoot(i).oddball_signal=lwdataFoot(i).bl_avg(:,idx_oddball);
    lwdataFoot(i).mean_base_signal=mean(lwdataFoot(i).base_signal,2);
    lwdataFoot(i).mean_oddball_signal=mean(lwdataFoot(i).oddball_signal,2);
end;

%average (baseline-corrected) signal across the relevant harmonics for each
%participant
for i=1:8;
    lwdataFoot(i).base_signal_individual=lwdataFoot(i).chanavg_bl(:,idx_base);
    lwdataFoot(i).oddball_signal_individual=lwdataFoot(i).chanavg_bl(:,idx_oddball);
    lwdataFoot(i).mean_base_signal_individual=mean(lwdataFoot(i).base_signal_individual,2);
    lwdataFoot(i).mean_oddball_signal_individual=mean(lwdataFoot(i).oddball_signal_individual,2);
end;

for i = 1:8
mean_oddball_foot(i).mean_oddball_signal_individual =  lwdataFoot(i).mean_oddball_signal_individual;
end

save('mean_oddball_foot_1509.mat','mean_oddball_foot');

%% NEXT, without WN2
mean_oddball_wn2excl_foot(1).mean_oddball_foot = ...
    mean_oddball_foot(1).mean_oddball_signal_individual;

mean_oddball_wn2excl_foot(2).mean_oddball_foot = ...
    mean_oddball_foot(3).mean_oddball_signal_individual;

mean_oddball_wn2excl_foot(3).mean_oddball_foot = ...
    mean_oddball_foot(4).mean_oddball_signal_individual;

mean_oddball_wn2excl_foot(4).mean_oddball_foot = ...
    mean_oddball_foot(5).mean_oddball_signal_individual;

mean_oddball_wn2excl_foot(5).mean_oddball_foot = ...
    mean_oddball_foot(6).mean_oddball_signal_individual;

mean_oddball_wn2excl_foot(6).mean_oddball_foot = ...
    mean_oddball_foot(7).mean_oddball_signal_individual;

mean_oddball_wn2excl_foot(7).mean_oddball_foot = ...
    mean_oddball_foot(8).mean_oddball_signal_individual;

save('mean_oddball_wn2excl_foot_1509.mat','mean_oddball_wn2excl_foot');
%%
for n=1:8;
    [Ah Ad]=CLW_load(['hilTranA' num2str(n)]);
    [Bh Bd]=CLW_load(['hilTranB' num2str(n)]);
    
    A=squeeze(Ad);
    B=squeeze(Bd);
    
    D(n,:)=abs(B-A);
    Dsum(n)=sum(D(n,:));
    Dsum_onset(n)=sum(D(n,1:1837));
end;
% 
Dsum = repmat(Dsum,17,1);
Dsum_onset= repmat(Dsum_onset,17,1);

wn1Foot = mean_oddball_foot(1).mean_oddball_signal_individual;
wn2Foot = mean_oddball_foot(2).mean_oddball_signal_individual;
wn3Foot = mean_oddball_foot(3).mean_oddball_signal_individual;
wn4Foot = mean_oddball_foot(4).mean_oddball_signal_individual;
wn5Foot = mean_oddball_foot(5).mean_oddball_signal_individual;
wn6Foot = mean_oddball_foot(6).mean_oddball_signal_individual;
wn7Foot = mean_oddball_foot(7).mean_oddball_signal_individual;
wn8Foot = mean_oddball_foot(8).mean_oddball_signal_individual;
oddball_allFoot = [wn1Foot wn2Foot wn3Foot wn4Foot wn5Foot wn6Foot wn7Foot wn8Foot];
oddball_foot_excl_Wn2 = [wn1Foot wn3Foot wn4Foot wn5Foot wn6Foot wn7Foot wn8Foot];

Dsum_excl_wn2 = Dsum(1,[1 3:8]);
Dsum_excl_wn2  = repmat(Dsum_excl_wn2, 17,1);
Dsum_onset_excl_wn2 = Dsum_onset(1,[1 3:8]);
Dsum_onset_excl_wn2  = repmat(Dsum_onset_excl_wn2, 17,1);
%% correlation with wn2
[RFootwn2,PFootwn2]= corrcoef(Dsum(:,:),oddball_allFoot(:,:));
[RFootwn2_onset,PFootwn2_onset]= corrcoef(Dsum_onset(:,:),oddball_allFoot(:,:));

%% correlation without wn2
[RFoot_excl_wn2,PFoot_excl_wn2]= corrcoef(Dsum(:,[1 3:8]),oddball_allFoot(:,[1 3:8]));
[RFoot_onset_excl_wn2,PFoot_onset_excl_wn2]= corrcoef(Dsum_onset(:,[1 3:8]),oddball_allFoot(:,[1 3:8]));
% 
% 
 %% save workspace
 
filename = 'analyisis_individual_wn_trials_foot_1509.mat';
save (filename, '-v7.3');
