%% Load the merged datasets
cd('C:\data\gesposito\xpVibroAnalysisFinal\FullAnalysis_Wn2Included');
[lwdata(1).header,lwdata(1).data]=CLW_load('merged_epochs_hand_s4');
[lwdata(2).header,lwdata(2).data]=CLW_load('merged_epochs_foot_s4');
[lwdata(3).header,lwdata(3).data]=CLW_load('merged_epochs_hand_s2');
[lwdata(4).header,lwdata(4).data]=CLW_load('merged_epochs_foot_s2');

for i=1:4;
    lwdata(i).data=squeeze(lwdata(i).data);
end;

%%
%assign conditions
%stim_location : 1=hand 2=foot
lwdata(1).stim_location=1;
lwdata(2).stim_location=2;
lwdata(3).stim_location=1;
lwdata(4).stim_location=2;
%stim_type : 1=frequency/intensity 2=noise
lwdata(1).stim_type=1;
lwdata(2).stim_type=1;
lwdata(3).stim_type=2;
lwdata(4).stim_type=2;

%% Frequency axis bins>frequency
tpx=1:size(lwdata(1).data,3);
tpx=(tpx-1)*lwdata(1).header.xstep;

%% Baseline correction and Z-score computation
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:4;
    tp=lwdata(i).data;
    lwdata(i).bl_data=tp;
    lwdata(i).Z_data=tp;
    for epochpos=1:size(tp,1);
        for channelpos=1:size(tp,2);
            for dx=1:size(tp,3);
                dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
                dxbins(dxbins<1)=[];
                dxbins(dxbins>size(tp,3))=[];
                tpy=tp(epochpos,channelpos,dxbins);
                mean_tpy=mean(tpy);
                SD_tpy=std(tpy);
                lwdata(i).bl_data(epochpos,channelpos,dx)=lwdata(i).data(epochpos,channelpos,dx)-mean_tpy;
                lwdata(i).Z_data(epochpos,channelpos,dx)=(lwdata(i).data(epochpos,channelpos,dx)-mean_tpy)/SD_tpy;
            end;
        end;
    end;
end;

%% group-level average amplitude + average across scalp channels

for i=1:4
    lwdata(i).avg=squeeze(mean(lwdata(i).data,1));
    lwdata(i).avg_bl=squeeze(mean(lwdata(i).bl_data,1));
    lwdata(i).avg_Z=squeeze(mean(lwdata(i).Z_data,1));
    lwdata(i).chanavg_avg=squeeze(mean(lwdata(i).avg,1));
    lwdata(i).chanavg_avg_bl=squeeze(mean(lwdata(i).avg_bl,1));
    lwdata(i).chanavg_bl=squeeze(mean(lwdata(i).bl_data,2));
    lwdata(i).chanavg_Z=squeeze(mean(lwdata(i).Z_data,2));   
    lwdata(i).chanavg=squeeze(mean(lwdata(i).data,2));
end;

%% zscore of average amplitude 
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:4;
    num_channels=size(lwdata(i).data,2);
    num_bins=size(lwdata(i).data,3);
    %avg
    for channelpos=1:num_channels;
        for dx=1:num_bins;
            dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
            dxbins(dxbins<1)=[];
            dxbins(dxbins>num_bins)=[];
            tpy=lwdata(i).avg(channelpos,dxbins);
            mean_tpy=mean(tpy);
            SD_tpy=std(tpy);
            lwdata(i).bl_avg(channelpos,dx)=lwdata(i).avg(channelpos,dx)-mean_tpy;
           lwdata(i).Z_avg(channelpos,dx)=(lwdata(i).avg(channelpos,dx)-mean_tpy)/SD_tpy;
        end;
    end;
    lwdata(i).chanvg_Z_avg=mean(lwdata(i).Z_avg,1);
    %chanavg_avg
    for dx=1:num_bins;
        dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
        dxbins(dxbins<1)=[];
        dxbins(dxbins>num_bins)=[];
        tpy=lwdata(i).chanavg_avg(dxbins);
        mean_tpy=mean(tpy);
        SD_tpy=std(tpy);
        lwdata(i).bl_chanavg_avg(dx)=lwdata(i).chanavg_avg(dx)-mean_tpy;
        lwdata(i).Z_chanavg_avg(dx)=(lwdata(i).chanavg_avg(dx)-mean_tpy)/SD_tpy;
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
    [a,b]=min(abs(tpx-base_freq(i)));
    base_freq_dx(i)=b;
end;
for i=1:length(oddball_freq);
    [a,b]=min(abs(tpx-oddball_freq(i)));
    oddball_freq_dx(i)=b;
end;

results.base_freq=base_freq;
results.oddball_freq=oddball_freq;
results.base_freq_dx=base_freq_dx;
results.oddball_freq_dx=oddball_freq_dx;


%% ttests to determine whether or not there is a periodic response

for i=1:4;
    %baseline_corrected amplitude at all base frequency harmonics, averaged across all channels
    %base frequencies
    tp_base=lwdata(i).chanavg_bl(:,base_freq_dx);
    mean_tp_base=mean(tp_base,2);
    [a,b]=ttest(mean_tp_base,0,'Tail','right');
    lwdata(i).ttest.base_chan_avg_bl=tp_base;
    lwdata(i).ttest.mean_base_chan_avg_bl=mean_tp_base;
    lwdata(i).ttest.base_significant=a;
    lwdata(i).ttest.base_pvalue=b;
    %oddball frequencies
    tp_oddball=lwdata(i).chanavg_bl(:,oddball_freq_dx);
    mean_tp_oddball=mean(tp_oddball,2);
    [a,b]=ttest(mean_tp_oddball,0,'Tail','right');
    lwdata(i).ttest.oddball_chan_avg_bl=tp_oddball;
    lwdata(i).ttest.mean_oddball_chan_avg_bl=mean_tp_oddball;
    lwdata(i).ttest.oddball_significant=a;
    lwdata(i).ttest.oddball_pvalue=b;
 
end;

%% ttests at each frequency

for i=1:4;
    lwdata(i).ttest_base_freq=[];
    lwdata(i).ttest_oddball_freq=[];
    tp=lwdata(i).chanavg_bl;
    %loop through base_frequencies
    for j=1:length(results.base_freq_dx);
        tp2=tp(:,results.base_freq_dx(j));
        lwdata(i).ttest_base_freq(j)=ttest(tp2,0,'Tail','right');
    end;
    %loop through oddball_frequencies
    for j=1:length(results.oddball_freq_dx);
        tp2=tp(:,results.oddball_freq_dx(j));
        lwdata(i).ttest_oddball_freq(j)=ttest(tp2,0,'Tail','right');
    end;
end;

%% extra plot
% for i=1:4;
%     tpBase=[];
%     tpOdd=[];
%     for epochpos= 1:17
%         for channelpos=1:64
%             %baseline_corrected amplitude at all base frequency harmonics, averaged across all channels
%             %base frequencies
%             for j=1:length(base_freq_dx)
%                 tpBase(epochpos,channelpos,j,1:2)=squeeze(lwdata(i).bl_data(epochpos,channelpos,base_freq_dx(j)));
%             end;
%             %oddball frequencies
%             for j=1:length(oddball_freq_dx)
%                 tpOdd(epochpos,channelpos,j,1:2)=squeeze(lwdata(i).bl_data(epochpos,channelpos,oddball_freq_dx(j)));
%             end;
%         end
%     end
%     lwdata(i).bl_data_avg_harm_base=squeeze(mean(tpBase,3));
%     lwdata(i).bl_data_avg_harm_oddball=squeeze(mean(tpOdd,3));
% end
% 
% base_hand_s2=lwdata(1).bl_data_avg_harm_base;
% base_foot_s2=lwdata(2).bl_data_avg_harm_base;
% base_hand_s4=lwdata(3).bl_data_avg_harm_base;
% base_foot_s4=lwdata(4).bl_data_avg_harm_base;
% oddball_hand_s2=lwdata(1).bl_data_avg_harm_oddball;
% oddball_foot_s2=lwdata(2).bl_data_avg_harm_oddball;
% oddball_hand_s4=lwdata(3).bl_data_avg_harm_oddball;
% oddball_foot_s4=lwdata(4).bl_data_avg_harm_oddball;



%% find frequencies having Z score > threshold
Z_threshold=1.64;
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins
% 
% 
% %compute average amplitude spectrum across all conditions
tp=[];
for i=1:4;
    tp(:,i)=lwdata(i).chanavg_avg;
end;
avg_tp=mean(tp,2);
% 
% %compute Z score using neighbouring bins
num_bins=length(avg_tp);
for dx=1:num_bins;
    dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
    dxbins(dxbins<1)=[];
    dxbins(dxbins>num_bins)=[];
    tpy=avg_tp(dxbins);
    mean_tpy=mean(tpy);
    SD_tpy=std(tpy);
    bl_avg_tp(dx)=avg_tp(dx)-mean_tpy;
    Z_avg_tp(dx)=(avg_tp(dx)-mean_tpy)/SD_tpy;
end;
results.avgALL_chanavg_avg=avg_tp;
results.bl_avgALL_chanavg_avg=bl_avg_tp;
results.Z_avgALL_chanavg_avg=z_avg_tp;

%% plot spectra and topogtaphies (averaged across subjects and channels with results of ttest
titles= {'Frequency/intensity contrast (hand)' 'Frequency/intensity contrast (foot)' 'Spectrotemporal contrast (hand)' 'Spectrotemporal contrast (foot)'};
for i = 1:4
maxYamp(i) = ceil(max(lwdata(i).chanavg_avg_bl(1,12:1921))*100)/100; % to get ylim for each subplot
end

positions = zeros(4,4); % to position each subplot 
positions(1,2) = 0.78;
for i = 1:4
    positions(i,1) = 0.1;
    positions(i,3) = 0.85;
    positions(i,4) = 0.2;
    if i>1
    positions(i,2) = positions(i-1,2)-0.23;
    end
end

%%topoplot parameters
idx_base=[];
idx_oddball=[];
for i=1:4;
    idx_base=[idx_base results.base_freq_dx(find(lwdata(i).ttest_base_freq))];
    idx_oddball=[idx_oddball results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))];
end;
idx_base=sort(unique(idx_base));
idx_oddball=sort(unique(idx_oddball));

results.idx_base=idx_base;
results.idx_oddball=idx_oddball;

%average (baseline-corrected) signal across the relevant harmonics
for i=1:4;
    lwdata(i).base_signal=lwdata(i).bl_avg(:,idx_base);
    lwdata(i).oddball_signal=lwdata(i).bl_avg(:,idx_oddball);
    lwdata(i).mean_base_signal=mean(lwdata(i).base_signal,2);
    lwdata(i).mean_oddball_signal=mean(lwdata(i).oddball_signal,2);
end;

%set limits for colorbars
for i = 1:4
maxTopoAmpBase(i) = ceil(max(lwdata(i).mean_base_signal)*100)/100; % to get ylim for each subplot
end

for i = 1:4
minTopoAmpBase(i) = min(lwdata(i).mean_base_signal)*100/100; % to get ylim for each subplot
end

for i = 1:4
maxTopoAmpOddball(i) = ceil(max(lwdata(i).mean_oddball_signal)*100)/100; % to get ylim for each subplot
end

for i = 1:4
minTopoAmpOddball(i) = min(lwdata(i).mean_oddball_signal); % to get ylim for each subplot
end


%topoplot of base response
positionsTopoBase = zeros(4,4); % to position each subplot 
positionsTopoBase(1,2) = 0.8;
for i = 1:4
    positionsTopoBase(i,1) = 0.72;
    positionsTopoBase(i,3) = 0.13;
    positionsTopoBase(i,4) = 0.13;
    if i>1
    positionsTopoBase(i,2) = positionsTopoBase(i-1,2)-0.23;
    end
end

positionsTopoOddball = zeros(4,4); % to position each subplot 
positionsTopoOddball(1,2) = 0.8;
for i = 1:4
    positionsTopoOddball(i,1) = 0.85;
    positionsTopoOddball(i,3) = 0.13;
    positionsTopoOddball(i,4) = 0.13;
    if i>1
    positionsTopoOddball(i,2) = positionsTopoOddball(i-1,2)-0.23;
    end
end

%%make the plot
figure;
for i=1:4;
    axes('Position',positions(i,:));
    h(i)= plot(tpx,lwdata(i).chanavg_avg_bl,'k-');
    text(20,maxYamp(i),titles{i},'FontSize',18,'FontWeight','bold','HorizontalAlignment','center');
    box off
    xlim([0 40]);
    hold on
    %highlight significant base frequencies
    sig_base_freqs=results.base_freq(find(lwdata(i).ttest_base_freq));
    sig_base_amps=lwdata(i).chanavg_avg_bl(results.base_freq_dx(find(lwdata(i).ttest_base_freq)));
    scatter(sig_base_freqs,sig_base_amps,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
    stem(tpx(results.base_freq_dx(find(lwdata(i).ttest_base_freq))),...
        lwdata(i).chanavg_avg_bl(results.base_freq_dx(find(lwdata(i).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
    
    %highlight significant oddball frequencies
    sig_oddball_freqs=results.oddball_freq(find(lwdata(i).ttest_oddball_freq));
    sig_oddball_amps=lwdata(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq)));
    
    scatter(sig_oddball_freqs,sig_oddball_amps,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
    
    stem(tpx(results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))),...
        lwdata(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    
    %scatter(sig_oddball_freqs,sig_oddball_amps,'r','filled');
    %hold off;
    
    set(h(i),'LineWidth',1);
    set(gca,'xlim',[0.25 40]); %change back to 40
    set(gca,'ylim',[0 maxYamp(i)+0.01]);
    set(gca,'TickLength',[0 0]);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',16)
    set(gca,'YTick',[0 maxYamp(i)]);
%     if i <4
%         set(gca, 'XTickLabel', [])
%     end
end;
axes('visible','off')
xlabel('Frequency (Hz)','FontSize',16','FontWeight','bold','Position',[0.5,-0.09,1])
ylabel('Amplitude (\muV)','FontSize',16,'FontWeight','bold','Position',[-0.1,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gcf, 'Color', 'White');
%hold off;

%%Topographical plot (baseline-corrected data)

for i=1:4;
    %subplot(2,2,i);
    axes('Position',positionsTopoBase(i,:));
   % if i == 1 || i == 2
    title('Base response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(i).header,lwdata(i).mean_base_signal,'maplimits',[minTopoAmpBase(i) maxTopoAmpBase(i)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpBase(i) maxTopoAmpBase(i)])
    set(cB,'FontSize',14);
    ylabel(cB,'\muV','rotation',0, 'Position',[7 (minTopoAmpBase(i)+minTopoAmpBase(i)*0.1) 1]);
  
end;
hold on;
%topoplot of oddball response
%figure;
for i=1:4;
    axes('Position',positionsTopoOddball(i,:));
    title('Oddball response','fontweight','bold','fontsize',14);

    CLW_topoplot_vector(lwdata(i).header,lwdata(i).mean_oddball_signal,'maplimits',[minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
    cO = colorbar;
    set(cO, 'ylim', [minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
   % set(cO,'Location','WestOutside');
    set(cO,'FontSize',14);
    ylabel(cO,'\muV','rotation',0, 'Position',[7 (minTopoAmpOddball(i)+minTopoAmpOddball(i)*0.1) 1]);
end;

%% Topographical plot (Z-score data)
% 
%find harmonics that were significant (ttest) in at least one condition
idx_base=[];
idx_oddball=[];
for i=1:4;
    idx_base=[idx_base results.base_freq_dx(find(lwdata(i).ttest_base_freq))];
    idx_oddball=[idx_oddball results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))];
end;
idx_base=sort(unique(idx_base));
idx_oddball=sort(unique(idx_oddball));

results.sig_idx_base=idx_base;
results.sig_idx_oddball=idx_oddball;

%average (baseline-corrected) signal across the relevant harmonics
% for i=1:4;
%     lwdata(i).Z_base_signal=lwdata(i).Z_avg(:,idx_base);
%     lwdata(i).Z_oddball_signal=lwdata(i).Z_avg(:,idx_oddball);
%     lwdata(i).Z_mean_base_signal=mean(lwdata(i).Z_base_signal,2);
%     lwdata(i).Z_mean_oddball_signal=mean(lwdata(i).Z_oddball_signal,2);
% end;

%% Create template using base response of intensity/frequency contrast

for i=1:4;
    stim_locations(i)=lwdata(i).stim_location;
    stim_type(i)=lwdata(i).stim_type;
end;

%hand_template and foot_template
%hand
a=find(stim_locations==1);
b=find(stim_type==1); %frequency/intensity
c_hand=intersect(a,b); %the dataset for hand stimulation with freq/intensity contrast
idx_hand=results.base_freq_dx(find(lwdata(c_hand).ttest_base_freq));
%foot
a=find(stim_locations==2); %foot
b=find(stim_type==1); %frequency/intensity
c_foot=intersect(a,b); %the dataset for foot stimulation with freq/intensity contrast
idx_foot=results.base_freq_dx(find(lwdata(c_foot).ttest_base_freq));

results.idx_hand=idx_hand;
results.idx_foot=idx_foot;


%take topography of mean of baseline_corrected amplitudes
foot_template=mean(lwdata(c_foot).avg_bl(:,idx_foot),2);
hand_template=mean(lwdata(c_hand).avg_bl(:,idx_hand),2);

%set to zero negative baseline_corrected amplitudes
foot_template(find(foot_template<0))=0;
hand_template(find(hand_template<0))=0;

%divide by sum - in this way, the sum of the template = 1
foot_template=foot_template/sum(foot_template);
hand_template=hand_template/sum(hand_template);


maxTopoTemplateHand = ceil(max(hand_template)*100)/100; % to get ylim for each subplot
minTopoTemplateHand = min(hand_template);

maxTopoTemplateFoot = ceil(max(foot_template)*100)/100; % to get ylim for each subplot
minTopoTemplateFoot = min(foot_template);

%% Apply and test the templates
%extract values for the other conditions
%hand_noise
a=find(stim_locations==1); %hand
b=find(stim_type==2); %noise
cB=intersect(a,b); 
hand_noise_base=mean(lwdata(cB).bl_data(:,:,results.sig_idx_base),3);
hand_noise_oddball=mean(lwdata(cB).bl_data(:,:,results.sig_idx_oddball),3);
%foot_noise
a=find(stim_locations==2); %foot
b=find(stim_type==2); %noise
cB=intersect(a,b); 
foot_noise_base=mean(lwdata(cB).bl_data(:,:,results.sig_idx_base),3);
foot_noise_oddball=mean(lwdata(cB).bl_data(:,:,results.sig_idx_oddball),3);

%hand_noise X hand and foot templates
hand_noise_base_hand_template=[];
hand_noise_base_foot_template=[];
foot_noise_base_hand_template=[];
foot_noise_base_foot_template=[];
hand_noise_oddball_hand_template=[];
hand_noise_oddball_foot_template=[];
foot_noise_oddball_hand_template=[];
foot_noise_oddball_foot_template=[];

for i=1:size(hand_noise_base,1);
    hand_noise_base_hand_template(i)=sum(hand_noise_base(i,:)'.*hand_template);
    hand_noise_base_foot_template(i)=sum(hand_noise_base(i,:)'.*foot_template);
    foot_noise_base_hand_template(i)=sum(foot_noise_base(i,:)'.*hand_template);
    foot_noise_base_foot_template(i)=sum(foot_noise_base(i,:)'.*foot_template);
end;

for i=1:size(hand_noise_oddball,1);
    hand_noise_oddball_hand_template(i)=sum(hand_noise_oddball(i,:)'.*hand_template);
    hand_noise_oddball_foot_template(i)=sum(hand_noise_oddball(i,:)'.*foot_template);
    foot_noise_oddball_hand_template(i)=sum(foot_noise_oddball(i,:)'.*hand_template);
    foot_noise_oddball_foot_template(i)=sum(foot_noise_oddball(i,:)'.*foot_template);
end;

%% plot templates and test
subplot(2,2,1);
hold on
for i=1:length(hand_noise_base_hand_template)
    a(1)=hand_noise_base_hand_template(i);
    a(2)=hand_noise_base_foot_template(i);
    b(1)=1;
    b(2)=2;
    plot(b,a,'Color',[24, 116, 205]/255);
    set(gca,'Position',[0.1 0.55 0.3 0.4]);
    set(gca,'XTick',[1 2]);
    set(gca,'XTickLabel',{'Homotopic', 'Heterotopic'});
%     set(gca,XAxis
end;
a(1)=mean(hand_noise_base_hand_template);
a(2)=mean(hand_noise_base_foot_template);
plot(b,a,'k','LineWidth',4);
set(gca,'Position',[0.1 0.55 0.3 0.4]);
title(' Hand base response','FontSize',18,'FontWeight','bold');
% text(2.1,0.03,'Foot template','FontSize',18);
% text(0.1,0.03,'Hand template','FontSize',18);
xlim([0 3]);
set(gca,'FontSize',14)
hold off

%%
subplot(2,2,2);
hold on
for i=1:length(hand_noise_oddball_hand_template)
    a(1)=hand_noise_oddball_hand_template(i);
    a(2)=hand_noise_oddball_foot_template(i);
    b(1)=1;
    b(2)=2;
    plot(b,a,'Color',[24, 116, 205]/255);
    %set(gca,'Position',[0.13 0.056 0.3 0.3])
    set(gca,'Position',[0.45 0.55 0.3, 0.4]);
    set(gca,'XTick',[1 2]);
    set(gca,'XTickLabel',{'Homotopic', 'Heterotopic'});
end;
a(1)=mean(hand_noise_oddball_hand_template);
a(2)=mean(hand_noise_oddball_foot_template);
plot(b,a,'k','LineWidth',4);
%set(gca,'Position',[0.13 0.056 0.3 0.3])
set(gca,'Position',[0.45 0.55 0.3, 0.4]);
title('Hand oddball response','FontSize',18,'FontWeight','bold');
% text(2.05,0.004,'Foot template','FontSize',18);
% text(0.1,0.004,'Hand template','FontSize',18);
xlim([0 3]);
set(gca,'FontSize',14)
%hold off
%%
subplot(2,2,3);
hold on
for i=1:length(hand_noise_base_hand_template)
    a(1)=foot_noise_base_foot_template(i);
    a(2)=foot_noise_base_hand_template(i);
    %a(2)=foot_noise_base_foot_template(i);
    b(1)=1;
    b(2)=2;
    plot(b,a,'Color',[24, 116, 205]/255);
    set(gca,'Position',[0.55 0.46 0.3 0.3]);
  set(gca,'Position',[0.1 0.056 0.3 0.4])
    set(gca,'XTick',[1 2]);
    set(gca,'XTickLabel',{'Homotopic', 'Heterotopic'});
end;
a(1)=mean(foot_noise_base_foot_template);
a(2)=mean(foot_noise_base_hand_template);
%a(2)=mean(foot_noise_base_foot_template);
plot(b,a,'k','LineWidth',4);
%set(gca,'Position',[0.57 0.46 0.3, 0.3]);
set(gca,'Position',[0.1 0.056 0.3 0.4])
title('Foot base response','FontSize',18,'FontWeight','bold');
% text(2.1,0.035/2,'Foot template','FontSize',18);
% text(0.1,0.035/2,'Hand template','FontSize',18);
xlim([0 3]);
set(gca,'FontSize',14)
hold off
%%
subplot(2,2,4);
hold on
for i=1:length(hand_noise_oddball_hand_template)
    a(1)=foot_noise_oddball_foot_template(i);
    a(2)=foot_noise_oddball_hand_template(i);
    %a(2)=foot_noise_oddball_foot_template(i);
    b(1)=1;
    b(2)=2;
    plot(b,a,'Color',[24, 116, 205]/255);
    set(gca,'Position',[0.45 0.056 0.3 0.4]);
    set(gca,'XTick',[1 2]);
    set(gca,'XTickLabel',{'Homotopic', 'Heterotopic'});
end;
a(1)=mean(foot_noise_oddball_foot_template);
a(2)=mean(foot_noise_oddball_hand_template);
% a(2)=mean(foot_noise_oddball_foot_template);
plot(b,a,'k','LineWidth',4);
set(gca,'Position',[0.45 0.056 0.3 0.4]);
title('Foot oddball response','FontSize',18,'FontWeight','bold');
% text(2.05,1*10^-3,'Foot template','FontSize',18);
% text(0.1,1*10^-3,'Hand template','FontSize',18);
xlim([0 3]);
set(gca,'FontSize',14)
hold off

%% plot templates
%figure;
%subplot(1,2,1);
axes('Position',[0.75,0.65,0.2,0.2]);
CLW_topoplot_vector(lwdata(1).header,hand_template,'maplimits',[minTopoTemplateHand maxTopoTemplateHand]);
title('Hand Template','FontSize',14,'FontWeight','Bold')
cH = colorbar;
    set(cH, 'ylim', [minTopoTemplateHand maxTopoTemplateHand])
    set(cH,'FontSize',18);
    ylabel(cH,'\muV','rotation',0, 'Position',[7 minTopoTemplateHand 1], 'FontSize',14);
%subplot(1,2,2);
hold on
axes('Position',[0.75,0.156,0.2,0.2]);
CLW_topoplot_vector(lwdata(1).header,foot_template,'maplimits',[minTopoTemplateFoot maxTopoTemplateFoot]);
title('Foot Template','FontSize',14,'FontWeight','Bold')
cF = colorbar;
    set(cF, 'ylim', [minTopoTemplateFoot maxTopoTemplateFoot])
    set(cF,'FontSize',18);
    ylabel(cF,'\muV','rotation',0, 'Position',[7 minTopoTemplateFoot 1], 'FontSize',14);
    
 axes('visible','off')
xlabel([]);
ylabel('Amplitude (\muV)','FontSize',24,'FontWeight','bold','Position',[-0.1,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gcf, 'Color', 'White');


%% ttest
[results.hand_noise_base_result,results.hand_noise_base_pvalue]=ttest(hand_noise_base_hand_template,hand_noise_base_foot_template, 'tail', 'right');
[results.foot_noise_base_result,results.foot_noise_base_pvalue]=ttest(foot_noise_base_hand_template,foot_noise_base_foot_template, 'tail', 'left');
[results.hand_noise_oddball_result,results.hand_noise_oddball_pvalue]=ttest(hand_noise_oddball_hand_template,hand_noise_oddball_foot_template, 'tail', 'right');
[results.foot_noise_oddball_result,results.foot_noise_oddball_pvalue]=ttest(foot_noise_oddball_hand_template,foot_noise_oddball_foot_template, 'tail', 'left');

%% Create ''mismatch tables''
for i = 1:64
    chanlocs= lwdata(1).header.chanlocs;
end
for i = 1:64
chanlocs_all{i} = chanlocs(i).labels;
end
chanlocs_all = chanlocs_all';
chanlocs_all = cell2table(chanlocs_all);

hand_template_table = array2table(hand_template);
hand_template_table=[hand_template_table chanlocs_all];
hand_table_sorted = sortrows(hand_template_table,'hand_template','descend');
hand_table_sorted.Properties.VariableNames ={'template','chanlocs'};


foot_template_table = array2table(foot_template);
foot_template_table=[foot_template_table chanlocs_all];
foot_table_sorted = sortrows(foot_template_table,'foot_template','descend');
foot_table_sorted.Properties.VariableNames ={'template','chanlocs'};

%ismember(hand_table_sorted.chanlocs,foot_table_sorted.chanlocs,'chanlocs')

hand_table_MM = hand_table_sorted(1:7,:);
hand_template_MM =table2array(hand_table_MM(:,1));
foot_table_MM = foot_table_sorted(1:7,:);
foot_template_MM =table2array(foot_table_MM(:,1));
% apply mismatch templates
chanlocs_all.Properties.VariableNames ={'chanlocs'};
[channelsH, iah ibh]= intersect(hand_table_sorted(:,2),chanlocs_all);
[channelsF, iaf ibf]= intersect(foot_table_sorted(:,2),chanlocs_all);

for i=1:size(hand_noise_base,1);
    hand_noise_base_hand_template_MM(i)=sum(hand_noise_base(i,ibh)'.*hand_template_MM(iah));
    hand_noise_base_foot_template_MM(i)=sum(hand_noise_base(i,ibf)'.*foot_template_MM(iaf));
    foot_noise_base_hand_template_MM(i)=sum(foot_noise_base(i,ibh)'.*hand_template_MM(iah));
    foot_noise_base_foot_template_MM(i)=sum(foot_noise_base(i,ibf)'.*foot_template_MM(iaf));
end;

for i=1:size(hand_noise_oddball,1);
    hand_noise_oddball_hand_template_MM(i)=sum(hand_noise_oddball(i,ibh)'.*hand_template_MM(iah));
    hand_noise_oddball_foot_template_MM(i)=sum(hand_noise_oddball(i,ibf)'.*foot_template_MM(iaf));
    foot_noise_oddball_hand_template_MM(i)=sum(foot_noise_oddball(i,ibh)'.*hand_template_MM(iah));
    foot_noise_oddball_foot_template_MM(i)=sum(foot_noise_oddball(i,ibf)'.*foot_template_MM(iaf));
end;


%% New plot to compare FFT of EEG and of signal
%spectrotemporal hand
sig_base_freqs_3=results.base_freq(find(lwdata(3).ttest_base_freq));
sig_oddball_freqs_3=results.oddball_freq(find(lwdata(3).ttest_oddball_freq));
sig_base_amps_3=lwdata(3).chanavg_avg_bl(results.base_freq_dx(find(lwdata(3).ttest_base_freq)));
sig_oddball_amps_3= lwdata(3).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(3).ttest_oddball_freq)));

subplot(2,1,1);
hold on
plot(tpx,lwdata(3).chanavg_avg_bl,'k-');
xlim([0.25 40]);
scatter(sig_base_freqs_3,sig_base_amps_3,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
stem(tpx(results.base_freq_dx(find(lwdata(3).ttest_base_freq))),...
        lwdata(3).chanavg_avg_bl(results.base_freq_dx(find(lwdata(3).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
    
scatter(sig_oddball_freqs_3,sig_oddball_amps_3,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
stem(tpx(results.oddball_freq_dx(find(lwdata(3).ttest_oddball_freq))),...
        lwdata(3).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(3).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    ylim([0,0.06]);
title('Spectrotemporal contrast (Hand)','FontWeight','bold','FontSize',20);
%spectrotemporal foot
sig_base_freqs_4=results.base_freq(find(lwdata(4).ttest_base_freq));
sig_oddball_freqs_4=results.oddball_freq(find(lwdata(4).ttest_oddball_freq));
sig_base_amps_4=lwdata(4).chanavg_avg_bl(results.base_freq_dx(find(lwdata(4).ttest_base_freq)));
sig_oddball_amps_4= lwdata(4).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(4).ttest_oddball_freq)));

subplot(2,1,2);
plot(tpx,lwdata(4).chanavg_avg_bl,'k-');
box off;
xlim([0.25 40]);
hold on
scatter(sig_base_freqs_4,sig_base_amps_4,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
hold on;
stem(tpx(results.base_freq_dx(find(lwdata(4).ttest_base_freq))),...
        lwdata(4).chanavg_avg_bl(results.base_freq_dx(find(lwdata(4).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
hold on    
scatter(sig_oddball_freqs_4,sig_oddball_amps_4,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
stem(tpx(results.oddball_freq_dx(find(lwdata(4).ttest_oddball_freq))),...
        lwdata(4).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(4).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
ylim([0,0.06]);
title('Spectrotemporal contrast (Foot)','FontWeight','bold','FontSize',20);
    
%topoplots foot
axes('Position',[0.7 0.25 0.14 0.14]);   

title('Base response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(4).header,lwdata(4).mean_base_signal,'maplimits',[minTopoAmpBase(4) maxTopoAmpBase(4)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpBase(4) maxTopoAmpBase(4)])
    set(cB,'FontSize',14);
    
axes('Position',[0.85 0.25 0.14 0.14]);    
title('Oddball response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(4).header,lwdata(4).mean_oddball_signal,'maplimits',[minTopoAmpOddball(4) maxTopoAmpOddball(4)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpOddball(4) maxTopoAmpOddball(4)])
    set(cB,'FontSize',14);
% signal
    
% subplot(3,1,3);
% plot(tpxSignal,signalFFT.data,'k-');
% box off;
% hold on;
% xlim([0.25 40]);
% stem(sig_base_freqs_3,signalFFT.data(base_freq_dx),...
%         'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
%     hold on
%     
% scatter(base_freq,signalFFT.data(base_freq_dx(find(lwdata(3).ttest_base_freq))),'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
% 
% hold on
% stem(oddball_freq,signalFFT.data(oddball_freq_dx),...
%         'Color',[238, 44, 44]/255,'marker','none','LineWidth',2.5);
% scatter(oddball_freq,signalFFT.data(oddball_freq_dx),'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
% box off
% ylim([0,0.06]);
% title('Spectrum of averaged spectrotemporal sequences envelopes','FontWeight','bold','FontSize',20);
% topoplots
axes('Position',[0.7 0.65 0.14 0.14]);   

title('Base response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(3).header,lwdata(3).mean_base_signal,'maplimits',[minTopoAmpBase(3) maxTopoAmpBase(3)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpBase(3) maxTopoAmpBase(3)])
    set(cB,'FontSize',14);
    
axes('Position',[0.85 0.65 0.14 0.14]);    
title('Oddball response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(3).header,lwdata(3).mean_oddball_signal,'maplimits',[minTopoAmpOddball(3) maxTopoAmpOddball(3)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpOddball(3) maxTopoAmpOddball(3)])
    set(cB,'FontSize',14);
hold on         
axes('visible','off')
xlabel('Frequency (Hz)','FontSize',16','FontWeight','bold','Position',[0.5,-0.04,1])
ylabel('Amplitude (\muV)','FontSize',16,'FontWeight','bold','Position',[-0.04,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gcf, 'Color', 'White');

%%
%% New plot to compare FFT of EEG and of signal
%spectrotemporal hand
sig_base_freqs_1=results.base_freq(find(lwdata(1).ttest_base_freq));
sig_oddball_freqs_1=results.oddball_freq(find(lwdata(1).ttest_oddball_freq));
sig_base_amps_1=lwdata(1).chanavg_avg_bl(results.base_freq_dx(find(lwdata(1).ttest_base_freq)));
sig_oddball_amps_1= lwdata(1).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(1).ttest_oddball_freq)));
figure;
subplot(2,1,1);
hold on
plot(tpx,lwdata(1).chanavg_avg_bl,'k-');
xlim([0.25 40]);
scatter(sig_base_freqs_1,sig_base_amps_1,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
stem(tpx(results.base_freq_dx(find(lwdata(1).ttest_base_freq))),...
        lwdata(1).chanavg_avg_bl(results.base_freq_dx(find(lwdata(1).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
    
scatter(sig_oddball_freqs_1,sig_oddball_amps_1,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
stem(tpx(results.oddball_freq_dx(find(lwdata(1).ttest_oddball_freq))),...
        lwdata(1).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(1).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    ylim([0,0.13]);
title('Frequency/intensity contrast (Hand)','FontWeight','bold','FontSize',20);
%spectrotemporal foot
sig_base_freqs_2=results.base_freq(find(lwdata(2).ttest_base_freq));
sig_oddball_freqs_2=results.oddball_freq(find(lwdata(2).ttest_oddball_freq));
sig_base_amps_2=lwdata(2).chanavg_avg_bl(results.base_freq_dx(find(lwdata(2).ttest_base_freq)));
sig_oddball_amps_2= lwdata(2).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(2).ttest_oddball_freq)));

subplot(2,1,2);
plot(tpx,lwdata(2).chanavg_avg_bl,'k-');
box off;
xlim([0.25 40]);
hold on
scatter(sig_base_freqs_2,sig_base_amps_2,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
hold on;
stem(tpx(results.base_freq_dx(find(lwdata(2).ttest_base_freq))),...
        lwdata(2).chanavg_avg_bl(results.base_freq_dx(find(lwdata(2).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
hold on    
scatter(sig_oddball_freqs_2,sig_oddball_amps_2,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
stem(tpx(results.oddball_freq_dx(find(lwdata(2).ttest_oddball_freq))),...
        lwdata(2).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(2).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
ylim([0,0.06]);
title('Frequncy/intensity contrast (Foot)','FontWeight','bold','FontSize',20);
    
%topoplots foot
axes('Position',[0.7 0.25 0.14 0.14]);   

title('Base response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(2).header,lwdata(2).mean_base_signal,'maplimits',[minTopoAmpBase(2) maxTopoAmpBase(2)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpBase(2) maxTopoAmpBase(2)])
    set(cB,'FontSize',14);
    
axes('Position',[0.85 0.25 0.14 0.14]);    
title('Oddball response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(4).header,lwdata(2).mean_oddball_signal,'maplimits',[minTopoAmpOddball(2) maxTopoAmpOddball(2)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpOddball(2) maxTopoAmpOddball(2)])
    set(cB,'FontSize',14);

axes('Position',[0.7 0.65 0.14 0.14]);   

title('Base response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(1).header,lwdata(1).mean_base_signal,'maplimits',[minTopoAmpBase(1) maxTopoAmpBase(1)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpBase(1) maxTopoAmpBase(1)])
    set(cB,'FontSize',14);
    
axes('Position',[0.85 0.65 0.14 0.14]);    
title('Oddball response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(1).header,lwdata(1).mean_oddball_signal,'maplimits',[minTopoAmpOddball(1) maxTopoAmpOddball(1)]);
    cB = colorbar;
    %set(cB,'Location','WestOutside');
    set(cB, 'ylim', [minTopoAmpOddball(1) maxTopoAmpOddball(1)])
    set(cB,'FontSize',14);
hold on         
axes('visible','off')
xlabel('Frequency (Hz)','FontSize',16','FontWeight','bold','Position',[0.5,-0.04,1])
ylabel('Amplitude (\muV)','FontSize',16,'FontWeight','bold','Position',[-0.04,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gcf, 'Color', 'White');
%%

base_signal_sign= signalFFT.data(base_freq_dx);%(find(lwdata(3).ttest_base_freq)));
oddball_signal_sign=signalFFT.data(oddball_freq_dx);%(find(lwdata(3).ttest_base_freq)));
ratio_signal = mean(oddball_signal_sign)/mean(base_signal_sign);
ratio_signal = mean(oddball_signal_sign)/mean(base_signal_sign);
ratio_hand= mean(sig_oddball_amps_3)/mean(sig_base_amps_3);
ratio_foot= mean(sig_oddball_amps_4)/mean(sig_base_amps_4);

%%
filename = 'analysis_group_level_wn1_included_1809.mat';
save (filename, '-v7.3');


%%
base_data_anova_incl_wn2 = [hand_noise_base_hand_template' hand_noise_base_foot_template' foot_noise_base_foot_template' foot_noise_base_hand_template'];
oddball_data_anova_incl_wn2 = [hand_noise_oddball_hand_template' hand_noise_oddball_foot_template' foot_noise_oddball_foot_template' foot_noise_oddball_hand_template'];
base_data_anova_incl_wn2=array2table(base_data_anova_incl_wn2);
oddball_data_anova_incl_wn2=array2table(oddball_data_anova_incl_wn2);

writetable(base_data_anova_incl_wn2,'base_data_anova_wn2_incl_1509.csv')
writetable(oddball_data_anova_incl_wn2,'oddball_data_anova_wn2_incl_1509.csv')
