%% Load the merged datasets
cd('C:\data\gesposito\xpVibroAnalysisFinal\FullAnalysis_Wn2Included');
[lwdata(1).header,lwdata(1).data]=CLW_load('merged_epochs_hand_s4');
[lwdata(2).header,lwdata(2).data]=CLW_load('merged_epochs_foot_s4');
[lwdata(3).header,lwdata(3).data]=CLW_load('merged_epochs_hand_s2');
[lwdata(4).header,lwdata(4).data]=CLW_load('merged_epochs_foot_s2');

for i=1:4;
    lwdata(i).data=squeeze(lwdata(i).data);
end;

%% assign conditions
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

%% Baseline correction
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:4;
    tp=lwdata(i).data;
    lwdata(i).bl_data=tp;
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
            end;
        end;
    end;
end;

%% group-level average amplitude + average across scalp channels
for i=1:4
    lwdata(i).avg=squeeze(mean(lwdata(i).data,1)); % average across subjects (dim 1: epochs) - FFT data
    lwdata(i).avg_bl=squeeze(mean(lwdata(i).bl_data,1)); % average across subjects (dim 1: epochs) BL FFT data
    lwdata(i).chanavg_avg=squeeze(mean(lwdata(i).avg,1)); % average across channels (on grand averaged FFT data)
    lwdata(i).chanavg_avg_bl=squeeze(mean(lwdata(i).avg_bl,1)); % average across channels (on grand averaged BL FFT data)
    lwdata(i).chanavg_bl=squeeze(mean(lwdata(i).bl_data,2)); % average across channels (dim 2: channels) across BL FFT data
    lwdata(i).chanavg=squeeze(mean(lwdata(i).data,2)); % average across channels (dim 2: channels) - FFT data
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

%% find harmonics that were significant (ttest) in at least one condition
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

%% plot spectra and topogtaphies (averaged across subjects and channels with results of ttest
%average (baseline-corrected) signal across the relevant harmonics
for i=3:4;
    lwdata(i).base_signal=lwdata(i).avg_bl(:,idx_base);
    lwdata(i).oddball_signal=lwdata(i).avg_bl(:,idx_oddball);
    lwdata(i).mean_base_signal=mean(lwdata(i).base_signal,2);
    lwdata(i).mean_oddball_signal=mean(lwdata(i).oddball_signal,2);
end;

titles= {'Frequency/intensity contrast (hand)' 'Frequency/intensity contrast (foot)' 'Spectrotemporal contrast (hand)' 'Spectrotemporal contrast (foot)'};
for i = 3:4
    maxYamp(i) = ceil(max(lwdata(i).chanavg_avg_bl(1,12:1921))*100)/100; % to get ylim for each subplot
end

%  to position each subplot
positions = zeros(2,4); % to position each subplot 
positions(1,2) = 0.55;
for i = 1:2
    positions(i,1) = 0.1;
    positions(i,3) = 0.85;
    positions(i,4) = 0.4;
    if i==2
    positions(i,2) = positions(i-1,2)-0.46;
    end
end

positions= repmat(positions, 2,1);

%%topoplot parameters
idx_base=[];
idx_oddball=[];
for i=3:4;
    idx_base=[idx_base results.base_freq_dx(find(lwdata(i).ttest_base_freq))];
    idx_oddball=[idx_oddball results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))];
end;
idx_base=sort(unique(idx_base));
idx_oddball=sort(unique(idx_oddball));

results.idx_base=idx_base;
results.idx_oddball=idx_oddball;

% set limits for colorbars
for i = 3:4
    maxTopoAmpBase(i) = ceil(max(lwdata(i).mean_base_signal)*100)/100; % to get ylim for each subplot
end

for i = 3:4
    minTopoAmpBase(i) = min(lwdata(i).mean_base_signal)*100/100; % to get ylim for each subplot
end

for i = 3:4
    maxTopoAmpOddball(i) = ceil(max(lwdata(i).mean_oddball_signal)*100)/100; % to get ylim for each subplot
end

for i = 3:4
    minTopoAmpOddball(i) = min(lwdata(i).mean_oddball_signal); % to get ylim for each subplot
end

%topoplot of base response
positionsTopoBase = zeros(2,4); % to position each subplot 
positionsTopoBase(1,2) = 0.65;
for i = 1:2
    positionsTopoBase(i,1) = 0.7;
    positionsTopoBase(i,3) = 0.14;
    positionsTopoBase(i,4) = 0.14;
    if i>1
    positionsTopoBase(i,2) = positionsTopoBase(i-1,2)-0.46;
    end
end
positionsTopoBase = repmat(positionsTopoBase,2,1);

positionsTopoOddball = zeros(2,4); % to position each subplot 
positionsTopoOddball(1,2) = 0.65;
for i = 1:2
    positionsTopoOddball(i,1) = 0.85;
    positionsTopoOddball(i,3) = 0.14;
    positionsTopoOddball(i,4) = 0.14;
    if i==2
    positionsTopoOddball(i,2) = positionsTopoOddball(i-1,2)-0.46;
    end
end
positionsTopoOddball= repmat(positionsTopoOddball,2,1);%%make the plot
figure;
for i=3:4;
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

% Topographical plot (baseline-corrected data)
for i=3:4;
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
for i=3:4;
    axes('Position',positionsTopoOddball(i,:));
    title('Oddball response','fontweight','bold','fontsize',14);
    
    CLW_topoplot_vector(lwdata(i).header,lwdata(i).mean_oddball_signal,'maplimits',[minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
    cO = colorbar;
    set(cO, 'ylim', [minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
    % set(cO,'Location','WestOutside');
    set(cO,'FontSize',14);
    ylabel(cO,'\muV','rotation',0, 'Position',[7 (minTopoAmpOddball(i)+minTopoAmpOddball(i)*0.1) 1]);
end;

%% average hand spectrotemporal responses across clusters of harmonics to
% look at topography
lwdata(3).bl_avg_cluster1=mean(lwdata(3).avg_bl(:,oddball_freq_dx([1:4])),2);
lwdata(3).bl_avg_cluster2=mean(lwdata(3).avg_bl(:,oddball_freq_dx([5:8])),2);
lwdata(3).bl_avg_cluster3=mean(lwdata(3).avg_bl(:,oddball_freq_dx([9:12])),2);
lwdata(3).bl_avg_cluster4=mean(lwdata(3).avg_bl(:,oddball_freq_dx([13:16])),2);
lwdata(3).bl_avg_cluster5=mean(lwdata(3).avg_bl(:,oddball_freq_dx([17:20])),2);

harm_clusters_hand = [lwdata(3).bl_avg_cluster1 lwdata(3).bl_avg_cluster2 lwdata(3).bl_avg_cluster3 lwdata(3).bl_avg_cluster4 lwdata(3).bl_avg_cluster5];
min_clust_hand = min(harm_clusters_hand);
max_clust_hand = max(harm_clusters_hand);

% Define subplot properties
subplotWidth = 0.18;  % Adjust the width of each subplot
subplotHeight = 0.7;  % Adjust the height of each subplot
subplotSpacing = 0.02;  % Adjust the spacing between subplots

% Loop through the clusters and create subplots - hand
for i = 1:5
    % Calculate the position of the current subplot
    subplotX = (i - 1) * (subplotWidth + subplotSpacing) + 0.05;
    
    % Create subplot
    subplot('Position', [subplotX, 0.1, subplotWidth, subplotHeight]);
    
    % Plot the topoplot for the current cluster
    CLW_topoplot_vector(lwdata(3).header, lwdata(3).(['bl_avg_cluster' num2str(i)]), 'maplimits', [min_clust_hand(i) max_clust_hand(i)]);
    
    % Add colorbar
    cB = colorbar;
    set(cB, 'ylim', [min_clust_hand(i) max_clust_hand(i)]);
    set(cB, 'FontSize', 16);
      ylabel(cB,'\muV','rotation',0);
    
    % Clear colorbar variable
    clear cB;
end

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
xlim([0 3]);
set(gca,'FontSize',14)
hold off

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
xlim([0 3]);
set(gca,'FontSize',14)
%hold off

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
xlim([0 3]);
set(gca,'FontSize',14)
hold off

% plot templates
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


%% ttest topographical test template
[results.hand_noise_base_result,results.hand_noise_base_pvalue]=ttest(hand_noise_base_hand_template,hand_noise_base_foot_template, 'tail', 'right');
[results.foot_noise_base_result,results.foot_noise_base_pvalue]=ttest(foot_noise_base_hand_template,foot_noise_base_foot_template, 'tail', 'left');
[results.hand_noise_oddball_result,results.hand_noise_oddball_pvalue]=ttest(hand_noise_oddball_hand_template,hand_noise_oddball_foot_template, 'tail', 'right');
[results.foot_noise_oddball_result,results.foot_noise_oddball_pvalue]=ttest(foot_noise_oddball_hand_template,foot_noise_oddball_foot_template, 'tail', 'left');

%% save workspace
filename = 'analysis_group_level_wn2_included_240524.mat';
save (filename, '-v7.3');

%% save csv files for topographical test RM Anova in JASP
base_data_anova_incl_wn2 = [hand_noise_base_hand_template' hand_noise_base_foot_template' foot_noise_base_foot_template' foot_noise_base_hand_template'];
oddball_data_anova_incl_wn2 = [hand_noise_oddball_hand_template' hand_noise_oddball_foot_template' foot_noise_oddball_foot_template' foot_noise_oddball_hand_template'];
base_data_anova_incl_wn2=array2table(base_data_anova_incl_wn2);
oddball_data_anova_incl_wn2=array2table(oddball_data_anova_incl_wn2);

writetable(base_data_anova_incl_wn2,'base_data_anova_wn2_incl_3011.csv')
writetable(oddball_data_anova_incl_wn2,'oddball_data_anova_wn2_incl_3011.csv')
