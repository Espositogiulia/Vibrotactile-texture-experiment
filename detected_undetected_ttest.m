%% Load the merged datasets
cd('C:\data\gesposito\xpVibroAnalysisFinal\IndividualTrialsAnalysisWN\preprocessing');
[lwdata(1).header,lwdata(1).data]=CLW_load('fft merged_epochs DETECTED');
[lwdata(2).header,lwdata(2).data]=CLW_load('fft merged_epochs NOTDETECTED');


for i=1:2;
    lwdata(i).data=squeeze(lwdata(i).data);
end;
%% assign conditions

lwdata(1).detected=1;
lwdata(2).notdetected=2;


%% Frequency axis bins>frequency
tpx=1:size(lwdata(1).data,3);
tpx=(tpx-1)*lwdata(1).header.xstep;

%% Baseline correction and Z-score computation
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:2;
    tp=lwdata(i).data;
    lwdata(i).bl_data=tp;
    % lwdata(i).Z_data=tp;
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
                %lwdata(i).Z_data(epochpos,channelpos,dx)=(lwdata(i).data(epochpos,channelpos,dx)-mean_tpy)/SD_tpy;
            end;
        end;
    end;
end;

%% group-level average amplitude + average across scalp channels

for i=1:2
    lwdata(i).avg=squeeze(mean(lwdata(i).data,1));
    lwdata(i).avg_bl=squeeze(mean(lwdata(i).bl_data,1));
    % lwdata(i).avg_Z=squeeze(mean(lwdata(i).Z_data,1));
    lwdata(i).chanavg_avg=squeeze(mean(lwdata(i).avg,1));
    lwdata(i).chanavg_avg_bl=squeeze(mean(lwdata(i).avg_bl,1));
    lwdata(i).chanavg_bl=squeeze(mean(lwdata(i).bl_data,2));
    % lwdata(i).chanavg_Z=squeeze(mean(lwdata(i).Z_data,2));
    lwdata(i).chanavg=squeeze(mean(lwdata(i).data,2));
end;

%% zscore of average amplitude
start_bin=3; %to exclude the first two neighbouring bins
end_bin=start_bin+12; %to include the 12 neighbouring bins

for i=1:2;
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
            %  lwdata(i).Z_avg(channelpos,dx)=(lwdata(i).avg(channelpos,dx)-mean_tpy)/SD_tpy;
        end;
    end;
    %  lwdata(i).chanvg_Z_avg=mean(lwdata(i).Z_avg,1);
    %chanavg_avg
    for dx=1:num_bins;
        dxbins=[dx-(end_bin-1):dx-start_bin dx+start_bin:dx+(end_bin-1)];
        dxbins(dxbins<1)=[];
        dxbins(dxbins>num_bins)=[];
        tpy=lwdata(i).chanavg_avg(dxbins);
        mean_tpy=mean(tpy);
        SD_tpy=std(tpy);
        lwdata(i).bl_chanavg_avg(dx)=lwdata(i).chanavg_avg(dx)-mean_tpy;
        % lwdata(i).Z_chanavg_avg(dx)=(lwdata(i).chanavg_avg(dx)-mean_tpy)/SD_tpy;
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

for i=1:2;
    %baseline_corrected amplitude at all base frequency harmonics, averaged across all channels
    %base frequencies
    tp_base=lwdata(i).chanavg_bl(:,base_freq_dx);
    mean_tp_base=mean(tp_base,2);
    %[a,b]=ttest(mean_tp_base);
    [a,b]=ttest(mean_tp_base,0,'Tail','right');
    lwdata(i).ttest.base_chan_avg_bl=tp_base;
    lwdata(i).ttest.mean_base_chan_avg_bl=mean_tp_base;
    lwdata(i).ttest.grand_avg_mean_base_chan_avg_bl = mean( lwdata(i).ttest.mean_base_chan_avg_bl);
    lwdata(i).ttest.stdev_mean_base_chan_avg_bl = std(lwdata(i).ttest.mean_base_chan_avg_bl);
    lwdata(i).ttest.base_significant=a;
    lwdata(i).ttest.base_pvalue=b;
    %oddball frequencies
    tp_oddball=lwdata(i).chanavg_bl(:,oddball_freq_dx);
    mean_tp_oddball=mean(tp_oddball,2);
    [a,b]=ttest(mean_tp_oddball,0,'Tail','right');
    %     [a,b]=ttest(mean_tp_oddball);
    lwdata(i).ttest.oddball_chan_avg_bl=tp_oddball;
    lwdata(i).ttest.mean_oddball_chan_avg_bl=mean_tp_oddball;
    lwdata(i).ttest.grand_avg_mean_oddball_chan_avg_bl = mean( lwdata(i).ttest.mean_oddball_chan_avg_bl);
    lwdata(i).ttest.stdev_mean_oddball_chan_avg_bl = std( lwdata(i).ttest.mean_oddball_chan_avg_bl);
    lwdata(i).ttest.oddball_significant=a;
    lwdata(i).ttest.oddball_pvalue=b;
    
end;
%% ttest at each frequency
for i=1:2;
    lwdata(i).ttest_base_freq=[];
    lwdata(i).ttest_oddball_freq=[];
    tp=lwdata(i).chanavg_bl;
    %loop through base_frequencies
    for j=1:length(results.base_freq_dx);
        tp2=tp(:,results.base_freq_dx(j));
        lwdata(i).ttest_base_freq(j)=ttest(tp2,0,'Tail','right');
        % lwdata(i).ttest_base_freq(j)=ttest(tp2);
    end;
    %loop through oddball_frequencies
    for j=1:length(results.oddball_freq_dx);
        tp2=tp(:,results.oddball_freq_dx(j));
        lwdata(i).ttest_oddball_freq(j)=ttest(tp2,0,'Tail','right');
        % lwdata(i).ttest_oddball_freq(j)=ttest(tp2);
    end;
end;

%% plot spectra and topogtaphies (averaged across subjects and channels with results of ttest - only for white noise condition
titles= {'Detected' 'Not detected'};
for i = 1:length(lwdata)
    maxYamp(i) = ceil(max(lwdata(i).chanavg_avg_bl(1,12:1921))*100)/100; % to get ylim for each subplot
end

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

% load idx_base and idx_oddball (original harmonics)
%idx_oddball=sort(unique(idx_oddball));
results.idx_base=idx_base;
results.idx_oddball=idx_oddball;

%average (baseline-corrected) signal across the relevant harmonics
for i=1:length(lwdata);
    lwdata(i).base_signal=lwdata(i).bl_avg(:,idx_base);
    lwdata(i).base_signal_ALL= lwdata(i).bl_data(:,:,idx_base);% to keep all participants separate
    lwdata(i).oddball_signal=lwdata(i).bl_avg(:,idx_oddball);
    lwdata(i).oddball_signal_ALL= lwdata(i).bl_data(:,:,idx_oddball);
    lwdata(i).mean_base_signal=mean(lwdata(i).base_signal,2);
    lwdata(i).mean_base_signal_ALL=mean(lwdata(i).base_signal_ALL,3);% to keep all participants separate
    lwdata(i).mean_oddball_signal=mean(lwdata(i).oddball_signal,2);
    lwdata(i).mean_oddball_signal_ALL=mean(lwdata(i).oddball_signal_ALL,3);
end;

%set limits for colorbars
for i = 1:length(lwdata)
    maxTopoAmpBase(i) = ceil(max(lwdata(i).mean_base_signal)*100)/100; % to get ylim for each subplot
end

for i = 1:length(lwdata)
    minTopoAmpBase(i) = min(lwdata(i).mean_base_signal); % to get ylim for each subplot
end

for i = 1:length(lwdata)
    maxTopoAmpOddball(i) = ceil(max(lwdata(i).mean_oddball_signal)*100)/100; % to get ylim for each subplot
end

for i = 1:length(lwdata)
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
positionsTopoOddball= repmat(positionsTopoOddball,2,1);
%%make the plot
figure;
for i=1:length(lwdata);
    axes('Position',positions(i,:));
    h(i)= plot(tpx,lwdata(i).chanavg_avg_bl,'k-');
    %text(20,maxYamp(i),titles{i},'FontSize',18,'FontWeight','bold','HorizontalAlignment','center');
    box off
    xlim([0 40]);
    hold on
    %highlight significant base frequencies
    sig_base_freqs=results.base_freq(find(lwdata(i).ttest_base_freq));
    sig_base_amps=lwdata(i).chanavg_avg_bl(results.base_freq_dx(find(lwdata(i).ttest_base_freq)));
    % non significant base frequencies
   base_amps = lwdata(i).chanavg_avg_bl(results.base_freq_dx);
    
   scatter(base_freq,base_amps,'MarkerEdgeColor',[24, 116, 205]/255,'LineWidth',2.5);
    
    scatter(sig_base_freqs,sig_base_amps,'MarkerEdgeColor',[24, 116, 205]/255,'MarkerFaceColor',[24, 116, 205]/255);
    stem(tpx(results.base_freq_dx(find(lwdata(i).ttest_base_freq))),...
        lwdata(i).chanavg_avg_bl(results.base_freq_dx(find(lwdata(i).ttest_base_freq))),...
        'Color',[24, 116, 205]/255,'marker','none','LineWidth',2.5);
    
    %highlight significant oddball frequencies
    sig_oddball_freqs=results.oddball_freq(find(lwdata(i).ttest_oddball_freq));
    sig_oddball_amps=lwdata(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq)));
    % non significant frequecies
    oddball_amps = lwdata(i).chanavg_avg_bl(results.oddball_freq_dx);
    
    scatter(oddball_freq,oddball_amps,'MarkerEdgeColor',[238, 44, 44]/255,'LineWidth',2.5);
    scatter(sig_oddball_freqs,sig_oddball_amps,'MarkerEdgeColor',[238, 44, 44]/255,'MarkerFaceColor',[238, 44, 44]/255);
    
    stem(tpx(results.oddball_freq_dx(find(results.oddball_freq))),...
        lwdata(i).chanavg_avg_bl(results.oddball_freq_dx(find(results.oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    %
    stem(tpx(results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))),...
        lwdata(i).chanavg_avg_bl(results.oddball_freq_dx(find(lwdata(i).ttest_oddball_freq))),...
        'Color', [238, 44, 44]/255,'marker','none','LineWidth',2.5);
    
    %scatter(sig_oddball_freqs,sig_oddball_amps,'r','filled');
    %hold off;
    
    set(h(i),'LineWidth',1);
    set(gca,'xlim',[0.25 40]);
    set(gca,'ylim',[-0.01 maxYamp(i)+0.01]);
    set(gca,'TickLength',[0 0]);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',16)
    set(gca,'YTick',[0 maxYamp(i)]);
    if i <4
        set(gca, 'XTickLabel', [])
    end
    % Add title to each subplot
    title(titles{i}, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end;


axes('visible','off')
xlabel('Frequency (Hz)','FontSize',16','FontWeight','bold','Position',[0.5,-0.09,1])
ylabel('Amplitude (\muV)','FontSize',16,'FontWeight','bold','Position',[-0.1,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gcf, 'Color', 'White');
%hold off;

%%Topographical plot (baseline-corrected data)

for i=1:length(lwdata);
    %subplot(2,2,i);
    axes('Position',positionsTopoBase(i,:));
    % if i == 1 || i == 2
    title('Base response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    CLW_topoplot_vector(lwdata(i).header,lwdata(i).mean_base_signal,'maplimits',[minTopoAmpBase(i) maxTopoAmpBase(i)]);
    cB = colorbar;
    set(cB, 'ylim', [minTopoAmpBase(i) maxTopoAmpBase(i)])
    ylabel(cB,'\muV','rotation',0, 'Position',[7 minTopoAmpBase(i) 1], 'FontSize',14);
    set(cB,'FontSize',14);
end;
hold on;
%topoplot of oddball response
%figure;
for i=1:length(lwdata);
    axes('Position',positionsTopoOddball(i,:));
    title('Oddball response','fontweight','bold','fontsize',14,'Position',[0 0.713 5]);
    
    CLW_topoplot_vector(lwdata(i).header,lwdata(i).mean_oddball_signal,'maplimits',[minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
    cO = colorbar;
    set(cO, 'ylim', [minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
    %set(cO, 'YTickLabel', [minTopoAmpOddball(i) maxTopoAmpOddball(i)]);
    ylabel(cO,'\muV','rotation',0, 'Position',[7 minTopoAmpOddball(i) 1],'FontSize',14);
    set(cO,'FontSize',14);
end;

%% test detected vs undetected
for i = 1:2
    lwdata(i).chanavg_mean_oddball_signal_ALL = mean(lwdata(i).mean_oddball_signal_ALL,2);
end

[hPvsNP pPvsNP]= ttest(lwdata(1).chanavg_mean_oddball_signal_ALL, lwdata(2).chanavg_mean_oddball_signal_ALL([3 5 6 7 8 9 ...
    10 11 12 13 14 15 16])); % remove ppts


filename = 'analysis_perceived_not_perceived_2405.mat';
save (filename, '-v7.3');

