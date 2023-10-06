% to run this, you need to load the foot and hand individual wn results.
clear;
clc;
%% Load the fft of band pass filtered sequences
[signalFFT(1).header,signalFFT(1).data]=CLW_load('fft but hilbert wnTrial1');
[signalFFT(2).header,signalFFT(2).data]=CLW_load('fft but hilbert wnTrial2');
[signalFFT(3).header,signalFFT(3).data]=CLW_load('fft but hilbert wnTrial3');
[signalFFT(4).header,signalFFT(4).data]=CLW_load('fft but hilbert wnTrial4');
[signalFFT(5).header,signalFFT(5).data]=CLW_load('fft but hilbert wnTrial5');
[signalFFT(6).header,signalFFT(6).data]=CLW_load('fft but hilbert wnTrial6');
[signalFFT(7).header,signalFFT(7).data]=CLW_load('fft but hilbert wnTrial7');
[signalFFT(8).header,signalFFT(8).data]=CLW_load('fft but hilbert wnTrial8');

for i=1:8;
    signalFFT(i).data=squeeze(signalFFT(i).data);
end;

for i=1:8;
signalFFT(i).oddball_signal=signalFFT(i).data(idx_oddball);
end;


for i=1:8;
signalFFT(i).base_signal=signalFFT(i).data(idx_base);
end;

for i = 1:8
    signalFFT(i).mean_oddball= mean(signalFFT(i).oddball_signal);
end
wn1 = signalFFT(1).mean_oddball;
wn2 = signalFFT(2).mean_oddball;
wn3 = signalFFT(3).mean_oddball;
wn4 = signalFFT(4).mean_oddball;
wn5 = signalFFT(5).mean_oddball;
wn6 = signalFFT(6).mean_oddball;
wn7 = signalFFT(7).mean_oddball;
wn8 = signalFFT(8).mean_oddball;

for i = 1:8
    signalFFT(i).mean_base= mean(signalFFT(i).base_signal);
end

for i = 1:8
signalFFT(i).ratio= signalFFT(i).mean_oddball/signalFFT(i).mean_base;
end

wn1_ratio = signalFFT(1).ratio;
wn2_ratio = signalFFT(2).ratio;
wn3_ratio = signalFFT(3).ratio;
wn4_ratio = signalFFT(4).ratio;
wn5_ratio = signalFFT(5).ratio;
wn6_ratio = signalFFT(6).ratio;
wn7_ratio = signalFFT(7).ratio;
wn8_ratio = signalFFT(8).ratio;

%% ratio hand
    for i = 1:8
mean_base_hand(i).mean_base_hand =  lwdataHand(i).mean_base_signal_individual;
    end
wn1hand_base = mean_base_hand(1).mean_base_hand;
wn2hand_base = mean_base_hand(2).mean_base_hand;
wn3hand_base = mean_base_hand(3).mean_base_hand;
wn4hand_base = mean_base_hand(4).mean_base_hand;
wn5hand_base = mean_base_hand(5).mean_base_hand;
wn6hand_base = mean_base_hand(6).mean_base_hand;
wn7hand_base = mean_base_hand(7).mean_base_hand;
wn8hand_base = mean_base_hand(8).mean_base_hand;
base_allHand = [wn1hand_base wn2hand_base wn3hand_base wn4hand_base wn5hand_base wn6hand_base...
    wn7hand_base wn8hand_base];

    for i = 1:8
mean_oddball_hand(i).mean_oddball_hand =  lwdataHand(i).mean_oddball_signal_individual;
    end
    
wn1hand_oddball = mean_oddball_hand(1).mean_oddball_hand;
wn2hand_oddball = mean_oddball_hand(2).mean_oddball_hand;
wn3hand_oddball = mean_oddball_hand(3).mean_oddball_hand;
wn4hand_oddball = mean_oddball_hand(4).mean_oddball_hand;
wn5hand_oddball = mean_oddball_hand(5).mean_oddball_hand;
wn6hand_oddball = mean_oddball_hand(6).mean_oddball_hand;
wn7hand_oddball = mean_oddball_hand(7).mean_oddball_hand;
wn8hand_oddball = mean_oddball_hand(8).mean_oddball_hand;
oddball_allHand = [wn1hand_oddball wn2hand_oddball wn3hand_oddball wn4hand_oddball wn5hand_oddball wn6hand_oddball...
    wn7hand_oddball wn8hand_oddball];

oddball_hand_excl_Wn2= [wn1hand_oddball wn3hand_oddball wn4hand_oddball wn5hand_oddball ...
    wn6hand_oddball wn7hand_oddball wn8hand_oddball];
ratio_hand = oddball_allHand./base_allHand;
ratio_hand_excl_wn2 = ratio_hand(:,[1 3:8]);
%% ratio foot
    for i = 1:8
mean_base_foot(i).mean_base_foot =  lwdataFoot(i).mean_base_signal_individual;
    end
wn1foot_base = mean_base_foot(1).mean_base_foot;
wn2foot_base = mean_base_foot(2).mean_base_foot;
wn3foot_base = mean_base_foot(3).mean_base_foot;
wn4foot_base = mean_base_foot(4).mean_base_foot;
wn5foot_base = mean_base_foot(5).mean_base_foot;
wn6foot_base = mean_base_foot(6).mean_base_foot;
wn7foot_base = mean_base_foot(7).mean_base_foot;
wn8foot_base = mean_base_foot(8).mean_base_foot;
base_allFoot = [wn1foot_base wn2foot_base wn3foot_base wn4foot_base wn5foot_base wn6foot_base...
    wn7foot_base wn8foot_base];

    for i = 1:8
mean_oddball_foot(i).mean_oddball_foot =  lwdataFoot(i).mean_oddball_signal_individual;
    end
    
wn1foot_oddball = mean_oddball_foot(1).mean_oddball_foot;
wn2foot_oddball = mean_oddball_foot(2).mean_oddball_foot;
wn3foot_oddball = mean_oddball_foot(3).mean_oddball_foot;
wn4foot_oddball = mean_oddball_foot(4).mean_oddball_foot;
wn5foot_oddball = mean_oddball_foot(5).mean_oddball_foot;
wn6foot_oddball = mean_oddball_foot(6).mean_oddball_foot;
wn7foot_oddball = mean_oddball_foot(7).mean_oddball_foot;
wn8foot_oddball = mean_oddball_foot(8).mean_oddball_foot;
oddball_allFoot = [wn1foot_oddball wn2foot_oddball wn3foot_oddball wn4foot_oddball wn5foot_oddball wn6foot_oddball...
    wn7foot_oddball wn8foot_oddball];


ratio_foot = oddball_allFoot./base_allFoot;
ratio_foot_excl_wn2 = ratio_foot(:,[1 3:8]);
%%
oddball_wnAll = [wn1 wn2 wn3 wn4 wn5 wn6 wn7 wn8];
oddball_wnAll = repmat(oddball_wnAll,17,1);
oddball_wnAll_excl_wn2 = [wn1 wn3 wn4 wn5 wn6 wn7 wn8];
oddball_wnAll_excl_wn2 = repmat(oddball_wnAll_excl_wn2,17,1);
oddball_wnAll_ratio = [wn1_ratio wn2_ratio wn3_ratio wn4_ratio wn5_ratio wn6_ratio wn7_ratio wn8_ratio];
oddball_wnAll_excl_wn2_ratio = [wn1_ratio wn3_ratio wn4_ratio wn5_ratio wn6_ratio wn7_ratio wn8_ratio];
oddball_wnAll_excl_wn2_ratio = repmat(oddball_wnAll_excl_wn2_ratio,17,1);

oddball_wnAll_ratio = repmat(oddball_wnAll_ratio,17,1);

[Rhandwn2_FFT_ratio,Phandwn2_FFT_ratio]= corrcoef(oddball_wnAll_ratio(:,:),ratio_hand(:,:));
[Rhand_excl_wn2_FFT_ratio,Phand_excl_wn2_FFT_ratio]= corrcoef(oddball_wnAll_ratio(:,[1 3:8]),ratio_hand(:,[1 3:8]));

% foot
[Rfootwn2_FFT_ratio,Pfootwn2_FFT_ratio]= corrcoef(oddball_wnAll_ratio(:,:),ratio_foot(:,:));
[Rfoot_excl_wn2_FFT_ratio,Pfoot_excl_wn2_FFT_ratio]= corrcoef(oddball_wnAll_ratio(:,[1 3:8]),ratio_foot(:,[1 3:8]));

%% fit 1 hand wn2

[xData1, yData1] = prepareCurveData( oddball_wnAll_ratio, ratio_hand );

%[xData1, yData1] = prepareCurveData( oddball_wnAll, oddball_allHand );

% Set up fittype and options.
ft1 = fittype( 'poly1' );
% Fit model to data.
[fitresult1, gof1] = fit( xData1, yData1, ft1 );

%% Fit: 'untitled fit 2'. hand no wn2
[xData2, yData2] = prepareCurveData( oddball_wnAll_excl_wn2_ratio, ratio_hand_excl_wn2 );
%[xData2, yData2] = prepareCurveData( oddball_wnAll_excl_wn2, oddball_hand_excl_Wn2);

% Set up fittype and options.
ft2 = fittype( 'poly1' );
% Fit model to data.
[fitresult2, gof2] = fit( xData2, yData2, ft2 );

%% fit 3 (foot wn2
[xData3, yData3] = prepareCurveData( oddball_wnAll_ratio, ratio_foot );
%[xData3, yData3] = prepareCurveData( oddball_wnAll, oddball_allFoot);

% Set up fittype and options.
ft3 = fittype( 'poly1' );

% Fit model to data.
[fitresult3, gof3] = fit( xData3, yData3, ft3 );

%% Fit: 'untitled fit 4'. foot no wn2
[xData4, yData4] = prepareCurveData( oddball_wnAll_excl_wn2_ratio, ratio_foot_excl_wn2 );
%[xData4, yData4] = prepareCurveData( oddball_wnAll_excl_wn2,oddball_foot_excl_Wn2);

% Set up fittype and options.
ft4 = fittype( 'poly1' );
ft4_onset = fittype( 'poly1' );
% Fit model to data.
[fitresult4, gof4] = fit( xData4, yData4, ft4 );

%% plot
fig = figure;
f.WindowState = 'maximized';
axes('Position',[0.2,0.57,0.6,0.4]);
x1= plot(fitresult1, xData1, yData1,'ko');
box off
set(x1,'LineWidth',1.5);
set(x1,'MarkerSize',3);
set(x1, {'MarkerFaceColor'}, get(x1,'Color'));

hold on

% x2 = plot(fitresult2,xData2,yData2);
% set(x2,'LineWidth',1.5);
% set(x2,'Color',[24, 116, 205]/255);
% set(x2,'Marker','none');
 xlabel('')
 ylabel('');
 set (gca,'XTick', [0.15  0.2 0.25 0.3  0.35 0.4]);
% set (gca,'Tickdir', 'out');
% a = get(gca,'XTickLabel');  
 set(gca,'XTickLabel',a,'fontsize',18)
 set(gca,'XTickLabelMode','auto')
 legend('off')

axes('Position',[0.2,0.1,0.6,0.4]);
%axes('Position',[0.2,0.57,0.6,0.4]);
x3= plot(fitresult3, xData3, yData3,'ko');
box off
set(x3,'LineWidth',1.5);
set(x3,'MarkerSize',3);
set(x3, {'MarkerFaceColor'}, get(x3,'Color'));


hold on

% x4 = plot(fitresult4,xData4,yData4);
% set(x4,'LineWidth',1.5);
% set(x4,'Color',[24, 116, 205]/255);
% set(x4,'Marker','none');
 xlabel('')
 ylabel('');
set (gca,'XTick', [0.15  0.2 0.25 0.3  0.35 0.4]);
% set (gca,'Tickdir', 'out');
% a = get(gca,'XTickLabel');  
 set(gca,'XTickLabel',a,'fontsize',18)
 set(gca,'XTickLabelMode','auto')
 legend('off')


% title('(b)','FontSize',14);
legend('off')

hold on
axes('visible','off')
xlabel('Oddball/Base ratio (signal) (AU)','FontSize',24','FontWeight','bold','Position',[0.5,-0.08,1])
ylabel('Oddball/Base ratio (EEG) (AU)','FontSize',24,'FontWeight','bold','Position',[0.03,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')
