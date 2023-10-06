%% Fit: 'untitled fit 1'. HAND
[xData1, yData1] = prepareCurveData( Dsum, oddball_allHand );
[xData1_onset, yData1_onset] = prepareCurveData( Dsum_onset, oddball_allHand );

% Set up fittype and options.
ft1 = fittype( 'poly1' );
ft1_onset = fittype( 'poly1' );

% Fit model to data.
[fitresult1, gof1] = fit( xData1, yData1, ft1 );
[fitresult1_onset, gof1_onset] = fit( xData1_onset, yData1_onset, ft1_onset );

%% Fit: 'untitled fit 1'.
[xData2, yData2] = prepareCurveData( Dsum_excl_wn2, oddball_hand_excl_Wn2);
%[xData2_onset, yData2_onset] = prepareCurveData( Dsum_onset_excl_wn2, oddball_hand_excl_Wn2);


% Set up fittype and options.
ft2 = fittype( 'poly1' );
ft2_onset = fittype( 'poly1' );


% Fit model to data.
[fitresult2, gof2] = fit( xData2, yData2, ft2 );
%[fitresult2_onset, gof2_onset] = fit( xData2_onset, yData2_onset, ft2_onset );

%% Fit: 'untitled fit 1'. FOOT
[xData3, yData3] = prepareCurveData( Dsum, oddball_allFoot );
[xData3_onset, yData3_onset] = prepareCurveData( Dsum_onset, oddball_allFoot );

% Set up fittype and options.
ft3 = fittype( 'poly1' );
ft3_onset = fittype( 'poly1' );

% Fit model to data.
[fitresult3, gof3] = fit( xData3, yData3, ft3 );
[fitresult3_onset, gof3_onset] = fit( xData3_onset, yData3_onset, ft3_onset );
%% Fit: 'untitled fit 1'.
[xData4, yData4] = prepareCurveData( Dsum_excl_wn2, oddball_foot_excl_Wn2);
[xData4_onset, yData4_onset] = prepareCurveData( Dsum_onset_excl_wn2, oddball_foot_excl_Wn2);


% Set up fittype and options.
ft4 = fittype( 'poly1' );
ft4_onset = fittype( 'poly1' );

% Fit model to data.
[fitresult4, gof4] = fit( xData4, yData4, ft4 );
[fitresult4_onset, gof4_onset] = fit( xData4_onset, yData4_onset, ft4_onset );

%% figure Dsum
fig = figure;
f.WindowState = 'maximized';
axes('Position',[0.2,0.57,0.6,0.4]);
x1= plot(fitresult1, xData1, yData1,'ko');
box off
set(x1,'LineWidth',1.5);
set(x1,'MarkerSize',3);
set(x1, {'MarkerFaceColor'}, get(x1,'Color'));
% text(Dsum(1,1),-0.005, '1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,2),-0.005, '2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,3),-0.007, '3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,4),-0.005, '4', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,5),-0.005, '5', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,6),-0.005, '6', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,7),-0.005, '7', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,8),-0.005, '8', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

hold on

x2 = plot(fitresult2,xData2,yData2);
set(x2,'LineWidth',1.5);
set(x2,'Color',[24, 116, 205]/255);
set(x2,'Marker','none');
xlabel('')
ylabel('');
set (gca,'YTick', [0  0.01  0.02  0.03]);
set (gca,'Tickdir', 'out');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',18)
set(gca,'XTickLabelMode','auto')
legend('off')

%%
%h(1)= subplot(2,2,1);
axes('Position',[0.2,0.1,0.6,0.4]);
%axes('Position',[0.2,0.57,0.6,0.4]);
x3= plot(fitresult3, xData3, yData3,'ko');
box off
set(x3,'LineWidth',1.5);
set(x3,'MarkerSize',3);
set(x3, {'MarkerFaceColor'}, get(x3,'Color'));

% text(Dsum(1,1),-0.007, '1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,2),-0.007, '2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,3),-0.009, '3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,4),-0.007, '4', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,5),-0.007, '5', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,6),-0.007, '6', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,7),-0.007, '7', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum(1,8),-0.007, '8', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold on

x4 = plot(fitresult4,xData4,yData4);
set(x4,'LineWidth',1.5);
set(x4,'Color',[24, 116, 205]/255);
set(x4,'Marker','none');
xlabel('')
ylabel('');
set (gca,'YTick', [-0.005 0 0.005 0.01]);
set (gca,'Tickdir', 'out');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',18)
set(gca,'XTickLabelMode','auto')
legend('off')


% title('(b)','FontSize',14);
legend('off')

hold on
axes('visible','off')
xlabel('Sequence envelope dissimilarity (AU)','FontSize',24','FontWeight','bold','Position',[0.5,-0.08,1])
ylabel('Amplitude (\muV)','FontSize',24,'FontWeight','bold','Position',[0.03,0.5,1])
set(gca, 'visible', 'off')
set(findall(gca, 'type', 'text'), 'visible', 'on')

% %% figure Dsum_onset
% fig_onser = figure;
% f.WindowState = 'maximized';
% axes('Position',[0.2,0.57,0.6,0.4]);
% x1_onset= plot(fitresult1_onset, xData1_onset, yData1_onset,'ko');
% box off
% set(x1_onset,'LineWidth',1.5);
% set(x1_onset,'MarkerSize',3);
% set(x1_onset, {'MarkerFaceColor'}, get(x1_onset,'Color'));
% text(Dsum_onset(1,1),-0.007, '1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,2),-0.007, '2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,3),-0.007, '3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,4),-0.007, '4', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,5),-0.007, '5', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,6),-0.007, '6', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,7),-0.007, '7', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,8),-0.007, '8', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% hold on
% 
% x2_onset = plot(fitresult2_onset,xData2_onset,yData2_onset);
% set(x2_onset,'LineWidth',1.5);
% set(x2_onset,'Color',[24, 116, 205]/255);
% set(x2_onset,'Marker','none');
% xlabel('')
% ylabel('');
% set (gca,'YTick', [0  0.01  0.02  0.03]);
% set (gca,'Tickdir', 'out');
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',18)
% set(gca,'XTickLabelMode','auto')
% legend('off')
% 
% %%
% %h(1)= subplot(2,2,1);
% axes('Position',[0.2,0.1,0.6,0.4]);
% %axes('Position',[0.2,0.57,0.6,0.4]);
% x3_onset= plot(fitresult3_onset, xData3_onset, yData3_onset,'ko');
% box off
% set(x3_onset,'LineWidth',1.5);
% set(x3_onset,'MarkerSize',3);
% set(x3_onset, {'MarkerFaceColor'}, get(x3_onset,'Color'));
% text(Dsum_onset(1,1),-0.007, '1', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,2),-0.007, '2', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,3),-0.007, '3', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,4),-0.007, '4', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,5),-0.007, '5', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,6),-0.007, '6', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,7),-0.007, '7', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(Dsum_onset(1,8),-0.007, '8', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% hold on
% 
% x4_onset = plot(fitresult4_onset,xData4_onset,yData4_onset);
% set(x4_onset,'LineWidth',1.5);
% set(x4_onset,'Color',[24, 116, 205]/255);
% set(x4_onset,'Marker','none');
% xlabel('')
% ylabel('');
% set (gca,'YTick', [-0.005 0 0.005 0.01]);
% set (gca,'Tickdir', 'out');
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',18)
% set(gca,'XTickLabelMode','auto')
% legend('off')
% 
% 
% % title('(b)','FontSize',14);
% legend('off')
% 
% hold on
% axes('visible','off')
% xlabel('Sequence envelope dissimilarity (AU)','FontSize',24','FontWeight','bold','Position',[0.5,-0.08,1])
% ylabel('Amplitude (\muV)','FontSize',24,'FontWeight','bold','Position',[0.03,0.5,1])
% set(gca, 'visible', 'off')
% set(findall(gca, 'type', 'text'), 'visible', 'on')
