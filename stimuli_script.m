% stimuli script 0302
clear;
clc; 
% cd('C:\Users\Nocions\Desktop\vibroExpGiuliaE\new');
rng('shuffle');
%% Sequence condition 1 = large contrast vibrations
% define frequencies and characteristics of vibrations
%freq1 = 47;
freq1 = 200;
%freq2 = 197; 
freq2 = 300; 
amp = 1;
samFreq = 44100;
durSound = 0.125; % duration of A and B
durSeq = durSound*5; % AAAB duration
rTime = 0.015; %ramp time
t = [0:1/samFreq:durSound-1/samFreq];
lt = length(t);
tr = [0:1/samFreq:rTime-1/samFreq];
lr = size(tr,2);
rampUp = ( (cos(2*pi*tr/rTime/2+pi)+1)/2 ).^2;
rampDown = ( (cos(2*pi*tr/rTime/2)+1)/2 ).^2;

%% generate sounds for cond 1
% v1
v1 = amp*sin(2*pi*freq1*t); % generate sine for 47Hz tone
v1(1:lr) = rampUp.*v1(1:lr); % ramp up
v1(lt-lr+1:lt) = rampDown.*v1(lt-lr+1:lt); % ramp down
normV1 = v1/rms(v1); % rms norm v1
%v2  
v2 = amp*sin(2*pi*freq2*t); % generate sine for 197Hz tone
v2(1:lr) = rampUp.*v2(1:lr); % ramp up
v2(lt-lr+1:lt) = rampDown.*v2(lt-lr+1:lt); % ramp down
normV2 = v2/rms(v2); % rms norm v2

%% put together as AAAAB for cond1
cond1smallContrastSeq = [normV1 normV1 normV1 normV1 normV2];

%% concatenate 64 AAAAB sequences to have final 40 sequence for cond 1
durSeq = durSound*5; %duration of a single AAAAB sequence
durTotSeq = durSound*5*64; %duration for 64 repetitions
nseq = durTotSeq/durSeq;%number of repetitions

% cond1smallContrastSeqTot = [];
% for i = 1:nseq
%     cond1smallContrastSeqTot = [cond1smallContrastSeqTot cond1smallContrastSeq];
% end
% % normalise amplitude between -1 and 1 to center around 0
% cond1smallContrastSeqTot=cond1smallContrastSeqTot/max(abs(minmax(cond1smallContrastSeqTot)));
% save('condition1smallContr.mat','cond1smallContrastSeqTot');


%% put together as AAAAB cond2
cond2smallContrastSeq = [normV2 normV2 normV2 normV2 normV1];

%% concatenate 64 AAAAB sequences to have final 40 sequence cond2
durSeq = durSound*5; %duration of a single AAAAB sequence
durTotSeq = durSound*5*64; %duration for 64 repetitions
nseq = durTotSeq/durSeq;%number of repetitions

cond2smallContrastSeqTot = [];
for i = 1:nseq
    cond2smallContrastSeqTot = [cond2smallContrastSeq cond2smallContrastSeqTot];
end
% normalise amplitude between -1 and 1 to center around 0
cond2smallContrastSeqTot=cond2smallContrastSeqTot/max(abs(minmax(cond2smallContrastSeqTot)));
cond2smallContrastSeqTot = cond2smallContrastSeqTot *0.50 ; % reduce amplitude of sequence to match wn 
% need to check how it is with broadband wn
%save('condition2smallContr.mat', 'cond2smallContrastSeqTot');

%filename = 'freqContrast.wav';
audiowrite(filename, cond2smallContrastSeqTot, 44100);
[y,Fs] = audioread(filename) ;

clear filename;

%% create an array of 20 x length of tone A for white noise 
wnTableA = zeros(8, length(v1)); 
%% generate 0 white noise sequences and put in array white noise 
for i = 1:size(wnTableA) %loop through each row
wnTableA(i,:) = rand(1,length(v1)); %this way you get WN amplitude range of 0 to 1
wnTableA(i,:) = 2*(wnTableA(i,:)-0.5); % this way you get amplitude -0.99 to 0.99 
end

% ranomise A to get B

[M,N] = size(wnTableA);
% Preserve the row indices
rowIndex = repmat((1:M)',[1 N]);
% Get randomized column indices by sorting a second random array
[~,randomizedColIndex] = sort(rand(M,N),2);
% Need to use linear indexing to create B
newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
wnTableB = wnTableA(newLinearIndex);
clear i;


%% design and apply bandpass filter to each white noise white noise 

% % julien's way using filtfilt
 cutoff_threshold = 40;
% 
% high_cutoff_threshold = 400;
% %high_cutoff_threshold = 1000;
% 
 d = designfilt('highpassiir','FilterOrder',4,'PassbandFrequency',40, 'SampleRate',samFreq);
% 
 filtWnTableA = zeros(size(wnTableA));
 for k = 1:size(filtWnTableA)
 filtWnTableA(k,:) = filtfilt(d,wnTableA(k,:));
 end
 clear k;
% 
 filtWnTableB = zeros(size(wnTableB));
 for k = 1:size(filtWnTableB)
 filtWnTableB(k,:) = filtfilt(d,wnTableB(k,:));
 end
 clear k;
 
 %% lowpass
 d = designfilt('lowpassiir','FilterOrder',4,'PassbandFrequency',200, 'SampleRate',samFreq);
% 

 for i = 1:size(filtWnTableA)
 filtWnTableA(i,:) = filtfilt(d,filtWnTableA(i,:));
 end

% 

 for i = 1:size(filtWnTableB)
 filtWnTableB(i,:) = filtfilt(d,filtWnTableB(i,:));
 end
 clear k;

 %% ramp up and down each row (ie wn sequence in the array) white noise 
for j = 1:size(filtWnTableA)
filtWnTableA(j,1:lr) = rampUp.*filtWnTableA(j,1:lr); % ramp up wn1
filtWnTableA(j,lt-lr+1:lt) = rampDown.*filtWnTableA(j,lt-lr+1:lt); % ramp down wn1
end 
clear j;

for j = 1:size(filtWnTableB)
filtWnTableB(j,1:lr) = rampUp.*filtWnTableB(j,1:lr); % ramp up wn1
filtWnTableB(j,lt-lr+1:lt) = rampDown.*filtWnTableB(j,lt-lr+1:lt); % ramp down wn1
end 
clear j;

%% rms normalisation white noise 
normWnTableA = zeros(size(filtWnTableA));
for l = 1:size(normWnTableA)
normWnTableA(l,:) = filtWnTableA(l,:)/rms(filtWnTableA(l,:));
end
clear l;

normWnTableB = zeros(size(filtWnTableB));
for l = 1:size(normWnTableB)
normWnTableB(l,:) = filtWnTableB(l,:)/rms(filtWnTableB(l,:));
end
clear l;


wnA = normWnTableA;
wnB = normWnTableB;

%% now create an AAAAB pattern 10 times
wn1 = [repmat(wnA(1,:),1,4) wnB(1,:)];%repmat creates a 1 by length of a single row of wnA repeated 4 times
wn2 = [repmat(wnA(2,:),1,4) wnB(2,:)];
wn3 = [repmat(wnA(3,:),1,4) wnB(3,:)];
wn4 = [repmat(wnA(4,:),1,4) wnB(4,:)];
wn5 = [repmat(wnA(5,:),1,4) wnB(5,:)];
wn6 = [repmat(wnA(6,:),1,4) wnB(6,:)];
wn7 = [repmat(wnA(7,:),1,4) wnB(7,:)];
wn8 = [repmat(wnA(8,:),1,4) wnB(8,:)];

totWn = [wn1;wn2;wn3;wn4;wn5;wn6;wn7;wn8]; %make an array with all AAAAB sequences

% next I need to make each 40 seconds
lengthTot = [zeros(1,length(cond1smallContrastSeq)) ]; % 
wnTrial = zeros(8,length(lengthTot)); % make an array for all white noise seq (40 sec)

%% concatenate 64 AAAAB sequences to have final 40 + 5 sec pause sequence for cond 3
durTotSeq = durSound*5*64; %duration for 64 repetitions
nseq = durTotSeq/durSeq;%number of repetitions

% this should give me an array where each row is a complete 64xAAAAB(white
% noise) lasting 40 seconds 
wnTrial = [];
for i = 1:nseq
    wnTrial = [wnTrial totWn(:,:)]; % in this way I should concatenate everything                                                        
end
clear i;

%% normalise amplotude between -1 and 1 white noise 
maxMinWnTable = zeros(size(wnTrial));
for m = 1:size(maxMinWnTable)
maxMinWnTable(m,:)=wnTrial(m,:)/max(abs(minmax(wnTrial(m,:))));
end
clear m;
%% make individual sequences = can do this in for loop
        
wnTrial1 = maxMinWnTable(1,:);

save('wnTrial1.mat','wnTrial1');
wnTrial2 = maxMinWnTable(2,:);

save('wnTrial2.mat','wnTrial2');
wnTrial3 = maxMinWnTable(3,:);

save('wnTrial3.mat','wnTrial3');
wnTrial4 = maxMinWnTable(4,:);

save('wnTrial4.mat','wnTrial4');
wnTrial5 = maxMinWnTable(5,:);

save('wnTrial5.mat','wnTrial5');
wnTrial6 = maxMinWnTable(6,:);

save('wnTrial6.mat','wnTrial6');
wnTrial7 = maxMinWnTable(7,:);

save('wnTrial7.mat','wnTrial7');    
wnTrial8 = maxMinWnTable(8,:);

save('wnTrial8.mat','wnTrial8');

%% check envelope of white noise
Fs = 44100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(wnTrial1);             % Length of signal
t = (0:L-1)/L;        % Time vector

hilTranSignal = abs(hilbert(wnTrial1));

NFFT = 2^nextpow2(length(hilTranSignal));

S1=fft(hilTranSignal,NFFT)/L;
S2=fft(abs(hilTranSignal),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);


figure;plot(f,abs(S2(1:NFFT/2+1)));
xlim([0 30]);
ylim([0 0.2]);


%% spectrogram
% 
freqSp = cond2smallContrastSeqTot(1:55120);
wnSp = wnTrial1(1:55120);
% this is to bring it to the signal analyser app to get the time/frequency analysis
SFreq = timetable(seconds(0:length(freqSp)-1)'/Fs,freqSp');
Swn = timetable(seconds(0:length(wnSp)-1)'/Fs,wnSp');
% 
% % Parameters
% timeLimits = seconds([0 1.249864]); % seconds
% frequencyLimits = [0 500]; % Hz
% leakage = 0.85;
% timeResolution = 25; % seconds
% overlapPercent = 50;
% 
% pspectrum([cond1smallContrastSeq cond1smallContrastSeq],Fs,'spectrogram', 'Leakage',0.85,'OverlapPercent',0, ...
%     'MinThreshold',-10,'FrequencyLimits',[0, 500]), 'TimeLimits', [0, 1.2], 'TimeResolution', 0.025;

filename = 'wn1.wav';
audiowrite(filename, wnTrial1, 44100);
[y,Fs] = audioread(filename) ;
