 % generate matrices to present blocks and white noise trials
clear;
clc;
% cd('C:\Users\Nocions\Desktop\vibroExpGiuliaE\new');

% cond1smallContrInt= load('condition1smallContr.mat');% all condition matrices are saved in the stimuli script
% cond1smallContrInt = cond1smallContrInt.cond1smallContrastSeqTot;

cond2smallContrInt = load('condition2SmallContr.mat');
cond2smallContrInt = cond2smallContrInt.cond2smallContrastSeqTot;

SR=44100;
%% triggers
%triggers for frequency contrast
trig1freq = zeros(1,length(cond2smallContrInt)); % make it the same length as the sound sequence
trig2freq = zeros(1,length(cond2smallContrInt));
trig1freq(1:10*SR/1000)=5;
cond2smallContrTrig = [cond2smallContrInt' trig1freq' trig2freq'];
%save('cond2smallContrTrig.mat','cond2smallContrTrig');

%triggers for white noise contrast
trig1wn = zeros(1,length(cond2smallContrInt));
trig2wn = zeros(1,length(cond2smallContrInt));
trig2wn(1:10*SR/1000) = 5;

%% white noise trials
cond3wn1 = load('wnTrial1.mat');% i did this one by one because extracting the whole struct ...
%                                 I saved in the other script for some reason put ...
%                                  everything into a single row...
cond3wn1 = cond3wn1.wnTrial1;

 cond3wn2 = load('wnTrial2.mat');
 cond3wn2 = cond3wn2.wnTrial2;

cond3wn3 = load('wnTrial3.mat');
cond3wn3 = cond3wn3.wnTrial3;

cond3wn4 = load('wnTrial4.mat');
cond3wn4 = cond3wn4.wnTrial4;

cond3wn5 = load('wnTrial5.mat');
cond3wn5 = cond3wn5.wnTrial5;

 cond3wn6 = load('wnTrial6.mat');
 cond3wn6 = cond3wn6.wnTrial6;

cond3wn7 = load('wnTrial7.mat');
cond3wn7 = cond3wn7.wnTrial7;

cond3wn8 = load('wnTrial8.mat');
cond3wn8 = cond3wn8.wnTrial8;


s.whiteNoise = {cond3wn1 cond3wn2 cond3wn3 cond3wn4 ...
    cond3wn5 cond3wn6 cond3wn7 cond3wn8}; %structure containing all white noise trials
whiteNoiseTrials = s.whiteNoise; % whiteNoiseTrials is a 1x10 cell array. Each cell array is a vector containing the white noise sequence for a single trial
%save('whiteNoiseTrialsNew.mat', 'whiteNoiseTrials');% this saves it as a struct to be loaded to the presentation script

blocks = [1 2]; % one block per condition

% %% this is for adding eeg triggers when using the daq

 wnTrig1 = [cond3wn1' trig1wn' trig2wn']; % create array where the sequence and the triggers are on 2 columns
 wnTrig2 = [cond3wn2' trig1wn' trig2wn'];
 wnTrig3= [cond3wn3' trig1wn' trig2wn'];
 wnTrig4= [cond3wn4' trig1wn' trig2wn'];
 wnTrig5= [cond3wn5' trig1wn' trig2wn'];
 wnTrig6= [cond3wn6' trig1wn' trig2wn'];
 wnTrig7= [cond3wn7' trig1wn' trig2wn'];
 wnTrig8= [cond3wn8' trig1wn' trig2wn'];


 s.trials = {cond2smallContrTrig cond2smallContrTrig cond2smallContrTrig cond2smallContrTrig cond2smallContrTrig ...
     cond2smallContrTrig cond2smallContrTrig cond2smallContrTrig wnTrig1 wnTrig2 wnTrig3 wnTrig4 ...
     wnTrig5 wnTrig6 wnTrig7 wnTrig8};
 
 trials = s.trials;
 %save('trials.mat','trials');
%  s.whiteNoiseTrig = {wnTrig1 wnTrig3 wnTrig4 wnTrig5 ...
%      wnTrig7 wnTrig8 wnTrig9 wnTrig10}; %structure containing all white noise trials with eeg triggers
% whiteNoiseTrialsTrig = s.whiteNoiseTrig; % extract cell array from struct

%% save rand white noise matrix
%block1 hand
randTrialsHand2 = zeros(1,length(trials)); % 10 white noise trials (10 columns)

for i = 1:18 % since the main array has 12 rows (participants) after the counterbalancing, make 12 different random combination of trials for white noise
             % this will then be matched with the ppt number in the
             % presentation script
randTrialsHand2(i,:) = randperm(numel(trials)); % for each row added and all columns, shuffle numbers 1-10 to have random order of white noise sequence presentation
end

save('randTrialsHand2.mat','randTrialsHand2'); % save a matrix to be loaded to presentation script

% %block2 hand
% randTrialsFoot = zeros(1,length(trials)); % 10 white noise trials (10 columns)
% 
% for i = 1:36 % since the main array has 12 rows (participants) after the counterbalancing, make 12 different random combination of trials for white noise
%              % this will then be matched with the ppt number in the
%              % presentation script
% randTrialsFoot(i,:) = randperm(numel(trials)); % for each row added and all columns, shuffle numbers 1-10 to have random order of white noise sequence presentation
% end
% 
% save('randTrialsFoot.mat','randTrialsFoot'); % save a matrix to be loaded to presentation script
% 
% % %block1 foot
% % randTrialsBlock1Foot = zeros(1,length(trials)); % 10 white noise trials (10 columns)
% % 
% % for i = 1:15 % since the main array has 12 rows (participants) after the counterbalancing, make 12 different random combination of trials for white noise
% %              % this will then be matched with the ppt number in the
% %              % presentation script
% % randTrialsBlock1Foot(i,:) = randperm(numel(trials)); % for each row added and all columns, shuffle numbers 1-10 to have random order of white noise sequence presentation
% % end
% % 
% % save('randTrialsBlock1Foot.mat','randTrialsBlock1Foot'); % save a matrix to be loaded to presentation script
% % 
% % %block2 foot
% % for i = 1:15 % since the main array has 12 rows (participants) after the counterbalancing, make 12 different random combination of trials for white noise
% %              % this will then be matched with the ppt number in the
% %              % presentation script
% % randTrialsBlock2Foot(i,:) = randperm(numel(trials)); % for each row added and all columns, shuffle numbers 1-10 to have random order of white noise sequence presentation
% % end
% % 
% % save('randTrialsBlock2Foot.mat','randTrialsBlock2Foot'); % save a matrix to be loaded to presentation script
% %  

randTrialsHand = [randTrialsHand; randTrialsHand2];
save('randTrialsHand.mat','randTrialsHand');
