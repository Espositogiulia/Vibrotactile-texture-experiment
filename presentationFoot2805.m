%% presentation script hand - 19/05
clear;
clc;
%cd('C:\Users\Nocions\Desktop\GIULIAE\vibroExpGiuliaE\newDesign1905');

%% load stimuli and counterbalancing/randomisation matrices 
load('trials.mat');
load('randTrialsFoot.mat');

 %% start daq session
NI.NI= daq.createSession('ni');
NI.NI=NI_initializeEEG;
%
% %NI.ID contains the ID of the connected NI device
a=daq.getDevices;
%  NI.ID=a(1).ID;
% % %set the DAQ samplingrate to 44100 Hz
%  NI.NI.Rate=44100;
% 
% %
% 
% addAnalogOutputChannel(NI.NI,NI.ID,'ao0','Voltage');
% addAnalogOutputChannel(NI.NI,NI.ID,'ao1','Voltage');
% addAnalogOutputChannel(NI.NI,NI.ID,'ao2','Voltage');
% %

%% get participant number to find order of blocks 
prompt = 'Enter ppt nummber '; % check switch case
ppt = input(prompt); % input number betwee  n 1 and 24 to find row to play condition order

%% number of trials
nTrials = 8;

%% present experiment for 16 trials, 2 blocks (order depends on condCountTot row ie ppt ID)
for thisTrial = 1:size(randTrialsFoot,2) % 16 trials
        for thisPpt = 1:size(randTrialsFoot,1) %15 ppt
            if ppt ==thisPpt
                data = trials{randTrialsFoot(thisPpt,thisTrial)};
                queueOutputData (NI.NI,data);
                prepare(NI.NI);   %prepare the NI
                display(['queueing data ' num2str(randTrialsFoot(thisPpt,thisTrial))]); %present white noise according to participant and trial number
                pause(5);
                NI.NI.startForeground();
                prompt = 'Did you feel one or more types of vibrations within the sequence? Press 0 for no and 1 for more';
                deviantsFoot(thisTrial) = input(prompt);
            end
        end
    end

save(['deviantsFoot' num2str(ppt) '.mat'], 'deviantsFoot');