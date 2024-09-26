%% behavioural script
% import all behavioural results .mat file (4/participant - 2xfoot and 2xhand)
clear;
clc;
cd('C:\data\xvibro\behavioural');
files = dir('*.mat'); % get all mat files in folder
files = natsortfiles(files); %need to add natsortfiles to path for this to work

%% foot
for i = 1:32
    contFoot(i) = load(files(i).name); % the first 18 files are the foot trials
end
clear i;

dataFoot = struct2cell(contFoot); % convert to cell to access the responses
dataFoot = squeeze(dataFoot); % get 2D instead of 3D structure
dataFoot = cell2mat(dataFoot); % convert to double to create array of foot responses (same order as trials)

randTrialsFoot = struct2cell(load(files(65).name)); % this is not in the same order as dataFoot!! - sorted for now - the problem won't be there once I have all ppts
                                                    % row 37 of structure
                                                    % contains foot trials
                                                    % for now!! 
randTrialsFoot = cell2mat(randTrialsFoot);  % conver to double to create array of foot trials

trialsFoot1 = randTrialsFoot(1:2,:); % all this will be unnecessary once I have the whole dataset, some data are on cookie and I need to copy it
trialsFoot2 = randTrialsFoot(3:4,:); % I am missing ppt 16(15) i.e. trials 31-32
trialsFoot3 = randTrialsFoot(5:6,:);
trialsFoot4 = randTrialsFoot(7:8,:);
trialsFoot5 = randTrialsFoot(9:10,:);
trialsFoot6 = randTrialsFoot(11:12,:);
trialsFoot7 = randTrialsFoot(13:14,:);
trialsFoot8 = randTrialsFoot(17:18,:);
trialsFoot9 = randTrialsFoot(19:20,:);
trialsFoot10 = randTrialsFoot(21:22,:);
trialsFoot11 = randTrialsFoot(23:24,:);
trialsFoot12 = randTrialsFoot(25:26,:);
trialsFoot13 = randTrialsFoot(27:28,:);
trialsFoot14 = randTrialsFoot(29:30,:);
trialsFoot16 = randTrialsFoot(35:36,:);
trialsFoot17 = randTrialsFoot(15:16,:);

trialsFootSoFar = [trialsFoot1; trialsFoot2; trialsFoot3; trialsFoot4; ...
    trialsFoot5; trialsFoot6; trialsFoot7; trialsFoot8; trialsFoot9; ...
    trialsFoot10; trialsFoot11; trialsFoot12; trialsFoot13; trialsFoot14; trialsFoot16; trialsFoot17]; % this will also be unnecessary

for i = 1:size(trialsFootSoFar)
    for j = 1:size(trialsFootSoFar,2)
        if trialsFootSoFar(i,j) <= 8 % trials 1-8 are frequency
            dataFreqFoot(i,j) = dataFoot(i,j); % get only freq reponses from the main array that has all freq and wn responses
        else
            dataFreqFoot(i,j) = NaN; %trials 9-16 are white noise, wn responses set to NaN
        end
    end
end
clear i;
clear j;

for i = 1:size(trialsFootSoFar)
    for j = 1:size(trialsFootSoFar,2)
        if trialsFootSoFar(i,j) > 8 % trials 9-16 are white noise, exclude wn2
            dataWnFoot(i,j) = dataFoot(i,j); % get wn responses
        else
            dataWnFoot(i,j) = NaN; %trials 1-8 are frequency, freq reponses set to NaN
        end
        
    end
end
clear i;
clear j;
%% 
dataFreqFootOnly = zeros(32,8);  % create an array to store 18 (2blocks x ppt - 18 rows, 9 ppts) x 8 (8 trials x block - 8 columns)

for i = 1:size(dataFreqFoot,1)
dataFreqFootOnly(i,:) = rmmissing(dataFreqFoot(i,:)); % remove NaN values to only get array with freq reponses
end
clear i;

ppt1FreqFoot = dataFreqFootOnly(1:2,:); % I need to put 2 rows (ie 2 blocks for 1 ppt) together in one single row
%ppt1FreqFoot = reshape(ppt1FreqFoot,1,16);
ppt2FreqFoot = dataFreqFootOnly(3:4,:);
%ppt2FreqFoot = reshape(ppt2FreqFoot,1,16);
ppt3FreqFoot = dataFreqFootOnly(5:6,:);
%ppt3FreqFoot = reshape(ppt3FreqFoot,1,16);
 ppt4FreqFoot = dataFreqFootOnly(7:8,:);
 %ppt4FreqFoot = reshape(ppt4FreqFoot,1,16);
ppt5FreqFoot = dataFreqFootOnly(9:10,:);
%ppt5FreqFoot = reshape(ppt5FreqFoot,1,16);
ppt6FreqFoot = dataFreqFootOnly(11:12,:);
%ppt6FreqFoot = reshape(ppt6FreqFoot,1,16);
ppt7FreqFoot = dataFreqFootOnly(13:14,:);
%ppt7FreqFoot = reshape(ppt7FreqFoot,1,16);
ppt8FreqFoot = dataFreqFootOnly(15:16,:);
%ppt8FreqFoot = reshape(ppt8FreqFoot,1,16);
ppt9FreqFoot = dataFreqFootOnly(17:18,:);
%ppt9FreqFoot = reshape(ppt9FreqFoot,1,16);
ppt10FreqFoot = dataFreqFootOnly(19:20,:);
%ppt10FreqFoot = reshape(ppt10FreqFoot,1,16);
ppt11FreqFoot = dataFreqFootOnly(21:22,:);
%ppt11FreqFoot = reshape(ppt11FreqFoot,1,16);
ppt12FreqFoot = dataFreqFootOnly(23:24,:);
%ppt12FreqFoot = reshape(ppt12FreqFoot,1,16);
ppt13FreqFoot = dataFreqFootOnly(25:26,:);
%ppt13FreqFoot = reshape(ppt13FreqFoot,1,16);
ppt14FreqFoot = dataFreqFootOnly(27:28,:);
%ppt14FreqFoot = reshape(ppt14FreqFoot,1,16);
ppt16FreqFoot = dataFreqFootOnly(29:30,:);
%ppt16FreqFoot = reshape(ppt16FreqFoot,1,16);
ppt17FreqFoot = dataFreqFootOnly(31:32,:);
%ppt17FreqFoot = reshape(ppt17FreqFoot,1,16);

dataAllFootFreq = [ppt1FreqFoot; ppt2FreqFoot; ppt3FreqFoot; ppt4FreqFoot; ...
    ppt5FreqFoot; ...
    ppt6FreqFoot; ppt7FreqFoot; ppt8FreqFoot; ppt9FreqFoot; ppt10FreqFoot; ...
    ppt11FreqFoot; ppt12FreqFoot; ppt13FreqFoot;ppt14FreqFoot;ppt16FreqFoot; ppt17FreqFoot]; % put data for all ppts together in one array

for i = 1:size(dataAllFootFreq,1)
        percFeltFootFreq(i) = sum(dataAllFootFreq(i,:)== 1)/... % how many 1 responses there are (sum)
            numel(dataAllFootFreq(i,:))*100; % divide by number of total responses (16) and get % of yes responses
end
clear i;

avgFreqFoot = mean(percFeltFootFreq,2); % total average
sDevFreqFoot = std(percFeltFootFreq);
%% white noise foot
dataWnFootOnly = zeros(32,8);  % create an array to store 18 (2blocks x ppt - 18 rows, 9 ppts) x 8 (8 trials x block - 8 columns)

for i = 1:size(dataWnFoot,1)
dataWnFootOnly(i,:) = rmmissing(dataWnFoot(i,:)); % remove NaN values to only get array with freq reponses
end
clear i;

ppt1WnFoot = dataWnFootOnly(1:2,:); % I need to put 2 rows (ie 2 blocks for 1 ppt) together in one single row
%ppt1WnFoot = reshape(ppt1WnFoot,1,16); % there is a mistake here!!!!!!
ppt2WnFoot = dataWnFootOnly(3:4,:);
%ppt2WnFoot = reshape(ppt2WnFoot,1,16);
ppt3WnFoot = dataWnFootOnly(5:6,:);
%ppt3WnFoot = reshape(ppt3WnFoot,1,16);
 ppt4WnFoot = dataWnFootOnly(7:8,:);
% ppt4WnFoot = reshape(ppt4WnFoot,1,16);
ppt5WnFoot = dataWnFootOnly(9:10,:);
%ppt5WnFoot = reshape(ppt5WnFoot,1,16);
ppt6WnFoot = dataWnFootOnly(11:12,:);
%ppt6WnFoot = reshape(ppt6WnFoot,1,16);
ppt7WnFoot = dataWnFootOnly(13:14,:);
%ppt7WnFoot = reshape(ppt7WnFoot,1,16);
ppt8WnFoot = dataWnFootOnly(15:16,:);
%ppt8WnFoot = reshape(ppt8WnFoot,1,16);
ppt9WnFoot = dataWnFootOnly(17:18,:);
%ppt9WnFoot = reshape(ppt9WnFoot,1,16);
ppt10WnFoot = dataWnFootOnly(19:20,:);
%ppt10WnFoot = reshape(ppt10WnFoot,1,16);
ppt11WnFoot = dataWnFootOnly(21:22,:);
%ppt11WnFoot = reshape(ppt11WnFoot,1,16);
ppt12WnFoot = dataWnFootOnly(23:24,:);
%ppt12WnFoot = reshape(ppt12WnFoot,1,16);
ppt13WnFoot = dataWnFootOnly(25:26,:);
%ppt13WnFoot = reshape(ppt13WnFoot,1,16);
ppt14WnFoot = dataWnFootOnly(27:28,:);
%ppt14WnFoot = reshape(ppt14WnFoot,1,16);
ppt16WnFoot = dataWnFootOnly(29:30,:); % ppt 17
%ppt16WnFoot = reshape(ppt16WnFoot,1,16);
ppt17WnFoot = dataWnFootOnly(31:32,:); % ppt 18
%ppt17WnFoot = reshape(ppt17WnFoot,1,16);

dataAllFootWn = [ppt1WnFoot; ppt2WnFoot; ppt3WnFoot; ppt4WnFoot; ...
    ppt5WnFoot; ...
    ppt6WnFoot; ppt7WnFoot; ppt8WnFoot; ppt9WnFoot; ppt10WnFoot; ...
    ppt11WnFoot; ppt12WnFoot; ppt13WnFoot;ppt14WnFoot;ppt16WnFoot; ppt17WnFoot]; % put data for all ppts together in one array

for i = 1:size(dataAllFootWn,1)
        percFeltFootWn(i) = sum(dataAllFootWn(i,:)== 1)/... % how many 1 responses there are (sum)
            numel(dataAllFootWn(i,:))*100; % divide by number of total responses (16) and get % of yes responses
end
clear i;

avgWnFoot = mean(percFeltFootWn,2); % total average
sDevWnFoot = std(percFeltFootWn);

%% Hand
for i = 33:64
    contHand(i) = load(files(i).name); % the first 18 files are the Hand trials
end
clear i;

dataHand = struct2cell(contHand); % convert to cell to access the responses
dataHand = squeeze(dataHand); % get 2D instead of 3D structure
dataHand = cell2mat(dataHand); % convert to double to create array of Hand responses (same order as trials)

randTrialsHand = struct2cell(load(files(66).name)); % this is not in the same order as dataHand!! - sorted for now - the problem won't be there once I have all ppts
                                                    % row 37 of structure
                                                    % contains Hand trials
                                                    % for now!! 
randTrialsHand = cell2mat(randTrialsHand);  % conver to double to create array of Hand trials

trialsHand1 = randTrialsHand(1:2,:); % all this will be unnecessary once I have the whole dataset, some data are on cookie and I need to copy it
trialsHand2 = randTrialsHand(3:4,:);
trialsHand3 = randTrialsHand(5:6,:);
trialsHand4 = randTrialsHand(7:8,:);
trialsHand5 = randTrialsHand(9:10,:);
trialsHand6 = randTrialsHand(11:12,:);
trialsHand7 = randTrialsHand(13:14,:);
trialsHand8 = randTrialsHand(17:18,:);
trialsHand9 = randTrialsHand(19:20,:);
trialsHand10 = randTrialsHand(21:22,:);
trialsHand11 = randTrialsHand(23:24,:);
trialsHand12 = randTrialsHand(25:26,:);
trialsHand13 = randTrialsHand(27:28,:);
trialsHand14 = randTrialsHand(29:30,:);
trialsHand16 = randTrialsHand(35:36,:);
trialsHand17 = randTrialsHand(15:16,:);

trialsHandSoFar = [trialsHand1; trialsHand2; trialsHand3; trialsHand4; ...
    trialsHand5; trialsHand6; trialsHand7; trialsHand8; trialsHand9; ...
    trialsHand10; trialsHand11; trialsHand12; trialsHand13; trialsHand14; trialsHand16; trialsHand17]; % this will also be unnecessary

for i = 1:size(trialsHandSoFar)
    for j = 1:size(trialsHandSoFar,2)
        if trialsHandSoFar(i,j) <= 8 % trials 1-8 are frequency
            dataFreqHand(i,j) = dataHand(i,j); % get only freq reponses from the main array that has all freq and wn responses
        else
            dataFreqHand(i,j) = NaN; %trials 9-16 are white noise, wn responses set to NaN
        end
    end
end
clear i;
clear j;

for i = 1:size(trialsHandSoFar)
    for j = 1:size(trialsHandSoFar,2)
        if trialsHandSoFar(i,j) > 8 % trials 9-16 are white noise
            dataWnHand(i,j) = dataHand(i,j); % get wn responses
        else
            dataWnHand(i,j) = NaN; %trials 1-8 are frequency, freq reponses set to NaN
        end
    end
end
clear i;
clear j;
%% 
dataFreqHandOnly = zeros(32,8);  % create an array to store 18 (2blocks x ppt - 18 rows, 9 ppts) x 8 (8 trials x block - 8 columns)

for i = 1:size(dataFreqHand,1)
dataFreqHandOnly(i,:) = rmmissing(dataFreqHand(i,:)); % remove NaN values to only get array with freq reponses
end
clear i;

ppt1FreqHand = dataFreqHandOnly(1:2,:); % I need to put 2 rows (ie 2 blocks for 1 ppt) together in one single row
%ppt1FreqHand = reshape(ppt1FreqHand,1,16);
ppt2FreqHand = dataFreqHandOnly(3:4,:);
%ppt2FreqHand = reshape(ppt2FreqHand,1,16);
ppt3FreqHand = dataFreqHandOnly(5:6,:);
%ppt3FreqHand = reshape(ppt3FreqHand,1,16);
 ppt4FreqHand = dataFreqHandOnly(7:8,:);
% ppt4FreqHand = reshape(ppt4FreqHand,1,16);
ppt5FreqHand = dataFreqHandOnly(9:10,:);
%ppt5FreqHand = reshape(ppt5FreqHand,1,16);
ppt6FreqHand = dataFreqHandOnly(11:12,:);
%ppt6FreqHand = reshape(ppt6FreqHand,1,16);
ppt7FreqHand = dataFreqHandOnly(13:14,:);
%ppt7FreqHand = reshape(ppt7FreqHand,1,16);
ppt8FreqHand = dataFreqHandOnly(15:16,:);
%ppt8FreqHand = reshape(ppt8FreqHand,1,16);
ppt9FreqHand = dataFreqHandOnly(17:18,:);
%ppt9FreqHand = reshape(ppt9FreqHand,1,16);
ppt10FreqHand = dataFreqHandOnly(19:20,:);
%ppt10FreqHand = reshape(ppt10FreqHand,1,16);
ppt11FreqHand = dataFreqHandOnly(21:22,:);
%ppt11FreqHand = reshape(ppt11FreqHand,1,16);
ppt12FreqHand = dataFreqHandOnly(23:24,:);
%ppt12FreqHand = reshape(ppt12FreqHand,1,16);
ppt13FreqHand = dataFreqHandOnly(25:26,:);
%ppt13FreqHand = reshape(ppt13FreqHand,1,16);
ppt14FreqHand = dataFreqHandOnly(27:28,:);
%ppt14FreqHand = reshape(ppt14FreqHand,1,16);
ppt16FreqHand = dataFreqHandOnly(29:30,:);
%ppt16FreqHand = reshape(ppt16FreqHand,1,16);
ppt17FreqHand = dataFreqHandOnly(31:32,:);
%ppt17FreqHand = reshape(ppt17FreqHand,1,16);

dataAllHandFreq = [ppt1FreqHand; ppt2FreqHand; ppt3FreqHand; ppt4FreqHand; ...
    ppt5FreqHand; ...
    ppt6FreqHand; ppt7FreqHand; ppt8FreqHand; ppt9FreqHand; ppt10FreqHand; ...
    ppt11FreqHand; ppt12FreqHand; ppt13FreqHand;ppt14FreqHand;ppt16FreqHand; ppt17FreqHand]; % put data for all ppts together in one array

for i = 1:size(dataAllHandFreq,1)
        percFeltHandFreq(i) = sum(dataAllHandFreq(i,:)== 1)/... % how many 1 responses there are (sum)
            numel(dataAllHandFreq(i,:))*100; % divide by number of total responses (16) and get % of yes responses
end
clear i;

avgFreqHand = mean(percFeltHandFreq,2); % total average
sDevFreqHand = std(percFeltHandFreq);

%% white noise Hand
dataWnHandOnly = zeros(32,8);  % create an array to store 18 (2blocks x ppt - 18 rows, 9 ppts) x 8 (8 trials x block - 8 columns)

for i = 1:size(dataWnHand,1)
dataWnHandOnly(i,:) = rmmissing(dataWnHand(i,:)); % remove NaN values to only get array with freq reponses
end
clear i;

ppt1WnHand = dataWnHandOnly(1:2,:); % I need to put 2 rows (ie 2 blocks for 1 ppt) together in one single row
ppt1WnHand = [ppt1WnHand(1,:) ppt1WnHand(2,:)];
%ppt1WnHand = reshape(ppt1WnHand,1,16);
ppt2WnHand = dataWnHandOnly(3:4,:);
ppt2WnHand = [ppt2WnHand(1,:) ppt2WnHand(2,:)];
%ppt2WnHand = reshape(ppt2WnHand,1,16);
ppt3WnHand = dataWnHandOnly(5:6,:);
ppt3WnHand = [ppt3WnHand(1,:) ppt3WnHand(2,:)];
%ppt3WnHand = reshape(ppt3WnHand,1,16);
 ppt4WnHand = dataWnHandOnly(7:8,:);
 ppt4WnHand = [ppt4WnHand(1,:) ppt4WnHand(2,:)];
% ppt4WnHand = reshape(ppt4WnHand,1,16);
ppt5WnHand = dataWnHandOnly(9:10,:);
ppt5WnHand = [ppt5WnHand(1,:) ppt5WnHand(2,:)];
%ppt5WnHand = reshape(ppt5WnHand,1,16);
ppt6WnHand = dataWnHandOnly(11:12,:);
ppt6WnHand = [ppt6WnHand(1,:) ppt6WnHand(2,:)];
%ppt6WnHand = reshape(ppt6WnHand,1,16);
ppt7WnHand = dataWnHandOnly(13:14,:);
ppt7WnHand = [ppt7WnHand(1,:) ppt7WnHand(2,:)];
%ppt7WnHand = reshape(ppt7WnHand,1,16);
ppt8WnHand = dataWnHandOnly(15:16,:);
ppt8WnHand = [ppt8WnHand(1,:) ppt8WnHand(2,:)];
%ppt8WnHand = reshape(ppt8WnHand,1,16);
ppt9WnHand = dataWnHandOnly(17:18,:);
ppt9WnHand = [ppt9WnHand(1,:) ppt9WnHand(2,:)];
%ppt9WnHand = reshape(ppt9WnHand,1,16);
ppt10WnHand = dataWnHandOnly(19:20,:);
ppt10WnHand = [ppt10WnHand(1,:) ppt10WnHand(2,:)];
%ppt10WnHand = reshape(ppt10WnHand,1,16);
ppt11WnHand = dataWnHandOnly(21:22,:);
ppt11WnHand = [ppt11WnHand(1,:) ppt11WnHand(2,:)];
%ppt11WnHand = reshape(ppt11WnHand,1,16);
ppt12WnHand = dataWnHandOnly(23:24,:);
ppt12WnHand = [ppt12WnHand(1,:) ppt12WnHand(2,:)];
%ppt12WnHand = reshape(ppt12WnHand,1,16);
ppt13WnHand = dataWnHandOnly(25:26,:);
ppt13WnHand = [ppt13WnHand(1,:) ppt13WnHand(2,:)];
%ppt13WnHand = reshape(ppt13WnHand,1,16);
ppt14WnHand = dataWnHandOnly(27:28,:);
ppt14WnHand = [ppt14WnHand(1,:) ppt14WnHand(2,:)];
%ppt14WnHand = reshape(ppt14WnHand,1,16);
ppt16WnHand = dataWnHandOnly(29:30,:);
ppt16WnHand = [ppt16WnHand(1,:) ppt16WnHand(2,:)];
%ppt16WnHand = reshape(ppt16WnHand,1,16);
ppt17WnHand = dataWnHandOnly(31:32,:);
ppt17WnHand = [ppt17WnHand(1,:) ppt17WnHand(2,:)];
%ppt17WnHand = reshape(ppt17WnHand,1,16);

dataAllHandWn = [ppt1WnHand; ppt2WnHand; ppt3WnHand; ppt4WnHand; 
    ppt5WnHand; ...
    ppt6WnHand; ppt7WnHand; ppt8WnHand; ppt9WnHand; ppt10WnHand; ...
    ppt11WnHand; ppt12WnHand; ppt13WnHand;ppt14WnHand;ppt16WnHand; ppt17WnHand]; % put data for all ppts together in one array

for i = 1:size(dataAllHandWn,1)
        percFeltHandWn(i) = sum(dataAllHandWn(i,:)== 1)/... % how many 1 responses there are (sum)
            16*100; % divide by number of total responses (16) and get % of yes responses
end
clear i;

avgWnHand = mean(percFeltHandWn,2); % total average
sDevWnHand = std(percFeltHandWn);
%% check if seq 2 was perceived more often
% hand
for i = 1:8
    sequence(i).resp=NaN(size(trialsHandSoFar));
end

for i = 1:size(trialsHandSoFar,1)
    for j = 1:size(trialsHandSoFar,2)
        if trialsHandSoFar(i,j)==9
        sequence(1).resp(i,j)= dataHand(i,j);
            elseif trialsHandSoFar(i,j)==10
        sequence(2).resp(i,j)= dataHand(i,j);
            elseif trialsHandSoFar(i,j)==11
        sequence(3).resp(i,j)= dataHand(i,j);
        elseif trialsHandSoFar(i,j)==12
        sequence(4).resp(i,j)= dataHand(i,j);
        elseif trialsHandSoFar(i,j)==13
        sequence(5).resp(i,j)= dataHand(i,j);
        elseif trialsHandSoFar(i,j)==14
        sequence(6).resp(i,j)= dataHand(i,j);
        elseif trialsHandSoFar(i,j)==15
        sequence(7).resp(i,j)= dataHand(i,j);
        elseif trialsHandSoFar(i,j)==16
        sequence(8).resp(i,j)= dataHand(i,j);
        end
    end
end

% find responses to specific wn sequences (each 2 rows = 2 trials for the
% same participant
for i = 1:size(sequence(1).resp,1)
    respHand1(i,:)= rmmissing(sequence(1).resp(i,:));
    respHand2(i,:)= rmmissing(sequence(2).resp(i,:));
    respHand3(i,:)= rmmissing(sequence(3).resp(i,:));
    respHand4(i,:)= rmmissing(sequence(4).resp(i,:));
    respHand5(i,:)= rmmissing(sequence(5).resp(i,:));
    respHand6(i,:)= rmmissing(sequence(6).resp(i,:));
    respHand7(i,:)= rmmissing(sequence(7).resp(i,:));
    respHand8(i,:)= rmmissing(sequence(8).resp(i,:));
end

respHand1 = reshape(respHand1, 16, 2);
respHand2 = reshape(respHand2, 16, 2);
respHand3 = reshape(respHand3, 16, 2);
respHand4 = reshape(respHand4, 16, 2);
respHand5 = reshape(respHand5, 16, 2);
respHand6 = reshape(respHand6, 16, 2);
respHand7 = reshape(respHand7, 16, 2);
respHand8= reshape(respHand8, 16, 2);


respHand1Tot= mean(sum(respHand1,2)/2*100);
respHand2Tot= mean(sum(respHand2,2)/2*100);
respHand3Tot= mean(sum(respHand3,2)/2*100);
respHand4Tot= mean(sum(respHand4,2)/2*100);
respHand5Tot= mean(sum(respHand5,2)/2*100);
respHand6Tot= mean(sum(respHand6,2)/2*100);
respHand7Tot= mean(sum(respHand7,2)/2*100);
respHand8Tot= mean(sum(respHand8,2)/2*100);

respHand= [respHand1Tot respHand2Tot respHand3Tot respHand4Tot respHand5Tot respHand6Tot respHand7Tot respHand8Tot];
% foot

for i = 1:8
    sequence(i).resp=NaN(size(trialsFootSoFar));
end

for i = 1:size(trialsFootSoFar)
    for j = 1:size(trialsFootSoFar,2)
        if trialsFootSoFar(i,j)==9
        sequence(1).resp(i,j)= dataFoot(i,j);
            elseif trialsFootSoFar(i,j)==10
        sequence(2).resp(i,j)= dataFoot(i,j);
            elseif trialsFootSoFar(i,j)==11
        sequence(3).resp(i,j)= dataFoot(i,j);
        elseif trialsFootSoFar(i,j)==12
        sequence(4).resp(i,j)= dataFoot(i,j);
        elseif trialsFootSoFar(i,j)==13
        sequence(5).resp(i,j)= dataFoot(i,j);
        elseif trialsFootSoFar(i,j)==14
        sequence(6).resp(i,j)= dataFoot(i,j);
        elseif trialsFootSoFar(i,j)==15
        sequence(7).resp(i,j)= dataFoot(i,j);
        elseif trialsFootSoFar(i,j)==16
        sequence(8).resp(i,j)= dataFoot(i,j);
        end
    end
end

% find responses to specific wn sequences (each 2 rows = 2 trials for the
% same participant
for i = 1:size(sequence(1).resp,1)
    respFoot1(i,:)= rmmissing(sequence(1).resp(i,:));
    respFoot2(i,:)= rmmissing(sequence(2).resp(i,:));
    respFoot3(i,:)= rmmissing(sequence(3).resp(i,:));
    respFoot4(i,:)= rmmissing(sequence(4).resp(i,:));
    respFoot5(i,:)= rmmissing(sequence(5).resp(i,:));
    respFoot6(i,:)= rmmissing(sequence(6).resp(i,:));
    respFoot7(i,:)= rmmissing(sequence(7).resp(i,:));
    respFoot8(i,:)= rmmissing(sequence(8).resp(i,:));
end

respFoot1 = reshape(respFoot1, 16, 2);
respFoot2 = reshape(respFoot2, 16, 2);
respFoot3 = reshape(respFoot3, 16, 2);
respFoot4 = reshape(respFoot4, 16, 2);
respFoot5 = reshape(respFoot5, 16, 2);
respFoot6 = reshape(respFoot6, 16, 2);
respFoot7 = reshape(respFoot7, 16, 2);
respFoot8= reshape(respFoot8, 16, 2);


respFoot1Tot= mean(sum(respFoot1,2)/2*100);
respFoot2Tot= mean(sum(respFoot2,2)/2*100);
respFoot3Tot= mean(sum(respFoot3,2)/2*100);
respFoot4Tot= mean(sum(respFoot4,2)/2*100);
respFoot5Tot= mean(sum(respFoot5,2)/2*100);
respFoot6Tot= mean(sum(respFoot6,2)/2*100);
respFoot7Tot= mean(sum(respFoot7,2)/2*100);
respFoot8Tot= mean(sum(respFoot8,2)/2*100);

respFoot= [respFoot1Tot respFoot2Tot respFoot3Tot respFoot4Tot respFoot5Tot respFoot6Tot respFoot7Tot respFoot8Tot];
