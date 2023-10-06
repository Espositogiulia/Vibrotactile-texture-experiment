% script to create new events to analyse EEG data individually for each
% spectrotemporal sequence
for i = 1:size(randTrialsHand,1)
    for j = 1:size(randTrialsHand,2)
    if randTrialsHand(i,j) < 9
        wnHand(i,j) = 0;
    else
        wnHand(i,j) = randTrialsHand(i,j)';
    end
end
end


for i = 1:size(randTrialsFoot,1)
    for j = 1:size(randTrialsFoot,2)
    if randTrialsFoot(i,j) < 9
        wnFoot(i,j) = 0;
    else
        wnFoot(i,j) = randTrialsFoot(i,j)';
    end
end
end


vFoot = nonzeros(wnFoot');
wnFoot = reshape(vFoot,8,36);


vHand = nonzeros(wnHand');
wnHand = reshape(vHand,8,36);

for i = 1:size(wnFoot,1)
    for j = 1:size(wnFoot,2)
        if wnFoot(i,j) > 8
            wnFoot(i,j) = wnFoot(i,j) - 8;
        end
    end
end

for i = 1:size(wnHand,1)
    for j = 1:size(wnHand,2)
        if wnHand(i,j) > 8
            wnHand(i,j) = wnHand(i,j) - 8;
        end
    end
end

xF = arrayfun(@num2str, wnFoot, 'UniformOutput', 0);      
xH = arrayfun(@num2str, wnHand, 'UniformOutput', 0) ;    
% wnFoot= randTrialsFoot(~isnan(randTrialsFoot));
% wnFoot = reshape(wnFoot,[36,8]);
% wnHand= randTrialsHand(~isnan(randTrialsHand));

