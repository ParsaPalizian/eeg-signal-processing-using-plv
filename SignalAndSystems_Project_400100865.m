close all,clear,clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%please notice to comment and break points%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%preprocessing
%step 1,2
load Subject1.mat
load Subject2.mat

new_subject1=table2array(subject1)';
new_subject2=table2array(subject2)';

new_subject1=new_subject1(1:19,:);
new_subject2=new_subject2(1:19,:);

%You should break point this line
%step 3,4,5
%%step 3
% size(EEG.data)
ns3 = []
remove_time = 14*200;
ns2 = EEG.data(:,remove_time:end );
for i=1:120
    a = ns2(:,i*600 + 4200:(i+1)*600-1 +  4200);
    ns3(: , : , i) = a;
end
%%step 4
[channels, samples, trials] = size(ns3);
noisyTrials = [];
for ch = 1 : channels
    p = nan(trials, samples); 
    for trial = 1:trials
        trialData = squeeze(ns3(ch, :, trial)); % Extract
        %data for the current trial
        power = fft(trialData).^2; % power spectrum for the current channel
        p(trial, :) = power;
    end    
    % suggested commands
    variances = sum(nanstd(p, [], 2).^2, 2);
    noisyTrialsCh = find(abs(zscore(variances)) > 3.5);   
    noisyTrials = union(noisyTrials, noisyTrialsCh);
end
cleanEpochData = ns3(:, :, setdiff(1:trials, noisyTrials));
save('cleand_data.mat', 'cleanEpochData');
%%step 5
channelsToKeep = [1, 5, 10, 15];
numChannels = 19;
channelsIdx = zeros(1, numChannels);
for i = channelsToKeep
    channelsIdx(i) = 1;
end
selectedEpochData = ns3(logical(channelsIdx),:,:);
save('final.mat', 'selectedEpochData');
Subject1_final = selectedEpochData;
Normal = load('Normal.mat','normal');
odor = Normal.normal.odor;
noisy = Normal.normal.noisy;
myStruct = struct(...
    'Subject1_final', Subject1_final, ...
    'odor', odor, ...
    'noisy',noisy);
save('subject1_struct.mat', 'myStruct');

%You should break point this line
%4.1 and 4.2
close all , clc , clear
load Normal.mat
load AD.mat
normalFrequentPLV = [];
for np = normal
    normalFrequentPLV=[normalFrequentPLV PLV2(np.epoch , 200 , [35 40] , 2 , 3 , np.odor , 0)];
end
subplot(2,2,1);
boxplot(normalFrequentPLV);
title("NORMAL FREQUENT");
normalRarePLV  = [];
for np = normal
    normalRarePLV  =[normalRarePLV   PLV2(np.epoch , 200 , [35 40] , 2 , 3 , np.odor , 1)];
end
subplot(2,2,2);
boxplot(normalRarePLV);
title("NORMAL RARE");
ADFrequentPLV= [];
for np = AD
    ADFrequentPLV=[ADFrequentPLV PLV2(np.epoch , 200 , [35 40] , 2 , 3 , np.odor , 0)];
end
subplot(2,2 , 3);
boxplot(ADFrequentPLV);
title("AD FREQUENT");
ADRarePLV = [];
for np = AD
    ADRarePLV =[ADRarePLV PLV2(np.epoch , 200 , [35 40] , 2 , 3 , np.odor , 1)];
end
subplot(2,2,4);
boxplot(ADRarePLV);
title("AD RARE");
figure;
nfPLVm = mean(normalFrequentPLV)
nfPLVstd = std(normalFrequentPLV)
subplot(2,2,1);
% hist(normrnd(nfPLVm,nfPLVstd , 15))
pd = fitdist(normalFrequentPLV','Normal')
x = -2:0.01:2
y = pdf(pd,x);
plot(x,y)
title("NORMAL FREQUENT")
nrPLVm = mean(normalRarePLV);
nrPLVstd = std(normalRarePLV);
subplot(2,2,2);
pd = fitdist(normalRarePLV','Normal')
x = -2:0.01:2
y = pdf(pd,x);
plot(x,y)
title("NORMAL RARE")
ADfPLVm = mean(ADFrequentPLV);
ADfPLVstd = std(ADFrequentPLV);
subplot(2,2,3);
pd = fitdist(ADFrequentPLV','Normal')
x = -2:0.01:2
y = pdf(pd,x);
plot(x,y)
title("AD FREQUENT")
ADrPLVm = mean(ADRarePLV);
ADrPLVstd = std(ADRarePLV);
subplot(2,2,4);
pd = fitdist(ADRarePLV','Normal')
x = -2:0.01:2
y = pdf(pd,x);
plot(x,y)
title("AD RARE")
figure 
[h , PvalNormalFreqADFreq] = ttest2(normalFrequentPLV , ADFrequentPLV);
PvalNormalFreqADFreq
[h , PvalNormalRareADRare] = ttest2(normalRarePLV , ADRarePLV);
PvalNormalRareADRare
subplot(1,2,1);
polarhistogram(normalFrequentPLV, 2)
title("NORMAL FREQUENT")
subplot(1,2,2);
polarhistogram(ADFrequentPLV, 2)
title("AD FREQUENT")

%You should break point this line
%4.4
clear, clc, close all

load Normal.mat
load AD.mat

%%Mean Of Phase Difference of Any Subject
mms = []
for fr = 1:1:15
    odor = normal(fr ).odor;
    fzf = normal(fr ).epoch(2 , :, odor == 0) ;
    czf = normal(fr ).epoch(3 , :, odor == 0) ;
    phis = [];
    for i = 1:1:size(fzf, 3)
        x = fzf(1, : , i);
        y = czf(1, : , i);
        PhDiff = phdiffmeasure(x, y);
        phis = [phis  PhDiff ] ;
    end
    mm = mean(phis,2)
    mms = [mms  mm]
end
subplot(2,2,1)
polarhistogram(mms)
title("Mean Of Phase Difference - POLAR")
subplot(2,2,2)
plot(1:1:15 , mms)
title("Mean Of Phase Difference - PLOT")
%%all Phase Differences 
phis = [];
for fr = 1:1:15
    odor = normal(fr ).odor;
    fzf = normal(fr ).epoch(2 , :, odor == 0) ;
    czf = normal(fr ).epoch(3 , :, odor == 0) ;
    for i = 1:1:size(fzf, 3)
        x = fzf(1, : , i);
        y = czf(1, : , i);
        PhDiff = phdiffmeasure(x, y);
        phis = [phis  PhDiff ] ;
    end
end
subplot(2,2,3)
polarhistogram(phis)
title("all Phase Differences - POLAR")
subplot(2,2,4)
plot(1:1:size(phis, 2) , phis)
title("all Phase Differences -  PLOT")

%You should break point this line
%4.5
close all , clc , clear
load Normal.mat;
load AD.mat;
normalFreq = [];
normalRare= [];
i = 1;
for c1=1:4
    for c2=1:4
        nf = [];
        rf = [];
        for np = normal
            plv = PLV2(np.epoch , 200 , [35 40] , c1 , c2 , np.odor , 0);
            nf = [nf  plv ] ;
            plv = PLV2(np.epoch , 200 , [35 40] , c1 , c2 , np.odor , 1);
            rf = [rf  plv ] ;
        end
        normalFreq(c1 , c2)  = mean(nf,2);
        normalRare(c1 , c2)= mean(rf,2);;
    end
end
subplot(2,2,1);
heatmap(normalFreq ,  'XData', ["Fp1" "Fz" "Cz", "Pz"] , "YData",["Fp1" "Fz" "Cz", "Pz"] );
title("NORMAL FREQUENT");
subplot(2,2,2);
heatmap(normalRare,  'XData', ["Fp1" "Fz" "Cz", "Pz"] , "YData",["Fp1" "Fz" "Cz", "Pz"] );
title("NORMAL RARE");
ADFreq = [];
ADRare= [];
i = 1;
for c1=1:4
    for c2=1:4
        nf = [];
        rf = [];
        for np = AD
            plv = PLV2(np.epoch , 200 , [35 40] , c1 , c2 , np.odor , 0);
            nf = [nf  plv ] ;
            plv = PLV2(np.epoch , 200 , [35 40] , c1 , c2 , np.odor , 1);
            rf = [rf  plv ] ;
        end
        ADFreq (c1 , c2)  = mean(nf,2);
        ADRare(c1 , c2)= mean(rf,2);;
    end
end
subplot(2,2,3);
heatmap(ADFreq ,  'XData', ["Fp1" "Fz" "Cz", "Pz"] , "YData",["Fp1" "Fz" "Cz", "Pz"] );
title("AD FREQUENT");
subplot(2,2,4);
heatmap(ADRare,  'XData', ["Fp1" "Fz" "Cz", "Pz"] , "YData",["Fp1" "Fz" "Cz", "Pz"] );
title("AD RARE");

%You should break point this line
%5.1.2

mci = load('MCI.mat');
freqRange= [35, 40]; 
fs = 200
mplv = zeros(7, 2); 
%MCI  person
for i = 1:7
    [~ , ~ , trials] = size(mci.MCI(i).epoch);  
    epoch = mci.MCI(i).epoch; 
    odors = mci.MCI(i).odor; 
    % calc PLV
    for trial = 1 : trials
        freqOdor = sum(odors(:) == 0);
        rareOdor = sum(odors(:) == 1);
        odor = odors(trial,1)
        sfz = epoch(2,:,trial);
        scz = epoch(3,:,trial);
        mplv(i,odor+1) = mplv(i,odor+1) + calcPLV(sfz, scz, fs , freqRange);
    end
    % mean
    mplv(i,1) = mplv(i,1)/ freqOdor;
    mplv(i,2) = mplv(i,2)/ rareOdor;
end 
disp('MCI Group PLV:');
disp(mplv);
save('mplv.mat', 'mplv');

function plv = calcPLV(sfz, fcz, fs, range)
    sfzf= bandpass(sfz, range, fs);
    fczf= bandpass(fcz, range, fs);
    ph1= angle(hilbert(sfzf));
    ph2= angle(hilbert(fczf));
    dff= ph2- ph1;
    plv = abs(mean(exp(1i * dff)));
end

function plv = PLV2(epoch, samrate, filterRange, channelA, channelB, odors, type)
    samples = size(epoch,2);
    trials =  size(epoch,3);
    filterData = filter(fir1(50, 2/samrate*filterRange),1, epoch, [], 2);
    z = zeros(2, size(epoch,2), size(epoch,3));
    z(1,:,:) = filterData(channelA, :, :)+1i*hilbert(filterData(channelA, :, :));
    z(2,:,:) = filterData(channelB, :, :)+1i*hilbert(filterData(channelB, :, :));
    plv_hat = zeros(1,trials);
    num = 0;
    for trial = 1:trials
        if odors(trial)==type
            phi = angle(z(1,:,trial).*conj(z(2,:,trial)));
            plv_hat(trial) = abs(sum(exp(1i*phi))/samples);
            num = num+1;
        end
    end
    plv=sum(plv_hat)/num;
%     'done'
end