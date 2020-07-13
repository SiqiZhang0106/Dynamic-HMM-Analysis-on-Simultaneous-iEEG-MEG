% 
config = utils.get_studydetails;

% Define colours to use in state plots
set1_cols = utils.set1_cols;

% Define sample rate
sample_rate = 250;

%% Load in results from envelope data

% Meta data
method = 'embedded';
K = 5;

hmm_outfile1 = 'C:\home\Data\new_embedded_HMM_results\embedded_hmm_data_11SUB.mat';
load( hmm_outfile1 ,'R')
hmm_outfile2 = 'C:\home\Data\new_embedded_HMM_results\embedded_HMM_11SUB_K5.mat';
load( hmm_outfile2 ,'Gamma','T')

% account for delay embedding in state gammas
pad_options.embeddedlags = -7:7;
Gamma = padGamma(Gamma, T, pad_options);

% Temporal statistics

%Here we compute the global temporal statistics from the Time-Delay-Embedded
%HMM. These are computed per subject and visualised as violin plots

scan_T = [R(1,2) diff(R(:,2))']; % Indexing individual scan sessions
subj_T = sum(reshape(scan_T,5,[])); % Indexing individal subjects

%Compute temporal stats
subj_T=T;
% Fractional Occupancy is the proportion of time spent in each state
FO = getFractionalOccupancy( Gamma, subj_T, 2);
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma, subj_T, []);
ITmerged = cellfun(@mean,IT);clear IT
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma, subj_T, []);
LTmerged = cellfun(@mean,LT); clear LT


% Plot temporal stats
fontsize = 18;

figure;subplot(111);
distributionPlot(FO,'showMM',2,'color',{set1_cols{1:size(FO,2)}});
set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
%print([savebase '_temporalstats_FO'],'-depsc')

figure;subplot(111);
distributionPlot(LTmerged ./ sample_rate * 1000,'showMM',2,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;
set(gca,'YLim',[0 300],'FontSize',fontsize,'FontSize',fontsize)
%print([savebase '_temporalstats_LT'],'-depsc')

figure;subplot(111);
distributionPlot(ITmerged ./ sample_rate,'showMM',2,'color',{set1_cols{1:size(FO,2)}})
title('Interval Times');xlabel('State');ylabel('Time (secs)');grid on
set(gca,'YLim',[0 3],'FontSize',fontsize)
%print([savebase '_temporalstats_IT'],'-depsc')