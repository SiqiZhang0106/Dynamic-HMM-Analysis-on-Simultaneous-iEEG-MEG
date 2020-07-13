function dbs_meg_hmm_TDE(subjects)

data = [];            % HMM ready dataset
T = [];               % Length of continuous good segments
R = [];               % Indices of single run within data
P = [];
runlen = [];          % Length of run per good segment
B = {};      % Indices of bad samples per session
index_g = {};%Indices of good samples

%Main loop through subjects and sessions
ind = 0;

for s = 1:length(subjects)
    
    initials = subjects{s};
    
    cd('C:\home\Code\Shared');
    [~, seq, root, details] = dbs_subjects(initials, 0);
    
    cd(fullfile(root, 'SPMhmm'));
    
    files = cellstr(spm_select('FPList','.',['^d' initials '.*\.mat$']));
    
    for f = 1:numel(files)
        
        ind = ind+1;
        
        D = spm_eeg_load(char(files{f}));
        
        D = D.montage('switch', 0);save(D);

        S = [];
        S.D = D;
        S.band  = 'bandpass';
        S.freq  = [1 40];
        D = spm_eeg_filter(S);
        
        D = D.montage('switch', 5);save(D);
        
        % get data and orthogonalise
        dat = D(:,:,1);
        datorig = dat;
        
        figure(1)
        plot(zscore(std(datorig,[],1)));
        bs = ~good_samples( D ) | any(isnan(dat));
        
        %remove the samples as outliers in the distribution
        inds1 = zscore(std(datorig,[],1))>6;
        bs2 = any( cat(1,inds1,bs) );

        bs3 = movmax(bs2,250);
        dat(:,bs3) = 0;
        
        % plot(zscore(std(dat,[],1)));hold on
        
        dat = ROInets.remove_source_leakage(dat, 'symmetric');

        %-------------------------------
        % Get badsamples
        runlen(end+1) = size(dat,2);
        
        bs = ~good_samples( D ) | any(isnan(dat)) | bs3 ;
        
        nsamples = size(dat,2);
        
        %remove the samples in the begining and in the end
        bs(:,nsamples-124:nsamples)=1;
        bs(:,1:125)=1;
        
        % find single good samples - bug when we have consecutive bad segments
        xx = find(diff(diff(bs)) == 2)+1;
        if ~isempty(xx)
            bs(xx) = 1;
        end
        
        % store bad samples
        B{end+1} = find(bs);
        
        % indices of good samples
        good_inds = setdiff(1:runlen(end),B{end});
        
        % remove bad samples,
        % replace with >> a = zeros(44,size(D,2))*nan;a(:,inds) = dat;
        
        if any(bs)
            
            t_good = ~bs;
            db = find(diff([0; t_good(:); 0]));
            onset = db(1:2:end);
            offset = db(2:2:end);
            t = offset-onset;
            
            %for each sequence, if its length is less than 125, it will be
            %removed
            for i=1:size(t,1)
                if t(i)<125
                    t_good(:,onset(i):offset(i)-1)=0;
                end
            end
            t(t<125)=[];
            onset(t<125)=[];
            offset(t<125)=[];
            dat = dat(:,t_good);
            % sanity check
            if size(dat,2) ~= sum(t)
                disp('Mismatch between Data and T!!');
            end
        else
            t = size(dat,2);
        end
        %index_g{ind} = t_good;
        %--------------------------------
        % Store info
        
        
        
        offset = sum(T);
        
        R = cat(1,R,[offset+1 offset+size(dat,2)]);
        
        
        T = cat(1,T,t);
        
        
        P = cat(1, P, [0*t+s 0*t+f]);
        
        
        data = cat(2, data, dat);
        %check the data after processing if you want
        %          figure(1)
        %          plot( std(dat,[],1))

    end
    
end
%save index_g_10sub index_g

outdir = fullfile(root, '..', 'new_embedded_HMM_results');
%spm_mkdir(outdir);
outfile = fullfile(outdir, 'embedded_hmm_data_11SUB.mat' );
save( outfile, 'data', 'R', 'T', 'B', 'runlen', '-v7.3' );
%%
run_flip = true;
if run_flip
    % need to compute erp above
%     x = squeeze(nanmean(erf,3));
    x = data; %39*times sries of all subs
    T2 = repmat(size(x,2),1,1);

    options_sf = struct();
    options_sf.maxlag = 5;
    options_sf.noruns = 20;
    options_sf.nbatch = 3;
    options_sf.verbose = 1;
    flips = findflip(x',T2',options_sf);

    T3 = [R(1,2) sum(diff(R(:,2)))]; % Single scan sessions
    data = flipdata(data',T3',flips);
end



%%
% Prepare options structure
options = struct();
options.verbose = 1;

% These options specify the data and preprocessing that hmmmar might perform. Further options are discussed here
options.onpower = 0;
options.standardise = 0;
options.Fs = 250;

% Here we specify the HMM parameters
options.K = 5;  	         % The number of states to infer
options.order = 0; 	         % The lag used, this is only relevant when using MAR observations
options.zeromean = 1; 	     % We do not want to model the mean, so zeromean is set on
options.covtype = 'full';    % We want to model the full covariance matrix
options.embeddedlags = -7:7; % 15 lags are used from -7 to 7
options.pca = 39*4;          % The PCA dimensionality reduction is 4 times the number of ROIs

% These options specify parameters relevant for the Stochastic inference. They
% may be omitted to run a standard inference, but this will greatly increase
% the memory and CPU demands during processing. A detailed description of the
% Stochastic options and their usage can be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#stochastic
options.BIGNinitbatch = 15;
options.BIGNbatch = 15;
options.BIGtol = 1e-7;
options.BIGcyc = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay = 5;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;

% The following loop performs the main HMM inference. We start by
% estimating a 6 state HMM as used in the manuscript.
states_to_infer = [5];

% Optionally, we can explore a wider range of values for K by looping through
% several values. This can be done by uncommenting the line below.
% Warning: This is likely to be extremely time-consuming to infer!

%states_to_infer = 2:2:12; % uncomment this line to explore different numbers of states

% The HMM inference is repeated a number of times and the results based on
% the iteration with the lowest free energy. Note that this can be
% extremely time-consuming for large datasets. For a quick exploration of
% results, nrepeats can be set to a smaller value or even 1. The full inference
% is run over 10 repeats.
nrepeats = 1;
%data=data';
for kk = states_to_infer
    best_freeenergy = nan;
    options.K = kk;
    
    for irep = 1:nrepeats
        % Run the HMM, note we only store a subset of the outputs
        % more details can be found here: https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimation
        [hmm_iter, Gamma_iter, ~, vpath_iter, ~, ~, ~, ~, fehist] = hmmmar (data,T',options);
        
        if isnan(best_freeenergy) || fehist(end) < best_freeenergy
            hmm = hmm_iter;
            Gamma = Gamma_iter;
            vpath = vpath_iter;
        end
    end
    % Save the HMM outputs
    hmm_outfile = fullfile(outdir, sprintf('embedded_HMM_11SUB_K%d',options.K));
    save( hmm_outfile ,'hmm','Gamma','vpath','T', 'P', 'subjects');
end

%% Load 6 state HMM
% The state-wise spectra will be computed for 6 states
%data=data';
nstates = 5;
hmm_outfile = 'C:\home\Data\new_embedded_HMM_results\embedded_HMM_11SUB_K5.mat';
load( hmm_outfile ,'hmm','Gamma','vpath','T')

%% Statewise-Spectra
% Next we estimate the state-wise multitaper for each subject and each parcel.

% account for delay embedding in state gammas
pad_options = struct;
pad_options.embeddedlags = -7:7;
Gamma = padGamma(Gamma, T, pad_options);

if size(Gamma,1) ~= size(data,1)
    warning('The size of data and Gamma do not match');
end

% These options specify the how the spectra will be computed. Full details can
% be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#spectra

spec_options = struct();
spec_options.fpass = [1 40];
spec_options.p = 0; % no confidence intervals
spec_options.to_do = [1 0]; % no pdc
spec_options.win = 256;
spec_options.embeddedlags = -7:7;
spec_options.Fs = 250;

N = size(R,1);
psd = zeros(N,nstates,39,39,39);%(Nf x ndim x ndim)
coh = zeros(N,nstates,39,39,39);
for ind = 1:N
    disp(ind);
    subj_data = data(R(ind,1):R(ind,2),:);
    
    fit = hmmspectramt(subj_data,R(ind,2)-R(ind,1),Gamma(R(ind,1):R(ind,2),:),spec_options);
    for jj = 1:nstates
        psd(ind,jj,:,:,:) = fit.state(jj).psd;
        coh(ind,jj,:,:,:) = fit.state(jj).coh;
    end
    clear fit subj_data
end

% Save the MT outputs
mt_outfile = fullfile( outdir, sprintf('embedded_HMM_11SUB_K%d_spectra',nstates));
save( mt_outfile ,'psd','coh')

