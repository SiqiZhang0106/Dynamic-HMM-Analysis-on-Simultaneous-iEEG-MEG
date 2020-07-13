clc
clear

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
    fitsub{ind} = fit;
    clear fit subj_data
end

options.Ncomp = 3;
options.Method = 'NNMF';
options.Base = 'psd';

[sp_fit,sp_fit_group,sp_profiles] = spectdecompose(fitsub,options)
figure(2)
plot(sp_profiles)
for ind = 1:N
    for jj = 1:nstates
        psd(ind,jj,:,:,:) = sp_fit{1,ind}.state(jj).psd;
        coh(ind,jj,:,:,:) = sp_fit{1,ind}.state(jj).coh;
    end
end
save sp sp_fit sp_fit_group sp_profiles

for jj = 1:5
        psd_group(jj,:,:,:) = sp_fit_group.state(jj).psd;
        coh_group(jj,:,:,:) = sp_fit_group.state(jj).coh;
    end
save psd_group psd_group