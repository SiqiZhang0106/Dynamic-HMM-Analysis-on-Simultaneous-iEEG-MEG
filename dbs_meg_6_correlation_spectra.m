clc
clear
load('C:\Users\szhang\Desktop\Inside_Cluster_14sub\timelock_11sub.mat');%iEEG time series of all subjects
hmm_outfile1 = 'C:\home\Data\new_embedded_HMM_results\embedded_hmm_data_11SUB.mat';
load( hmm_outfile1 ,'R')
hmm_outfile2 = 'C:\home\Data\new_embedded_HMM_results\embedded_HMM_11SUB_K5.mat';
load( hmm_outfile2 ,'Gamma','T')
% account for delay embedding in state gammas
pad_options = struct;
pad_options.embeddedlags = -7:7;
Gamma = padGamma(Gamma, T, pad_options);
clear T
Fs = 250;
N = size(R,1);
for ind = 1:N
    disp(ind);
    
    gamma = Gamma(R(ind,1):R(ind,2),:);
    TL = timelock{ind};
    TL=TL';
    if size(gamma,1) ~= size(TL,1)
        warning('The size of data and gamma do not match');
    end
    if mod(size(gamma,1),Fs)~=0;
        tlen = fix(size(gamma,1)/Fs);
        gamma = gamma(1:Fs*tlen,:);
        TL = TL(1:Fs*tlen,:);
    end
    
    K=5;
    ndim = size(TL,2);
    for e = 1:ndim  % electrodes number
        x = TL(:,e);
        [sp,f,t] = pspectrum(x,Fs,'spectrogram','FrequencyLimits',[1,125],'TimeResolution',1);
        %         images(t,f,P)
        %         pspectrum(x,Fs,'spectrogram','FrequencyLimits',[1,40],'FrequencyResolution',25);
        t = single(t.*250);
        cor = zeros(size(sp,1),K);
        for k = 1:K
            g = gamma(t(:,1),k);
            for i = 1:size(sp,1)
                [r,p] = corrcoef(sp(i,:)',g);
                cor(i,k) = r(1,2);
                p(i,k) = p(1,2);
            end
        end
        plot(f,cor)
        %clear cor
        xlabel('frequency','FontSize',12);
        ylabel('correlation','FontSize',12);
        set(gca,'YLim',[-0.2 0.2]);
        legend('state1','state2','state3','state4','state5');
    end
end
