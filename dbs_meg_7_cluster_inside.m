clc
clear

subjects={'OXSH_EP11','OXSH_EP14','OXSH_EP3_RT','OXSH_EP3_LT','OXSH_EP3_RP','OXSH_EP3_LP','OXSH_EP5_LH','OXSH_EP5_RH','OXSH_EP6','OXSH_EP7','OXSH_EP8','OXSH_EP9','OXSH_EP21','OXSH_EP22','OXSH_EP23'};
%load( 'C:\home\Data\new_embedded_HMM_results\embedded_HMM_11SUB_K5_spectra.mat' )

K = 5;
%Chan_Locs = [];
%net_mean = zeros(39,K);

%MNI coordinates corresponding
parc = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');
mask = parc.template_mask;
dat2 = parc.template_coordinates;

mask2 = find(reshape(mask,prod(size(mask)),1)>0);
dat4 = zeros(prod(size(mask)),size(dat2,2));
try
    dat4(mask2,:) = dat2;
catch
    error('mask and nvoxels in the data incompatible');
end;
dat4 = reshape(dat4,size(mask,1),size(mask,2),size(mask,3),size(dat2,2));%size(dat4):23*27*23*3:mask*coordinates
dat4 = reshape(dat4,23*27,23,3);

p = 0;
for s=1:length(subjects)
    initials = subjects{s};
    cd('C:\home\Code\Shared');
    [~, sequence, root, details] = dbs_subjects_shanghai(initials, 0);
    Chan_Loc = [];
    Chan_Loc = details.chanpos;

    cd(fullfile(root, 'SPMhmm'));
    files = cellstr(spm_select('FPList','.',['^d' initials '.*\.mat$']));

    for f = 1:numel(files)
        p = p+1;
        for kk = 1:K
            %tmp = squeeze( mean(mean(abs(psd(:,kk,1:39,:,:)),3),1) );
            tmp = squeeze( mean(abs(psd(p,kk,1:39,:,:)),3) );
            net_mean = zscore(diag(tmp));


            sampleThr = 0.95;
            sampleThr = quantile(net_mean(:),sampleThr);
            candidates = (net_mean>=sampleThr);
            net_mean = net_mean.*candidates;

            data = parc.to_vol(net_mean);

            data = reshape(data,23*27,23);

            [m1,m2] = find(data~=0);

            Clus_Loc = [];%(size(m1,1),3);

            for j = 1:size(m1,1)

                Clus_Loc = cat(1,Clus_Loc,reshape(dat4(m1(j),m2(j),:),1,3));
                for i = 1:size(Chan_Loc,1)
                    Dis_Chans_Clus(i,j)=sqrt(sum((Chan_Loc(i,:)-Clus_Loc(j,:)).^2));
                end
            end
            %Dis_Chans_Clus(find(Dis_Chans_Clus<=8))
            Dis_Chans_Clus_K{p,kk} = Dis_Chans_Clus;
            clear Dis_Chans_Clus;
        end
    end
end
save Dis_Chans_Clus_K Dis_Chans_Clus_K

load Dis_Chans_Clus_K.mat;
%load SPEC_m1_files;
[F,K] = size(Dis_Chans_Clus_K);
%if the location of electrode f was close to the cluster k, the value of
%M(k,m) could be 1, otherwise, could be 0.
for f = 1:F
    %SPEC = SPEC_m1{f};
    M = zeros(K,size(Dis_Chans_Clus_K{f,1},1));
    for k = 1:K
        DI = Dis_Chans_Clus_K{f,k};
        [m,n] = find(DI<=8);
        m = unique(m);
%         SP = SPEC(k,:,:);
%         SP = reshape(SP,size(SP,2),size(SP,3));%size(SP,2):spectrum resolution;size(SP,3):electrodes mumber
        if length(m)~=0
            M(k,m) = 1;
%             for i = 1:length(m)
%                 ins_spec (:,i) = SP(:,m(i));
%             end
%             clu_ins_spec{f,k} = ins_spec;
%             clear ins_spec
        end
    end
    MM{f} = M;
    clear M
end

save MM_index MM

