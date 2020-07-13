%<<<<<<< .mine
clc
clear
load ('H:\FIL\home\Data\new_embedded_HMM_results\embedded_HMM_11SUB_K5.mat');
% =======
% load envelope_HMM_K6.mat
%>>>>>>> .r174

options = [];
options.Fs = 250;
%FO = getFractionalOccupancy( Gamma, subj_T, 2);
FO = getFractionalOccupancy(Gamma,T,2);

subind = unique(P(:, 1));

figure;
for i = 1:length(subind)
    subplot(length(subind), 1, i);
    cFO = sum(FO(P(:, 1)==i, :).*repmat(T(P(:, 1)==i)./sum(T(P(:, 1)==i)), 1, size(FO, 2)), 1);
    bar(cFO);
    title(subjects{subind(i)}, 'Interpreter', 'none');
end

[nsample nstate]=size(Gamma);
vmatrix = zeros(nsample,nstate);
for ii = 1:nstate
    vmatrix(vpath==ii,ii) = ii;
end
figure;
for i = 1:size(vmatrix, 2)
    subplot(size(vmatrix, 2), 1, i);
    imagesc(vmatrix(:,i)');figure(gcf);
    temp=[0 6];
    caxis(temp)
end
%%
% %% Mean Activation Maps

parc = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');

%% Broadband power plots
load( 'H:\FIL\home\Data\new_embedded_HMM_results\embedded_HMM_11SUB_K5_spectra.mat' )

net_mean = zeros(39,size(psd,2));
for f = size(psd,1)
    for kk = 1:size(psd,2)
        tmp = squeeze( mean(mean(abs(psd(:,kk,1:39,:,:)),3),1) );
        net_mean(:,kk) = zscore(diag(tmp));
    end
end


% visualise state in OSLEYES
parc.osleyes(net_mean);

% Optionally save a nifti of the results, these are used to generate the
% figures in the paper via HCP Workbench
parc.savenii( net_mean, 'res_meanactivations');
%%
gunzip('*.gz');
%%
file = 'res_meanactivations.nii'; %'res_meanactivations.nii';%

v = spm_vol(file);
M = gifti( fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
figure;
for i = 1:numel(v)
    AX = subplot(numel(v), 3, 3*(i-1)+1);
    spm_mesh_render('Disp',M, 'parent', AX);
    spm_mesh_render('Overlay',AX,spm_file(file, 'number', i));
    spm_mesh_render('ColourMap',AX,jet);
    view(0, 90);
    light('Position', [0 1 1]);
    
    AX = subplot(numel(v), 3, 3*(i-1)+2);
    spm_mesh_render('Disp',M, 'parent', AX);
    spm_mesh_render('Overlay',AX,spm_file(file, 'number', i));
    spm_mesh_render('ColourMap',AX,jet)
    view(90, 0);
    light('Position', [1 0 1]);
    
    AX = subplot(numel(v), 3, 3*i);
    spm_mesh_render('Disp',M, 'parent', AX);
    spm_mesh_render('Overlay',AX,spm_file(file, 'number', i));
    spm_mesh_render('ColourMap',AX,jet);
    view(-90, 0);
    light('Position', [1 1 0]);
    colorbar
end
