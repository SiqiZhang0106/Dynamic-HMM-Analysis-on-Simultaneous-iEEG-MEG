function dbs_meg_hmm_prepare(subjects)

p = parcellation( 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );

%Main loop through subjects and sessions
fD = {};
TS=[];
for s = 1:length(subjects)
    
    initials = subjects{s};
 
    cd('C:\home\Code\Shared');
    [files, seq, root, details] = dbs_subjects(initials, 0);
    
    cd(fullfile(root, 'SPMhmm'));
    
    files = cellstr(spm_select('FPList','.',['^' initials '.*\.mat$']));
    
    for f = 1:numel(files)
        D = spm_eeg_load(char(files{f}));
        
        
        
        % Downsample and copy, or just copy
        if D.fsample > 250
            D = spm_eeg_downsample(struct('D',D,'fsample_new',250));
        end
        
        % Apply a 1-45Hz passband filter
        D = osl_filter(D,[1 45],'prefix','');
                
        % Apply automatric bad segment detection
        D = osl_detect_artefacts(D,'badchannels',false);
        
        % Run ICA artefact detection. This will automatically reject components
        % which have correlations larger than .5 with either of the artefact
        % channels.
        D = osl_africa(D, 'precompute_topos', false,'used_maxfilter',true);
       % D = osl_africa(D,'used_maxfilter',true);
        % Though the automatic correlations generally work well, we should be
        % careful to check for unsusual datasets and possibly manually correct the
        % automatic assessment. This is particularly important for relatively noisy
        % data or when analysing a new dataset for the first time.
        %
        
        
        % This is where manual Africa would go
        %D = D.montage('remove',1:D.montage('getnumber'));
        %D = osl_africa(D,'do_ident','manual');

        % Run LCMV Beamformer
        D = osl_inverse_model(D,p.template_coordinates,'pca_order',50);
        
        % Do parcellation
        D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');
        
        save(D);
        
        fD{end+1} = D;
        
        D = D.montage('switch', 0);save(D);
        %%
        S = [];
        S.D = D;
        S.band  = 'bandpass';
        S.freq  = [7 13];
        D = spm_eeg_filter(S);


        %% Beamformer sanity check - does BF data contain expected power distribution over space?
        
        
        D = D.montage('switch', 4);save(D);
        
        env_data = osl_envelope(D, 'downsample_env',1,'orthogonalise',false); %'filter',[8 13],
        netmat = ROInets.aec( env_data );
        
        mean_env = nanmean(env_data,2); % dim 1 might not be correct
        
        p.savenii(mean_env, spm_file(fullfile(D), 'ext', '.nii.gz')); % to view externally
        
        gunzip(spm_file(fullfile(D), 'ext', '.nii.gz'));
        
        spm_check_registration(spm_file(fullfile(D), 'ext', '.nii'));
        colormap(jet);
        
        print('-dtiff', '-r600', spm_file(files{f}, 'suffix', '_alphatopo', 'ext', '.tiff'));
        
        p.plot_network( netmat );
        view(0, 90);
        print('-dtiff', '-r600', spm_file(files{f}, 'suffix', '_alphanet', 'ext', '.tiff'));              
        
        delete(spm_file(fullfile(D), 'ext', '.nii.gz'));
        delete(D);    
        rmdir('osl_bf*', 's');
    end
end

