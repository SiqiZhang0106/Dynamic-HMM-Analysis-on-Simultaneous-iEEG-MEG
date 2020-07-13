function fD = dbs_meg_epilepsy_prepare_spm12(initials, drug, task)

druglbl = {'off', 'on'};

keep = 0;
tsss = 0;

if nargin<3
    task = 'R';
end

if tsss
    prefix = 'tsss_';
else
    prefix = '';
end

try
    [files, seq, root, details] = dbs_subjects(initials, drug);
catch
    D = [];
    return
end
%%
spm_mkdir(fullfile(root, 'hmm'))

cd(fullfile(root, 'hmm'));


if details.cont_head_loc
    [alldewar, alldist, allsens, allfid] = dbs_meg_headloc(files);
    
    [alldewar, alldist, allsens, allfid] = dbs_meg_headloc(files, alldewar);
else
    alldewar = [];
    alldist = [];
end

fD = {};
%%
for f = 1:numel(files)
    if ~any(f == strmatch(task, seq));
        continue;
    end
    
    
    S = [];
    S.dataset = files{f};
    S.channels = details.chanset;
    S.checkboundary = 0;
    S.saveorigheader = 1;
    S.conditionlabels = seq{f};
    S.mode = 'continuous';
    S.outfile = [task num2str(f) '_' spm_file(S.dataset, 'basename')];
    
    if details.brainamp
        if details.hamburg
            S.ref2 = 'UADC004';
        elseif isfield(details, 'megnoise')
            S.ref2 = details.megnoise;
        end
        
        D = dbs_meg_brainamp_preproc(S);
    else
        D = spm_eeg_convert(S);
    end
    
    D = sensors(D, 'EEG', []);
        
    if details.berlin        
        D = fiducials(D,  ft_convert_units(ft_read_headshape(spm_file(S.dataset, 'filename', details.markers(f).files{1}))));
        
        hdr = ft_read_header(S.dataset);
        
        event     = ft_read_event(S.dataset, 'detectflank', 'down', 'trigindx', ...
            spm_match_str({'EEG157', 'EEG159', 'EEG160'}, hdr.label), 'threshold', 5e-3);
        
        save(D);
    else
        event    = ft_read_event(fullfile(D), 'detectflank', 'both');
    end
    
    
    eventdata = zeros(1, D.nsamples);
    trigchanind = D.indchannel(details.eventchannel);
    
    if ~isempty(event)
        trigind  = find(strcmp(details.eventtype, {event.type}));
        eventdata([event(trigind).sample]) = 1;
    end
    
    D = chanlabels(D, trigchanind, 'event');
    
    D(D.indchannel('event'), :) = eventdata;        
    
    save(D);
    
    if details.cont_head_loc
        D = sensors(D, 'MEG', allsens(f));
        D = fiducials(D, allfid(f));
        
        S = [];
        S.D = D;
        if ~isempty(alldewar)
            S.valid_fid = squeeze(alldewar(:, :, f));
        end
        D = spm_eeg_fix_ctf_headloc(S);
    end
    
    if ~isempty(details.montage)
        montage = details.montage;
        
        S = [];
        S.D = D;
        
        if details.oxford
            ind = [];
            
            neeg = length(strmatch('EEG', D.chanlabels));
            if neeg < 8
                ind = [ind; strmatch('Cz', montage.labelnew, 'exact')];
            end
            
            if neeg < 5
                ind = [ind; strmatch('HEOG', montage.labelnew, 'exact')];
                ind = [ind; strmatch('VEOG', montage.labelnew, 'exact')];
            end
            
        elseif details.berlin
            montage.tra(end+1, end+1) = 1;
            montage.labelorg{end+1} = 'event';
            montage.labelnew{end+1} = 'event';
        elseif details.hamburg
        else
            ind = [];
            if isempty(strmatch('HLC', D.chanlabels))
                ind = [ind; strmatch('HLC', montage.labelnew)];
            end
            
            neeg = length(strmatch('EEG', D.chanlabels));
            if neeg < 16
                ind = [ind; strmatch('Cz', montage.labelnew, 'exact')];
            end
            
            if neeg < 13
                ind = [ind; strmatch('HEOG', montage.labelnew, 'exact')];
                ind = [ind; strmatch('VEOG', montage.labelnew, 'exact')];
            end
            
            if neeg == 0
                ind = [ind; strmatch('EMG', montage.labelnew)];
                ind = [ind; strmatch('LFP', montage.labelnew)];
            end
            
            if ~isempty(ind)
                montage.labelnew(ind) = [];
                montage.tra(ind, :) = [];
            end
        end
        
        S.montage = montage;
        S.keepothers = 0;
        D = spm_eeg_montage(S);
        
        if ~keep, delete(S.D);  end
        
    end
    
    D = chantype(D, D.indchannel(details.chan), 'LFP');
    
    D = chantype(D, D.indchannel(ft_senslabel('eeg1005')), 'PHYS');
    
    if details.hamburg
        D = chantype(D, D.indchannel(details.eegchan), 'EEG');
    end
    
    
    save(D);
    %%
    %D = chantype(D, D.indchannel('LFP_L0R0'), 'PHYS');
    %%
    
    S = [];
    S.D = D;
    S.mode = 'mark';
    S.badchanthresh = 0.8;
    if details.oxford
        S.methods(1).channels = {'MEGMAG'};
    else
        S.methods(1).channels = {'MEGGRAD'};
    end
    S.methods(1).fun = 'flat';
    if details.berlin
        S.methods(1).settings.threshold = flatthresh;
    else
        S.methods(1).settings.threshold = 1e-010;
    end
    S.methods(1).settings.seqlength = 10;
    S.methods(2).channels = {'MEGPLANAR'};
    S.methods(2).fun = 'flat';
    S.methods(2).settings.threshold = 0.1;
    S.methods(2).settings.seqlength = 10;
    if details.oxford
        S.methods(3).channels = {'MEGMAG'};
    else
        S.methods(3).channels = {'MEGGRAD'};
    end
    S.methods(3).fun = 'jump';
    if details.oxford
        S.methods(3).settings.threshold = 50000;
    elseif details.berlin
        S.methods(3).settings.threshold = 1e4;
    else
        S.methods(3).settings.threshold = 20000;
    end
    S.methods(3).settings.excwin = 2000;
    S.methods(4).channels = {'MEGPLANAR'};
    S.methods(4).fun = 'jump';
    S.methods(4).settings.threshold = 5000;
    S.methods(4).settings.excwin = 2000;
    
    if details.berlin
        S.methods(5).channels = {'MEG'};
        S.methods(5).fun = 'threshchan';
        S.methods(5).settings.threshold = ampthresh;
        S.methods(5).settings.excwin = 1000;
    end
    D = spm_eeg_artefact(S);        
    
    if ~keep
        delete(S.D);
    end
    
    spikes_file= spm_file(files{f}, 'path', fullfile(spm_file(files{f}, 'path'), '..', 'spikes'), 'ext', '.evt');
    if exist(spikes_file, 'file')
        fid = fopen(spikes_file);
        textscan(fid, '%s\t%s', 1);
        ep = textscan(fid, '%u%s');
        fclose(fid);
        
        times = 1e-6*double(ep{1});
        lbl   = ep{2};
        
        ev = struct([]);
        evind = 1;
        for i = 1:numel(lbl)       
            if isequal(lbl{i}, 'spike')
                ev(evind).type     = 'artefact_manual';
                ev(evind).value    = 'all';
                ev(evind).time     = max(D.timeonset, times(i)-1);
                ev(evind).duration = 2;
                evind = evind+1;
            elseif isequal(lbl{i}, 'start')
                ev(evind).type     = 'artefact_manual';
                ev(evind).value    = 'all';
                ev(evind).time     = times(i);
                
                if i==numel(lbl) || ~isequal(lbl{i+1}, 'end')
                    error('Error parsing the spike specification file');
                end
                
                ev(evind).duration = times(i+1)-times(i);
                evind = evind+1;
            end
        end
        
        D = events(D, 1, spm_cat_struct(events(D, 1), ev));
    end
    
    save(D);
    
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 1;
    S.dir = 'twopass';
    S.order = 5;
    D = spm_eeg_filter(S);
    
    
    if ~keep, delete(S.D);  end
    
    if D.fsample > 250
        
        S = [];
        S.D = D;
        S.fsample_new = 250;
        
        D = spm_eeg_downsample(S);
        
        if ~keep, delete(S.D);  end
    end
    
    if details.berlin
        [ampthresh, flatthresh] = berlin_gain(spm_file(files{f}, 'ext', 'gain.txt'));
    end
    
    
    if  tsss
        S = [];
        S.D = D;
        S.magscale   = 100;
        S.tsss       = 1;
        S.t_window   = 5;
        S.xspace     = 0;
        %S.Dref       = 'C:\home\Data\Rafal\asens.mat';
        D = tsss_spm_enm(S);
        
        D = badchannels(D, D.indchantype('MEGANY'), 0);save(D);
        
        if ~keep, delete(S.D);  end
        
        prefix = 'tsss';
    end
    
    
    D.initials = initials;
    
    if details.neuromag
        fids = D.fiducials;
        [~, sel] = spm_match_str({'Nasion', 'LPA', 'RPA'}, fids.fid.label);
        fids.fid.pnt = fids.fid.pnt(sel, :);
        fids.fid.label = {'nas'; 'lpa'; 'rpa'};
        
        D = fiducials(D, fids);
    end
    
    save(D);
            
    cd('C:\home\Code\Shared');
    D = dbs_meg_headmodelling(D);
    D = D.move([prefix initials '_' druglbl{drug+1} '_' task num2str(f)]);
    
    fD{f} = D;
end

%%


