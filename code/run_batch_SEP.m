function run_batch_SEP(INDIR, OUTDIR, epoch_win, baseline_win)
% RUN_BATCH_SEP
% Minimal working example for generating SEPs (somatosensory/ERP) across
% subjects
%
% INPUTS:
% data - samples x channels
% fs - sampling rate (Hz)
% event_samples - event indices (stim onset) in samples
% epoch_win - [t_pre t_post] window (s)
% baseline_win - [t_start t_end] baseline window (s)
%
% OUTPUTS:
% time - time vector (s)
% SEP - averaged ERP/SEP (samples x channels)
% epochs - 3D array (trials x samples x channels)

    % Convert time windows to sample offsets
    sample_offsets = round(epoch_win * fs); % e.g., [-50, 150]
    nSamples = sample_offsets(2) - sample_offsets(1) + 1;
    time = (sample_offsets(1):sample_offsets(2)) / fs;

    baseline_offsets = round(baseline_win * fs); % e.g., [-50, 0]

    nEvents = numel(event_samples);
    nChannels = size(data, 2);

    % Preallocate trials x samples x channels
    epochs = nan(nEvents, nSamples, nChannels);

    for t = 1:nEvents

        ev = event_samples(t);

        % Index range for this epoch
        idx = ev + (sample_offsets(1):sample_offsets(2));

        % Skip trials too close to edges
        if idx(1) < 1 || idx(end) > size(data, 1)
            continue;
        end

        % Extract epoch
        e = data(idz, :); % samples x channels

        % Baseline correction
        base_idx = (baseline_offsets(1):baseline_offsets(2)) - sample_offsets(1) + 1;
        base_idx(base_idx < 1) = 1;
        base_idx(base_idx > nSamples) = nSamples;

        baseline = mean(e(base_idx, :), 1);
        e_bc = e - baseline;

        epochs(t,:,:) = e_bc;
    end

    % Drop empty (NaN-Only) epochs
    good = -all(all(isnan(epochs),3),2);

    % Average across trials + ERP/SEP
end