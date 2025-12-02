function run_batch_PAF(INDIR, OUTDIR, alpha_range)

% RUN_BATCH_PAF
% Minimal working version; batch compute Peak Alpha Frequency (PAF)
% using Welch PSD + Center of Gravity (CoG) for each participant.
% 
% INPUTS:
% INDIR - folder containing subject .mat files
%           each file should have variables (data, fs)
% OUTDIR - folder where PAF results will be saved
% alpha_range - [low high] alpha band in Hz, e.g., [8 13]
%
% EXAMPLE USAGE:
%   run_batch_PAF('D:\EEG\Cleaned\', 'D:\EEG\PAF\', [8 13])
%
% NOTE: 
% Later, can replace INDIR/OUTDIR with paths from a CONFIG.m file

% -------
% 0. Handle Inputs
% -------
    if nargin < 3
        alpha_range = [8 13]; % default alpha band
    end

    if ~exist(INDIR, 'dir')
        error('Input directory does not exist: %s', INDIR);
    end

    if ~exist(OUTDIR, 'dir')
        mkdir(OUTDIR);
    end

% -------
% 1. Find subject files
% -------
    files = dir(fullfile(INDIR, 'sub-*.mat'));

    if isempty(files)
        error('No files matching "sub-*.mat" found in %s.', INDIR);
    end

    % Preallocate a simple results struct for group-level summary
    results = struct('subjID', {}, ...
        'PAF_mean', {}, ...
        'PAF_channels', {}, ...
        'alpha_range', {});

% -------
% 2. Loop over subjects
% -------
    for iSub = 1:numel(files)

        inFile = fullfile(files(iSub).folder, files(iSub).name);
        [~, subjID, ~] = fileparts(inFile); % e.g., 'sub-01'

        fprintf('n\=== Processing %s ===\n', subjID);

        % -------
        % 2.1 Load Data
        % -------
        % We assume each .mat file contains
        % data: samples x channels
        % fs: sampling rate in Hz
        S = load(inFIle);

        if ~isfield(S, 'data') || ~isfield(S, 'fs')
            error('File %s must contain variables "data" and "fs".', inFile);
        end

        data = S.data;
        fs = S.fs;

        % Sanity check: data dimensions
        if size(data, 1) > size(data, 2)
            % If you data is channels x samples, flip it
            warning('%s: data appears to be channels x samples, transposing...', subjID);
            data = data.';
        end

        % -------
        % 2.2 Compute PAF per channel
        % -------
        % This helper does:
        % - Welch PSD per channel
        % - Restrict to alpha band
        % - Compute CoG-based PAF per channel
        [PAF_chan, f, pxx] = compute_PAF_CoG_from_data(data, fs, alpha_range);

        % For now, just take the mean across channels as a simple subject
        % PAF
        % Later, replace this with the actual method
        PAF_mean = mean(PAF_chan, 'omitnan');

        fprintf('  Mean PAF for %s: %2.f Hz\n', subjID, PAF_mean);

        % -------
        % 2.3 Save per-subject result
        % -------
        outFile = fullfile(OUTDIR, [subjID '_PAT.mat']);

        % Save the essentials: pxx/f are optional but nice for debugging
        save(outFile, 'subjID', 'PAF_chan', 'PAF_mean',"pxx", ...
            'alpha_range', 'f', 'pxx', '-v7.3');

        % -------
        % 2.4 Store in group results struct
        % -------
        results(iSub).subjID = subjID;
        results(iSub).PAF_mean = PAF_mean;
        results(iSub).PAF_channels = PAF_chan;
        results(iSub).alpha_range = alpha_range;

    end

    % -------
    % 3. Save group level results
    % -------
    groupFile = fullfile(OUTDIR, 'PAF_group_results.mat');
    save(groupFile, 'results', '-v7.3');

    fprintf('\n=== Done. Group results saved to: %s ===\n', groupFile);

end

%% ================================================
% Helper function: Welch PSD + CoG-based PAF per channel
% =================================================
function [PAF, f, pxx] = compute_PAF_CoG_from_data(data, fs, alpha_range)
% COMPUTE_PAF_COG_FROM_DATA
% Compute channel-wise PAF using:
% 1) Welch power spectral density (PSD)
% 2) Center of Gravity (CoG) within alpha_range
%
% INPUTS:
% data - samples x channels (double)
% fs - sampling rate in Hz
% alpha_range - [low high] alpha band in Hz (e.g., [8 13])
%
% OUTPUTS:
% PAF - 1 x channels CoG-based PAF (Hz)
% f - frequency vector from pwelch (Hz)
% pxx - PSD power, freqBins x channels

    if nargin < 3
        alpha_range = [8 12];
    end

    % -------
    % 1. Welch PSD parameters
    % -------
    % win: length of each segment in samples (here: 2 s)
    % nover: ovelap between segments (50%)
    % nfft: number of FFT points (controls frequency resolution)
    win = fs * 2; % 2 second windows
    nover = win / 2; % 50% overlap
    nfft = max(1024, 2^nextpow2(win)); % at least 1024

    % pwelch expects data as columns = channels, so we transpose
    % (our data is samples x channels)
    [pxx, f] = pwelch(data, win, nover, nfft, fs);

    % pxx is freqBins x channels

    % -------
    % 2. Restrict to alpha band
    % -------
    alpha_idx = f >= alpha_range(1) & f <= alpha_range(2);

    f_alpha = f(alpha_idx); % freq values within alpha band
    pxx_alpha = pxx(alpha_idx, :); % power within alpha band

    % -------
    % 3. CoG-based PAF per channel
    % -------
    % CoG formula:
    % PAF = sum(f * P) / sum(P)
    % 
    % Here:
    % f_alpha: nAlphaBins x 1
    % pxx_alpha: nAlphaBins x nChannels

    % We want a 1 x nChannels result:
    numerator = f_alpha' * pxx_alpha; % 1 x nChannels (sum(f * P))
    denominator = sum(pxx_alpha, 1); % 1 x nChannels (sum(P))

    PAF = numerator ./ denominator; % 1 x nChannels

end