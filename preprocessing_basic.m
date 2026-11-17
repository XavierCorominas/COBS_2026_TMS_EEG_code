%% Preprocessing pipeline for early iTEPs.

% DRCMR 2026 - YFS XCT - BrainStim Methods group

% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide a simple and fast eeg
% preprocessing focusing on the essential steps. The preprocessing might be
% sufficient to explore iTEPS in datasets non contaminated with large
% muscular (or other type of) artefacts.
% NOTE: Not every step here is needed for every dataset. Choose the steps,
% their order and their parameters to suit your own data.

% ---> STEPS:
% 0.   Load data
% 1.   Epoch
% 2.   Baseline Correct
% 3.   Remove bad trials and channels 
% 3.1  Identify bad trials and remove
% 3.2  Plot and identify bad channels and remove
% 3.3  (Optional) Interpolate removed channels
% 4.   (Optional) Remove TMS pulse artifact
% 5.   (Optional) Interpolate missing data
% 6.   Low-pass filtering
% 7.   Re-Baseline correction
% 8.   (Optional) Re-reference 
% 9.   Save

%% ADDPATHS 
clearvars -except cobsroot
if ~exist('cobsroot', 'var')
    cobsroot = '/mnt/projects/iTEPs/COBS_2026_data_example';
end

data_path = char(fullfile(cobsroot, 'data', 'examples'));
addpath(fullfile(cobsroot, 'dependencies', 'eeglab2025.0.0', 'FastICA_25'));

close all;
clc

fprintf('\n>> ADDPATHS (workspace cleared, data path added) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 0. LOAD DATA

% Open EEGLAB
eeglab;

% choose which file to load: 

% --- single pulse TMS ---
% name_dataset = 'X99880_sp_110rmt_block1.vhdr';


% --- sici ---
% name_dataset = 'X99880_sici_90amt_110rmt_block1.vhdr';

% --- lici ---
% name_dataset = 'X99880_lici_110rmt_110rmt_100msisi_block1.vhdr';

% 
[~, save_name] = fileparts(name_dataset);

% Load file
EEG = pop_loadbv(data_path, name_dataset);

% % --- right shift marker timing for 'sici' and 'lici' files---
% This makes TS pulse appear at at t=0 ms.
isi = [];
if contains(name_dataset, "sici",'IgnoreCase',true)
    isi = 2.1;  % sici
elseif contains(name_dataset, "100ms",'IgnoreCase',true)
    isi = 100;  % lici
end

if ~isempty(isi)
    allEvent = string({EEG.event.type}');
    markerLabelNorm = 'R  8';
    markerIdx = find(allEvent == markerLabelNorm);
    lat = [EEG.event(markerIdx).latency]';  % latency in samples 
    markershift = isi / 1000 * EEG.srate;
    changedLatency = num2cell(lat + markershift);
    [EEG.event(markerIdx).latency] = changedLatency{:};
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

eeg_setname = 'continuous eeg';
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
        'overwrite', 'on', 'setname', eeg_setname);  %[give each dataset a name]  
eeglab redraw


fprintf('\n>> Step 0 (Load data: %s) ended\n', name_dataset);
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 1. EPOCH
eeg_setname = 'epoch';
% Epoch around TMS pulses 'R  8' or 'S 11' markers
EEG = pop_epoch(EEG, {'R  8', 'S 11'}, [-1 1]); % sec, x-axis range
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname);
eeglab redraw;

fprintf('\n>> Step 1 (Epoch around TMS pulses, -1 to 1 s: %d epochs) ended\n', EEG.trials);
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 2. Baseline Correct
eeg_setname = 'baseline correction';
baselineW = [-210 -10];
if contains(name_dataset, "100ms",'IgnoreCase',true)
 baselineW = [-310 -110];
end

EEG = pop_rmbase( EEG, baselineW ,[]); 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

fprintf('\n>> Step 2 (Baseline correction) ended\n'); %[not Demeam anymore]
fprintf('%s\n\n', repmat('-', 1, 70));
%% Step 3. Remove bad trials, remove and interpolate bad channels
%% trials rejection - 3.1.1 (identify bad trials) --> Flag the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

fprintf('\n>> Step 3.1.1 (Identify bad trials: data displayed, flag trials and press REJECT) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));
 
%% trials rejection - 3.1.2 (remove bad trials)
eeg_setname = 'trials rejection';
if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

fprintf('\n>> Step 3.1.2 (Remove bad trials: %d removed, %d remaining) ended\n', numel(EEG.BadTr), EEG.trials);
fprintf('%s\n\n', repmat('-', 1, 70));
%% 3.1.3 Plot the data to inspect and identify bad channels

% Plot channels: hover over a line to see the channel number
plotTEP(EEG, [2.1 3.3 4.5 9 13 17 35 48 70 150 240], [-55 55]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-10  20], [2.1 3.3 4.5], 'iTEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 3.1.3 (Plot data to inspect and identify bad channels) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));


%% 3.4 Manual rejection and interpolation of channels.
eeg_setname = 'channel rejection and interpolate';
% Select bad channel numbers
badChannels = [22]; % example: channel 2

% Save original channel locations before removal
originalChanLocs = EEG.chanlocs;

% Remove bad channels
EEG = pop_select(EEG, 'nochannel', badChannels);

% Interpolate the removed channels 
EEG = pop_interp(EEG, originalChanLocs, 'spherical'); % or use 'invdist' or 'nearest'

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot channels: hover over a line to see the channel number
plotTEP(EEG, [2.1 3.3 4.5 9 13 17 35 48 70 150 240], [-55 55]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-10  20], [2.1 3.3 4.5], 'iTEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 3.4 (Manual rejection/interpolation of channels: removed [%s]) ended\n', num2str(badChannels));
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 4. (Optional) Remove TMS pulse artifact
eeg_setname = 'remove TMS pulse artifact';
% TMS artefact lengths in ms
TMSremoval_range = [-1 2]; % in ms
TMSremoval_range2 = [ ];

% sici and lici: TMSremoval_range needs to cover TS+CS pulse
if contains(name_dataset, "sici",'IgnoreCase',true)
    TMSremoval_range = [-3.1 2];  
elseif contains(name_dataset, "100ms",'IgnoreCase',true)
    TMSremoval_range2 = [-101 -98];  
end
EEG = pop_tesa_removedata( EEG, TMSremoval_range );


if ~isempty(TMSremoval_range2)
 EEG = pop_tesa_removedata( EEG, TMSremoval_range2 ); % lici
end 

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot channels: hover over a line to see the channel number
plotTEP(EEG, [2.1 3.3 4.5 9 13 17 35 48 70 150 240], [-55 55]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-10  20], [2.1 3.3 4.5], 'iTEP and scalp maps','maplimits', [-30 30]);
%}


fprintf('\n>> Step 4 (Remove TMS artifact: %g to %g ms) ended\n', TMSremoval_range(1), TMSremoval_range(2));
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 5. (Optional) Interpolate missing data
eeg_setname = 'interpolate missing data';
% Interpolate
EEG = pop_tesa_interpdata( EEG, 'cubic', [1 1] ); 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot channels: hover over a line to see the channel number
plotTEP(EEG, [2.1 3.3 4.5 9 13 17 35 48 70 150 240], [-55 55]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-10  20], [2.1 3.3 4.5], 'iTEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 5 (Interpolate missing data) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 6. Low pass filtering 
eeg_setname = 'low pass filter';
EEG = tesa_filtbutter(EEG, [], 2000, 2, 'lowpass');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

%Plot the data 
plotTEP(EEG, [2.1 3.3 4.5 9 13 17 35 48 70 150 240], [-55 55]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-10  20], [2.1 3.3 4.5], 'iTEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 6 (Filtering) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 7.Baseline correction
eeg_setname = 're baseline correction';
EEG = pop_rmbase( EEG, baselineW ,[]);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

%Plot the data 
plotTEP(EEG, [2.1 3.3 4.5 10 13 17 35 48 70 150 240], [-55 55]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-10  20], [2.1 3.3 4.5], 'iTEP and scalp maps','maplimits', [-30 30]);
%}
fprintf('\n>> Step 7 (Baseline correction, -210 to -10 ms) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 8. (Optional) REREFERENCING If necessary
eeg_setname = 're-ref';
EEG = pop_reref( EEG, []);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

fprintf('\n>> Rereferencing step ended (currently commented out: no change applied)\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% Step 9. SAVE DATASETS
eeg_setname = [save_name,'_clean.set'];
EEG = pop_saveset( EEG, 'filename',eeg_setname,'filepath',data_path);
 
fprintf('\n>> Step 8 (Save dataset: %s_clean_no_pulse_basic.set) ended\n', save_name);
fprintf('%s\n\n', repmat('-', 1, 70));
%% END


%% ******************* CODE PLOT FUNCTIONS
function fig = plotTEP(EEG, topoLat, mapLimits)
% PLOTTEP  Butterfly TEP plot (full + early zoom) with topographies.
%   Hover a trace to see the channel; click it to open channel details.
%
%   plotTEP(EEG)
%   plotTEP(EEG, [2.1 3.3 4.5 9 13 17 35 48 70 150 240])
%   plotTEP(EEG, topoLat, [-45 55])
%
%   EEG       - EEGLAB EEG structure (channels x time x trials)
%   topoLat   - latencies (ms) for the topographies (optional)
%   mapLimits - colour limits for the topographies (optional; default
%               symmetric around the max absolute value at those latencies)

if nargin < 2 || isempty(topoLat)
    topoLat = [2.1 3.3 4.5 9 13 17 35 48 70 150 240];
end

erp   = mean(EEG.data, 3);            % average over trials: channels x time
nTopo = numel(topoLat);

% Sample closest to each requested latency
topoIdx = zeros(1, nTopo);
for k = 1:nTopo
    [~, topoIdx(k)] = min(abs(EEG.times - topoLat(k)));
end

if nargin < 3 || isempty(mapLimits)
    m = max(max(abs(erp(:, topoIdx))));
    mapLimits = [-m m];
end

fig = figure('Color', 'w', 'Position', [10 5 1000 600]);

% Row 1: GMFP of the averaged TEP (SD across channels of the ERP).
% Computed AFTER averaging trials, so non-phase-locked noise cancels and the
% baseline sits close to 0. (GMFP is an SD, so it is always >= 0.)
gmfp = std(erp, 1, 1);                    % 1 x time
ax0 = subplot(4, nTopo, 1:nTopo);
plot(EEG.times, gmfp, 'b', 'LineWidth', 2); hold(ax0, 'on');
xline(ax0, topoLat, '--k', 'LineWidth', 1, 'HitTest', 'off');
xlim([-100 600]);
ylim([-1 50]);
xlabel('Time (ms)'); ylabel('GMFP (\muV)');
title(sprintf('GMFP of the average TEP (%d trials)', size(EEG.data, 3)));
set(ax0, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);

% Row 2: full time window
ax1 = subplot(4, nTopo, nTopo+1:2*nTopo);
h1 = plot(EEG.times, erp');
xlim([-100 600]); ylim([-80 80]);
xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
title('Hover to see channel, click to open channel details');
set(ax1, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);
hold(ax1, 'on');
xline(ax1, topoLat, '--k', 'LineWidth', 1, 'HitTest', 'off');

% Middle: zoom on -5 to 20 ms
ax2 = subplot(4, nTopo, 2*nTopo+1:3*nTopo);
h2 = plot(EEG.times, erp');
xlim([-5 20]); ylim([-80 80]);
xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
set(ax2, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);
hold(ax2, 'on');
xline(ax2, topoLat(topoLat >= -5 & topoLat <= 20), '--k', 'LineWidth', 1, 'HitTest', 'off');

% Bottom: topographies
for k = 1:nTopo
    axT = subplot(4, nTopo, 3*nTopo + k);
    topoplot(erp(:, topoIdx(k)), EEG.chanlocs, ...
        'maplimits', mapLimits, 'electrodes', 'on', 'style', 'map');
    title(sprintf('%g ms', EEG.times(topoIdx(k))), 'FontSize', 14, 'FontWeight', 'bold');
end

cb = colorbar(axT, 'Position', [0.93 0.08 0.012 0.15]);
cb.Label.String = '\muV';
cb.FontSize = 12;

% Hover labels + click action on both ERP plots
for ch = 1:EEG.nbchan
    lbl = sprintf('Channel %d (%s)', ch, EEG.chanlocs(ch).labels);
    h1(ch).DataTipTemplate.DataTipRows(1).Label = lbl;
    h2(ch).DataTipTemplate.DataTipRows(1).Label = lbl;
    set([h1(ch) h2(ch)], 'ButtonDownFcn', @(~, ~) showChannel(EEG, erp, ch));
end
end


function showChannel(EEG, erp, ch)
% Detail figure for one channel: single trials + channel ERP vs all channels
tr  = squeeze(EEG.data(ch, :, :));
lbl = EEG.chanlocs(ch).labels;

figure('Color', 'w', 'Name', sprintf('Channel %d (%s)', ch, lbl));

axA = subplot(2, 1, 1);
plot(EEG.times, tr, 'Color', [0.75 0.75 0.75]); hold on;
plot(EEG.times, mean(tr, 2), 'k', 'LineWidth', 2.5); hold off;
xlim([-100 600]); ylim([-150 150]);
xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
title(sprintf('Channel %d (%s): single trials (grey), ERP (black)', ch, lbl));
set(axA, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);

axB = subplot(2, 1, 2);
m_all = mean(erp, 1); s_all = std(erp, 0, 1);
fill([EEG.times fliplr(EEG.times)], [m_all + s_all fliplr(m_all - s_all)], ...
    [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on;
plot(EEG.times, m_all, 'k', 'LineWidth', 1.5);
plot(EEG.times, erp(ch, :), 'r', 'LineWidth', 2.5); hold off;
xlim([-100 600]); ylim([-80 80]);
xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
title('Channel ERP (red) vs. mean \pm SD across all channels');
legend('Mean \pm SD all channels', 'Mean all channels', sprintf('Channel %d', ch), 'Location', 'best');
set(axB, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);
end