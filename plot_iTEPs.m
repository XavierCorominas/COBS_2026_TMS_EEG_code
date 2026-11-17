%% Plot early iTEPs
% DRCMR 2026 - YFS XCT - BrainStim Methods group
%
% Run this after the preprocessing pipeline. It plots TEPs from the
% single-pulse (SP), SICI and LICI paradigms.

%% Load the three cleaned datasets
clearvars -except cobsroot
if ~exist('cobsroot', 'var')
    cobsroot = '/mnt/projects/iTEPs/COBS_2026_data_example';
end

data_path = char(fullfile(cobsroot, 'data', 'examples'));
addpath(fullfile(cobsroot, 'dependencies', 'eeglab2025.0.0', 'FastICA_25'));

close all;
clc

% Open EEGLAB
eeglab;

filenames = { ...
    'X99880_sp_110rmt_block1_clean.set', ...
    'X99880_sici_90amt_110rmt_block1_clean.set', ...
    'X99880_lici_110rmt_110rmt_100msisi_block1_clean.set'};
condNames = {'SP', 'SICI', 'LICI'};   % same order as filenames

for k = 1:numel(filenames)
    EEG = pop_loadset('filename', filenames{k}, 'filepath', data_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, k);
end
eeglab redraw

%% 0. Plot settings (shared by all sections)
ch2use  = '17';          % channel label to highlight
chIdx = find(strcmp({ALLEEG(1).chanlocs.labels}, ch2use));

sp_xlimt   = [-10 20];     % ms
xlimt   = [-110 40];
ylimt   = [-60 60];     

% Line colours 
cSP      = [0   114 178] / 255;   % blue
cSICI    = [213 94  0  ] / 255;   % red (vermillion)
cLICI    = [0   158 115] / 255;   % green
cLICIsub = [0   80  60 ] / 255;   % dark green (LICI - CS)
condColors = {cSP, cSICI, cLICI}; 

%% 1. SP dataset: butterfly plot + topoplot
EEG   = ALLEEG(1);                                  % single pulse
tep   = mean(EEG.data, 3);                          % average over trials
chIdx = find(strcmp({EEG.chanlocs.labels}, ch2use));

t2use = 3.5;                                          % ms, topoplot time
tIdx  = dsearchn(EEG.times', t2use);
clim  = [-15 15];                                   

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.3 0.6 0.35]);

% Butterfly
subplot(1, 3, [1 2]); hold on
plot(EEG.times, tep, 'Color', [0.7 0.7 0.7],  'HandleVisibility', 'off');
plot(EEG.times, tep(chIdx, :), 'Color',cSP, 'LineWidth', 2);
xline(0, '--'); yline(0, '--');
set(gca, 'XLim', sp_xlimt, 'YLim', ylimt);
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title('SP - TEP butterfly');
legend(['Ch ' ch2use],'Location','southwest');

% Topoplot
subplot(1, 3, 3);
topoplot(tep(:, tIdx), EEG.chanlocs, 'maplimits', clim, 'electrodes', 'on');
title(sprintf('%.1f ms', EEG.times(tIdx)));
colormap('parula');
cb = colorbar;
cb.Label.String = 'Amplitude (µV)';

%% 2. One channel: SP vs SICI vs LICI overlaid
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2 0.1 0.4 0.3]);
hold on

for k = 1:numel(ALLEEG)
    plot(ALLEEG(k).times, mean(ALLEEG(k).data(chIdx, :, :), 3),...
         'Color', condColors{k},'LineWidth', 2);
end
xline(0, '--', 'HandleVisibility', 'off'); yline(0, '--', 'HandleVisibility', 'off');
set(gca, 'XLim', xlimt, 'YLim', ylimt);
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title(['Ch ' ch2use ' - TEP by condition']);
legend(condNames,'Location','southwest');

%% 3. LICI: subtract the conditioning-pulse (CS) response (for visualisation)
% In LICI the test pulse is at 0 ms and the CS is at -100 ms. The CS was
% given at the same intensity as the SP (110% RMT), so the SP TEP shifted
% back by 100 ms estimates what the CS alone adds to the LICI trace.
% Subtracting it leaves the response to the test pulse.

isi = 100;                                          % ms, interstimulus interval

sp    = ALLEEG(1);
lici  = ALLEEG(3);

chIdx = find(strcmp({sp.chanlocs.labels}, ch2use));
tepSP   = mean(sp.data(chIdx, :, :), 3);
tepLICI = mean(lici.data(chIdx, :, :), 3);

% Shift the SP TEP 100 ms earlier: LICI time t gets the SP value at t + 100 ms
shift = round(isi * sp.srate / 1000);               % ms -> samples
tepCS = nan(size(tepSP));                         
tepCS(1:end-shift) = tepSP(1+shift:end);

tepLICIsub = tepLICI - tepCS;

% Plot
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2 0.1 0.4 0.3]);
hold on

plot(sp.times, tepSP,  'Color', cSP,     'LineWidth', 2);
plot(sp.times, tepLICI,  'Color', cLICI,   'LineWidth', 2);
plot(sp.times, tepLICIsub, 'Color', cLICIsub, 'LineWidth', 2);

xline(0, '--', 'HandleVisibility', 'off'); yline(0, '--', 'HandleVisibility', 'off');
set(gca, 'XLim', xlimt, 'YLim', ylimt);
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title(['Ch ' ch2use ' - LICI with CS response subtracted']);
legend({'SP', 'LICI', 'LICI - CS'},'Location','southwest');