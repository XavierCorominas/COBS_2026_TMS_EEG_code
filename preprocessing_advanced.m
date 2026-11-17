%% Preprocessing pipeline example for TEPs.
% DRCMR 2026 - YFS XCT - BrainStim Methods group

% The preprocessing is inspired in the following references:
% Bracco, M., Mutanen, T. P., Veniero, D., Thut, G., & Robertson, E. M. (2025). Protocol to assess changes in brain network resistance to perturbation during offline processing using TMS-EEG. STAR protocols, 6(1), 103622. https://doi.org/10.1016/j.xpro.2025.103622
% Bertazzoli, G., Esposito, R., Mutanen, T. P., Ferrari, C., Ilmoniemi, R. J., Miniussi, C., & Bortoletto, M. (2021). The impact of artifact removal approaches on TMS-EEG signal. NeuroImage, 239, 118272. https://doi.org/10.1016/j.neuroimage.2021.118272
% TESA tmseeg preprocessing steps: 
% https://nigelrogasch.gitbook.io/tesa-user-manual/an_overview_of_tms-eeg_analysis

% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide an elaborated and meticulous eeg
% preprocessing focusing on removing backgrond noise, TMS induced artifacts, muscular artifacts and eye movements. The preprocessing might be
% sufficient to explore iTEPS/TEPS in datasets contaminated with muscular (or other type of) artifacts.
% NOTE: Not every step here is needed for every dataset. Choose the steps,
% their order and their parameters to suit your own data.

% ---> STEPS:
% 0.   LOAD DATA
% 1.   EPOCH
% 2.   Baseline Correct
% 3.   Remove bad trials and channels 
% 3.1  Identify bad trials and remove
% 3.2  Plot and identify bad channels and remove
% 4.   Remove TMS pulse artifact
% 5.   Interpolate missing data
% 5.1  (Optional) Remove and interpolate recharging artefact if necessary
% 6.   First round of ICA --> remove blinks or eye movements
% 6.1    Check maximum number of independent components (SVD)
% 6.2    MAIN ICA - First round
% 6.3    Plot the independent components and reject components associated to eye movement/blink
% 7.   Baseline correction
% 8.  (Optional) Deal with noise if necessary with SOUND
% 9.  (Optional) SSP-SIR for muscular evoked activity
% 10.  Second round of ICA --> identify residual muscular and persisting artifacts
% 10.1   Check maximum number of independent components (SVD)
% 10.2   MAIN ICA - second round
% 10.3   Plot the independent components and reject components associated to persistent muscular activity
% 11.  Lowpass filter
% 12.  Baseline correction
% 13.  Rereferencing (if necessary)
% 14.  Save dataset


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

%% 0. LOAD DATA
% Open EEGLAB
eeglab;

% Define file to load
name_dataset = 'X99880_sp_110rmt_block1.vhdr';
[~, save_name] = fileparts(name_dataset);

% Load file
eeg_setname = 'continuous eeg';
EEG = pop_loadbv(data_path, name_dataset);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
        'overwrite', 'on', 'setname', eeg_setname);  
eeglab redraw

fprintf('\n>> Step 0 (Load data: %s) ended\n', name_dataset);
fprintf('%s\n\n', repmat('-', 1, 70));

%% 1. Epoch
eeg_setname = 'epoch';
% Epoch around TMS pulses 'R  8' or 'S 11' markers
EEG = pop_epoch(EEG, {'R  8', 'S 11'}, [-1 1]); % sec, x-axis range
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname);
eeglab redraw;

fprintf('\n>> Step 1 (Epoch around TMS pulses, -1 to 1 s: %d epochs) ended\n', EEG.trials);
fprintf('%s\n\n', repmat('-', 1, 70));

%% 2. Baseline Correct
eeg_setname = 'baseline correction';

EEG = pop_rmbase( EEG, [-210 -10] ,[]); 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

fprintf('\n>> Step 2 (Baseline correction) ended\n'); %[not Demeam anymore]
fprintf('%s\n\n', repmat('-', 1, 70));

%% 3. Remove bad trials, remove and interpolate bad channels
%% 3.1.1 (identify bad trials) --> Flag the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

fprintf('\n>> Step 3.1.1 (Identify bad trials: data displayed, flag trials and press REJECT) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));
 
%% 3.1.2 (remove bad trials)
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
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}


fprintf('\n>> Step 3.1.3 (Plot data to inspect and identify bad channels) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% 3.4 Manual rejection and interpolation of channels.
eeg_setname = 'channel rejection';
% Select bad channel numbers
badChannels = [22]; % example: channel 2

% Save original channel locations before removal
originalChanLocs = EEG.chanlocs;

% Remove bad channels
EEG = pop_select(EEG, 'nochannel', badChannels);

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

fprintf('\n>> Step 3.4 (Manual rejection/interpolation of channels: removed [%s]) ended\n', num2str(badChannels));
fprintf('%s\n\n', repmat('-', 1, 70));

%% 4. Remove TMS pulse artifact(s)
eeg_setname = 'remove TMS pulse artifact';
% TMS artefact lengths in ms
TMSremoval_range =[-2 6]; % in ms
EEG = pop_tesa_removedata( EEG, TMSremoval_range );

% Remove TMS pulse
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot channels: hover over a line to see the channel number
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 4 (Remove TMS artifact: %g to %g ms) ended\n', TMSremoval_range(1), TMSremoval_range(2));
fprintf('%s\n\n', repmat('-', 1, 70));


%% 5. Interpolate removed window
eeg_setname = 'interpolate missing data';
% Interpolate
EEG = pop_tesa_interpdata( EEG, 'cubic', [1 1] ); 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot channels: hover over a line to see the channel number
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 5 (Interpolate removed window) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%{
%% 5.1 (Optional) Remove and interpolate recharging artefact if necessary

EEG = EEGLAB_remove_and_interpolate_recharging_artifact(EEG, 1); % 1 + single pulse, 2 = repetitive pulse

fprintf('\n>> Step 5.1 (Remove and interpolate recharging artifact) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));
%}

%% 6. First round of ICA --> remove blinks or eye movements
%{
% Resample first for ica:
eeg_setname = 'resample';
EEG = pop_resample( EEG, 1000); %or 5000
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw
%}

% === load resampled data ===
EEG = pop_loadset('filepath', data_path,...
    'filename','X99880_sp_110rmt_block1_resample.set');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', 'resample'); 
eeglab redraw


% Prepare the data for ICA: First remove for a while the very worst channels!
EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


% Labeling the very worst channels to not affect the ICA run
badC = find(sigmas > (median(sigmas) + 3*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);
EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC

fprintf('\n>> Step 6 (Prepare first ICA: resampled to %d Hz, %d worst channel(s) excluded [%s]) ended\n', EEG.srate, numel(badC), num2str(badC(:)'));
fprintf('%s\n\n', repmat('-', 1, 70));


%% 6.1 Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition

tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca

fprintf('\n>> Step 6.1 (SVD check: max number of independent components = %d) ended\n', Cpca);
fprintf('%s\n\n', repmat('-', 1, 70));

%% 6.2 MAIN ICA - First round

% Set the number of dimentions accordingly
% You might need to go one below the previous value
% Here we only remove eye-blinks and horizontal movements: to make your
% choise, if unsure check: 
% https://link.springer.com/article/10.1007/s11517-011-0748-9
% https://www.sciencedirect.com/science/article/abs/pii/S105381191400620X?via%3Dihub
% https://www.sciencedirect.com/science/article/pii/S1053811916305845?via%3Dihub

EEG2ICA = pop_tesa_pcacompress( EEG2ICA, 'compVal', 30, 'plot','on' ); % Set the number of dimentions here
EEG2ICA = pop_tesa_fastica( EEG2ICA, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

fprintf('\n>> Step 6.2 (MAIN ICA - first round: PCA compression + FastICA) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));
 

%% 6.3 Plot the independent components following ICA. Here the user can verify the results of ICA
% as well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
eeg_setname = '1st round ICA';
% Alternative
%{
figure; 
for i=1:30
    subplot(5,6,i)
    topoplot(EEG2ICA.icawinv(:,i), EEG2ICA.chanlocs);
end
%}

EEG2ICA = pop_tesa_compplot( EEG2ICA,'figSize','medium','plotTimeX',[-500 499],'plotFreqX',[1 100], 'freqScale','log', 'saveWeights','off' );


% Update corrected data and plot
EEG = EEG2ICA;
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw


plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 6.3 (Plot/classify independent components, first round) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% 7. Baseline correction
% === load data after 1st round ica ===
EEG = pop_loadset('filepath', data_path,...
    'filename','X99880_sp_110rmt_block1_1st round ICA.set');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', '1st round ICA'); 
eeglab redraw

eeg_setname = 're baseline correction';
EEG = pop_rmbase( EEG, [-210 -10] ,[]);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 10 (Baseline correction, -210 to -10 ms) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% 8. Deal with noise with SOUND
% Mutanen, T. P., Metsomaa, J., Liljander, S., & Ilmoniemi, R. J. (2018). Automatic and robust noise suppression in EEG and MEG: The SOUND algorithm. NeuroImage, 166, 135–151. https://doi.org/10.1016/j.neuroimage.2017.10.021
% TESA tmseeg preprocessing steps: 
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/sound

% Note1: after SOUND, the dataset is in the average reference by default'
% Note2: if a EEGLAB dataset with all the original channels is provided,
% SOUND automatically interpolate removed channels

eeg_setname = 'SOUND';
% make a copy of the data before
data_before_sound = EEG;
% re-ref to align with data after SOUND
data_before_sound = pop_reref( data_before_sound, []);

% Modify channel type name for runing sound
for i = 1:length(EEG.chanlocs)
EEG.chanlocs(i).type = 'EEG';
end 

allchannels = [data_path,filesep, 'eeglab_dataset_with_allchannels.set'];

lambda = 0.05; % try 0.01, 0.05, 0.1 to see the smoothness effect
% USE SOUND
EEG = pop_tesa_sound(EEG, 'lambdaValue', lambda, 'iter', 10, 'replaceChans', allchannels);


[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

%
%Plot data
figure; pop_timtopo(EEG, [-10  200], [10 30  45], 'TEP and scalp maps');

% Compare before and after SOUND
figure('color','w', 'Position', [100 100 900 492]);

hBefore = plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'r');
hold on
hAfter  = plot(EEG.times, mean(EEG.data, 3)', 'b');
hold off

xlim([-10 200]);   % ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
legend([hBefore(1) hAfter(1)], {'Before SOUND', 'After SOUND'}, ...
    'FontSize', 20, 'Location', 'northeast');

fprintf('\n>> Step 11 (SOUND noise suppression, lambda 0.1, 10 iterations) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%{
%% 9. (Optional) SSP-SIR for muscular evoked activity
% Mutanen, T. P., Kukkonen, M., Nieminen, J. O., Stenroos, M., Sarvas, J., & Ilmoniemi, R. J. (2016). Recovering TMS-evoked EEG responses masked by muscle artifacts. NeuroImage, 139, 157–166. https://doi.org/10.1016/j.neuroimage.2016.05.028
% TESA tmseeg preprocessing steps: 
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/ssp-sir

eeg_setname = 'SSP-SIR';
% Make a copy before
data_before_sspsir = EEG;

% Plot
[EEG] = pop_tesa_sspsir(EEG, 'artScale', 'manual','timeRange',[5,35],'PC',[]);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw
%
%Plot data
figure; pop_timtopo(EEG, [-10  200], [10  30  45], 'TEP data and scalp maps');

% Compare before and after SSP-SIR
figure('Position', [100 100 900 492]);  

hBefore = plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k');
hold on
hAfter  = plot(EEG.times, mean(EEG.data, 3)', 'b');
hold off

xlim([-10 200]);   % ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
legend([hBefore(1) hAfter(1)], {'Before SSP-SIR', 'After SSP-SIR'}, ...
    'FontSize', 20, 'Location', 'northeast');

fprintf('\n>> Step 12 (SSP-SIR for muscular evoked activity, 5-35 ms) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));
%}

%% 10. Second round of ICA --> indentify residual muscular and persisting artifacts
% Prepare the data for ICA: First remove for a while the very worst channels!
EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')

% Labeling the very worst channels to not affect the ICA run
badC = find(sigmas > (median(sigmas) + 3*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);
EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC

fprintf('\n>> Step 13 (Prepare second ICA: %d worst channel(s) flagged [%s]) ended\n', numel(badC), num2str(badC(:)'));
fprintf('%s\n\n', repmat('-', 1, 70));


%% 10.1 Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition
tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca

fprintf('\n>> Step 13.1 (SVD check: max number of independent components = %d) ended\n', Cpca);
fprintf('%s\n\n', repmat('-', 1, 70));


%% 10.2 MAIN ICA - second round

% Set the number of dimentions accordingly
% Here we remove persistent muscular artifacts and noise
% choise, if unsure check: 
% https://link.springer.com/article/10.1007/s11517-011-0748-9
% https://www.sciencedirect.com/science/article/abs/pii/S105381191400620X?via%3Dihub
% https://www.sciencedirect.com/science/article/pii/S1053811916305845?via%3Dihub

EEG2ICA = pop_tesa_pcacompress( EEG2ICA, 'compVal', 25, 'plot','on' ); % Set the number of dimentions here
EEG2ICA = pop_tesa_fastica( EEG2ICA, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

fprintf('\n>> Step 13.2 (MAIN ICA - second round: PCA compression + FastICA) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));
 

%% 10.3 Plot the independent components following ICA. Here the user can verify the results of ICA
% as well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
eeg_setname = '2nd round ICA';

% Alternative
%{
figure; 
for i=1:20
    subplot(5,4,i)
    topoplot(EEG2ICA.icawinv(:,i), EEG2ICA.chanlocs);
end
%}

EEG2ICA = pop_tesa_compplot( EEG2ICA,'figSize','medium','plotTimeX',...
    [-500 499],'plotFreqX',[1 100], 'freqScale','log', 'saveWeights','off' );

EEG = EEG2ICA;
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 13.3 (Plot/classify independent components, second round) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));


%% 11. Lowpass Filter
% === load data after 2nd round ica ===
EEG = pop_loadset('filepath', data_path,...
    'filename','X99880_sp_110rmt_block1_2nd round ICA.set');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', '2nd round ICA'); 
eeglab redraw

eeg_setname = 'filter data';
EEG = tesa_filtbutter(EEG, [], 80, 2, 'lowpass');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw

% Plot
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-10 10]);
%}
fprintf('\n>> Step 14 (Low-pass filter 2000 Hz, Butterworth order 2) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));


%%  12. Baseline Correction
eeg_setname = 're baseline correction';
EEG = pop_rmbase( EEG, [-210 -10] ,[]);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','off', 'setname', eeg_setname); 
eeglab redraw


%Plot the data 
plotTEP(EEG, [10 30 45 60 90 150 240], [-20 20]);

% Alternative plot
%{
figure; 
pop_timtopo(EEG, [-50  200], [10 30 45], 'TEP and scalp maps','maplimits', [-30 30]);
%}

fprintf('\n>> Step 15 (Baseline correction, -210 to -10 ms) ended\n');
fprintf('%s\n\n', repmat('-', 1, 70));

%% 13. (Optional) REREFERENCING If necessary

% EEG = pop_reref( EEG, []);

% fprintf('\n>> Rereferencing step ended (currently commented out: no change applied)\n');
% fprintf('%s\n\n', repmat('-', 1, 70));


%% 14. SAVE DATASETS
eeg_setname = [save_name,'_clean_no_pulse_advanced.set'];
EEG = pop_saveset( EEG, 'filename',eeg_setname,'filepath',data_path);

fprintf('\n>> Step 16 (Save dataset: %s_clean_no_pulse_advanced.set) ended\n', save_name);
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
ylim([-1 10]);
xlabel('Time (ms)'); ylabel('GMFP (\muV)');
title(sprintf('GMFP of the average TEP (%d trials)', size(EEG.data, 3)));
set(ax0, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);

% Row 2: full time window
ax1 = subplot(4, nTopo, nTopo+1:2*nTopo);
h1 = plot(EEG.times, erp');
xlim([-600 600]); ylim([-30 30]);
xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
title('Hover to see channel, click to open channel details');
set(ax1, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5);
hold(ax1, 'on');
xline(ax1, topoLat, '--k', 'LineWidth', 1, 'HitTest', 'off');

% Middle: zoom on -5 to 20 ms
ax2 = subplot(4, nTopo, 2*nTopo+1:3*nTopo);
h2 = plot(EEG.times, erp');
xlim([-5 20]); ylim([-40 40]);
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