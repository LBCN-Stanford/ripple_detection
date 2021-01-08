% SWR detection main script:
% detecting ripples after exclusion of electrical, muscular and IED artifacts.
%
% Author: Yitzhak Norman, 2019, Malach Lab @ Weizmann Institute
close all;
clear all;
clc; warning('off');
% set the relevant path:
path_to_toolboxes = 'D:\MATLAB_ToolBoxes\';

% Set Path:
rmpath(fullfile(path_to_toolboxes,'eeglab13_4_4b'));
addpath(fullfile(path_to_toolboxes,'eeglab14_1_2b'));
rmpath(genpath(fullfile(path_to_toolboxes,'chronux_2_12')));
addpath(genpath(fullfile('D:\ECoG\MMR\data_sharing\Ripple_detection','matlab_scripts')));

% load eeglab:
[ALLEEG, EEG, CURRENTSET] = eeglab;

subjects={'SUB01','SUB02','SUB03'};

% Load data:
for subjid=subjects(1)
    subjid=cell2mat(subjid);
    
    clearvars -except subjects subjid run ALLEEG EEG HFOepoches HFOch
    ALLEEG=[]; EEG=[]; CURRENTSET=1;
    HFOepoches = {}; HFOch = {};
    RIPPLES = [];
    % Set Path:
    warning('off','all')
    scenario = 'task';
    maindir = ['D:\ECoG\MMR\data_sharing\Ripple_detection\' subjid]; % adjust single-subject data directory
    CURRENTSET=1;
    
    % Set relevant filenames and folders:
    ref_flag=2; % 1 = Common Ref; 2 = Bipolar montage;
    switch ref_flag
        case 1
            datadir=[maindir '\EEGLAB_datasets'];
            file_ending='_preprocessed.set';
        case 2
            datadir=[maindir '\EEGLAB_datasets_BP'];
            file_ending='_preprocessed_BP_montage.set';
    end
    
    % Load EEGLAB datasets:
    cd(datadir);
    filenames=dir(['*' scenario '*.set']);
    for SET=1:numel(filenames)        
        filename=filenames(SET).name;
        [EEG] = pop_loadset('filename', filename);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
    end
    
    
    %% Defining Hippocampus:
    set_figure_colors;
    [hippocampus,hippocampus_all_channels,WM_ref] = define_hippocampal_channels(subjid,ref_flag); % hippocampal CA1/Sub channel for ripple detection
    
    ripples = [];
    for current_hippocampus = hippocampus_all_channels
        
        cnum1=find(strcmpi({EEG.chanlocs.labels},current_hippocampus)); % CA1 channel, identified anatomically
        cnum2=find(strcmpi({EEG.chanlocs.labels},'CREF')); % common ref. channel
        channel=EEG.chanlocs(cnum1).labels;           
        RIPPLES(end+1).channel = channel;
               
        % we use the averaged LFP across all good
        % channels (common ref) to perform control detection of ripples,
        % which allows to discard muscular/electrical artifacts.
        
        if isempty(current_hippocampus)
            error('Missing hippocampal channel!')
        end
        
        %% Data normalization and epoching:
        
        %%%%% Ripple detection PARAMETERS %%%%%
        minDistance=0.030; % in sec
        minRippleDuration=0.020; % in sec [0.015]
        maxRippleDuration=0.200; % in sec
        % =================================================================
        th=[2 4]; % ripple detection thresholds (in std from the mean)
        ripple_band=[70 180]; % in Hz
        % note: correct for a transition bandwith of 5Hz (assuming forder is 166)
        % (expand the frequency band by +/- 5 Hz)
        IED_frequency_band=[25 60]; % in Hz
        LPcutoff=round(mean(ripple_band)/pi); % in Hz (see Stark et al. 2014)
        %LPcutoff=40; % in Hz (fixed size LP filter)
        Fs=EEG.srate; % in Hz
        forder_BP=330; % in samples
        wtype = 'hamming';
        warg = 5.653; % beta value, in case you choose Kaiser window
        % =================================================================
        
        % set output folder:
        outdir=fullfile('D:\ECoG\MMR\data_sharing\Ripple_detection\swrs',...
            sprintf('Ripple_times_%s_%d-%dstd_%d-%dHz_%dHz_%dms_%dms_ref%d',...
            wtype,th(1),th(2),ripple_band(1),ripple_band(2),LPcutoff,...
            minRippleDuration*1000,minDistance*1000,ref_flag),subjid); % ADJUST OUTDIR
        
        if ~exist(outdir,'dir')
            mkdir(outdir);
            disp('Creating Output Directory...')
        end
              
       
        for SET = 1:numel(ALLEEG)
            
            clear ripples avg stdev
            close all
            
            % Compute mean and std of ripple band amplitude across the entire run: 
            [EEG_bandpass, ~,b] = pop_firws(ALLEEG(SET), 'fcutoff', ripple_band, 'ftype', 'bandpass', 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            rawSignal = double(ALLEEG(SET).data(cnum1,:));
            BANDPASS = EEG_bandpass.data(cnum1,:);
            absSignal = double(abs(hilbert(BANDPASS))); % hilbert envelope
            HFOonsets = zeros(size(rawSignal));         % prepare timeseries for HFOs
            % Try to load HFOs onsets (optional): 
            HFOsdir = fullfile('D:\ECoG\MMR\Results\data\HFOs',subjid);
            HFOfile = fullfile(HFOsdir,[subjid ' run ' num2str(SET) ' HFOs ' cell2mat(current_hippocampus) '.mat']);
            HFOs = [];
            HFOflag = 0;
            if exist(HFOfile,'file') && HFOflag, load(HFOfile); end
            if ~isempty(HFOs)
                for ii = 1:length(HFOs)
                    [~,mn] = min(abs((ALLEEG(SET).times/ALLEEG(SET).srate)-HFOs(ii)));
                    HFOonsets(mn)=1;
                end
                [ep,eplim] = epoch(rawSignal,HFOs,[-1.5 1.5],'srate',ALLEEG(SET).srate);
                if ~isempty(ep) 
                        ep = squeeze(ep);
                        if size(ep,1)==1, ep = ep'; end
                        HFOepoches = cat(1,HFOepoches,{ep}); 
                        HFOch = cat(1,HFOch,[subjid '_' cell2mat(current_hippocampus)]);
                end
            end
            % Defining thr for clipping:
            [robustAvg,robustStdev] = robustMean(absSignal,2,th(2));     % option 1: robust mean and std used for clipping            
            topLim=robustAvg+th(2)*robustStdev;
            %topLim = nanmean(absSignal)+th(2)*nanstd(absSignal);        % option 2: mean & std
            %topLim =  nanmedian(absSignal)+(th(2)/1.35)*iqr(absSignal); % option 3:median & IQR 
            
            absSignal(absSignal>topLim) = topLim; % clipping the signal
            squaredSignal = absSignal.^2;         % squared amplitude
            % Lowpass fliter:            
            lpFilt = designfilt('lowpassfir','PassbandFrequency',LPcutoff, ...
                'StopbandFrequency',LPcutoff+10,'PassbandRipple',0.001, ...
                'StopbandAttenuation',60,'DesignMethod','kaiserwin','SampleRate',Fs);% (FIR kaiserwin lowpass filter)
            squaredSignal = filtfilt(lpFilt,double(squaredSignal)); 
            % Compute mean and std of the clipped & squared time series:
            avg = nanmean(squaredSignal);
            stdev = nanstd(squaredSignal,[],2);
            clear squaredSignal BANDPASS rawSignal robustStdev robustAvg
            close all
                        
            EEG=ALLEEG(SET);
            EEG_BANDPASS_ripple = pop_firws(ALLEEG(SET), 'fcutoff', ripple_band, 'ftype', 'bandpass', 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            EEG_BANDPASS_IED = pop_firws(ALLEEG(SET), 'fcutoff', IED_frequency_band, 'ftype', 'bandpass', 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
            
            
            BANDPASS_A = EEG_BANDPASS_ripple.data(cnum1,:);
            BANDPASS_B = EEG_BANDPASS_ripple.data(cnum2,:);
            BANDPASS_C = EEG_BANDPASS_IED.data(cnum1,:);
            
            strind = round(EEG.event(1).latency); % finds the first event (starting point)
            tail_correction = EEG.times(strind)./1000;                        
            %tail_correction = 5; % in Sec
            T=double(EEG.times./1000) - tail_correction; % in Sec
            
            % Hilbert envelope (rectification):
            absSignalA=double(abs(hilbert(BANDPASS_A)));
            absSignalB=double(abs(hilbert(BANDPASS_B)));
            absSignalC=double(abs(hilbert(BANDPASS_C)));
            
            % Squaring the signal:
            squaredSignalA = absSignalA.^2;
            squaredSignalB = absSignalB.^2;
            squaredSignalC = absSignalC.^2;
            
            
            % FIR filter (after squaring):
            squaredSignalA = filtfilt(lpFilt,double(squaredSignalA));
            squaredSignalB = filtfilt(lpFilt,double(squaredSignalB));
            squaredSignalC = filtfilt(lpFilt,double(squaredSignalC));
                        
            % ZSCORE:
            squaredSignalNormA = (squaredSignalA-avg)/stdev; % z-score hippocampal ripple band amplitude
            squaredSignalNormB = (squaredSignalB-nanmean(squaredSignalB))/nanstd(squaredSignalB); % z-score "emg" signal (extracted from common average)
            squaredSignalNormC = (squaredSignalC-nanmean(squaredSignalC))/nanstd(squaredSignalC); % z-score hippocmapal channel - IED
            
            % Find ripples:
            [ripples,ripples_stat,rnoise] = ripples_detection_excluding_IED(squaredSignalNormA,BANDPASS_A,T,Fs,th,minDistance,...
                minRippleDuration,maxRippleDuration,squaredSignalNormB,squaredSignalNormC,HFOonsets); 
            outfilename = [subjid ' run ' num2str(SET) ' ripples ' cell2mat(current_hippocampus)];
            eval([sprintf('t_run%d',SET) '=T;'])
            save(fullfile(outdir,outfilename),'ripples','ripples_stat',sprintf('t_run%d',SET))
            
            % Store the ripples in a data structure:
            RIPPLES(end).swramp = mean(ripples.amplitude);
            RIPPLES(end).blavg = avg;
            RIPPLES(end).blstdev = stdev;
            RIPPLES(end).count = size(ripples,1);
            RIPPLES(end).rnoise = rnoise;
            
            
            %% Plot ripple detection figure:
            close all            
            scenarioLabel=capitalize(sprintf('%s run %d',scenario,SET)); scenarioLabel(scenarioLabel=='_')=' ';
            H1=figure('Name',[subjid ' ' scenarioLabel ' ripples detection Ch ' channel],'position',[0 100 1000 300],'Color','w');
            hold on;
            title(['Ripples Detection - ' scenarioLabel ' (ripple band envelope)']);
       
            h1 = plot(T,zscore(BANDPASS_A)+13,'linewidth',0.5,'color',[1,0.7,0.7]); hold on;
            h2 = plot(T,squaredSignalNormA,'k','Linesmoothing','on');
            h0 = plot(T,ones(size(T))*th(2),'b--','Linewidth',0.5);
            if ~isempty(ripples)
                h3=scatter(ripples.peak,ones(size(ripples,1),1)*-2,30,'yo','fill'); hold on                
            end
            ylabel(sprintf('%d-%dHz Amplitude^2 (Z-score)',ripple_band(1),ripple_band(2)))
            xlim([T(1),T(end)])
            ylim([-2,18])
            set(gca,'ytick',[0:2:10])          
            legend([h2 h1 h3 h0],{'Ripple Band Amplitude^2 (zscore)','Ripple Band (A.U)','Ripple','Detection thr.'}); legend boxoff;
            xlim([190 210])
            set(gca,'xtick',190:5:210,'xticklabel',{0:5:20})
                     
            
            %% save figures;
            saveflag = 1;
            if saveflag
                figdir=fullfile(outdir,'detection_figures');
                if ~exist(figdir,'dir')
                    mkdir(figdir);
                    disp('Creating Output Directory...')
                end
                set(0,'DefaultAxesFontName', 'Arial')
                set(gcf,'renderer','painters')
                saveas(H1,fullfile(figdir,get(gcf,'name')),'fig');
                export_fig(fullfile(figdir,get(gcf,'name')),'-nocrop','-jpg','-r150','-transparent')
                export_fig(fullfile(figdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap','-transparent')

            end
            
        end
    end
    
    %% Find the channel with the strongest SNR: (this part is still in development...) 
    % SNR = ripple power divided by the background ripple-band activity 
    channels= {RIPPLES.channel};
    swramp  = [RIPPLES.swramp];
    blavg   = [RIPPLES.blavg];
    blstdev   = [RIPPLES.blstdev];
    swrsnr = 10*log10(swramp./blavg); % SNR calculation
    %swrsnr = (swramp-blavg)./blstdev; % Cohen's d alternative
    %swrsnr = swramp./blstdev; %  alternative SNR calculation
    figure('color','w','name',['compare ripples between electrodes ref ' num2str(ref_flag)],'position',[0 0 600 200]);
    subplot(1,3,1); hold on; superbar(1:numel(RIPPLES),swrsnr); axis square
    set(gca,'xtick',[1:numel(RIPPLES)],'xticklabel',{RIPPLES.channel})
    rotateXLabels(gca,60)
    ylabel('Ripple SNR (dB)')
    subplot(1,3,2); hold on; superbar(1:numel(RIPPLES),[RIPPLES.count]); axis square
    set(gca,'xtick',[1:numel(RIPPLES)],'xticklabel',{RIPPLES.channel})
    rotateXLabels(gca,60)
    ylabel('Ripple Count')
    subplot(1,3,3); hold on; superbar(1:numel(RIPPLES),[RIPPLES.rnoise].^2); axis square
    set(gca,'xtick',[1:numel(RIPPLES)],'xticklabel',{RIPPLES.channel})
    rotateXLabels(gca,60)
    ylabel('noise corr (R^2)')
    
    [ampmax, indmax] = max(swrsnr);
    fprintf(' \n BEST RIPPLES SNR WAS FOUND IN CHANNEL: %s \n ', channels{indmax});
    suptitle(sprintf('BEST RIPPLES SNR WAS FOUND IN CHANNEL: %s', channels{indmax}));
    set_font_size_and_type;
    saveas(gcf,fullfile(outdir,get(gcf,'name')),'fig');
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-jpg','-r150','-transparent')
    export_fig(fullfile(outdir,get(gcf,'name')),'-nocrop','-pdf','-painters','-nofontswap','-transparent')
    
    % Rejected HFO-related epochs:   
    if HFOflag
        if ~exist(HFOsdir,'dir'),mkdir(HFOsdir); disp('Creating HFOs Output Directory...'); end
        save(fullfile(HFOsdir,'rejected_HFO_events_epochs.mat'),'HFOch','HFOepoches'); disp('rejected HFO-related data - saved');
    end
end

%% OPTION 1: Compute stdev across the entire experiment - 
%         EEGmereged = pop_mergeset( ALLEEG, 1:numel(ALLEEG));
%         [EEG_bandpass, ~,b] = pop_firws(EEGmereged, 'fcutoff', ripple_band, 'ftype', 'bandpass', 'warg', warg, 'wtype', wtype, 'forder', forder_BP,  'minphase', 0);
%         figure; freqz(b);
%         eeglab redraw
%         rawSignal = double(EEGmereged.data(cnum1,:));
%         BANDPASS = EEG_bandpass.data(cnum1,:);
%         absSignal = double(abs(hilbert(BANDPASS))); % hilbert envelope
%         [robustAvg,robustStdev] = robustMean(absSignal,2,3);
%         % Clipping the signal:
%         topLim=robustAvg+th(2)*robustStdev;
%         absSignal(absSignal>topLim)=topLim;
%         squaredSignal = absSignal.^2;
%         % Lowpass and square:
%         % lpFilt = designfilt('lowpassfir','PassbandFrequency',LPcutoff, ...
%             'StopbandFrequency',LPcutoff+10,'PassbandRipple',0.001, ...
%             'StopbandAttenuation',60,'DesignMethod','kaiserwin','SampleRate',Fs);
%             % FIR kaiserwin lowpass filter%             
%         % fvtool(lpFilt,'OverlayedAnalysis','phase')
%         squaredSignal = filtfilt(lpFilt,double(squaredSignal));       
%         % Compute mean and std:
%         avg = nanmean(squaredSignal);
%         stdev = nanstd(squaredSignal,[],2);