% Ripples_analysis_load_raw_datasets:

cd(datadir);
for SET=1:numel(scenario)
    filename=dir(['*' scenario{SET} '_' file_ending]);
    filename=filename(1).name
    [EEG] = pop_loadset('filename', filename);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;    
end


