
function [ripples,ripples_stat,r]=ripples_detection_excluding_IED(signal,BP,t,Fs,th,minDistance,minRippleDuration,maxRippleDuration,noise_ch,IED_ch,HFOonsets)
% detecting ripples based on Strack et al. 2014 (Neuron) method
% signal = normalized squared bandpassed LFP signal from the hippocampus
%          (assuming sampling rate of 500Hz)
% BP = raw bandpass signal
% t = corresponding time in sec
% th = threshold in stdev [onset/offset peak]
% minDistance = in sec
% minRippleDuration = min ripple duration in Sec
% maxRippleDuration = max ripple duration in Sec
% noise_ch = channel to exclude EMG artefacts identifed as rippples
% IED_ch = channel to exclude IEDs that were identifed as rippples
% HFOonsets = timestamps of HFOs, detected using Su et al. (2016) method
%
% Author: Itzik Norman 10/03/20

if size(signal,2)<size(signal,1),signal=signal'; end
if size(noise_ch,2)<size(noise_ch,1),noise_ch=noise_ch'; end
if size(IED_ch,2)<size(IED_ch,1),IED_ch=IED_ch'; end

[pks,locs] = findpeaks(signal, 'MINPEAKHEIGHT', th(2));
ENV = abs(hilbert(BP)).^2;   % to calculate ripple amplitude (squared)
% ENV = abs(hilbert(BP)); ENV=ENV./nanmedian(ENV); ENV = 10*log10(ENV); % to calculate ripple amplitude in dB relative to median
% Noise rejection:
if ~isempty(noise_ch)
    [r,p]=corr(signal',noise_ch');
    fprintf('\n --- correlation with noise ch.: %.2f \n',r)
    rej_count=0;
    pre_rej_count=size(locs,2);    
    fprintf('\n => rejecting noise artefacts... \n')
    [~,noise_locs] = findpeaks(noise_ch, 'MINPEAKHEIGHT', th(2));
    % Ignore global electrical artifacts:
    for i=1:numel(noise_locs)
        tmp=find(abs(locs-noise_locs(i))<0.05*Fs);
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events based on noise channel correlation \n',rej_count,pre_rej_count)
    ripples_stat.noise_rejection = rej_count/pre_rej_count;
end

% IED rejection (Gelinas et al. 2016 Nat. Med.):
if ~isempty(IED_ch)
    rej_count=0;
    pre_rej_count=size(locs,2);    
    fprintf('\n => rejecting IED events... \n')
    [~,IED_locs] = findpeaks(IED_ch, 'MINPEAKHEIGHT', 4); % select 3 or 4 SD
    % Ignore IED-ripples events (event that coincide within 200 ms):
    for i=1:numel(IED_locs)
        tmp=find(abs(locs-IED_locs(i))<0.1*Fs);
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events based on IED channel correlation \n',rej_count,pre_rej_count)
    ripples_stat.IED_rejection = rej_count/pre_rej_count;
end

% stereotyped HFOs rejection: Su et al. (2016,2018)
if any(HFOonsets)
    pre_rej_count=size(locs,2); 
    % Ignore ripples that overlap with stereotyped HFO events (that coincide within 200 ms):
    HFO_locs = find(HFOonsets);
    for i=1:numel(HFO_locs)
        tmp=find(abs(locs-HFO_locs(i))<0.1*Fs);
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events that coincide with stereotyped HFOs \n',rej_count,pre_rej_count)
    ripples_stat.HFO_rejection = rej_count/pre_rej_count;
end

counter=1;
ripples=nan(1,4);
ripples=array2table(ripples,'VariableNames',{'str','peak','fin','amplitude'});
ripples(1,:)=[];
for k=locs
    % find the starting point of the peak:
    stop=0;
    str=k;
    while ~stop && ~str==0
        if str==1
            break
        end
        str=str-1;
        if signal(str)<th(1), stop=1; end
    end
    % find the ending point of the peak:
    stop=0;
    fin=k;
    while ~stop && ~(fin==numel(signal))
        fin=fin+1;
        if signal(fin)<th(1), stop=1; end
    end
    % =====================================================================
    % Alternative #1:
    % Detect negative peak position for each ripple (closest to ripple's power peak)
    minIndex = [];
    [~,minpos] = findpeaks(-double(BP(str:fin)));
    if isempty(minpos), [~,minpos] = min(BP(str:fin)); minpos = minpos(1); end
    [~,maxamp] = max(double(ENV(str:fin)));
    minpos=minpos-1; 
    maxamp=maxamp-1;
    %p=k-str;
    %[~,tmp] = min(abs(minpos-p));
    [~,tmp] = min(abs(minpos-maxamp));
    minIndex=minpos(tmp);
    peakPosition = min((str + minIndex),numel(signal));
    % =====================================================================
    % Alternative #2: (absolut min)
%     minIndex = [];
%     [~,minpos] = min(BP(str:fin));
%     minpos=minpos-1;
%     peakPosition = min((str + minpos),numel(signal));
    
    
    try
        M = [t(str), t(peakPosition), t(fin), ENV(peakPosition)];
        M = round(M.*1000)./1000;
        ripples(counter,:)=array2table(M);
    catch
        disp(ripples);
        fprintf('\n Error has occured in event # %d \n',counter);
    end
    counter=counter+1;
end
disp(['After detection by thresholding: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end


% Merge ripples if inter-ripple period is less than minDistance:
ripples_edit=ripples;
rej=zeros(size(ripples,1),1);
for k = 2:size(ripples,1)
    if (ripples.peak(k)-ripples.peak(k-1)) < minDistance,
        % Merge
        ripples_edit.fin(k-1) = ripples.fin(k);
        rej(k)=1;
    end
end
if any(rej), ripples_edit(find(rej),:)=[]; end
ripples=ripples_edit;
disp(['After ripple merge: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end

% duration test:
duration = ripples.fin-ripples.str;
ripples(duration<minRippleDuration,:) = [];
disp(['After min duration test: ' num2str(size(ripples,1)) ' events.']);
duration = ripples.fin-ripples.str;
ripples(duration>maxRippleDuration,:) = [];
disp(['After max duration test: ' num2str(size(ripples,1)) ' events.']);


end
