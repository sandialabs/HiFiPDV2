%% **************** High Fidelity PDV Analysis Version 2.0 ****************
%  ************************************************************************
%  ************************************************************************
%
%  This code investigates a broad array of reasonable spectrogram inputs 
%  and peak extraction methods, then analyzes the distribution of output
%  histories to determine the best fit & uncertainty.
%
%  Author:  Travis J. Voorhees, [email: tjvoorh@sandia.gov]
%  Most recent update:  October 2, 2023
%
%  Original HiFiPDV release was on December 12, 2020, embedded in PhD 
%  thesis and included in the associated copyright (Travis Voorhees 2020).
% =========================================================================

%% CODE DESCRIPTION & NOTES:
%  This code is intended to be run with parallel CPU computing only
%
%  Other Versions of this Code That Are Available:
%   - GPU computing version (deprecated)
%     -- Requires large amounts of GPU VRAM
%     -- Latest GPU code update is December 2020
%
%  Planned Updates/Upgrades (also see latest release notes):
%   - Ability to load in ROI bounds from previous run.
%   - Choose your favorite spectrogram to put in final plot.
% =========================================================================

%% INPUT VARIABLES INVESTIGATED:
%   - Duration Length (length of timebin)
%       -- 5 to 100 ns, in 5 ns steps
%       -- Duration Overlap is scaled to produce a constant 2.5 ns time-step
%   - Window function
%       -- Hamming
%       -- Hahn
%       -- Blackman
%   - Peak determination method
%       -- Max of power spectrum at time
%       -- Centroid (center of mass) of power spectrum at time
%       -- Peak of Gaussian fit to power spectrum at time
%
%  Total iterations to be performed:  [DurLength]*[WinFun]*[PeakEx]=180
%  NumSpec = 20*3;
%  NumHist = NumSpec*3; % number of histories per window function
% =========================================================================

%% INPUT VARIABLES THAT ARE HELD CONSTANT
%   - Number of frequency bins
%       -- 2^15 = 32,768 bins
%       -- Results in slighter smaller than 1 m/s velocity bin size for max
%          frequency on HSR lab high speed scope
%   - Time-step / point-to-point time resolution
%       -- ~2.5 ns skip between DFT calculations
%       -- Overlap is scaled with duration length accordingly
%       -- NOTE:  Time-step needs to be a factor of all window lengths, so
%                 the time-step may expand beyond 2.5ns to accomodate this.
%                 Code output will tell you if this happens, and explicitly
%                 state the new time-step.
% =========================================================================

%% ASSUMPTIONS MADE IN THIS CALCULATION
%  (1) A single velocity signal is present within the user's ROI bounds
%  (2) Signal spectral width does not collide with neighboring signals
%  (3) No window corrections are applied (can be added if requested)
%  (4) No PDV probe angle corrections are applied (assumed normal to vector)
% =========================================================================

%% LATEST RELEASE NOTES:
% ============================ (October 2, 2023) ===========================
%  Added:
%   - Advanced user options for colormap choice and baseline filtering
%   - More descriptive comments/notes
%   - More consistent output file naming
%
% ============================ (September 1, 2023) ===========================
%  Fixed:
%   - Gaussian fitting procedure for Heatmap's history extraction
%     -- Original method was throwing spurious results occasionally
%
% ============================ (August 16, 2023) ===========================
%  Changed:
%   - More efficient Gaussian fitting process -- thank you Dan Champion!
%   - Spectrum outside region of interest set to zero (instead of -60 dB)
%     -- Helps with new Gaussian fitting process
%   - Moved Cropping parameters outside of Spectrogram parfor loop
%     -- Makes cropping vectors permanent variables
%   - Removed negative velocities from Heatmap calculation
%   - Memory reduction via FP32 or FP16 data storage (laptop friendly).
%
%  Fixed:
%   - Calculation of baseline frequency
%     -- Previous version had an typo in conversion to real space for
%        Blackman and Hann windows
%
%  Beta Testing:
%   - Reduced spectrogram bit precision to FP16 instead of FP64
%     -- Significant RAM reduction, negligible precision reduction
%   - Faster Gaussian fitting program
%   - Better cropping of signal during peak fitting
%
%  Remaining Updates Before Full Version Release:
%   - Updates for down-shifted signals
%     -- Not enough urgency to add this functionality
%
% ============================ (March 6, 2023) ============================
%  Added:
%   - Ability to read in HDF5 (*.h5) files
%     -- Borrows "agilend_h5read.m" from SIRHEN; credit to Dan Dolan
%   - Baseline filtering/removal with a bandstop filter
%
%  Changed:
%   - Renaming "RequiredFunctions" to "private"
%     -- This should make all enclosed scripts automatically in path
%     -- Thanks to Dan Dolan for advice!
%
% ============================ (March 1, 2023) ============================
%  Added:
%   - Option to correct up-shifted signals
%     -- Additional yes/no option in first pop-up window
%     -- Early step to fit baseline with very long durations
%   - All prompts are now case insensitive
%
% ============================ (August 17, 2022) ==========================
%  Added:
%   - Method to read .dig files
%     -- Borrowed from SIRHEN; credit to Dan Dolan
%
% ============================ (December 9, 2021) =========================
%  Added:
%   - Better method to read in .dat and .csv files
%     -- any columns containing NaN data are removed
%
% ============================ (October 25, 2021) =========================
%  Added:
%   - Way to read *.isf files
%
% ============================ (December 12, 2020) ========================
%  HiFiPDV first released as an embedded tarball on Page 181 of PhD thesis.
%      Copyright Travis J Voorhees 2020
%      See license file in /private/associated_licenses/ for details
%
% ============================ (August 7, 2020) ===========================
%  Changed:
%   - Example spectrogram is now made with longest duration length (100ns)
%   - Example spectrogram bin-to-bin resolution is 10 ns
%   - Increased linewidth of output plots to 2
%  Added:
%   - Window correction factor options:  LiF, PMMA, more coming 'someday'
% ============================ (June 16, 2020) ============================
%  Added:
%   - Option to output Rise Times (10-50-90) of autodetected front
%   - You can now select whether or not to run Gaussian fitting process
%   - Improved method for extrapolation of user inputted ROI bounds
%     -- No longer extrapolates beyond greatest or lowest user input
%   - User ROI bounds are now shown in final spectrogram plot
%
%  Removed:
%   - Checking for "MonteCarlo" on the path (more user naming freedom)
%     -- Code still looks for "RequiredFunctions" directory
%
%  Considering working on:
%   - Moving user prompts from a pop-up window to being inside the actual
%     figure UI (not sure how this would affect performance).
%   - Lowering RAM usage (laptop friendly)
% =========================================================================





% *************************************************************************
% ****************************** BEGIN CODE *******************************
% *************************************************************************





%% Start fresh
close all; clear; clc


%% ADVANCED USER OPTIONS:
%  These are advanced options that users can vary if they wish.
%  Default values are:
%   >> SPECTROGRAM_COLORMAP = jet(256); % 'jet' colormap with 8-bit color
%   >>     HEATMAP_COLORMAP = hot(256); % 'hot' colormap with 8-bit color
%   >>             FILT_AMT = 0.01;     % ±1 percent baseline filter

SPECTROGRAM_COLORMAP = jet(256); % colormap for spectrograms
HEATMAP_COLORMAP     = hot(256); % colormap for heatmap
FILT_AMT = 0.01;                 % ± percentage of baseline filter, 0.01 = ±1 % filter


%% Initiate a parallel CPU pool if one hasn't already been setup
Pool = gcp; % Get Current Parallel pool


%% Read in the signal

% Use a GUI to select the PDV data file
%   Filetypes supported:  *.trc, *.csv, *.txt
[File,Path] = uigetfile({'*.trc;*.csv;*.txt;*.wfm;*.isf;*.dig;*.h5;*.H5'},'MultiSelect','off');
if File==0 % User did not select anything or hit cancel
    clear;
    error('User selected cancel or closed window. Process terminating now.');
else
    [~,~,FileFormat] = fileparts([Path File]); % updated from just '*.xxx'
    if ~FileFormat(1) == '.'
        error('File format could not be determined.')
    else
    end
    disp(['User selected file:  ' File]);
end

% Set up a GUI for user to input quick cropping parameters & Gauss select
prompt = {'Time Units (s,us,ns):','Approximate Start Time (μs):',...
    'Approximate End Time (μs):','Laser Wavelength (nm):',...
    'Run Gaussian Fit? (yes/no)', 'Calculate 10-50-90 Rise Times (yes/no)?',...
    'Window Material:', 'Correct for frequency-shifting? (up/down/no)',...
    'Filter out signal baseline? (yes/no)'};
dlgtitle = 'Data Import Settings';
dims = [1 45];
definput = {'s','0','300','1550','yes','no','none','no','no','no'};

% Launch a GUI for the user to input cropping parameters
C = inputdlg(prompt,dlgtitle,dims,definput);

% Determine Conversion Factor to Convert Time-Scale to Nanoseconds
switch lower(C{1})
    case {'s','seconds','sec','1'}
        ConvFactor = 1e9;
    case {'ns','nanoseconds','nano','1e-9'}
        ConvFactor = 1;
    case {'\mus','mus','microseconds','micro','1e-6'}
        ConvFactor = 1e3;
    case {'ms','mili','milli','miliseconds','milliseconds','1e-3'}
        FirstTimePoint = Data(1,1)*1e3;
        Time = Data(:,1).*1e3;
        disp('Are you sure your time format is milliseconds?');
    otherwise
        error('foo:bar',['========== USER DID NOT ENTER USEABLE TIMESCALE ==========' ...
            '\nACCEPTABLE TIMESCALES: s, ns, us, \mus, seconds, nanoseconds, microseconds, nano, micro' ...
            '\nNOTE:  TIME-SCALE ENTRIES ARE CASE-SENSITIVE, USE LOWER-CASE']);
end

% Determine if user wants to run Gaussian fitting or not
switch lower(C{5})
    case {'yes','y','1','yes I want to','yes please','true'}
        RunGauss = 1;
        disp('Gaussian fitting process selected.');
    case {'no','n','0','no thank you','no I do not want to','false'}
        RunGauss = 0;
        disp('Gaussian fitting process not selected.');
    case {''}
        RunGauss = 0;
        disp('User left Gaussian fitting question blank.');
        disp('Gaussian fitting process will not occur.');
    otherwise
        error('foo:bar',['========== GAUSSIAN FITTING ANSWER NOT ACCEPTED ==========' ...
            '\nACCEPTABLE ANSWERS: y, n, yes, no, Yes, No, YES, NO, etc.' ...
            '\nUser entered the following answer:  ' C{5}]);
end


% Determine the window correction factor (if needed)
switch lower(C{7})
    case {'yes','y','not sure'}
        WFdlg = inputdlg('What is the window correction factor?','Enter Window Correction Factor',[1 45],{'0.7895'});
        if WFdlg == ''
            error('User did not enter a window correction factor when asked');
        else
        end
        WF.a = str2double(WFdlg);
        WF.b = 1;
    case {'no','NO','none','No','N/A'}
        disp('No window correction factor applied')
        WF.a = 1;
        WF.b = 1;
    case {'LiF'}
        WF.a = 0.7895;
        WF.b = 0.9918;
    case {'PMMA'}
        disp('  ** PMMA requires a negligible window correction factor.')
        WF.a = 1; % negligible
        WF.b = 1; % negligible
    otherwise
        WFdlg = inputdlg('What is the window correction factor?','Enter Window Correction Factor',[1 45],'0.7895');
        if WFdlg == ''
            error('User did not enter a window correction factor when asked');
        else
        end
        WF.a = str2double(WF);
        WF.b = 1;
end


% Determine if the user wants to correct for frequency shifting or not:
switch lower(C{8})
    % The user might input 'yes', which is not descriptive enough
    % Ask them to specify if the signal is "up" or "down" shifted:
    case {'yes','y','ya'}
        FSdlg = inputdlg('Is the signal up-shifted or down-shifted?','Is the signal up-shifted or down-shifted?',[1 45],{'up'});
        if strcmpi(FSdlg,'')
            warning('User did not state if signal was up- or down-shifted');
            warning('Proceeding with analysis, assuming signal is not frequency shifted');
            warning('RESTART ANALYSIS IF THIS IS INCORRECT');
            C{8} = 'no';
        else
            C{8} = FSdlg{1};
        end
        if strcmpi(FSdlg,'yes')
            error('Please retry.  "Yes" is not a valid answer for direction of frequency shifting.')
        end
end
switch lower(C{8})
    case {'no','none','n/a'}
        disp('Signal is *NOT* frequency-shifted.')
    case {'up','up-shifted','UP','Up','upshifted'}
        SHIFT_SIGN = 1;  % will be used in later section
    case {'down','downshifted','down-shifted','DOWN','Down'}
        SHIFT_SIGN = -1;  % will be used in later section
        warning('Code is not well suited for down-shifted signals yet');
    otherwise
        warning(['Frequency-shifting answer could not be interpreted.  User input: ' C{8}])
        error('Please retry.')
end


% Determine if the user wants to filter out the baseline
switch lower(C{9})
    case {'yes','y','ya','true'}
        disp('User selected to filter baseline.')
        disp('Bandstop filter will be applied from 99-101% of baseline')
        BaseRemoval = 1;
    case {'no','n','false','none','n/a'}
        disp('User selected not to filter baseline from signal')
        BaseRemoval = 0;
    otherwise
        warning('Bad input for baseline filtering question')
        error(['When asked if signal should be filtered, user input: ' C{9}])
end

% Read in the user inputs
Start = str2double(C{2})*1e3; % convert to nanoseconds
End = str2double(C{3})*1e3; % convert to nanoseconds
Wavelength = str2double(C{4});

% Read in the data
switch FileFormat
    case '.trc'
        data = ReadLeCroyBinaryWaveform([Path File]);
        Signal = [data.x data.y];
    case '.txt'
        Signal = textscan([Path File],'%f%f');
        Signal(any(isnan(Signal), 2), :) = []; % remove any rows with NaN
    case {'.csv','.dat'}
        Signal = readmatrix([Path File]);
        Signal(:,isnan(mean(Signal))) = [];
        warning("CSV data read in such that column 1 is time, column 2 is voltage:");
        disp(Signal(1:10,:));
        disp("...")
        warning("If this is incorrect, please reformat input data")
        if width(Signal) == 1
            error("Only one column of data was read in.  Time and Voltage is needed.")
        else
        end 
    case '.wfm'
        [data.y,data.x] = wfm2read([Path File]);
        Signal = [data.x data.y];
    case '.isf'
        data = isfread([Path File]);
        Signal = [data.x data.y];
    case '.dig'
        [data.y,data.x] = digread([Path File]);
        Signal = [data.x data.y];
    case {'.h5','.H5'}
        % HDF5 can hold a lot of channels, have the user select one
        H5dlg = inputdlg('HDF5 Detected. Which Channel are you interested in?',...
            'HDF5 Detected, Choose Channel.',[1 55],{'1'});
        % If user does not input a selection, default to channel 1
        if strcmpi(H5dlg,'')
            warning('User did not select a channel. Default is 1.');
            warning('Proceeding analysis with Channel 1 selected.');
            warning('** RESTART ANALYSIS IF THIS IS INCORRECT **');
            Channel = 1;
        else
            % Try to evaluate user's channel selection
            % Throw an error if non-numeric value is entered
            try
                Channel = eval(H5dlg{1});
            catch
                error(['User input "' H5dlg{1} '" as channel selection']);
            end
        end
        try
            % Try to read the data
            [data.y,data.x] = agilent_h5read([Path File],Channel);
            Signal = [data.x data.y];
            warning('Reading HDF5 format is a *NEW* feature in HiFiPDV.');
            warning('Contact tjvoorh@sandia.gov if bugs arise.');
        catch
            % Throw an error if data cannot be read.
            error(['Cannot read Channel ' num2str(Channel) ' in ' Path File])
        end
    otherwise
        error('foo:bar',['===== USER DID NOT CHOOSE COMPATIBLE FILE FORMAT =====' ...
            '\nACCEPTABLE FORMATS: .trc, .csv, .txt, .wfm, .h5']);
end

% Convert data timescale to nanoseconds
Signal(:,1) = Signal(:,1).*ConvFactor;

% Determine Signal Data Sampling Rate (must be done prior to cropping)
SampleRate = 1e9/(mean(diff(Signal(:,1)))); %1e9 b/c time is in nanoseconds
% if Sample Rate is negative, then flip data across the X-axis
if sign(SampleRate)==-1
    Signal = flipud(Signal);
    SampleRate = 1e9/(mean(diff(Signal(:,1)))); %recalculate
else
end

% Spit out an error if sample rate is slower than 2 points/ns
if SampleRate<2.49e9
    error('foo:bar',['Signal Sample Rate is ' num2str(SampleRate./1e9,'%.2e') ...
        ' points/nanosecond\nSample Rate must be 2.00 points/nanosecond or' ...
        ' greater for this code to run correctly.'])
else
end

% Crop the data
[~,FirstInd] = min(abs(Signal(:,1)-Start)); % Finds closest value to start limit
[~,LastInd]=min(abs(Signal(:,1)-End)); % Finds closest value to end limit
Signal = Signal(FirstInd:LastInd,:);

% Time-shift signal
FirstTimePoint = Signal(1,1);
Signal(:,1) = Signal(:,1)-FirstTimePoint; % nanoseconds, starting at zero


%% Create STFT inputs and pre-allocate arrays
% REMEMBER:  time axis is now in nanoseconds, starting at zero
% First time point is "FirstTimePoint" (also in nanoseconds)

% Detector voltage is second column of signal, time is first column
DetectorVoltage = Signal(:,2);

% Determine the number of data points in a 1 nanosecond duration
Num1ns = SampleRate/1e9;
NumDataPointsFor1ns = round(Num1ns); % needs to be integer value

% NumDataPointsFor1ns needs to be even
% if this value is odd, then add 1 and let the user know the new time-resolution
if mod(NumDataPointsFor1ns,2) == 1 % odd case
    disp(['Number of data points per 1 ns is: ' num2str(Num1ns,'%.1f')]);
    NumDataPointsFor1ns = NumDataPointsFor1ns+1;
    disp(['Point-to-point time resolution set to: ' num2str(2.5*(NumDataPointsFor1ns...
        /Num1ns),'%.4f') 'ns']);
else % Even case
    disp(['Number of data points per 1 ns is: ' num2str(Num1ns,'%.1f')]);
    disp(['Point-to-point time resolution set to: ' num2str(2.5*(NumDataPointsFor1ns...
        /Num1ns),'%.4f') 'ns']);
end


% Create a range of duration lengths: ~5 ns to ~100 ns
RangeOfDurationLengths = [5:5:100].*NumDataPointsFor1ns;

% Define overlap as duration - 2 ns, this creates a 2 ns time-step
RangeOfOverlaps = RangeOfDurationLengths-round(2.5*NumDataPointsFor1ns);

% Create a variable for number of frequency bins
nFreqs = 2^15; % Reduce this value if RAM is being maxed out



%% Determine the reference baseline (for frequency shifted signals)
% If the signal is frequency-shifted, have the user draw lines around the
% baseline region (the region that exists before the region of interest)


% Only perform this process if the SHIFT_SIGN variable exists
if exist('SHIFT_SIGN')
    BaseClr = 'm';
    % Make temporary inputs for spectrogram
    TmpWin = hamming(RangeOfDurationLengths(end)); % longest window
    TmpOvr = RangeOfOverlaps(end); % corresponding overlap

    % Calculate a temporary spectrogram
    [zexample,f,texample] = spectrogram(DetectorVoltage,TmpWin,TmpOvr,nFreqs,SampleRate);
    zexample = real(zexample.*conj(zexample)); % get rid of imaginary numbers
    zexample = zexample./max(zexample(:)); % normalize globally
    zexample = 10.*log10(zexample); % convert to dB for plot
    texample = texample.*1e6 + FirstTimePoint/1e3; % convert time to microseconds and shift
    f = f./1e9; % convert to GHz
    
    % Show the spectrogram, get user input for cropping prior to lineout calcs
    figure
    hold on
    imagesc(texample,f.*(Wavelength/2),zexample)
    %colormap(jet(16)) % option for less color depth if X11 forwarding
    colormap(SPECTROGRAM_COLORMAP)
    colorbar
    caxis([-60 0])
    axis([-inf +inf -25 +inf]) % we don't usually see velocities over 5,000 m/s
    xlabel('Detector Time (\mus)')
    ylabel('Signal Velocity (m/s)')
    title('Reference/Baseline Selection')
    hold off
    
    % Have the user zoom into their region of interest
    zoom on
    Prompt{1} = 'ZOOM to your BASELINE/REFERENCE REGION';
    Prompt{2} = '';
    Prompt{3} = ' - LEFT-CLICK and DRAG a box around this region.';
    Prompt{4} = '';
    Prompt{5} = ' - PRESS ENTER when you are finished.';
    Prompt{6} = '';
    Prompt{7} = ' - DOUBLE LEFT-CLICK to undo/reset zoom level.';
    mb = msgbox(Prompt,'Zoom to REFERENCE region');
    pause
    zoom off
    delete(mb)
    clear mb
    
    % Upper Bound of Baseline
    PromptUPb{1} = 'Draw an UPPER BOUND for your BASELINE.';
    PromptUPb{2} = '';
    PromptUPb{3} = ' - Full screen viewing is recommended for higher accuracy.';
    PromptUPb{4} = '';
    PromptUPb{5} = ' - LEFT-CLICK to place points, RIGHT-CLICK when you are finished.';
    mb = msgbox(PromptUPb,'Draw an UPPER BOUND for your BASELINE');
    roiUPb = images.roi.Polyline('color',BaseClr);
    draw(roiUPb);
    ROIUPb = roiUPb.Position;
    ROIUPb = sortrows(ROIUPb); % must be in order for future interpolation
    % Clean up the ROI and replot it
    %   Don't allow multiple y-values per x-value
    ROIxstepUPb = ROIUPb(2:end,1)-ROIUPb(1:(end-1),1);
    if mean(ROIxstepUPb)>0 % moving rightwards
        ROIxstepUPb = [1;ROIxstepUPb];
        ROIUPb(ROIxstepUPb<=0,:)=[];
    else
        ROIxstepUPb = [-1;ROIxstepUPb];
        ROIUPb(ROIxstepUPb>=0,:)=[]; % moving leftwards
    end
    roiUPb.Position=ROIUPb; %update the drawing on the figure
    ROIUPxb = ROIUPb(:,1).*1e3 -FirstTimePoint; % convert to nanoseconds
    ROIUPyb = ROIUPb(:,2); % frequency in GHz
    delete(mb)
    clear mb

    % Lower Bound of Baseline
    PromptLOWb = PromptUPb;
    PromptLOWb{1} = 'Draw a LOWER BOUND for your BASELINE.';
    mb = msgbox(PromptLOWb,'Draw LOWER BOUND for BASELINE');
    roiLOWb = images.roi.Polyline('color',BaseClr);
    draw(roiLOWb);
    ROILOWb = roiLOWb.Position;
    ROILOWb = sortrows(ROILOWb); % must be in order for future interpolation
    % Clean up the ROI and replot it
    %   Don't allow multiple y-values per x-value
    ROIxstepLOWb = ROILOWb(2:end,1)-ROILOWb(1:(end-1),1);
    if mean(ROIxstepLOWb)>0 % moving rightwards
        ROIxstepLOWb = [1;ROIxstepLOWb];
        ROILOWb(ROIxstepLOWb<=0,:)=[];
    else
        ROIxstepLOWb = [-1;ROIxstepLOWb];
        ROILOWb(ROIxstepLOWb>=0,:)=[]; % moving leftwards
    end
    roiLOWb.Position=ROILOWb; %update the drawing on the figure
    ROILOWxb = ROILOWb(:,1).*1e3 - FirstTimePoint; % convert to nanoseconds, shift
    ROILOWyb = ROILOWb(:,2); % frequency in GHz
    delete(mb)
    clear mb
    
    
    % Check if user is okay with their frequency bounds before proceeding
    UserPrompt2 = questdlg('Are you satisfied with the current BASELINE frequency bounds?',...
        'Use These Bounds?','Yes','No','Cancel','No');
    switch UserPrompt2
        case 'Yes'
            disp('Baseline freqency bounds entered, beginning baseline extraction process.');
            disp(' ');
        case 'No'
            while strcmpi(UserPrompt2,'No')
                delete(roiUPb);
                roiUPb = images.roi.Polyline('color',BaseClr);
                delete(roiLOWb);
                roiLOWb = images.roi.Polyline('color',BaseClr);
                
                % Redo the upper bound
                mb = msgbox(PromptUPb,'Draw UPPER BOUND for your BASELINE');
                draw(roiUPb);
                ROIUPb = roiUPb.Position;
                ROIUPb = sortrows(ROIUPb); % must be in order for future interpolation
                % Clean up the ROI and replot it
                %   Don't allow multiple y-values per x-value
                ROIxstepUPb = ROIUPb(2:end,1)-ROIUPb(1:(end-1),1);
                if mean(ROIxstepUPb)>0 % moving rightwards
                    ROIxstepUPb = [1;ROIxstepUPb];
                    ROIUPb(ROIxstepUPb<=0,:)=[];
                else
                    ROIxstepUPb = [-1;ROIxstepUPb];
                    ROIUPb(ROIxstepUPb>=0,:)=[]; % moving leftwards
                end
                roiUPb.Position=ROIUPb; %update the drawing on the figure
                ROIUPxb = ROIUPb(:,1).*1e3 -FirstTimePoint; % convert to nanoseconds
                ROIUPyb = ROIUPb(:,2); % frequency in GHz
                delete(mb)
                clear mb

                % Redo the lower bound
                mb = msgbox(PromptLOWb,'Draw LOWER BOUND for your BASELINE');
                roiLOWb = images.roi.Polyline('color',BaseClr);
                draw(roiLOWb);
                ROILOWb = roiLOWb.Position;
                ROILOWb = sortrows(ROILOWb); % must be in order for future interpolation
                % Clean up the ROI and replot it
                %   Don't allow multiple y-values per x-value
                ROIxstepLOWb = ROILOWb(2:end,1)-ROILOWb(1:(end-1),1);
                if mean(ROIxstepLOWb)>0 % moving rightwards
                    ROIxstepLOWb = [1;ROIxstepLOWb];
                    ROILOWb(ROIxstepLOWb<=0,:)=[];
                else
                    ROIxstepLOWb = [-1;ROIxstepLOWb];
                    ROILOWb(ROIxstepLOWb>=0,:)=[]; % moving leftwards
                end
                roiLOWb.Position=ROILOWb; %update the drawing on the figure
                ROILOWxb = ROILOWb(:,1).*1e3 - FirstTimePoint; % convert to nanoseconds, shift
                ROILOWyb = ROILOWb(:,2); % frequency in GHz
                delete(mb)
                clear mb

                % Ask if the new bounds are satisfactory
                UserPrompt2 = questdlg('Are you satisfied with the new BASELINE frequency bounds?',...
                'Use These Bounds?','Yes','No','Cancel','No');
            end
        case 'Cancel'
            close all; clear; clc;
            error('User selected cancel. Terminating processes now.')
    end

    % Convert y-axis values to GHz
    ROIUPyb = ROIUPyb./(Wavelength/2);
    ROILOWyb = ROILOWyb./(Wavelength/2);

    close all % close the figure

    % Crop out the baseline for further analysis
    % Determine minimum and maximum time in drawn ROI limits
    BasexMin = min([ROIUPxb;ROILOWxb]);
    BasexMax = max([ROIUPxb;ROILOWxb]);
    BaseyMin = min([ROIUPyb;ROILOWyb]);
    BaseyMax = max([ROIUPyb;ROILOWyb]);

    % Crop signal to ROI left and right limits
    [~,BaseMinInd] = min(abs(Signal(:,1)-BasexMin));
    [~,BaseMaxInd] = min(abs(Signal(:,1)-BasexMax));
    BASELINE=Signal(BaseMinInd:BaseMaxInd,:);

    % Redefine DetectorVoltage with new time limits
    BaseTime = BASELINE(:,1);
    BaseVoltage = BASELINE(:,2);

    % Calculate a SINGLE velocity within the baseline region of interest
    nFreqsB = nFreqs*10; % 10 fold increase for the baseline determination
    BASEfft.x   = (1:nFreqs).*SampleRate./nFreqsB./1e9; % Convert to GHz
    BASEfft.ham = fft(BaseVoltage.*hamming(length(BaseVoltage)),nFreqsB);  % FFT with hamming window
    BASEfft.blk = fft(BaseVoltage.*blackman(length(BaseVoltage)),nFreqsB); % FFT with blackman window
    BASEfft.han = fft(BaseVoltage.*hann(length(BaseVoltage)),nFreqsB);     % FFT with hann window
    
    % Convert to real space
    BASEfft.ham = real(BASEfft.ham.*conj(BASEfft.ham)); % abs needs to be applied FIRST
    BASEfft.blk = real(BASEfft.blk.*conj(BASEfft.blk));
    BASEfft.han = real(BASEfft.han.*conj(BASEfft.han));
    
    % Crop each of the FFTs to just the region of interest (user input)
    [~,BASEfft.Low]  = min(abs(BASEfft.x-BaseyMin)); % find the low index
    [~,BASEfft.High] = min(abs(BASEfft.x-BaseyMax)); % find the high index
    BASEfft.ham = BASEfft.ham(BASEfft.Low:BASEfft.High); % crop
    BASEfft.blk = BASEfft.blk(BASEfft.Low:BASEfft.High); % crop
    BASEfft.han = BASEfft.han(BASEfft.Low:BASEfft.High); % crop
    BASEfft.x   = BASEfft.x(BASEfft.Low:BASEfft.High);   % crop
    
    % Normalize after cropping
    BASEfft.ham = BASEfft.ham/max(max(BASEfft.ham));
    BASEfft.blk = BASEfft.blk/max(max(BASEfft.blk));
    BASEfft.han = BASEfft.han/max(max(BASEfft.han));
    
    % Extract the peak location from each FFT
    BASEfft.MAX = zeros(1,3); % max extraction
    BASEfft.CEN = zeros(1,3); % centroid extraction
    BASEfft.GAS = zeros(1,3); % Gaussian extraction
    BASEfft.ALL = zeros(1,9); % 3 windows, 3 extractions, 9 answers
    %  MAXIMUM
    [~,tmp] = max(BASEfft.ham);
    BASEfft.MAX(1) = BASEfft.x(tmp); clear tmp
    [~,tmp] = max(BASEfft.blk);
    BASEfft.MAX(2) = BASEfft.x(tmp); clear tmp
    [~,tmp] = max(BASEfft.han);
    BASEfft.MAX(3) = BASEfft.x(tmp); clear tmp
    % CENTROID
    % f_centerofmass is equal to {int(f*Z)df / int(Z)df}
    % In other words:  integral of curve*x / area under curve
    BASEfft.CEN(1) = trapz(BASEfft.x,BASEfft.x.*BASEfft.ham')./...
        trapz(BASEfft.x,BASEfft.ham);
    BASEfft.CEN(2) = trapz(BASEfft.x,BASEfft.x.*BASEfft.blk')./...
        trapz(BASEfft.x,BASEfft.blk);
    BASEfft.CEN(3) = trapz(BASEfft.x,BASEfft.x.*BASEfft.han')./...
        trapz(BASEfft.x,BASEfft.han);
    % Gaussian fit
    % NOTE:  Uses a different Gaussian fit than Spectrogram and Heatmap extraction
    %  Why?  Old code that hasn't proven to be an issue (yet).  May update in next release.
    BaseFit.Xf = BASEfft.x(1):(0.01*2/Wavelength):BASEfft.x(end);
    BaseFit.Xv = BaseFit.Xf.*(Wavelength/2);
    tmp = fit(BASEfft.x',BASEfft.ham,'gauss1');
    BaseFit.ham = tmp(BaseFit.Xf); % do this for plots later
    BASEfft.GAS(1) = tmp.b1; clear tmp
    tmp = fit(BASEfft.x',BASEfft.blk,'gauss1');
    BaseFit.blk = tmp(BaseFit.Xf); % do this for plots later
    BASEfft.GAS(2) = tmp.b1; clear tmp
    tmp = fit(BASEfft.x',BASEfft.han,'gauss1'); 
    BaseFit.han = tmp(BaseFit.Xf); % do this for plots later
    BASEfft.GAS(3) = tmp.b1; clear tmp
    
    BASEfft.ALL = [BASEfft.MAX BASEfft.CEN BASEfft.GAS];
    SHIFT_AMT = SHIFT_SIGN*mean(BASEfft.ALL)*Wavelength/2; % converted to velocity
    SHIFT_UNC = 1.95*std(BASEfft.ALL)*Wavelength/2;
    
    Basefilter = [1-FILT_AMT 1+FILT_AMT].*(SHIFT_AMT*2/Wavelength).*1e9; % Bandstop filter limits
    
else
    disp(''); % signal is not shifted
    Basefilter = [0 0.5e8]; % bandstop filter limits (0 to 50 MHz/39 m/s)
    % NOTE TO USERS:
    %  - feel free to change this value, it is arbitrary
    %  - Remember:  Velocity = Signal_Frequency * (Laser_Wavelength/2)
end




%% Apply a bandstop filter to the filter out the baseline

if BaseRemoval == 1
    DetectorVoltage = bandstop(DetectorVoltage,Basefilter,SampleRate); % simple bandstop filter
else
end


%% User input Region of Interest for next step in calculations
% Also have them select the baseline if data is frequency shifted

% Define line colors for these steps
UpClr = 'w';   % line color for signal region of interest
LowClr = 'w';  % line color for signal region of interest

% Make temporary inputs for spectrogram
IndexForExample = 15; % default is end, but you can select a different spectrogram for the example
%  NOTE:  Consider allowing users to choose their own spectrogram for ROI
TmpWin = hamming(RangeOfDurationLengths(IndexForExample)); % Use Hamming for example
TmpOvr = RangeOfOverlaps(IndexForExample); % Corresponding overlap

% Calculate a temporary spectrogram
[zexample,f,texample] = spectrogram(DetectorVoltage,TmpWin,TmpOvr,nFreqs,SampleRate);
zexample = real(zexample.*conj(zexample)); % get rid of imaginary numbers
zexample = zexample./max(zexample(:)); % normalize globally
zexample = 10.*log10(zexample); % convert to dB for plot
texample = texample.*1e6 + FirstTimePoint/1e3; % convert time to microseconds and shift
f = f./1e9; % convert to GHz

% Show the spectrogram, get user input for cropping prior to lineout calcs
figure
hold on
imagesc(texample,f.*(Wavelength/2),zexample)
colormap(SPECTROGRAM_COLORMAP)
%colormap(jet(16)) % option for less color depth if X11 forwarding
colorbar
caxis([-60 0])
axis([-inf +inf -25 +inf]) % we don't usually see velocities over 5,000 m/s
xlabel('Detector Time (\mus)')
ylabel('Signal Velocity (m/s)')
hold off


% Have the user zoom into their region of interest
zoom on
Prompt{1} = 'ZOOM to your SIGNAL region of interest';
Prompt{2} = '';
Prompt{3} = ' - LEFT-CLICK and DRAG a box around your SIGNAL region of interest.';
Prompt{4} = '';
Prompt{5} = ' - PRESS ENTER when you are finished.';
Prompt{6} = '';
Prompt{7} = ' - DOUBLE LEFT-CLICK to undo/reset zoom level.';
Prompt{8} = '';
Prompt{9} = ' Note:  If your signal is frequency-shifted, include baseline';
mb = msgbox(Prompt,'Zoom to SIGNAL region of interest');
pause
zoom off
delete(mb)
clear mb



% =====
% User input for frequency limits/region of interest
% =====
% Upper Bound of Actual Signal
PromptUP{1} = 'Draw an UPPER BOUND for your SIGNAL.';
PromptUP{2} = '';
PromptUP{3} = ' - Full screen viewing is recommended for higher accuracy.';
PromptUP{4} = '';
PromptUP{5} = ' - LEFT-CLICK to place points, RIGHT-CLICK when you are finished.';
mb = msgbox(PromptUP,'Draw UPPER BOUND for Calculation');
roiUP = images.roi.Polyline('color',UpClr);
draw(roiUP);
ROIUP = roiUP.Position;
ROIUP = sortrows(ROIUP); % must be in order for future interpolation
% Clean up the ROI and replot it
%   Don't allow multiple y-values per x-value
ROIxstepUP = ROIUP(2:end,1)-ROIUP(1:(end-1),1);
if mean(ROIxstepUP)>0 % moving rightwards
    ROIxstepUP = [1;ROIxstepUP];
    ROIUP(ROIxstepUP<=0,:)=[];
else
    ROIxstepUP = [-1;ROIxstepUP];
    ROIUP(ROIxstepUP>=0,:)=[]; % moving leftwards
end
roiUP.Position=ROIUP; %update the drawing on the figure
ROIUPx = ROIUP(:,1).*1e3 -FirstTimePoint; % convert to nanoseconds
ROIUPy = ROIUP(:,2); % frequency in GHz
delete(mb)
clear mb

% Lower Bound of Actual Signal
PromptLOW = PromptUP;
PromptLOW{1} = 'Draw a LOWER BOUND for your SIGNAL.';
mb = msgbox(PromptLOW,'Draw LOWER BOUND for Calculation');
roiLOW = images.roi.Polyline('color',LowClr);
draw(roiLOW);
ROILOW = roiLOW.Position;
ROILOW = sortrows(ROILOW); % must be in order for future interpolation
% Clean up the ROI and replot it
%   Don't allow multiple y-values per x-value
ROIxstepLOW = ROILOW(2:end,1)-ROILOW(1:(end-1),1);
if mean(ROIxstepLOW)>0 % moving rightwards
    ROIxstepLOW = [1;ROIxstepLOW];
    ROILOW(ROIxstepLOW<=0,:)=[];
else
    ROIxstepLOW = [-1;ROIxstepLOW];
    ROILOW(ROIxstepLOW>=0,:)=[]; % moving leftwards
end
roiLOW.Position=ROILOW; %update the drawing on the figure
ROILOWx = ROILOW(:,1).*1e3 - FirstTimePoint; % convert to nanoseconds, shift
ROILOWy = ROILOW(:,2); % frequency in GHz
delete(mb)
clear mb


% Check if user is okay with their frequency bounds before proceeding
UserPrompt2 = questdlg('Are you satisfied with the current signal frequency bounds?',...
    'Use These Bounds?','Yes','No','Cancel','No');
switch UserPrompt2
    case 'Yes'
        disp('Freqency bounds entered, beginning extraction process.');
    case 'No'
        while strcmpi(UserPrompt2,'No')
            delete(roiUP);
            roiUP = images.roi.Polyline('color',UpClr);
            delete(roiLOW);
            roiLOW = images.roi.Polyline('color',LowClr);
            mb = msgbox(PromptUP,'Draw UPPER BOUND for your SIGNAL');
            
            % re-do the upper bound
            draw(roiUP);
            ROIUP = roiUP.Position;
            ROIUP = sortrows(ROIUP); % must be in order for future interpolation
            % Clean up the ROI and replot it
            %   Don't allow multiple y-values per x-value
            ROIxstepUP = ROIUP(2:end,1)-ROIUP(1:(end-1),1);
            if mean(ROIxstepUP)>0 % moving rightwards
                ROIxstepUP = [1;ROIxstepUP];
                ROIUP(ROIxstepUP<=0,:)=[];
            else
                ROIxstepUP = [-1;ROIxstepUP];
                ROIUP(ROIxstepUP>=0,:)=[]; % moving leftwards
            end
            roiUP.Position=ROIUP; %update the drawing on the figure
            ROIUPx = ROIUP(:,1).*1e3 - FirstTimePoint; % convert to nanoseconds, shift
            ROIUPy = ROIUP(:,2); % frequency in GHz
            delete(mb)
	    clear mb
            
            % re-do the lower bound
            mb = msgbox(PromptLOW,'Draw LOWER BOUND for your SIGNAL');
            draw(roiLOW);
            ROILOW = roiLOW.Position;
            ROILOW = sortrows(ROILOW); % must be in order for future interpolation
            % Clean up the ROI and replot it
            %   Don't allow multiple y-values per x-value
            ROIxstepLOW = ROILOW(2:end,1)-ROILOW(1:(end-1),1);
            if mean(ROIxstepLOW)>0 % moving rightwards
                ROIxstepLOW = [1;ROIxstepLOW];
                ROILOW(ROIxstepLOW<=0,:)=[];
            else
                ROIxstepLOW = [-1;ROIxstepLOW];
                ROILOW(ROIxstepLOW>=0,:)=[]; % moving leftwards
            end
            roiLOW.Position=ROILOW; %update the drawing on the figure
            ROILOWx = ROILOW(:,1).*1e3 - FirstTimePoint; % convert to nanoseconds, shift
            ROILOWy = ROILOW(:,2); % frequency in GHz
            delete(mb)
	    clear mb
            
            % Ask if the new bounds are satisfactory
            UserPrompt2 = questdlg('Are you satisfied with the new frequency bounds?',...
            'Use These Bounds?','Yes','No','Cancel','No');
        end
    case 'Cancel'
        close all; clear; clc;
        error('User selected cancel. Terminating processes now.')
end

% Convert ROIUPy and ROILOWy to GHz
ROIUPy = ROIUPy./(Wavelength/2);
ROILOWy = ROILOWy./(Wavelength/2);
    
close all % close the figure

% clear the temporary variables to free up memory
%   Don't clear f, it will be used in all future calcs
%   Don't clear texample or zexample, they will be used in an output plot
clear TmpWin TmpOvr roi Prompt Prompt2 ROI 



%% Crop the signal time according to the ROI drawn by the user
%   Doing this will reduce memory costs in spectrogram calculations

% Determine minimum and maximum time in drawn ROI limits
ROIxMin = min([ROIUPx;ROILOWx]);
ROIxMax = max([ROIUPx;ROILOWx]);

% Crop signal to ROI left and right limits
[~,LeftInd] = min(abs(Signal(:,1)-ROIxMin));
[~,RightInd] = min(abs(Signal(:,1)-ROIxMax));
Signal=Signal(LeftInd:RightInd,:);

% Redefine DetectorVoltage with new time limits
DetectorVoltage = Signal(:,2);

% Determine indices of frequency that are closest to bounds
% Preallocate variables populated in for loop:
fIndUP = zeros(1,length(ROIUPx));
fIndLOW = zeros(1,length(ROILOWx));
%   Upper ROI bound
parfor i = 1:length(ROIUPx)
    [~,fIndUP(i)] = min(abs(f-ROIUPy(i)));
end
%   Lower ROI bound
parfor i = 1:length(ROILOWx)
    [~,fIndLOW(i)] = min(abs(f-ROILOWy(i)));
end
% NOTE: frequency indices can be used directly in spectrograms
BottomInd = min([fIndUP fIndLOW]);
UpperInd = max([fIndUP fIndLOW])+1; % add 1 to be conservative

% Vertically crop f and zexample, this will be useful later
zexample = zexample(1:UpperInd,:);
f = f(1:UpperInd,:); % Keep f in double precision

% Re-normalize zexample:
%zexample = zexample./max(max(zexample));

% Redefine ROI Upper and Lower bound frequencies using f
ROIUPy = f(fIndUP);
ROILOWy = f(fIndLOW); % remove (-) numbers, will be helpful in future steps





%% Prep for calculation loop:  create progress bar and pre-allocate cells
% 60 varieties of spectrograms are to be calculated
% output variables will update after every calculation
% individual spectrograms are saved in cell format
% Code structure is optimized for "parfor" parallel computing

% Define loop variables
in = length(RangeOfDurationLengths);

% pre-allocate calculation cell-type variables
%   Note:  without pre-allocation, code run time is up to 5x longer
if RunGauss == 1
    FHAM = cell(3,in); % three extraction methods (max, centroid, gauss)
else
    FHAM = cell(2,in); % two extraction methods (max, centroid)
end
FHANN = FHAM;
FB = FHAM;
zHam = cell(1,in);
zB = zHam;
zHann = zHam;
t = zHam;




%% Begin the Calculation Loops - intended for parallel CPU
% Start a timer to say how long extraction took
ExtractStart = tic;

%% Calculate Spectrograms, Crop Data, Time-shift, & Account for ROI bounds
%   NOTE: f (y-axis) was already calculated in example problem
%   NOTE: t is the middle of the time-bin
%       Since longer durations start later and end earlier, an additional
%       step will be necessary after this loop to crop the left/right edges

% Make a progress bar for user feedback
%   NOTE: parfor makes doing this somewhat more complex
WB = waitbar(0,'Calculating spectrograms...','Name','Extracting Data from Signal');

disp('If the progress bar does not delete itself, please run the following commands:');
disp("First:  >> F = findall(0,'type','figure','tag','TMWWaitbar')");
disp("Then:   >> delete(F)");


SpecTimeStart = tic;
parfor i = 1:in
    %% Calculate the Spectrograms
    
    % Double precision:
    [zHam{i},~,t{i}] = spectrogram(DetectorVoltage,hamming(RangeOfDurationLengths(i)),...
        RangeOfOverlaps(i),nFreqs,SampleRate);
    zB{i} = spectrogram(DetectorVoltage,blackman(RangeOfDurationLengths(i)),...
        RangeOfOverlaps(i),nFreqs,SampleRate);
    zHann{i} = spectrogram(DetectorVoltage,hann(RangeOfDurationLengths(i)),...
        RangeOfOverlaps(i),nFreqs,SampleRate);
    
    % Convert to spectrograms from FP64 to FP32 or FP16 to reduce RAM:
    %  Also speeds up calculations a little?
    % FP32 Option:
    zHam{i}  = single(zHam{i});
    zB{i}    = single(zB{i});
    zHann{i} = single(zHann{i});
    % FP16 Option:  lead to gaps in heatmap
    %zHam{i}  = half(zHam{i});
    %zB{i}    = half(zB{i});
    %zHann{i} = half(zHann{i});
    
    
    %% Clean up the Spectrograms
    % Crop Z vertically to reduce memory cost
    zHam{i} = zHam{i}(1:UpperInd,:);
    zB{i} = zB{i}(1:UpperInd,:);
    zHann{i} = zHann{i}(1:UpperInd,:);

    % Convert complex spectrogram values to real
    zHam{i} = real(zHam{i}.*conj(zHam{i}));
    zB{i} = real(zB{i}.*conj(zB{i}));
    zHann{i} = real(zHann{i}.*conj(zHann{i}));

    % Normalize spectrogram values locally (for each time-bin)
    zHam{i} = zHam{i}./max(zHam{i});
    zB{i} = zB{i}./max(zB{i});
    zHann{i} = zHann{i}./max(zHann{i}); 
    
    
    %% Time-shift the spectrogram to match signal time
    t{i} = t{i}.*1e9 + ROIxMin + FirstTimePoint;
    
    
    %% Crop the spectrograms according to user-defined ROI bounds
    
    %   Interpolate the ROI bounds for later cropping
    UpCrop = interp1(ROIUPx+FirstTimePoint,ROIUPy,t{i},'linear','extrap');
    LowCrop = interp1(ROILOWx+FirstTimePoint,ROILOWy,t{i},'linear','extrap');

    % Remove any negative values (make zero)
    UpCrop(UpCrop<0) = 0;
    LowCrop(LowCrop<0) = 0;

    % Fix extrapolation to nonrealistic values (values outside input max/min)
    UpCrop(UpCrop>max(f))=max(f);
    LowCrop(LowCrop<min(f))=min(f);

    
    % Determine matrix indices to cut from Z
    %   min function is faster than find function
    [~,UpCropInd] = min(abs(repmat(UpCrop,length(f),1)-repmat(f,1,length(t{i}))));
    [~,LowCropInd] = min(abs(repmat(LowCrop,length(f),1)-repmat(f,1,length(t{i}))));

    %   LOOP n:  neglect all data outside of ROI bounds by setting to 0
    for n = 1:length(t{i})
        zHam{i}(UpCropInd(n):end,n)=0;
        zHam{i}(1:LowCropInd(n),n)=0;
        zB{i}(UpCropInd(n):end,n)=0;
        zB{i}(1:LowCropInd(n),n)=0;
        zHann{i}(UpCropInd(n):end,n)=0;
        zHann{i}(1:LowCropInd(n),n)=0;
        %{
        % Option to make -60 dB instead
        zHam{i}(UpCropInd(n):end,n)=1e-6;
        zHam{i}(1:LowCropInd(n),n)=1e-6;
        zB{i}(UpCropInd(n):end,n)=1e-6;
        zB{i}(1:LowCropInd(n),n)=1e-6;
        zHann{i}(UpCropInd(n):end,n)=1e-6;
        zHann{i}(1:LowCropInd(n),n)=1e-6;
        %}
        %{
        % Option to make NaN instead
        zHam{i}(UpCropInd(n):end,n)=NaN;
        zHam{i}(1:LowCropInd(n),n)=NaN;
        zB{i}(UpCropInd(n):end,n)=NaN;
        zB{i}(1:LowCropInd(n),n)=NaN;
        zHann{i}(UpCropInd(n):end,n)=NaN;
        zHann{i}(1:LowCropInd(n),n)=NaN;
        %}
        
    end
end

% Update waitbar
if RunGauss == 1
    waitbar(0.1,WB,'Cropping to ROI bounds and time-bin edges...');
else
    waitbar(0.8,WB,'Cropping to ROI bounds and time-bin edges...');
end


%% Crop t and z* horizontally, according to 100 ns Duration edges
for i = 1:in
    % left edges, time must be last
    zHam{i}(:,t{i}<min(t{in})) = [];
    zB{i}(:,t{i}<min(t{in})) = [];
    zHann{i}(:,t{i}<min(t{in})) = [];
    t{i}(t{i}<min(t{in})) = [];
    
    % right edges, time must be last
    zHam{i}(:,t{i}>max(t{in})) = [];
    zB{i}(:,t{i}>max(t{in})) = [];
    zHann{i}(:,t{i}>max(t{in})) = [];
    t{i}(t{i}>max(t{in})) = []; % right edges
end

% REDEFINE UPCROP AND LOWCROP (values were temporary in the parfor loop)
%   Interpolate the ROI bounds for later cropping
UpCrop = interp1(ROIUPx+FirstTimePoint,ROIUPy,t{1},'linear','extrap');
LowCrop = interp1(ROILOWx+FirstTimePoint,ROILOWy,t{1},'linear','extrap');

% Remove any negative values (make zero)
UpCrop(UpCrop<0) = 0;
LowCrop(LowCrop<0) = 0;

% Fix extrapolation to nonrealistic values (values outside input max/min)
UpCrop(UpCrop>max(f))=max(f);
LowCrop(LowCrop<min(f))=min(f);

% Slightly shrink cropping edges:
ShrinkFactor = 0.001; % 0.1 percent
UpCrop = UpCrop.*(1-ShrinkFactor);
LowCrop = LowCrop.*(1+ShrinkFactor);

% Mark the time!
SpecTime = toc(SpecTimeStart); % time to calculate spectrograms
% NOTE:  At this point, t is a matrix of a repeating vector, where every 
%   column shows the exact same value.  This is useful for heatmap calcs.

% Update waitbar
if RunGauss == 1
    waitbar(0.12,WB,'Extracting histories: Max Method');
else
    waitbar(0.85,WB,'Extracting histories: Max Method');
end

%% CALCULATE TIME-FREQUENCY HISTORIES

% =========================================================================
% ============================== MAX  METHOD ==============================
% =========================================================================
MaxStart = tic;
% Max CAN handle NaN cells
parfor i = 1:in
    % CONSIDER:  FINDPEAKS.M INSTEAD OF MAX.M
    %   Find indices
    [~,fham_max] = max(zHam{i});
    [~,fb_max] = max(zB{i});
    [~,fhann_max] = max(zHann{i});
    
    % Convert indices to frequencies
    FHAM{1,i} = transpose(f(fham_max));
    FB{1,i} = transpose(f(fb_max));
    FHANN{1,i} = transpose(f(fhann_max));
    
    % Set any value outside ROI boundaries to NaN
    FHAM{1,i}(FHAM{1,i}<=LowCrop)=NaN;
    FB{1,i}(FB{1,i}<=LowCrop)=NaN;
    FHANN{1,i}(FHANN{1,i}<=LowCrop)=NaN;
    % upper extreme removal:
    FHAM{1,i}(FHAM{1,i}>=UpCrop)=NaN;
    FB{1,i}(FB{1,i}>=UpCrop)=NaN;
    FHANN{1,i}(FHANN{1,i}>=UpCrop)=NaN;
end
MaxTime = toc(MaxStart); % time to calculate max histories

% Update waitbar
if RunGauss == 1
    waitbar(0.13,WB,'Extracting histories: Centroid Method');
else
    waitbar(0.95,WB,'Extracting histories: Centroid Method');
end

% =========================================================================
% ======================== ROBUST CENTROID METHOD =========================
% =========================================================================
CentStart = tic;
parfor i = 1:in
    
    % f_centerofmass is equal to {int(f*Z)df / int(Z)df}
    % In other words:  integral of curve*x / area under curve
    fZham = trapz(f,f.*zHam{i});
    fZb = trapz(f,f.*zB{i});
    fZhann = trapz(f,f.*zHann{i});
    Aham = trapz(f,zHam{i});
    Ab = trapz(f,zB{i});
    Ahann = trapz(f,zHann{i});
    % calculate the centroid y-location
    FHAM{2,i} = fZham./Aham;
    FB{2,i} = fZb./Ab;
    FHANN{2,i} = fZhann./Ahann;
    
    % Set any value outside ROI boundaries to NaN
    FHAM{2,i}(FHAM{2,i}<=LowCrop)=NaN;
    FB{2,i}(FB{2,i}<=LowCrop)=NaN;
    FHANN{2,i}(FHANN{2,i}<=LowCrop)=NaN;
    % upper extreme removal:
    FHAM{2,i}(FHAM{2,i}>=UpCrop)=NaN;
    FB{2,i}(FB{2,i}>=UpCrop)=NaN;
    FHANN{2,i}(FHANN{2,i}>=UpCrop)=NaN;
end
CentTime = toc(CentStart); % time to calculate robust centroid histories

% Update waitbar
if RunGauss == 1
    waitbar(0.15,WB,'Extracting histories: Gaussian Method');
else
    waitbar(0.99,WB,'Finishing...')
end

% =========================================================================
% ============================ GAUSSIAN METHOD ============================
% =========================================================================
% NOTE: This fitting procedure takes 90-99% of total calculation time
%error('run stopped')
% Testing out faster Gaussian fitting from Dan Champion's recommendation:
%  1. locally normalize (get rid of A variable)
%  2. Take log of the spectrogram, and fit with quadratic instead of Gauss
%     -- original fit:     f(x)  = a*exp(-((x-b)/c)^2)
%     -- new fit:      log(f(x)) = -((x-b)/c)^2
%        --- fit with basic 'poly2' function for increased speed


if RunGauss ==1
    GaussStart = tic;
    for i = 1:in
        % Update the waitbar:
        waitbar(0.15+0.84*(i/in),WB,'Extracting histories: Gaussian Method');
        
        % preallocate outputs:
        gfham = zeros(1,length(t{i})); % gf means Gaussian Fit
        gfb = zeros(1,length(t{i}));
        gfhann = zeros(1,length(t{i}));
        
        % Begin the "Gaussian Extraction"
        parfor m = 1:length(t{i})
            % pre-process data before fitting:
            % pull out y-data and normalize
            yHam  = double( zHam{i}(:,m)./max(zHam{i}(:,m)));  % Hamming spectrum
            yB    = double(   zB{i}(:,m)./max(zB{i}(:,m)));     % Blackman spectrum
            yHann = double(zHann{i}(:,m)./max(zHann{i}(:,m))); % Hann spectrum
            % crop data to only spectral densities above 2% of max, "Gthresh"
            Gthresh = 0.02;
            xHam  = double(f(yHam>Gthresh));
            xB    = double(f(yB>Gthresh));
            xHann = double(f(yHann>Gthresh));
            yHam  = yHam(yHam>Gthresh);
            yB    = yB(yB>Gthresh);
            yHann = yHann(yHann>Gthresh);
            % Take the natural log of the Spectral density:
            yHam  = log(yHam);
            yB    = log(yB);
            yHann = log(yHann);
            
            % Try to do the Gaussian (quadratic) fitting:
            try
                gHam = fit(xHam,yHam,'poly2');
                gfham(m) = -(gHam.p2/(2*gHam.p1)); % center = -b/2a
                
                gB = fit(xB,yB,'poly2');
                gfb(m) = -(gB.p2/(2*gB.p1)); % center = -b/2a
                
                gHann = fit(xHann,yHann,'poly2');
                gfhann(m) = -(gHann.p2/(2*gHann.p1)); % center = -b/2a
                
            catch % if fit fails, then output NaN (neglected in Heatmap calc)
                gfham(m) = NaN;
                gfb(m) = NaN;
                gfhann(m) = NaN;
            end
        end
        % remove values outside of cropped boundaries
        % lower extreme removal:
        gfham(gfham<LowCrop)=NaN;
        gfb(gfb<LowCrop)=NaN;
        gfhann(gfhann<LowCrop)=NaN;
        % upper extreme removal:
        gfham(gfham>UpCrop)=NaN;
        gfb(gfb>UpCrop)=NaN;
        gfhann(gfhann>UpCrop)=NaN;
        
        % Output history at specified time-resolution
        % Interpolate using a cubic spline, do not include NaN cells
        FHAM{3,i} = gfham;
        FB{3,i} = gfb;
        FHANN{3,i} = gfhann;
    end
    GaussTime = toc(GaussStart); % time to calculate gaussian histories
else
end


% End the timer and tell the user how long extraction took
ExtractEnd = toc(ExtractStart);

% Delete the progress bar
delete(WB)

% Tell the user how long each process took
disp(' ');
disp('==============================================================');
disp('================ Extraction Process Finished! ================');
disp('==============================================================');
disp(['Total elapsed time: ' num2str(ExtractEnd,'%.2f') ' seconds']);
disp(['  - Extraction efficiency:  ' num2str(ExtractEnd/range(t{1}./1e3),'%.2f') ' seconds per microsecond of data'])
disp(['  - Spectrogram production time: ' num2str(SpecTime,'%.2f') ' seconds']);
disp(['  - Maximum extraction method time: ' num2str(MaxTime,'%.2f') ' seconds']);
disp(['  - Robust Centroid extraction method time: ' num2str(CentTime,'%.2f') ' seconds']);
if RunGauss ==1
    disp(['  - Gaussian Peak extraction method time: ' num2str(GaussTime,'%.2f') ' seconds']);
%    disp(['    ^ I know this is long.  Working on it.']);
else
end
disp(' ');


disp('Calculating final history and uncertainty bounds...');







%% Process the output history data

% Collect histories by window function
FHAM_all = double(vertcat(FHAM{:}));
FB_all = double(vertcat(FB{:}));
FHANN_all = double(vertcat(FHANN{:}));
% Collect histories by extraction method
Fmax = double(vertcat(FHAM{1,:},FB{1,:},FHANN{1,:}));
Fcentroid = double((vertcat(FHAM{2,:},FB{2,:},FHANN{2,:})));
if RunGauss ==1
    Fgaussian = double(vertcat(FHAM{3,:},FB{3,:},FHANN{3,:}));
else
end
% Collect all histories
F = vertcat(FHAM_all,FB_all,FHANN_all);
%F(isnan(F)) = 0; % not necessary because histcounts2 neglects NaN cells



%% Prepare for bi-variate histogram (heatmap) calculation

% Set a time point for each frequency point; Tmat will correspond with F
Tmat = repmat(vertcat(t{:}),size(F,1)/in,1);

% Define time bins for the heatmap at same time-resolution of calculation
%   To do this, set first bin 1/2 timestep before first data point
%   and set last bin 1/2 timestep after last data point
halftimestep = (t{1}(2)-t{1}(1))/2;
tedges = [t{1}-halftimestep t{1}(end)+halftimestep];
% By doing this, t{1} will reflect the heatmap time-resolution
tfinal = t{1};

% Do the same treatment for frequency (vertical bins)
onefreqstep = f(2)-f(1);
halffreqstep = (f(2)-f(1))/2;
fedges = [f'-halffreqstep f(end)+halffreqstep]; % apostrophe transposes

% Check number of frequency bins, create an addition of zeros if needed
% Bivariate heatmap function, hiscounts2, has max bins of 1024
if length(fedges)>1024
    InitSize = length(fedges);
    fedges = fedges(fedges>min(LowCrop));
    %if length(fedges)>1024
        %warning('Number of frequency bins greater than 1024')
        %warning('Heatmap may have horizontal black lines')
    %end
    Removed = InitSize-length(fedges);
    HeatMapAddition = zeros(Removed,length(tedges)-1);
end
    
%fedges = min(LowCrop):onefreqstep:max(UpCrop);
% By doing this, f will reflect the heatmap frequency-resolution


%% Calculate the bi-variate histogram (heatmap)
%   "edge" refers to the time or frequency position of bin edges
%   Because of this implementation, edges could be applied anisotropically
HeatMap = histcounts2(Tmat(:),F(:),tedges,fedges); % f gives about a 1 m/s vertical bin length
HeatMap = transpose(HeatMap); % histcounts2 outputs t by f, we want f by t
if exist('HeatMapAddition')
    HeatMap = [HeatMapAddition; HeatMap];
else
end

% Locally normalize heatmap in each time bin
parfor i = 1:(length(tedges)-1) % num_bins = num_edges-1
    HeatMap(:,i) = HeatMap(:,i)/max(HeatMap(:,i));
end


%% Calculate final histories
% Use the spectrogram extraction methods on the heatmap

% Max method:  Highest occurence
[~,HistoryMaxF] = max(HeatMap);
HistoryMaxF = f(HistoryMaxF);
HistoryMaxF(HistoryMaxF<0.025)=0;% remove any values below ~20 m/s

% Robust Centroid method: weighted average
HistoryCentroidF = transpose(trapz(f,f.*HeatMap)./trapz(f,HeatMap));
HistoryCentroidF(HistoryCentroidF<0.025)=0; % remove any values below ~20 m/s

% Gaussian fit method:
HistoryGaussianF = zeros(1,length(tfinal));
parfor i = 1:length(tfinal)
    yHM = HeatMap(:,i)./max(HeatMap(:,i));
    xHM = f(yHM>0.02);
    yHM = yHM(yHM>0.02);
    yHM = log(yHM);
    try
        gHM = fit(xHM,yHM,'poly2')
        HistoryGaussianF(i) = -(gHM.p2/(2*gHM.p1)); % center = -b/2a
    catch
        HistoryGaussianF(i) = NaN;
    end
end
% Remove values that are outside ROI (bad/wrong fits)
HistoryGaussianF(HistoryGaussianF>UpCrop)=NaN;
HistoryGaussianF(HistoryGaussianF<LowCrop)=NaN;

%HistoryGaussianF(HistoryGaussianF<0.025)=0; % remove any values below ~20 m/s

% Find the mean of these extractions, and use that value for the final
HistoryMeanF = zeros(1,length(tfinal));
parfor i = 1:length(tfinal)
    HistoryMeanF(i) = nanmean([HistoryMaxF(i) HistoryCentroidF(i) ...
        HistoryGaussianF(i)]);
end



%% Calculate 95% Confidence Interval of Systematic Uncertainty
% Find area and cumulative area of each heatmap column
AHM = trapz(f,HeatMap); % area of each heatmap column
AHMmat = cumtrapz(f,HeatMap); % cumulative area of each heatmap column
AHMmat = AHMmat./max(AHMmat); % normalize so that max area is 1.0


% Find the 95% confidence interval lower and upper limit
% To do this, find where cumulative area is 0.025 and 0.975
% interpolate using a spline
LowUncert = zeros(size(AHM));
UpUncert = LowUncert;
parfor k = 1:size(AHMmat,2)
    LowUncert(k) = mean(mean(fnzeros(spline(f,AHMmat(:,k)-0.025))));
    UpUncert(k) = mean(mean(fnzeros(spline(f,AHMmat(:,k)-0.975))));
end

% Apply a rolling average to the uncertainty, 10 points wide (25 ns)
LowUncert = movmean(LowUncert,10);
UpUncert = movmean(UpUncert,10);

% Remove Low Uncert that is higher than AvgV
%  B/C lower bound cannot be above answer
LowUncert(LowUncert>HistoryMeanF)=HistoryMeanF((LowUncert>HistoryMeanF));

% Remove High Uncert that is lower than AvgV
%  B/C upper bound cannot be below answer
UpUncert(UpUncert<HistoryMeanF)=HistoryMeanF((UpUncert<HistoryMeanF));

% set output time axes to microseconds
tfinal = tfinal./1e3;
tedges = tedges./1e3;


%% Convert frequency to velocity for output figures
Vfactor = Wavelength/2; % velocity is frequency * wavelength/2

v = f.*Vfactor; % vert axis for spectrogram and heatmap

VHAM = FHAM_all.*Vfactor;
VB = FB_all.*Vfactor;
VHANN = FHANN_all.*Vfactor;

Vmax = Fmax.*Vfactor;
Vcentroid = Fcentroid.*Vfactor;
if RunGauss ==1
    Vgaussian = Fgaussian.*Vfactor;
else
end

FinalHistory = HistoryMeanF.*Vfactor; % m/s, final output history
UpUncertV = UpUncert.*Vfactor; % m/s
LowUncertV = LowUncert.*Vfactor; % m/s

UpCrop = interp1((ROIUPx+FirstTimePoint)./1e3,ROIUPy,tfinal,'linear','extrap');
LowCrop = interp1((ROILOWx+FirstTimePoint)./1e3,ROILOWy,tfinal,'linear','extrap');
UpCropV = UpCrop.*Vfactor; % m/s
LowCropV = LowCrop.*Vfactor; % m/s


%% Apply a window correction factor
% format is Vreal = WF.a*Vapparent^WF.b

v = WF.a.*v.^WF.b;
VHAM = WF.a.*VHAM.^WF.b;
VB = WF.a.*VB.^WF.b;
VHANN = WF.a.*VHANN.^WF.b;
Vmax = WF.a.*Vmax.^WF.b;
Vcentroid = WF.a.*Vcentroid.^WF.b;
if RunGauss==1
    Vgaussian = WF.a.*Vgaussian.^WF.b;
else
end
FinalHistory = WF.a.*FinalHistory.^WF.b;
UpUncertV = WF.a.*UpUncertV.^WF.b;
LowUncertV = WF.a.*LowUncertV.^WF.b;
UpCropV = WF.a.*UpCropV.^WF.b;
LowCropV = WF.a.*LowCropV.^WF.b;

%% Apply the frequency shifting amount (up/down) if there is one:
%  Calculated in previous steps
%  Now we just need to subtract the baseline/reference (shift amount)

if exist('SHIFT_SIGN')
    v     = v     - SHIFT_AMT;
    VHAM  = VHAM  - SHIFT_AMT;
    VB    = VB    - SHIFT_AMT;
    VHANN = VHANN - SHIFT_AMT;
    Vmax  = Vmax  - SHIFT_AMT;
    Vcentroid = Vcentroid - SHIFT_AMT;
    if RunGauss==1
        Vgaussian = Vgaussian - SHIFT_AMT;
    else
    end
    FinalHistory = FinalHistory - SHIFT_AMT;
    UpUncertV  = UpUncertV  - SHIFT_AMT;
    LowUncertV = LowUncertV - SHIFT_AMT;
    UpCropV    = UpCropV    - SHIFT_AMT;
    LowCropV   = LowCropV   - SHIFT_AMT;
else
end

% Smoothed history
SmoothHistory=movmean(FinalHistory,10);

%% Detect Rise and Calculate the Rise Times (10-50-90)
% The function "CalculateRiseTimes" was created to auto-detect the steepest
% rise in velocity in the velocity profile, then determine the peak
% velocity, and 10-50-90 rise times of the shock front.
%   See:  .../RequiredFunctions/CalculateRiseTimes.m for more details

switch lower(C{6})
    case {'yes','y','true','yes please','please'}
        % Raw history:
        %[RiseT(1),RiseT(2),RiseT(3),PeakVelocity] = CalculateRiseTimes(tfinal,FinalHistory);
        % Smoothed history:
        [RiseT(1),RiseT(2),RiseT(3),PeakVelocity] = CalculateRiseTimes(tfinal,SmoothHistory);
        RiseV = PeakVelocity.*[.1 .5 .9];
    otherwise
        disp('User selected to NOT calculate 10-50-90 rise times');
end




%% CALCULATIONS FINISHED!  Now, plot the output for user feedback

disp('Creating plots for your review...');

% Plot output preparation:
FS = 10; % font size
IS = [0 0 12 6.75]; % image size, 16x9 ratio

% Axes limits:
LeftLim = min(tfinal);
RightLim = max(tfinal);
LowerLim = min(LowCropV); % cut lower frequencies for upshifted signals
UpperLim = max(UpCropV);
axlim = [LeftLim RightLim LowerLim UpperLim];

% =========================================================================
% Plot the example spectrogram with the final history & uncertainties
% =========================================================================
FinalOut = figure;
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
FinalSpec = nexttile;
%FinalSpec = subplot(2,3,1);
hold on
imagesc(texample,v,zexample) % Higher resolution version of example spectrogram
set(gca,'YDir','normal')
colormap(FinalSpec,SPECTROGRAM_COLORMAP)
caxis([-60 0]);
plot(tfinal,SmoothHistory,'color','k','linewidth',2); % smoothed version 25 ns
plot(tfinal,FinalHistory,'color','k','linewidth',2,'LineStyle','--');
%plot(tfinal,LowUncertV,'color','k','LineStyle',':')
%plot(tfinal,UpUncertV,'color','k','LineStyle',':')
plot(tfinal,UpCropV,'color','w','LineWidth',1','LineStyle','-')
plot(tfinal,LowCropV,'color','w','LineWidth',1','LineStyle','-')
hold off
box on
axis(axlim)
%legend({'Max','Centroid','95% Conf. Int.'},'Location','southeast','FontSize',FS);
title('Example Spectrogram','FontSize',FS,'fontweight','normal');
xlabel('Detector Time (\mus)','FontSize',FS);
ylabel('Velocity (m/s)','FontSize',FS);
set(gca,'fontsize',FS)
%set(gcf,'PaperUnits','inches');
%set(gcf,'PaperPosition',IS); % 3.5" wide, 3.0" tall



% =========================================================================
% Plot the final history with uncertainty bounds
% =========================================================================
%tArea = [tfinal tfinal(end) tfinal tfinal(end)];
%UncertArea = [LowUncertV 0 UpUncertV 0];

%HistoryFig = subplot(2,3,[2,3]);
HistoryFig = nexttile([1 2]);
hold on
plot(tfinal,FinalHistory,'color',[179 163 105]./255,'LineWidth',1.25) % GT Gold
LEG{1} = 'History, Raw';
plot(tfinal,movmean(FinalHistory,10),'color',[0 48 87]./255,'LineStyle','-','LineWidth',1.25); % smoothed data (50 ns rolling average), GT Blue
LEG{2} = 'History, Smoothed';
if exist('SHIFT_AMT','var')==1
    plot([tfinal(1) tfinal(end)],[0 0],'color',[0 0 0 0.25],'LineWidth',1.25,'linestyle',':')
    LEG{3} = 'Baseline';
end
if exist('RiseT','var')==1
    scatter(RiseT,RiseV,10,'MarkerEdgeColor','k','MarkerFaceColor','k');
    LEG{4} = '10-50-90 Rise Times';
else
end
plot(tfinal,UpUncertV,'color',[0.5 0.5 0.5],'linestyle','-.','LineWidth',1); % 95% Conf bound up
plot(tfinal,LowUncertV,'color',[0.5 0.5 0.5],'linestyle','-.','LineWidth',1); % 95% Conf bound low
LEG{5} = 'Systematic Uncertainty (95%)';
LEG = LEG(~cellfun('isempty',LEG)); % remove empty legend cells

axis(axlim)
hold off
box on
title('Final Output History and Uncertainties','FontSize',FS,'fontweight','normal');
%legend(LEG,'Location','southeast','FontSize',FS);
legend(LEG,'Location','best','FontSize',FS);
xlabel('Detector Time (\mus)','FontSize',FS);
%ylabel('Velocity (m/s)','FontSize',FS);
set(gca,'fontsize',FS)



% =========================================================================
% Plot all histories on the same plot
% =========================================================================

% Define differential colors - red, blue, black
%   Extraction Method Colors
C11 = [1 0 0]; % RED
C12 = [0 0 1]; % BLUE
C13 = [0 0 0]; % BLACK
%   Window Function Colors
C21 = [0 0.8 0]; % GREEN
C22 = [0 0 1]; % BLUE
C23 = [0 0 0]; % BLACK



% -----
% Color by window function, or show analysis of the frequency shift
% -----
if exist('SHIFT_AMT','var')
    FinalBase = nexttile;
    hold on
    scatter(BASEfft.x.*Wavelength/2,BASEfft.ham,10,'o','markeredgecolor','k',...
        'markerfacecolor',C21) % hamming
    scatter(BASEfft.x.*Wavelength/2,BASEfft.blk,10,'o','markeredgecolor','k',...
        'markerfacecolor',C22) % blackman
    scatter(BASEfft.x.*Wavelength/2,BASEfft.han,10,'o','markeredgecolor','k',...
    'markerfacecolor',C23) % hann
    %plot(Reference.Xv,Reference.Y,'color','k','linewidth',1.25) % removed March 23, 2023
    plot(BaseFit.Xv,BaseFit.ham,'color',C21)
    plot(BaseFit.Xv,BaseFit.blk,'color',C22)
    plot(BaseFit.Xv,BaseFit.han,'color',C23)
    
    hold off
    box on
    %legend({'Hamming','Blackman','Hann'},'Location','southeast','FontSize',FS)
    legend({'Hamming','Blackman','Hann'},'Location','best','FontSize',FS)
    xlim([SHIFT_AMT-SHIFT_UNC*60 SHIFT_AMT+SHIFT_UNC*90]) % slightly left, so legend doesn't obscure curve
    ylim([0 1.1])
    C{8}(1) = upper(C{8}(1)); % capitalize the first letter of shift direction
    C{8}(2:end) = lower(C{8}(2:end)); % lowercase for other letters in shift direction
    title([C{8} 'shifted Baseline: ' num2str(SHIFT_AMT,'%.2f') '± ' num2str(SHIFT_UNC,'%.2f') ' m/s'],'FontSize',FS,'fontweight','normal');
    xlabel('Baseline Velocity (m/s)','FontSize',FS);
    ylabel('Spectral Density, arb','FontSize',FS);
    set(gca,'fontsize',FS)
else
    %FinalWindow = subplot(2,3,4);
    FinalWindow = nexttile;
    hold on
    % Legend Entries
    plot([0 0], [0 0],'color',C21,'linestyle','-'); % Placeholder for legend entry
    plot([0 0], [0 0],'color',C22,'linestyle','-'); % Placeholder for legend entry
    plot([0 0], [0 0],'color',C23,'linestyle','-'); % Placeholder for legend entry
    % Hamming Window - RED
    plot(tfinal,VHAM,'color',[C21 0.01]);
    % Blackman Window - BLUE
    plot(tfinal,VB,'color',[C22 0.01]);
    % Hann Window - BLACK
    plot(tfinal,VHANN,'color',[C23 0.01]);
    hold off
    box on
    %legend({'Hamming','Blackman','Hann'},'Location','southeast','FontSize',FS)
    legend({'Hamming','Blackman','Hann'},'Location','best','FontSize',FS)
    axis(axlim);
    title('All Histories: Window Function','FontSize',FS,'fontweight','normal');
    xlabel('Detector Time (\mus)','FontSize',FS);
    ylabel('Velocity (m/s)','FontSize',FS);
    set(gca,'fontsize',FS)
    %set(gcf,'PaperUnits','inches');
    %set(gcf,'PaperPosition',IS); % 3.5" wide, 3.0" tall
end



% -----
% Color by extraction method
% -----
%FinalExtraction = subplot(2,3,5);
FinalExtraction = nexttile;
hold on
% Legend entries
plot([0 0], [0 0],'color',C11,'linestyle','-'); % Placeholder for legend entry
plot([0 0], [0 0],'color',C12,'linestyle','-'); % Placeholder for legend entry
plot([0 0], [0 0],'color',C13,'linestyle','-'); % Placeholder for legend entry
% max method - RED
plot(tfinal,Vmax,'color',[C11 0.01]);
% centroid method - BLUE
plot(tfinal,Vcentroid,'color',[C12 0.01]);
% gaussian method - BLACK
if RunGauss ==1
    plot(tfinal,Vgaussian,'color',[C13 0.01]);
else
end
hold off
box on
if RunGauss ==1
    %legend({'Maximum','Centroid','Gaussian'},'Location','southeast','FontSize',FS)
    legend({'Maximum','Centroid','Gaussian'},'Location','best','FontSize',FS)
else
    %legend({'Maximum','Centroid'},'Location','southeast','FontSize',FS)
    legend({'Maximum','Centroid'},'Location','best','FontSize',FS)
end
axis(axlim);
title('All Histories: Extraction Method','FontSize',FS,'fontweight','normal');
xlabel('Detector Time (\mus)','FontSize',FS);
ylabel('Velocity (m/s)','FontSize',FS);
set(gca,'fontsize',FS)
%set(gcf,'PaperUnits','inches');
%set(gcf,'PaperPosition',IS); % 3.5" wide, 3.0" tall



% =========================================================================
% Plot a heatmap of all of the output histories
% =========================================================================
%FinalHM = subplot(2,3,6);
FinalHM = nexttile;
hold on
imagesc(tedges,v,HeatMap)
plot(tfinal,FinalHistory,'color',[1 1 1 0.25]) % 75% transparent white line
set(gca,'YDir','normal')
colormap(FinalHM,HEATMAP_COLORMAP)
caxis([0 1]);
hold off
box on
axis(axlim);
title('Heatmap of All Histories','FontSize',FS,'fontweight','normal');
xlabel('Detector Time (\mus)','FontSize',FS);
%ylabel('Velocity (m/s)','FontSize',FS);
set(gca,'fontsize',FS)
%set(gcf,'PaperUnits','inches');
%set(gcf,'PaperPosition',IS); % 3.5" wide, 3.0" tall



set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',IS); % 16:9 ratio

sgtitle(['PDV Summary for: ' File(1:(end-3))],'FontWeight','bold','FontSize',11)




%% SAVE OUTPUTS
% Ask the user if they want to save their outputs
% If yes, save all outputs

SavePrompt = questdlg('Would you like to save your history, spectrogram, heatmap, and plots?',...
    'Save Data?','Yes','No','Cancel','Yes');
switch SavePrompt
    case 'No'
        warning('No data saved.  Data output to workspace in case this was a mistake.')
    case 'Cancel'
        close all; clear; clc
        error('No data saved.  All varaibles have been cleared.');
    case 'Yes'
        % make the figure inivisible while saving is occuring
        %FinalOut.Visible = 'off';
        SaveDir = uigetdir(Path,'SELECT A DIRECTORY TO SAVE PDV OUTPUTS TO:');
        if SaveDir==0 % user did not select anything
            SaveDir = Path; % Save output next to original file
        else
            SaveDir = [SaveDir '/']; % append a / for smaller save name
        end
        SaveNames = {'Shot ID:','Probe ID','Additional Description (optional):'};
        SNtitle = 'Save Data Options';
        SNdims = [1 50];
        SNdef = {'XXXX','X',''};
        SN = inputdlg(SaveNames,SNtitle,SNdims,SNdef);
        
        % Define Final output filenames
        if strcmp(SN{3},'')
            COMMON_NAME = ['Shot' SN{1} '_Probe' SN{2}];
        else
            COMMON_NAME = ['Shot' SN{1} '_Probe' SN{2} '_' SN{3}];
        end
        HistName  = [COMMON_NAME '_History.csv'];
        SpecName  = [COMMON_NAME '_Spectrogram.mat'];
        HMName    = [COMMON_NAME '_Heatmap.mat'];
        PlotName  = [COMMON_NAME '_Plots.png'];
        ShiftName = [COMMON_NAME '_FrequencyShiftAmount.csv'];
        RiseName  = [COMMON_NAME '_RiseTimes.csv'];
        
        % Print to workspace the filenames that will be saved
        disp(['All outputs saved to directory: ' SaveDir]);
        disp(['History saved as "'     HistName '"']);
        disp(['Spectrogram saved as "' SpecName '"']);
        disp(['Heatmap saved as "'     HMName '"']);
        disp(['Final plots saved as "' PlotName '"']);
        if exist('RiseT','var')
            disp(['Rise Times (10-50-90) saved as "' RiseName '"']);
        else
        end
        if exist('SHIFT_AMT','var')
            disp(['Frequency Shift Info saved as "' ShiftName '"']);
        else
        end
        disp('')
        disp('Saving your data now...');
        
        SaveStart = tic;
        
        % Write the History file with headers
        fid = fopen([SaveDir HistName],'w');
        fprintf(fid,'%s\r\n','Time(\mus),HistoryRaw(m/s),HistorySmoothed(m/s),SysUncUpBound(m/s),SysUncLowBound(m/s),ROIup(m/s),ROIlow(m/s)');
        fclose(fid);
        dlmwrite([SaveDir HistName],[tfinal',FinalHistory',movmean(FinalHistory',10),UpUncertV',LowUncertV',UpCropV',LowCropV'],'-append','delimiter',',','precision',8)
        
        % Save the Spectrogram to a MATLAB binary file
        SpectrogramT = texample;
        SpectrogramV = v;
        Spectrogram = zexample;
        save([SaveDir SpecName],'SpectrogramT','SpectrogramV','Spectrogram');
        
        % Save the Heatmap to a MATLAB binary file
        HeatMapT = tedges;
        HeatMapV = v;
        save([SaveDir HMName],'HeatMapT','HeatMapV','HeatMap');
        
        % Change title of the final figure to say shot ID and probe number
        %figure(FinalOut)
        if strcmp(SN{3},'') % if there is an additional description, add it to title
            sgtitle(['Shot ' SN{1} ', Probe ' SN{2}],'FontWeight','bold','FontSize',11);
        else
            sgtitle(['Shot ' SN{1} ', Probe ' SN{2} ': ' SN{3}],'FontWeight','bold','FontSize',11);
        end
            
        % Save the summary figure
        %  Remind the figure of its printing parameters, for good measure
        set(gca,'fontsize',FS)
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperPosition',IS); % 16:9 ratio
        print(FinalOut,[SaveDir PlotName],'-dpng','-r300'); % 300 dpi, high resolution
        SaveEnd = toc(SaveStart);
        disp(['Data saved, elapsed write time: ' num2str(SaveEnd) ' seconds'])
        
        % Write the Rise Time file (if RiseT exists)
         if exist('RiseT','var')
            fid = fopen([SaveDir RiseName],'w');
            fprintf(fid,'%s\r\n','RiseTime10(\mus),RiseTime50(\mus),RiseTime90(\mus),PeakVelocity(m/s)');
            fclose(fid);
            dlmwrite([SaveDir RiseName],[RiseT PeakVelocity],'-append','delimiter',',','precision',8)
         else
         end
         
         % Write out the baseline velocities (if SHIFT_AMT exists)
         if exist('SHIFT_AMT','var')
             fid = fopen([SaveDir ShiftName],'w');
             fprintf(fid,'%s\r\n','FrequencyShift(GHz),VelocityShift(ms/)');
             dlmwrite([SaveDir ShiftName],[mean(BASEfft.ALL) SHIFT_AMT],'-append','delimiter',',','precision',8)
         else
         end
         
end

disp(' ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('==============================================================');
disp('=================== PDV ANALYSIS COMPLETE! ===================');
disp('==============================================================');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%FinalOut.Visible = 'on';


% Write an output file with code runtimes and date
switch SavePrompt
    case 'Yes'
        StdoutName = [COMMON_NAME '_stdout.txt'];
        fid2 = fopen([SaveDir StdoutName],'w');
            fprintf(fid2,'%s\r\n',['Code executed on: ' date]);
            fprintf(fid2,'%s\r\n',['User selected file: ' Path File]);
            fprintf(fid2,'%s\r\n',['Number of data points per 1 ns is: ' num2str(NumDataPointsFor1ns)]);
            fprintf(fid2,'%s\r\n','Point-to-point time resolution set to: 2.5000 ns');
            fprintf(fid2,'%s\r\n',['Total extraction time: ' num2str(ExtractEnd,'%.2f') ' seconds']);
            fprintf(fid2,'%s\r\n',['  - Spectrogram production time: ' num2str(SpecTime,'%.2f') ' seconds']);
            fprintf(fid2,'%s\r\n',['  - Maximum extraction method time: ' num2str(MaxTime,'%.2f') ' seconds']);
            fprintf(fid2,'%s\r\n',['  - Robust Centroid extraction method time: ' num2str(CentTime,'%.2f') ' seconds']);
            if RunGauss ==1
                fprintf(fid2,'%s\r\n',['  - Gaussian Peak extraction method time: ' num2str(GaussTime,'%.2f') ' seconds']);
            else
                fprintf(fid2,'%s\r\n',['  - Gaussian Peak extraction method: NOT EXECUTED']);
            end
            fprintf(fid2,'%s\r\n',['All outputs saved to directory: ' SaveDir]);
            fprintf(fid2,'%s\r\n',['  - History saved as "' HistName '"']);
            fprintf(fid2,'%s\r\n',['  - Spectrogram saved as "' SpecName '"']);
            fprintf(fid2,'%s\r\n',['  - Heatmap saved as "' HMName '"']);
            fprintf(fid2,'%s\r\n',['  - Final plots saved as "' PlotName '"']);
            if exist('RiseT','var')
                fprintf(fid2,'%s\r\n',['Rise Times (10-50-90) saved as "' RiseName '"']);
            else
                fprintf(fid2,'%s\r\n',['Rise Times NOT CALCULATED']);
            end
            fprintf(fid2,'%s\r\n',['Saved data write time: ' num2str(SaveEnd) ' seconds']);
        fclose(fid2);
end


































