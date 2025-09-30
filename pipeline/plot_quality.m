function plot_quality(varargin)
if isempty(varargin)
    rootDir = '/project/TIBIR/Lega_lab/shared/lega_ansir/Single_Units_BR';
else
    rootDir = varargin{1};
end

%List directories and load list of recordings
dataDir = fullfile(rootDir,'micros_database');
qualityDir = fullfile(rootDir,'micros_pipeline/quality_metrics/plots');
sorted_file = fullfile(dataDir,'sorted_micros.mat');
load(sorted_file,'sorted_micros')

%%%%%%%%%%%%%%%%%%%%%%%% Spike Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Separate table to collect all available spike data for the recordings that
%were successfully clustered and timed.
timed_micros = sorted_micros(sorted_micros.spikes_timed,:);
[timed_micros,channels,neurons,classes] = process_micros_data(dataDir,timed_micros);


%%%%%%%%%%%%%%%%%% Neurons Quality Metrics Analyses %%%%%%%%%%%%%%%%%%%%%%%
    
%%% List of all variables to be statistically compared
neuron_vars = [{'n_spikes';'f_spikes';'n_groups';'n_classes';'stdNoise'},...
    {'p_spc';'p_1st';'p_2nd';'m_distances';'std_distances'},...
    {'mean_FR';'m_amp_uV';'std_amp_uV';'m_voltage';'std_voltage'},...
    {'mean_peak_uV';'peak_width';'p_sub3ms';'p_60Hz';'isi_SNR'},...
    {'dispersion';'m_mahal';'std_mahal';'isoD';'L_ratio'},...
    {'dunn_idx';'silhouette';'SNR_split';'SNR_reref';'SNR_clean'}];
%Corresponding plot titles
neuron_titles = [{'# of spikes in group';'fraction from all detections';'# of groups in channel';'# of classes in group';'STD of signal'},...
    {'proportion matched at SPC';'proportion matched at 1st template';'proportion matched at 2nd template';'mean distances';'STD distances'},...
    {'mean firing rate';'mean amplitude (uV)';'STD amplitude(uV)';'mean voltage(uV)';'STD voltage(uV)'},...
    {'mean waveform amplitude (uV)';'peak width (samples)';'proportion with ISI < 3ms';'proportion with ISIs of line noise';'ISI SNR'},...
    {'dispersion';'mean mahalanobis distance';'STD mahalanobis distance';'isolation distance';'L-ratio'},...
    {'Dunn index';'silhouette score';'SNR of raw signal';'SNR of rerefed signal';'SNR of denoised signal'}];

su_vs_noise_analysis(qualityDir,'neurons','all_neurons',neuron_vars,neuron_titles,neurons);
su_vs_noise_analysis(qualityDir,'neurons','digital',neuron_vars,neuron_titles,neurons(neurons.digital_headstage,:));
su_vs_noise_analysis(qualityDir,'neurons','analog',neuron_vars,neuron_titles,neurons(~neurons.digital_headstage,:));
digital_vs_analog_analysis(qualityDir,'neurons','all_neurons',neuron_vars,neuron_titles,neurons);
digital_vs_analog_analysis(qualityDir,'neurons','single_units',neuron_vars,neuron_titles,neurons(neurons.is_SUA,:));
digital_vs_analog_analysis(qualityDir,'neurons','multi_units',neuron_vars,neuron_titles,neurons(neurons.is_MUA,:));
digital_vs_analog_analysis(qualityDir,'neurons','noise_units',neuron_vars,neuron_titles,neurons(neurons.is_noise,:));


%%%%%%%%%%%%%%%%%% Classes Quality Metrics Analysis %%%%%%%%%%%%%%%%%%%%%%%
    
%%% List of all variables to be statistically compared
class_vars = [{'n_spikes';'f_spikes';'n_groups';'n_classes';'stdNoise'},...
    {'p_spc';'p_1st';'p_2nd';'m_distances';'std_distances'},...
    {'L_ratio';'m_amp_uV';'std_amp_uV';'m_voltage';'std_voltage'},...
    {'m_dispersion';'m_mahal';'std_mahal';'intragroup_isoD';'extragroup_isoD'},...
    {'p_sub3ms';'p_60Hz';'SNR_split';'SNR_reref';'SNR_clean'}];
%Corresponding plot titles
class_titles = [{'# of spikes in class';'fraction from all detections';'# of groups in channel';'# of classes in group';'STD of signal'},...
    {'proportion matched at SPC';'proportion matched at 1st template';'proportion matched at 2nd template';'mean distances';'STD distances'},...
    {'L-ratio';'mean amplitude (uV)';'STD amplitude(uV)';'mean voltage(uV)';'STD voltage(uV)'},...
    {'dispersion';'mean mahalanobis distance';'std mahalanobis distance';'intragroup isolation distance';'extragroup isolation distance'},...
    {'proportion with ISI < 3ms';'proportion with ISIs of line noise';'SNR of raw signal';'SNR of rerefed signal';'SNR of denoised signal'}];

su_vs_noise_analysis(qualityDir,'classes','all_classes',class_vars,class_titles,classes);
su_vs_noise_analysis(qualityDir,'classes','digital',class_vars,class_titles,classes(classes.digital_headstage,:));
su_vs_noise_analysis(qualityDir,'classes','analog',class_vars,class_titles,classes(~classes.digital_headstage,:));
digital_vs_analog_analysis(qualityDir,'classes','all_classes',class_vars,class_titles,classes);
digital_vs_analog_analysis(qualityDir,'classes','single_units',class_vars,class_titles,classes(classes.is_SUA,:));
digital_vs_analog_analysis(qualityDir,'classes','multi_units',class_vars,class_titles,classes(classes.is_MUA,:));
digital_vs_analog_analysis(qualityDir,'classes','noise_units',class_vars,class_titles,classes(classes.is_noise,:));

%%%%%%%%%%%%%%%%%% Channels Quality Metrics Analysis %%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%% Electrodes Quality Metrics Analysis %%%%%%%%%%%%%%%%%%%%%
figure('Units', 'pixels')
fig = gcf;

subplot(7,1,1)
hold on
for idx = 2:height(sorted_micros)
    SNR_split = sorted_micros(idx,:).SNR_split{:};
    if ~isempty(SNR_split)
        plot(repmat(idx,length(SNR_split),1)',SNR_split,'.')
    end
end
hold off
title('SNR split')

subplot(7,1,2)
hold on
for idx = 2:height(sorted_micros)
    SNR_reref = sorted_micros(idx,:).SNR_reref{:};
    if ~isempty(SNR_reref)
        plot(repmat(idx,length(SNR_reref),1),SNR_reref,'.')
    end
end
hold off
title('SNR reref')


subplot(7,1,3)
hold on
for idx = 2:height(sorted_micros)
    SNR_clean = sorted_micros(idx,:).SNR_clean{:};
    if ~isempty(SNR_clean)
        plot(repmat(idx,length(SNR_clean),1),SNR_clean,'.')
    end
end
hold off
title('SNR clean')

subplot(7,1,4)
hold on
for idx = 2:height(sorted_micros)
    isi_SNR = sorted_micros(idx,:).isi_SNR{:};
    if ~isempty(isi_SNR)
        plot(repmat(idx,length(isi_SNR),1),isi_SNR,'.')
    end
end
hold off
title('SNR ISI')

subplot(7,1,5)
hold on
for idx = 2:height(sorted_micros)
    su_counts = sorted_micros(idx,:).su_counts;
    if ~isnan(su_counts)
        plot(idx,su_counts,'.')
    end
end
hold off
title('su counts')

subplot(7,1,6)
hold on
for idx = 2:height(sorted_micros)
    mu_counts = sorted_micros(idx,:).mu_counts;
    if ~isnan(mu_counts)
        plot(idx,mu_counts,'.')
    end
end
hold off
title('mu counts')

subplot(7,1,7)
hold on
for idx = 2:height(sorted_micros)
    noise_counts = sorted_micros(idx,:).noise_counts;
    if ~isnan(noise_counts)
        plot(idx,noise_counts,'.')
    end
end
hold off
title('noise counts')

end

function [timed_micros,all_channels,all_neurons,all_classes] = process_micros_data(dataDir,timed_micros)
timed_micros.folder = repelem({''},height(timed_micros),1);
unique_channels = [];
bad_micros = false(height(timed_micros),1);
for idx = 1:height(timed_micros)
    % Get all folder names from date, task, and part
    subject = timed_micros(idx,:).subject{:};
    date = timed_micros(idx,:).date{:};
    task = timed_micros(idx,:).task{:};
    part = timed_micros(idx,:).part;
    
    folder = strcat(date,'_',task);
    if part > 0
        folder = strcat(folder,'_',sprintf('part%d',part));
    end
    timed_micros(idx,:).folder = {folder};
    
    %Get min and max channel numbers in bank
    channel_numbers = timed_micros(idx,:).channels{:};
    
    % Avoid cases where jacksheet labels did not allow for accurate identification of a bank's channel numbers
    % These would have to be corrected manually
    if ~(isempty(channel_numbers)) && sum(channel_numbers)>0
        %Get all channel numbers
        channel_numbers = channel_numbers(1):channel_numbers(2);
        channel_numbers = arrayfun(@(x) num2str(x),channel_numbers','UniformOutput',false);
        %Concatenate string to subject code to be able to identify all
        %unique channels
        these_channels = strcat(subject,'_',channel_numbers);
        unique_channels = [unique_channels;these_channels];
    else
        bad_micros(idx) = true;
    end
end

%Filter out those with innacurate channel number identification
timed_micros = timed_micros(~bad_micros,:);

%Get unique subjects, sessions, electrodes, and channels IDd
[~,~,subject_IDs] = unique(timed_micros.subject,'stable');
[~,~,electrode_IDs] = unique(strcat(timed_micros.subject,'_',timed_micros.bank),'stable');
[~,~,session_IDs] = unique(strcat(timed_micros.subject,'_',timed_micros.folder),'stable');
unique_channels = unique(unique_channels,'stable');
timed_micros.subject_ID = subject_IDs;
timed_micros.electrode_ID = electrode_IDs;
timed_micros.session_ID = session_IDs;

%Declare empty tables to collect all channels, neurons, and classes
all_channels = table;
all_neurons = table;
all_classes = table;

%Loop through timed micros table to load all channels, neurons, and
%classes; associate them to their respective unique IDs
for idx = 1:height(timed_micros)
    this_micro = timed_micros(idx,:);
    
    %Two sets of channel numbers:
    %One for indexing SNR of split and rerefed data (rerefed contains one
    %NaN)
    this_micro_channels = this_micro.channels{:};
    this_micro_channels = this_micro_channels(1):this_micro_channels(2);
    %One for indexing SNR of denoised reref data (no NaNs);
    this_micro_rerefs = this_micro_channels(this_micro_channels~=this_micro.reref_ch);
    
    subject = this_micro.subject{:};
    folder = this_micro.folder{:};
    bank = this_micro.bank{:};
    
    subject_ID = this_micro.subject_ID;
    session_ID = this_micro.session_ID;
    electrode_ID = this_micro.electrode_ID;
    
    %Signal to noise ratios of every channel in the electrode (bank)
    SNR_split = this_micro.SNR_split{:};
    SNR_reref = this_micro.SNR_reref{:};
    SNR_clean = this_micro.SNR_clean{:};
    
    %Whether digital or analog headstage was used
    digital_headstage = this_micro.digital_headstage;
    
    %Load channels from spike data directory
    spikeDir = fullfile(dataDir,subject,folder,'spike_data',sprintf('Bank%s',bank));
    
    channels_file = fullfile(spikeDir,'channels.mat');
    neurons_file = fullfile(spikeDir,'neurons.mat');
    classes_file = fullfile(spikeDir,'classes.mat');
    
    load(channels_file,'channels');
    load(neurons_file,'neurons');
    load(classes_file,'classes');
    
    %Avoid cases where no channels of the electrode yielded any spike
    %groups or classes (channels would be an empty array)
    if ~isempty(channels)
        %Convert channels, neurons, and classes structures to table for
        %concatenation and addition of upper level information
        channels = struct2table(channels);
        channel_numbers = channels.channel_number;
        
        channel_IDs = arrayfun(@(x) num2str(x),channel_numbers,'UniformOutput',false);
        channel_IDs = strcat(subject,'_',channel_IDs);
        channel_IDs = cell2mat(cellfun(@(x) find(strcmp(unique_channels,x),1,'first'),channel_IDs,'UniformOutput',false));
        
        channels.SNR_split = cell2mat(arrayfun(@(x) SNR_split(this_micro_channels == x),channel_numbers,'UniformOutput',false));
        channels.SNR_reref = cell2mat(arrayfun(@(x) SNR_reref(this_micro_channels == x),channel_numbers,'UniformOutput',false));
        channels.SNR_clean = cell2mat(arrayfun(@(x) SNR_clean(this_micro_rerefs == x),channel_numbers,'UniformOutput',false));
        channels.digital_headstage = repmat(digital_headstage,height(channels),1);
        channels.channel_ID = channel_IDs;
       
        neurons = struct2table(neurons);
        if iscell(neurons.peak_width)     %Dev: Need to fix variable for cases where it might be empty so that it does not become a cell array
            peak_widths = neurons.peak_width;
            peak_widths = cellfun(@(x) sum(x),peak_widths,'UniformOutput',false); %This takes care of empty cells, converting them to 0.
            neurons.peak_width = [peak_widths{:}]';
        end
        
        classes = struct2table(classes);
        
        %Add subject, session, electrode, and channel IDs to tables
        %Delete subject, folder, and bank columns
        channels.subject_ID = repmat(subject_ID,height(channels),1);
        channels.subject = [];
        channels.session_ID = repmat(session_ID,height(channels),1);
        channels.folder = [];
        channels.electrode_ID = repmat(electrode_ID,height(channels),1);
        channels.bank = [];
        
        neurons.subject_ID = repmat(subject_ID,height(neurons),1);
        neurons.subject = [];
        neurons.session_ID = repmat(session_ID,height(neurons),1);
        neurons.folder = [];
        neurons.electrode_ID = repmat(electrode_ID,height(neurons),1);
        neurons.bank = [];
        neurons.channel_ID = arrayfun(@(x) channel_IDs(find(channel_numbers == x,1,'first')),neurons.channel_number);     

        neurons.n_groups = NaN(height(neurons),1);
        neurons.SNR_split = NaN(height(neurons),1);
        neurons.SNR_reref = NaN(height(neurons),1);
        neurons.SNR_clean = NaN(height(neurons),1);
        neurons.digital_headstage = repmat(digital_headstage,height(neurons),1);
        
        for jdx = 1:height(neurons)
            channel_ID = neurons(jdx,:).channel_ID;
            sign = neurons(jdx,:).is_neg;
            option = neurons(jdx,:).option;
            
            has_channel = channels.channel_ID == channel_ID;
            has_sign = channels.is_neg == sign;
            has_option = channels.option == option;
            
            this_channel = channels(has_channel & has_sign & has_option,:);
            neurons(jdx,:).n_groups = this_channel.n_groups;
            neurons(jdx,:).SNR_split = this_channel.SNR_split;
            neurons(jdx,:).SNR_reref = this_channel.SNR_reref;
            neurons(jdx,:).SNR_clean = this_channel.SNR_clean;
        end
        
        classes.subject_ID = repmat(subject_ID,height(classes),1);
        classes.subject = [];
        classes.session_ID = repmat(session_ID,height(classes),1);
        classes.folder = [];
        classes.electrode_ID = repmat(electrode_ID,height(classes),1);
        classes.bank = [];
        classes.channel_ID = arrayfun(@(x) channel_IDs(find(channel_numbers == x,1,'first')),classes.channel_number);
        
        %Loop through all classes to add respective channel's n_groups and SNRs
        %(based on parent neuron) and also add parent neuron's n_classes, is_SUA,
        %is_MUA, and is_noise
        classes.n_groups = NaN(height(classes),1);
        classes.SNR_split = NaN(height(classes),1);
        classes.SNR_reref = NaN(height(classes),1);
        classes.SNR_clean = NaN(height(classes),1);
        classes.n_classes = NaN(height(classes),1);
        classes.is_SUA = false(height(classes),1);
        classes.is_MUA = false(height(classes),1);
        classes.is_noise = false(height(classes),1);
        classes.digital_headstage = repmat(digital_headstage,height(classes),1);
        
        for jdx = 1:height(classes)
            % Get neuron channel ID, session ID, sign and option to find
            % corresponding channel structure
            this_channel_ID = classes(jdx,:).channel_ID;
            this_sign = classes(jdx,:).is_neg;
            this_option = classes(jdx,:).option;
            this_group = classes(jdx,:).group;
            
            has_channel = neurons.channel_ID == this_channel_ID;
            has_sign = neurons.is_neg == this_sign;
            has_option = neurons.option == this_option;
            has_group = neurons.group == this_group;
            
            this_neuron = neurons(has_channel & has_sign & has_option & has_group,:);
            classes(jdx,:).n_groups = this_neuron.n_groups;
            classes(jdx,:).SNR_split = this_neuron.SNR_split;
            classes(jdx,:).SNR_reref = this_neuron.SNR_reref;
            classes(jdx,:).SNR_clean = this_neuron.SNR_clean;
            classes(jdx,:).n_classes = this_neuron.n_classes;
            classes(jdx,:).is_SUA = this_neuron.is_SUA;
            classes(jdx,:).is_MUA = this_neuron.is_MUA;
            classes(jdx,:).is_noise = this_neuron.is_noise;
        end
        
        all_channels = [all_channels;channels];
        all_neurons = [all_neurons;neurons];
        all_classes = [all_classes;classes];
    end
end
end

function su_vs_noise_analysis(qualityDir,observationDir,observation_type,observation_vars,observation_titles,all_observations)

%%% Plot parameters
figWidth = 4480;
figHeight = 2520; %%%Fit to screen
nRows = size(observation_vars,1);
nColumns = size(observation_vars,2);
subplotHeight = round(figHeight/nRows);
subplotWidth = round(figWidth/nColumns);

su_color = [17 138 178]/255;
mu_color = [252 158 79]/255;
noise_color = [222 108 131]/255;

yy1 = [0 4];
yy2 = 1:3;
yy3 = {'N','MU','SU'};

thisPlotDir = fullfile(qualityDir,observationDir);
if ~isfolder(thisPlotDir)
    mkdir(thisPlotDir);
end
plotName = fullfile(thisPlotDir,strcat(datestr(now,'yyyy-mm-dd'),'_',observation_type,'_sua_vs_noise'));

formula = 'var ~ is_SUA + (1|subject_ID) + (1|subject_ID:electrode_ID)';

figure('Units', 'pixels', 'Visible', 'off')
fig = gcf;
axesHandles = zeros(nRows,nColumns);
fig.Position(3) = figWidth;
fig.Position(4) = figHeight;

for idx = 1:nRows
    for jdx = 1:nColumns
        this_var = observation_vars{idx,jdx};
        this_title = observation_titles{idx,jdx};
        
        these_observations = table;
        these_observations.subject_ID = all_observations.subject_ID;
        these_observations.electrode_ID = all_observations.electrode_ID;
        these_observations.is_SUA = all_observations.is_SUA;
        these_observations.is_MUA = all_observations.is_MUA;
        these_observations.is_noise = all_observations.is_noise;
        these_observations.var = eval(sprintf('all_observations.%s',this_var));
        
        these_observations = these_observations(~isnan(these_observations.var),:);
        these_observations = these_observations(~isinf(these_observations.var),:);
        if ~isreal(these_observations.var)
            these_observations = these_observations(imag(these_observations.var)==0,:);
            these_observations.var = double(these_observations.var);
        end
        this_analysis = these_observations(~these_observations.is_MUA,:);
        
        this_analysis.subject_ID = categorical(this_analysis.subject_ID);
        this_analysis.electrode_ID = categorical(this_analysis.electrode_ID);
        
        lme_su_vs_noise = fitlme(this_analysis,formula);
        
        tstat = lme_su_vs_noise.Coefficients.tStat(2);
        pval = lme_su_vs_noise.Coefficients.pValue(2);
        
        lengths = zeros(3,1);
        
        single_units = these_observations(these_observations.is_SUA,:);
        single_units.y = (randn(height(single_units),1)*0.15) + 3;
        lengths(3) = height(single_units);
        median_sua = median(single_units.var);
        q1_sua = quantile(single_units.var, 0.25);
        q3_sua = quantile(single_units.var, 0.75);
        low_sua = q1_sua - (1.5*(q3_sua - q1_sua));
        high_sua = q3_sua + (1.5*(q3_sua - q1_sua));

        multi_units = these_observations(these_observations.is_MUA,:);
        multi_units.y = (randn(height(multi_units),1)*0.15) + 2;
        lengths(2) = height(multi_units);
        median_mua = median(multi_units.var);
        q1_mua = quantile(multi_units.var, 0.25);
        q3_mua = quantile(multi_units.var, 0.75);
        low_mua = q1_mua - (1.5*(q3_mua - q1_mua));
        high_mua = q3_mua + (1.5*(q3_mua - q1_mua));

        noise_units = these_observations(these_observations.is_noise,:);
        noise_units.y = (randn(height(noise_units),1)*0.15) + 1;
        lengths(1) = height(noise_units);
        median_noise = median(noise_units.var);
        q1_noise = quantile(noise_units.var, 0.25);
        q3_noise = quantile(noise_units.var, 0.75);
        low_noise = q1_noise - (1.5*(q3_noise - q1_noise));
        high_noise = q3_noise + (1.5*(q3_noise - q1_noise));
        
        diff0 = diff([min(these_observations.var),max(these_observations.var)])==0;
        if diff0
            xx1 = [min(these_observations.var)-0.5,min(these_observations.var)+0.5];
            xx2 = xx1(1):0.5:xx1(2);
        else
            xx1 = [min(these_observations.var),max(these_observations.var)];
            diffxx1 = diff(xx1);
            xx2 = xx1(1):diffxx1/6:xx1(2);
            xx1(1) = xx1(1) - (diffxx1*0.15);
            xx1(2) = xx1(2) + (diffxx1*0.15);
        end
        xx3 = arrayfun(@(x) sprintf('%.2f',x),xx2,'UniformOutput',false);
        
        axesHandles(idx,jdx) = axes('Parent', fig, 'Units', 'pixels', ...
        'Position', [(jdx-1)*subplotWidth, fig.Position(4)-(idx*subplotHeight), subplotWidth, subplotHeight]);
    
        hold on
        
        text(mean(xx1),3.85,this_title,'FontSize',16,'HorizontalAlignment','center')
        text(mean(xx1),3.7,sprintf('t:%.2f, p:%.3f',tstat,pval),'FontSize',14,'HorizontalAlignment','center','Color',[0.3 0.3 0.3])
        
        plot(single_units.var,single_units.y,'.','Color',su_color)
        plot([low_sua high_sua],[3 3],'-','LineWidth',1.5,'Color',[0.4 0.4 0.4])
        plot([q1_sua q3_sua],[3 3],'-','LineWidth',4.5,'Color',[0.4 0.4 0.4]);
        plot([median_sua median_sua],[2.9 3.1],'-w','LineWidth',2)
        plot(median_sua,3,'o','Color',[1 0 1])
        
        plot(multi_units.var,multi_units.y,'.','Color',mu_color)
        plot([low_mua high_mua],[2 2],'-','LineWidth',1.5,'Color',[0.4 0.4 0.4])
        plot([q1_mua q3_mua],[2 2],'-','LineWidth',4.5,'Color',[0.4 0.4 0.4]);
        plot([median_mua median_mua],[1.9 2.1],'-w','LineWidth',2)
        plot(median_mua,2,'o','Color',[1 0 1])
        
        plot(noise_units.var,noise_units.y,'.','Color',noise_color)
        plot([low_noise high_noise],[1 1],'-','LineWidth',1.5,'Color',[0.4 0.4 0.4])
        plot([q1_noise q3_noise],[1 1],'-','LineWidth',4.5,'Color',[0.4 0.4 0.4]);
        plot([median_noise median_noise],[0.9 1.1],'-w','LineWidth',2)
        plot(median_noise,1,'o','Color',[1 0 1])
        
        for kdx = 1:length(xx2)
            text(xx2(kdx),0.2,xx3{kdx},'FontSize',12,'HorizontalAlignment','center');
        end
        
        for kdx = 1:length(yy2)
            text(xx1(1)+(diff(xx1)*0.025),yy2(kdx),yy3{kdx},'FontSize',16);
            text(xx1(2)-(diff(xx1)*0.1),yy2(kdx),sprintf('n = %d',lengths(kdx)),'FontSize',12);
        end
        
        hold off
        
        xlim(xx1); xticks(xx2); xticklabels([]);
        ylim(yy1); yticks(yy2); yticklabels([]);
    end
end
print(plotName,'-dpng');
print(plotName,'-dsvg');
close all
end

function digital_vs_analog_analysis(qualityDir,observationDir,observation_type,observation_vars,observation_titles,observations)
%%% Plot parameters
figWidth = 4480;
figHeight = 2520; %%%Fit to screen
nRows = size(observation_vars,1);
nColumns = size(observation_vars,2);
subplotHeight = round(figHeight/nRows);
subplotWidth = round(figWidth/nColumns);

digital_color = [17 138 178]/255;
analog_color = [252 158 79]/255;

yy1 = [0 3];
yy2 = 1:2;
yy3 = {'A','D'};

thisPlotDir = fullfile(qualityDir,observationDir);
if ~isfolder(thisPlotDir)
    mkdir(thisPlotDir);
end
plotName = fullfile(thisPlotDir,strcat(datestr(now,'yyyy-mm-dd'),'_',observation_type,'_digital_vs_analog'));

formula = 'var ~ digital_headstage + (1|subject_ID) + (1|subject_ID:electrode_ID)';

figure('Units', 'pixels', 'Visible', 'off')
fig = gcf;
axesHandles = zeros(nRows,nColumns);
fig.Position(3) = figWidth;
fig.Position(4) = figHeight;

for idx = 1:nRows
    for jdx = 1:nColumns
        this_var = observation_vars{idx,jdx};
        this_title = observation_titles{idx,jdx};
        
        these_observations = table;
        these_observations.subject_ID = observations.subject_ID;
        these_observations.electrode_ID = observations.electrode_ID;
        these_observations.digital_headstage = observations.digital_headstage;
        these_observations.var = eval(sprintf('all_observations.%s',this_var));
        
        these_observations = these_observations(~isnan(these_observations.var),:);
        these_observations = these_observations(~isinf(these_observations.var),:);
        if ~isreal(these_observations.var)
            these_observations = these_observations(imag(these_observations.var)==0,:);
            these_observations.var = double(these_observations.var);
        end
        this_analysis = these_observations;
        
        this_analysis.subject_ID = categorical(this_analysis.subject_ID);
        this_analysis.electrode_ID = categorical(this_analysis.electrode_ID);
        
        lme_digital_vs_analog = fitlme(this_analysis,formula);
        
        tstat = lme_digital_vs_analog.Coefficients.tStat(2);
        pval = lme_digital_vs_analog.Coefficients.pValue(2);
        
        diff0 = diff([min(these_observations.var),max(these_observations.var)])==0;
        if diff0
            xx1 = [min(these_observations.var)-0.5,min(these_observations.var)+0.5];
            xx2 = xx1(1):0.5:xx1(2);
        else
            xx1 = [min(these_observations.var),max(these_observations.var)];
            diffxx1 = diff(xx1);
            xx2 = xx1(1):diffxx1/6:xx1(2);
            xx1(1) = xx1(1) - (diffxx1*0.15);
            xx1(2) = xx1(2) + (diffxx1*0.15);
        end
        xx3 = arrayfun(@(x) sprintf('%.2f',x),xx2,'UniformOutput',false);
        
        lengths = zeros(2,1);
        
        digital = these_observations(these_observations.digital_headstage,:);
        digital.y = (randn(height(digital),1)*0.15) + 2;
        lengths(2) = height(digital);
        median_digital = median(digital.var);
        q1_digital = quantile(digital.var, 0.25);
        q3_digital = quantile(digital.var, 0.75);
        low_digital = q1_digital - (1.5*(q3_digital - q1_digital));
        low_digital = max([low_digital,min(these_observations.var)]);
        high_digital = q3_digital + (1.5*(q3_digital - q1_digital));
        high_digital = min([high_digital,max(these_observations.var)]);

        analog = these_observations(~these_observations.digital_headstage,:);
        analog.y = (randn(height(analog),1)*0.15) + 1;
        lengths(1) = height(analog);
        median_analog = median(analog.var);
        q1_analog = quantile(analog.var, 0.25);
        q3_analog = quantile(analog.var, 0.75);
        low_analog = q1_analog - (1.5*(q3_analog - q1_analog));
        low_analog = max([low_analog,min(these_observations.var)]);
        high_analog = q3_analog + (1.5*(q3_analog - q1_analog));
        high_analog = min([high_analog,max(these_observations.var)]);
        
        axesHandles(idx,jdx) = axes('Parent', fig, 'Units', 'pixels', ...
        'Position', [(jdx-1)*subplotWidth, fig.Position(4)-(idx*subplotHeight), subplotWidth, subplotHeight]);
    
        hold on
        
        text(mean(xx1),2.85,this_title,'FontSize',16,'HorizontalAlignment','center')
        text(mean(xx1),2.7,sprintf('t:%.2f, p:%.3f',tstat,pval),'FontSize',14,'HorizontalAlignment','center','Color',[0.3 0.3 0.3])
        
        plot(digital.var,digital.y,'.','Color',digital_color)
        plot([low_digital high_digital],[2 2],'-','LineWidth',1.5,'Color',[0.4 0.4 0.4])
        plot([q1_digital q3_digital],[2 2],'-','LineWidth',4.5,'Color',[0.4 0.4 0.4]);
        plot([median_digital median_digital],[1.9 2.1],'-w','LineWidth',2)
        plot(median_digital,2,'o','Color',[1 0 1])
        
        plot(analog.var,analog.y,'.','Color',analog_color)
        plot([low_analog high_analog],[1 1],'-','LineWidth',1.5,'Color',[0.4 0.4 0.4])
        plot([q1_analog q3_analog],[1 1],'-','LineWidth',4.5,'Color',[0.4 0.4 0.4]);
        plot([median_analog median_analog],[0.9 1.1],'-w','LineWidth',2)
        plot(median_analog,1,'o','Color',[1 0 1])
        
        for kdx = 1:length(xx2)
            text(xx2(kdx),0.2,xx3{kdx},'FontSize',12,'HorizontalAlignment','center');
        end
        
        for kdx = 1:length(yy2)
            text(xx1(1)+(diff(xx1)*0.025),yy2(kdx),yy3{kdx},'FontSize',16);
            text(xx1(2)-(diff(xx1)*0.1),yy2(kdx),sprintf('n = %d',lengths(kdx)),'FontSize',12);
        end
        
        hold off
        
        xlim(xx1); xticks(xx2); xticklabels([]);
        ylim(yy1); yticks(yy2); yticklabels([]);
    end
end
print(plotName,'-dpng');
print(plotName,'-dsvg');
close all
end