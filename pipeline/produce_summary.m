function produce_summary(varargin)
if isempty(varargin)
    rootDir = '/project/TIBIR/Lega_lab/shared/lega_ansir/Single_Units_BR';
else
    rootDir = varargin{1};
end

dataDir = fullfile(rootDir,'micros_database');
sorted_file = fullfile(dataDir,'sorted_micros.mat');
load(sorted_file,'sorted_micros')

sorted_micros.manual = true(height(sorted_micros),1);

su_counts = zeros(height(sorted_micros),1);
mu_counts = zeros(height(sorted_micros),1);
noise_counts = zeros(height(sorted_micros),1);
error_number = zeros(height(sorted_micros),1);

for idx = 2:height(sorted_micros)
    
    this_micro = sorted_micros(idx,:);
    
    subject = this_micro.subject{:};
    date = this_micro.date{:};
    task = this_micro.task{:};
    part = this_micro.part;
    bank = this_micro.bank{:};
    spikes_timed = this_micro.spikes_timed;
    has_error = this_micro.has_error;
    
    folder = strcat(date,'_',task);
    if part > 0
        folder = strcat(folder,'_',sprintf('part%d',part));
    end
    if spikes_timed
        spikeDir = fullfile(dataDir,subject,folder,'spike_data',sprintf('Bank%s',bank));
        channels_file = fullfile(spikeDir,'channels.mat');
        load(channels_file,'channels')
        if ~isempty(channels)
            su_counts(idx) = sum([channels(~isnan([channels.n_sua]),:).n_sua]);
            mu_counts(idx) = sum([channels(~isnan([channels.n_mua]),:).n_mua]);
            noise_counts(idx) = sum([channels(~isnan([channels.n_noise]),:).n_noise]);
        end
        clear channels
    end
    if has_error
        progress = [this_micro{1,18:25}];
        if any(progress == 0)
            error_number(idx) = find(progress == 0,1,'first');
        end
    end
end

sorted_micros.su_counts = su_counts;
sorted_micros.mu_counts = mu_counts;
sorted_micros.noise_counts = noise_counts;
sorted_micros.error_number = error_number;

tasks = unique(sorted_micros.task);
tasks = tasks(~strcmp(tasks,'ZZ'));
tasks{end+1} = 'all';

summary_table = table;
for idx = 1:length(tasks)
    this_task = tasks{idx};
    
    if strcmp(this_task,'all')
        this_micros = sorted_micros(2:end,:);
    else
        this_micros = sorted_micros(strcmp(sorted_micros.task,this_task),:);
    end
    
    temp_summary = table;
    temp_summary.task = {this_task};
    
    subjects = this_micros.subject;
    unique_subjects = unique(subjects);    
    has_su = false(length(unique_subjects),1);
    avg_su = zeros(length(unique_subjects),1);
    has_mu = false(length(unique_subjects),1);
    avg_mu = zeros(length(unique_subjects),1);
    all_error = false(length(unique_subjects),1);
    
    for jdx = 1:length(unique_subjects)
        this_subject = unique_subjects{jdx};
        subject_micros = this_micros(strcmp(subjects,this_subject),:);
        has_su(jdx) = any(subject_micros.su_counts>0);
        avg_su(jdx) = mean(subject_micros(subject_micros.su_counts>0,:).su_counts,'omitnan');
        has_mu(jdx) = any(subject_micros.mu_counts>0);
        avg_mu(jdx) = mean(subject_micros(subject_micros.mu_counts>0,:).mu_counts,'omitnan');
        all_error(jdx) = ~any(subject_micros.error_number == 0);
    end
    
    temp_summary.n_subjects = length(unique_subjects);
    temp_summary.sub_has_su = sum(has_su);
    temp_summary.sub_avg_su = mean(avg_su(has_su),'omitnan');
    temp_summary.sub_has_mu = sum(has_mu);
    temp_summary.sub_avg_mu = mean(avg_mu(has_mu),'omitnan');
    temp_summary.sub_all_error = sum(all_error);
    
    sessions = strcat(this_micros.subject,this_micros.date,this_micros.task,cell2mat(arrayfun(@(x) num2str(x),this_micros.part,'UniformOutput',false)));
    unique_sessions = unique(sessions);
    has_su = false(length(unique_sessions),1);
    avg_su = zeros(length(unique_sessions),1);
    has_mu = false(length(unique_sessions),1);
    avg_mu = zeros(length(unique_sessions),1);
    all_error = false(length(unique_sessions),1);
    
    for jdx = 1:length(unique_sessions)
        this_session = unique_sessions{jdx};
        session_micros = this_micros(strcmp(sessions,this_session),:);
        has_su(jdx) = any(session_micros.su_counts>0);
        avg_su(jdx) = mean(session_micros(session_micros.su_counts>0,:).su_counts,'omitnan');
        has_mu(jdx) = any(session_micros.mu_counts>0);
        avg_mu(jdx) = mean(session_micros(session_micros.mu_counts>0,:).mu_counts,'omitnan');
        all_error(jdx) = ~any(session_micros.error_number == 0);
    end
    
    temp_summary.n_sessions = length(unique_sessions);
    temp_summary.sess_has_su = sum(has_su);
    temp_summary.sess_avg_su = mean(avg_su(has_su),'omitnan');
    temp_summary.sess_has_mu = sum(has_mu);
    temp_summary.sess_avg_mu = mean(avg_mu(has_mu),'omitnan');
    temp_summary.sess_all_error = sum(all_error);
    
    electrodes = strcat(this_micros.subject,this_micros.bank);
    unique_electrodes = unique(electrodes);
    has_su = false(length(unique_electrodes),1);
    avg_su = zeros(length(unique_electrodes),1);
    has_mu = false(length(unique_electrodes),1);
    avg_mu = zeros(length(unique_electrodes),1);
    all_error = false(length(unique_electrodes),1);
    
    for jdx = 1:length(unique_electrodes)
        this_electrode = unique_electrodes{jdx};
        electrode_micros = this_micros(strcmp(electrodes,this_electrode),:);
        has_su(jdx) = any(electrode_micros.su_counts>0);
        avg_su(jdx) = mean(electrode_micros(electrode_micros.su_counts>0,:).su_counts,'omitnan');
        has_mu(jdx) = any(electrode_micros.mu_counts>0);
        avg_mu(jdx) = mean(electrode_micros(electrode_micros.mu_counts>0,:).mu_counts,'omitnan');
        all_error(jdx) = ~any(electrode_micros.error_number == 0);
    end
    
    temp_summary.n_electrodes = length(unique_electrodes);
    temp_summary.elec_has_su = sum(has_su);
    temp_summary.elec_avg_su = mean(avg_su(has_su),'omitnan');
    temp_summary.elec_has_mu = sum(has_mu);
    temp_summary.elec_avg_mu = mean(avg_mu(has_mu),'omitnan');
    temp_summary.elec_all_error = sum(all_error);
    
    summary_table = [summary_table;temp_summary];
end
writetable(summary_table,'summary.xlsx')
end