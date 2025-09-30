function plot_option_metrics(varargin)
if isempty(varargin)
    %%% Enter manually
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% _directoryectory with pipeline and database
    subject = 'SC000';                                             %%% Subject code
    folder = 'yyyy-mm-dd_task_part1';                              %%% Folder of session in yyyy-mm-dd_task or yyyy-mm-dd_task_part1 format
    bank = 'A';                                                    %%% Character of bank
else
    %%% Or enter in function in this order
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    bank = varargin{4};
end

%%% List directories
data_directory = fullfile(root_directory, 'micros_database');
combinato_directory = fullfile(data_directory, subject, folder, 'combinato_files', sprintf('Bank%s', bank));
spikes_directory = strrep(combinato_directory, 'combinato_files', 'spike_data');
quality_directory = fullfile(data_directory, 'quality_metrics', subject, folder);
if ~isfolder(quality_directory)
    mkdir(quality_directory);
end

%%% Make list of all possible combinations of channels and signs
channels_file = fullfile(spikes_directory, 'channels.mat');
load(channels_file, 'channels');

channels = struct2table(channels);

channel_numbers = unique(channels.channel_number);
polarity_signs = [0; 1]; %1 being negative

[Ax, Bx] = ndgrid(1:numel(polarity_signs), 1:numel(channel_numbers));
plot_combos = table;
plot_combos.channel = channel_numbers(Bx(:));
plot_combos.polarity_sign = polarity_signs(Ax(:));

plot_variables = {'n_sua', 'n_groups', 'n_classes', 'CHI', 'n_mua', 'n_spikes', 'n_spc', 'DBI', 'n_noise', 'n_bad_spikes', 'n_2nd', 'm_silhouette'};
%%% Number of: 
%%% groups, classes, spikes, artifact or skipped (bad) spikes,
%%% spikes matched at SPC or 2nd template; and classified as single or multi unit activity or noise. 
%%% Clustering scores (considering SUs only):
%%% Calinski-Harabasz Index, Davies-Bouldin Index, Silhouette Score.

%%% Normalize plot parameters as a percentage of maxima (0-1) across channels, signs, and options. 
%%% CHI, DBI, and m_silhouette are NaNs in the case where there is only one group. Change to 0 for plotting. 
%%% m_silhouette NaNs changed to 0 and scaled based on minimum and maximum
%%% values, centered at 0.

%%% Replace NaNs with zeros
channels.CHI(find(isnan(channels.CHI))) = zeros(sum(isnan(channels.CHI)), 1);
channels.DBI(find(isnan(channels.DBI))) = zeros(sum(isnan(channels.DBI)), 1);
channels.m_silhouette(find(isnan(channels.m_silhouette))) = zeros(sum(isnan(channels.m_silhouette)), 1);

%%% Get maxima and control for possible 0s to avoid Inf values
max_n_groups = max(channels.n_groups);
max_n_classes = max(channels.n_classes);
max_n_spikes = max(channels.n_spikes);
max_n_bad_spikes = max(channels.n_bad_spikes);
max_n_spc = max(channels.n_spc);
max_n_2nd = max(channels.n_2nd);
max_n_sua = max(channels.n_sua);
max_n_mua = max(channels.n_mua);
max_n_noise = max(channels.n_noise);
max_CHI = max(channels.CHI);
max_DBI = max(channels.DBI);
max_m_silhouette = max([abs(min(channels.m_silhouette)), abs(max(channels.m_silhouette))]);
if max_n_groups == 0
    max_n_groups = 1;
end
if max_n_classes == 0
    max_n_classes = 1;
end
if max_n_spikes == 0
    max_n_spikes = 1;
end
if max_n_bad_spikes == 0
    max_n_bad_spikes = 1;
end
if max_n_spc == 0
    max_n_spc = 1;
end
if max_n_2nd == 0
    max_n_2nd = 1;
end
if max_n_sua == 0
    max_n_sua = 1;
end
if max_n_mua == 0
    max_n_mua = 1;
end
if max_n_noise == 0
    max_n_noise = 1;
end
if max_CHI == 0
    max_CHI = 1;
end
if max_DBI == 0
    max_DBI = 1;
end
if max_m_silhouette == 0
    max_m_silhouette = 0.1;
end

%%% Rescale data
channels.n_groups = channels.n_groups/max_n_groups;
channels.n_classes = channels.n_classes/max_n_classes;
channels.n_spikes = channels.n_spikes/max_n_spikes;
channels.n_bad_spikes = channels.n_bad_spikes/max_n_bad_spikes;
channels.n_spc = channels.n_spc/max_n_spc;
channels.n_2nd = channels.n_2nd/max_n_2nd;
channels.n_sua = channels.n_sua/max_n_sua;
channels.n_mua = channels.n_mua/max_n_mua;
channels.n_noise = channels.n_noise/max_n_noise;
channels.CHI = channels.CHI/max_CHI;
channels.DBI = channels.DBI/max_DBI;
channels.m_silhouette = channels.m_silhouette/max_m_silhouette;

%%% For each channel and sign combination, only if at least one option yielded
%%% at least one group, there will be 12 subplots (one for each metric),
%%% arranged in 4 rows and 3 columns. Each subplot consists of a grid that is
%%% 4 by 4 values (actually 2x2 grid of 2x2 grids), where top and left
%%% represent options of artifact rejection, grouping, clustering, and
%%% recursion that aim for more clusters, and bottom and right represent
%%% options aiming for less clusters.

%%% Figure parameters
figure_width = 4480;
figure_height = 2520; %%%Fit to screen

%%% Subplot color gradient
map_bad_counts = makecolormap_EF('single_gradient1');
map_good_counts = makecolormap_EF('single_gradient2');
map_neutral_counts = makecolormap_EF('single_gradient3');
map_divergent = makecolormap_EF('uniform4');

%%% Subplot parameters
n_rows = 4;          %%% number of rows in grid
n_columns = 4;       %%% number of columns in grid
xlim1 = [0.5 4.5];
xlabel1 = 'recursion(out) / clustering(in)';
xticks1 = [1 2 2.5 4.5];
xticklabels1 = {'O', 'U', 'O', 'U'};
ylim1 = [0.5 4.5];
ylabel1 = 'grouping(out) / artifacts(in)';
yticks1 = [0.5 2.5 3 4];
yticklabels1 = {'O', 'U', 'O', 'U'};


%%% Loop through channel and sign combination
for idx = 1:height(plot_combos)
    %%% Plot name is channel_number and sign
    channel_number = plot_combos(idx, :).channel;
    polarity_sign = plot_combos(idx, :).polarity_sign;
    channel_label = sprintf('ch%03d', channel_number);
    if polarity_sign == 1
        polarity_sign_label = 'neg';
    else
        polarity_sign_label = 'pos';
    end
    plot_name = strcat(channel_label, '_', polarity_sign_label);
    
    %%% Get all options: 0 (default), 1-16 (grid-search options)
    these_options = channels(channels.channel_number == channel_number & channels.is_neg == polarity_sign, :);
    
    %%% Only create plot for channel sign combinatons were options yielded at least one group
    if~isempty(these_options)
        figure('Units', 'pixels', 'Visible', 'off')
        figure_handle = gcf;
        figure_handle.Position(3) = figure_width;
        figure_handle.Position(4) = figure_height;
        
        %%% Separate default option value and grid-search option values
        default_option = these_options(these_options.option == 0, :);
        grid_options = these_options(these_options.option > 0, :);
        
        %%% Loop through variables to make subplots
        for jdx = 1:length(plot_variables)
            this_variable = plot_variables{jdx};
            
            %%% Get values of variable for default and grid options
            %%% Fill in values of options that yielded no groups with zeros
            if isempty(default_option)
                default_value = 0;
            else
                default_value = eval(sprintf('default_option.%s', this_variable));
            end
            
            if isempty(grid_options)
                grid_values = zeros(16, 1);
            else
                temp_values = eval(sprintf('grid_options.%s', this_variable));
                indices = grid_options.option;
                grid_values = zeros(16, 1);
                for kdx = 1:length(indices)
                    grid_values(indices(kdx)) = temp_values(kdx);
                end
            end            
            
            %%% Get maximum value for labeling
            max_value = eval(sprintf('max_%s', this_variable));
            
            %%% Define colorbar ticks and labels for default and max values
            if strcmp(this_variable, 'm_silhouette')
                if abs(default_value) == max_value
                    default_value = 0.95*max_value;
                elseif default_value == 0
                    default_value = 0.05;
                end
                c_lim = [-1, 1];
                c_ticks = sort(unique([-1, 0, 1, default_value]));
                c_ticklabels = arrayfun(@(x) num2str(x), c_ticks, 'UniformOutput', false);
                c_ticklabels{c_ticks == default_value} = '*';
                c_ticklabels{c_ticks == 1} = sprintf('%.3f', max_value);
                c_ticklabels{c_ticks == -1} = sprintf('%.3f', (max_value*-1));
            else
                if default_value == max_value
                    default_value = 0.95*max_value;
                end
                if default_value == 0
                    default_value = 0.05;
                end
                c_lim = [0, 1];
                c_ticks = sort(unique([0, 1, default_value]));
                c_ticklabels = arrayfun(@(x) num2str(x), c_ticks, 'UniformOutput', false);
                c_ticklabels{c_ticks == default_value} = '*';
                c_ticklabels{c_ticks == 1} = num2str(max_value);
            end
            
            %%% Choose colorbar color
            switch this_variable
                case {'n_noise', 'n_2nd', 'n_bad_spikes'}
                    mapx = map_bad_counts;
                case {'n_spikes', 'n_groups', 'n_classes', 'n_mua'}
                    mapx = map_neutral_counts;
                case {'DBI', 'CHI', 'n_sua', 'n_spc'}
                    mapx = map_good_counts;
                case {'m_silhouette'}
                    mapx = map_divergent;
            end
            
            %%% Reshape grid search values in 4 by 4 and plot with imagesc
            grid_values = reshape(grid_values, 4, 4);
            
            subplot(3, 4, jdx)
            imagesc(grid_values)
            
            %%% Plot lines that divide grid into inner and outter variables
            hold on
            plot(repmat(1.5:1:(n_columns-0.5), 2, 1), repmat([0.5;n_rows+0.5], 1, n_columns-1), '-k', 'Line_width', 1)
            plot(repmat([0.5;(n_columns+0.5)], 1, n_rows-1), repmat(1.5:1:(n_rows-0.5), 2, 1), '-k', 'Line_width', 1)
            plot(repmat(2.5, 2, 1), repmat([0.5;(n_rows+0.5)], 1), '-k', 'Line_width', 3)
            plot(repmat([0.5;n_columns+0.5], 1, 1), repmat(2.5, 2, 1), '-k', 'Line_width', 3)
            
            %%% Add labels
            xlabel(xlabel1);xlim(xlim1);xticks(xticks1);xticklabels(xticklabels1);
            ylabel(ylabel1);ylim(ylim1);yticks(yticks1);yticklabels(yticklabels1);
            
            %%% Add colorbar and set colorbar labels
            colormap(gca, mapx)
            caxis(c_lim)
            cb = colorbar;
            cb.Location = 'southoutside';
            set(cb, 'Ticks', c_ticks)
            set(cb, 'TickLabels', c_ticklabels)
            title(strrep(this_variable, '_', ' '))
            hold off
        end
        sgtitle({'* is default option value for this channel marked in respective colorbar.', ...
            'Colorbar gradients are based on results for all channels, signs, and options.'})
        print(fullfile(quality_directory, plot_name), '-dsvg')
        print(fullfile(quality_directory, plot_name), '-dpng')
        close all
    end
end
end