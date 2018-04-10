%% Threading Code  

% Aims to use data from State_Space_4 (clustered active and inactive bouts)
    % To study transition dynamics 
    
%% Required scripts 

%% Notes 
% The code for calculating the transition probabilities is adapted from 
 %https://uk.mathworks.com/matlabcentral/answers/57877-seeking-help-creating-a-transition-probability-matrix-for-a-markov-chain

%% Settings 
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 
set(0,'defaultfigurecolor',[1 1 1]); % white background

 %% Load data (Post State_Space_4)
    % Note that running State_Space_4 first is necessary to re-order the
        % Clusters by a suitable metric 
      
% Load  
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat');

% Load Delta_px_sq data 
for f = 1:size(filename,2) %For each file
    delta_px_sq{1,f} = load(strcat(pathname,filename{f}),'delta_px_sq'); % load delta_px_sq data  
end 

%% Threading With Multiple Shuffles of Data 
% Slow, Roughly 3 hours for 629 fish, each with 10 shuffles 
% Note that I'm now keeping only the real time windows
% Rather than storing multiple copies of these (which are redundant)

% Shuffles
shuffles = 10; % hard coded number of shuffles

% Pre-allocation
threads = cell(max(fish_tags{1,1}),3,(1+shuffles)); % fish x (clusters,times (start & stop),time windows) x (samples & shuffled controls)
idx_numComp_sorted{1,1} = idx_numComp_sorted{1,1} + max(idx_numComp_sorted{2,1}); % Assign higher numbers to the wake bouts
idx_numComp_sorted{2,1}(isnan(idx_numComp_sorted{2,1})) = 0; % replace nan-values with zero

tic
for f = 1:max(fish_tags{1,1}) % For each fish
    % Pre-allocate
    % Data
    threads{f,1,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % clusters
    threads{f,2,1} = nan(size(threads{f,1,1},1),2,'single'); % times (start & stop)
    threads{f,3,1} = nan(size(threads{f,1,1}),'single'); % time windows
    
    % Deterime starting state (a = active, b = inactive)
    if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active
        a = 1; b = 2; % start active
    else % If the fish starts inactive
        a = 2; b = 1; % start inactive
    end
    
    % Fill in data
    % Clusters
    threads{f,1,1}(a:2:end,1) = idx_numComp_sorted{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active clusters
    threads{f,1,1}(b:2:end,1) = idx_numComp_sorted{2,1}...
        (fish_tags{2,1} == f,1); % Fill in inactive clusters
    % Times
    threads{f,2,1}(a:2:end,1:2) = wake_cells...
        (fish_tags{1,1} == f,1:2); % Fill in active times
    threads{f,2,1}(b:2:end,1:2) = sleep_cells...
        (fish_tags{2,1} == f,1:2); % Fill in inactive times
    % Time Windows
    threads{f,3,1}(a:2:end,1) = parameter_indicies{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active time windows
    threads{f,3,1}(b:2:end,1) = parameter_indicies{2,1}...
        (fish_tags{2,1} == f,1); % Fill in inactive time windows
    
    % Generate Shuffled Control Data
    % Note that each time window (e.g. day 1 vs day 2 etc) is
    % shuffled individually within each fish
    % Preserving the developmental & day/night cluster frequencies
    
    % Clusters
    for tc = 2:size(threads,3) % for each shuffle
        
        % Determine Starting State
        % Note that it's necessary to do this after each shuffle to
        % re-set the starting state
        if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active
            a = 1; b = 2; % start active
        else % If the fish starts inactive
            a = 2; b = 1; % start inactive
        end
        
        % Shuffle data
        for tw = 1:max(parameter_indicies{1,1}(fish_tags{1,1} == f)) % for each time window
            data = []; scrap = []; % empty structures
            
            % take real clusters from this time window
            scrap{1,1} = idx_numComp_sorted{1,1}(fish_tags{1,1} == f & parameter_indicies{1,1} == tw,1); % active
            scrap{2,1} = idx_numComp_sorted{2,1}(fish_tags{2,1} == f & parameter_indicies{2,1} == tw,1); % inactive
            
            data = nan((size(scrap{1,1},1) + size(scrap{2,1},1)),1,'single'); % allocate (active + inactive clusters)
            data(a:2:end,1) = scrap{1,1}(randperm(length(scrap{1,1}))); % fill active clusters
            data(b:2:end,1) = scrap{2,1}(randperm(length(scrap{2,1}))); % fill inactive clusters
            
            threads{f,1,tc} = [threads{f,1,tc} ; data]; % store data
            
            % Deterime current final state (a = wake, b = sleep)
            if threads{f,1,tc}(end,1) <= numComp(2) % If the fish ends inactive
                a = 1; b = 2; % next will be active
            else % If the fish ends active
                a = 2; b = 1; % next will be inactive
            end
            
        end
    end
    
    % Report progress
    if mod(f,100) == 0 % every 100 fish
        disp(horzcat('Threaded Fish ',num2str(f),' of ',...
            num2str(max(fish_tags{1,1})))); % Report progress
    end
    
end
toc

clear f a b tc tw data scrap

%% Threading 

% % Shuffles 
% shuffles = 10; % hard coded number of shuffles 
% 
% % Pre-allocation
% threads = cell(max(fish_tags{1,1}),3,(1+shuffles)); % fish x (clusters,times (start & stop),time windows) x (samples vs controls) 
% idx_numComp_sorted{1,1} = idx_numComp_sorted{1,1} + max(idx_numComp_sorted{2,1}); % Assign higher numbers to the wake bouts 
% idx_numComp_sorted{2,1}(isnan(idx_numComp_sorted{2,1})) = 0; % replace nan-values with zero 
% 
% tic
% for f = 1:max(fish_tags{1,1}) % For each fish 
%     
%     % Pre-allocate
%     % Data
%     threads{f,1,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
%         size(find(fish_tags{2,1} == f),1),1,'single'); % clusters    
%     threads{f,2,1} = nan(size(threads{f,1,1},1),2,'single'); % times (start & stop) 
%     threads{f,3,1} = nan(size(threads{f,1,1}),'single'); % time windows 
%     
%     % Control 
%     threads{f,1,2} = nan(size(threads{f,1,1}),'single'); % clusters    
%     threads{f,3,2} = nan(size(threads{f,1,1}),'single'); % time windows
%     
%     % Deterime starting state (a = wake, b = sleep) 
%     if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active 
%         a = 1; b = 2; 
%     else % If the fish starts inactive 
%         a = 2; b = 1; 
%     end
%     
%     % Fill in data 
%     % Clusters
%     threads{f,1,1}(a:2:end,1) = idx_numComp_sorted{1,1}...
%         (fish_tags{1,1} == f,1); % Fill in active clusters
%     threads{f,1,1}(b:2:end,1) = idx_numComp_sorted{2,1}...
%         (fish_tags{2,1} == f,1); % Fill in inactive clusters
%     % Times
%     threads{f,2,1}(a:2:end,1:2) = wake_cells...
%         (fish_tags{1,1} == f,1:2); % Fill in active times
%     threads{f,2,1}(b:2:end,1:2) = sleep_cells...
%         (fish_tags{2,1} == f,1:2); % Fill in inactive times
%     % Time Windows
%     threads{f,3,1}(a:2:end,1) = parameter_indicies{1,1}...
%         (fish_tags{1,1} == f,1); % Fill in active time windows 
%     threads{f,3,1}(b:2:end,1) = parameter_indicies{2,1}...
%         (fish_tags{2,1} == f,1); % Fill in inactive time windows 
%    
%     % Generate Control data 
%     % Clusters
%     clear scrap; scrap = idx_numComp_sorted{1,1}...
%         (fish_tags{1,1} == f,1); % Take active bout data 
%     threads{f,1,2}(a:2:end,1) = scrap(randperm(length(scrap))); % Mix  
%     clear scrap; scrap = idx_numComp_sorted{2,1}...
%         (fish_tags{2,1} == f,1); % Take inactive bout data 
%     threads{f,1,2}(b:2:end,1) = scrap(randperm(length(scrap))); % Mix 
%     % Time windows 
%     threads{f,3,2} = threads{f,3,1}; % Break up using the same time windows  
%     
%     % Report progress 
%     if mod(f,100) == 0 % every 100 fish 
%         disp(horzcat('Threaded Fish ',num2str(f),' of ',...
%             num2str(max(fish_tags{1,1})))); % Report progress
%     end
% end 
% toc 
% 
% clear f a b scrap  

%% Start Here 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180220.mat')

clear cells wake_cells_norm 

%% Filling in Data 
tic
% Calculate time off-sets between repeats of experiments 
c = nan(max(experiment_reps),1,'single'); % experiment groups x 1 
for er = 1:max(experiment_reps) % for each group of experiments
    lb_merge{er,1} = []; % structure {experiment groups}(light boundaries x experiments)
    for e = find(experiment_reps == er) % for each experiment in this group 
        lb_merge{er,1} = [lb_merge{er,1} lb{1,e}]; % merge light boundaries 
        % Note that this assumes that each experiment in each group has the
        % same number of light boundaries 
    end 
    
    [~,c(er)] = find(lb_merge{er,1} == max(lb_merge{er,1}(2,:))); 
    % Find the column (c) with the longest starting window 
    
    for e = 1:size(lb_merge{er,1},2) % For each experiment in this group
        offset(er,e) = lb_merge{er,1}(2,c(er)) - lb_merge{er,1}(2,e); % Calculate the offset
    end
end

% Pre-allocate data
for er = 1:max(experiment_reps) % for each group of experiments  
    states{er,1} = nan(size(find(i_experiment_reps == er),1),...
        max(lb_merge{er,1}(end,:) + ...
        offset(er,1:size(find(experiment_reps == er),2))),'single'); % {experiment, groups} - fish x time
    raw_data{er,1} = nan(size(states{er,1}),'single'); % {experiment, groups} - fish x time
end

% Fill data in
for er = 1:max(experiment_reps) % for each group of experiments
    counter = 1;  % start counter (fish in each group of experiments) 
    experiments = find(experiment_reps == er); % find experiments in this group  
    
    for e = 1:size(experiments,2) % for each experiment in this group
        counter_2 = 1; % start counter (fish in each experiment) 
        for f = find(i_experiment_tags == experiments(e))' % for each fish
            for b = 1:size(threads{f,1,1},1) % for each bout
                states{er,1}(counter,(threads{f,2,1}(b,1)+offset(er,e)):...
                    (threads{f,2,1}(b,2)+offset(er,e))) = ...
                    threads{f,1,1}(b,1); % Fill in cluster number
            end
            
            raw_data{er,1}(counter,(offset(er,e)+1):...
                (size(delta_px_sq{1,experiments(e)}.delta_px_sq,1)+offset(er,e))) =...
                delta_px_sq{1,experiments(e)}.delta_px_sq(:,counter_2); % fill in raw data 
            
            counter = counter + 1; % add to counter (fish in each group of experiments)
            counter_2 = counter_2 + 1; % add to counter (fish in each experiment)
        end
    end
end

% Keep aligned light boundaries 
for er = 1:max(experiment_reps) % for each group of experiments
    lb_merge{er,1} = lb_merge{er,1}(:,c(er));
end

clear c er e counter experiments f b counter_2 delta_px_sq
toc 

%% WT States Matrix 
    % Note (170930 - would be cool to split into sleep/wake subplots) 
    % Note (180116) - to show state matricies with zero's in (water top ups) 
    % You'll need to set zeros to be transparent and adjust the colormap 

figure; 
ax = imagesc(states{1,1}(:,find(sum(isnan(states{1,1}))==0,1,'first'):...
    find(sum(isnan(states{1,1}))==0,1,'last'))); % plot times where all experiments have data 
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap  
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32); 
set(ax,'CDataMapping','direct'); 
c = colorbar; c.Label.String = 'Cluster'; 
x = 0:(fps{1}*60*60*12):size((find(sum(isnan(states{1,1}))==0,1,'first'):...
find(sum(isnan(states{1,1}))==0,1,'last')),2); % set x ticks every 12 hours 
set(gca,'XTick',x); 
set(gca,'XTickLabel',{(0:size(x,2)-1)*12}); 
xlabel('Time (Hours)','Fontsize',32); 
ylabel('Fish ID','Fontsize',32);
 
clear ax c x 

%% Smooth Cluster Frequencies Over Time  

% New Version
all_states = double(max(idx_numComp_sorted{1,1})); % possible states
% Note that subplot can't take single values
ai_states(1:all_states) = 2; % specify if states are active (1) or inactive (2)
ai_states(min(idx_numComp_sorted{1,1}):end) = 1; % specify active states (1)
time_bins = fps{1}*60*5; % set smoothing window (note assumes a constant frame rate across experiments)

% Smooth Data
for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings
    
    for s = 1:2 % for active/inactive clusters
        counter_2 = 1; % start counter (counts clusters)
        
        for c = find(ai_states == s) % for each active or inactive cluster 
            for g = 1:max(i_group_tags(i_experiment_reps == er)) % For each group
                clear data;
                
%                 % Normalise by Std 
%                 data = states{er,1}(i_group_tags(i_experiment_reps == er) == g,...
%                     lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1))==c;
%                 if sum(sum(data,2)) ~= 0 % if fish use this cluster 
%                     smoothed_clusters{er,s}(counter_2,:,g) = ...
%                         smooth((nanmean(data)./nanstd(nanmean(data))),time_bins);
%                 else % if fish don't 
%                     smoothed_clusters{er,s}(counter_2,:,g) = zeros(1,size(data,2)); % Fill in Zeros 
%                 end
                
                % Rescale from 0-1 
                data = states{er,1}(i_group_tags(i_experiment_reps == er) == g,...
                    lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1))==c;
                
                if sum(sum(data,2)) ~= 0 % if fish use this cluster
                    % Re-Scale to 0-1
                    data = smooth(nanmean(data),time_bins)'; % smooth 
                    smoothed_clusters{er,s}(counter_2,:,g) = ...
                        (data - nanmin(data(1,(1+time_bins):(end-time_bins))))/...
                        range(data(1,(1+time_bins):(end-time_bins))); % scale 
                else % if fish don't
                    smoothed_clusters{er,s}(counter_2,:,g) = zeros(1,size(data,2)); % Fill in Zeros
                end
                
            end
            counter_2 = counter_2 + 1; % add to counter
        end
        
    end
end

% Calculate y-axis offset 
for er = 1:max(experiment_reps) % for each group of experiments
    for s = 1:2 % for active & inactive clusters 
        smoothed_clusters_offset(er,s) = max(max(max(smoothed_clusters{er,s}(:,(1+time_bins):(end-time_bins),:)))); 
        smoothed_clusters_offset(er,s) = smoothed_clusters_offset(er,s) + ...
            (smoothed_clusters_offset(er,s)*0.05); 
    end
end

clear er set_token s counter_2 c g data 

%% Cluster Frequencies Figure (Run one @ A Time) 

for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings
    
    figure; % make a figure for this experiment
    for s = 1:2 % for active/inactive
        subplot(1,2,s); hold on; clear scrap; set(gca,'FontName','Calibri');
        title(horzcat(strings{s},' Clusters'));
        for c = 1:size(smoothed_clusters{er,s},1) % for each cluster
            for g = 1:size(smoothed_clusters{er,s},3) % For each group
                
                if er == 1 % for the WT Data
                    legend_lines(1,g) = plot(lb_merge{er,1}(time_window{set_token}(1)):...
                        lb_merge{er,1}(time_window{set_token}(2)+1),...
                        (smoothed_clusters{er,s}(c,:,g) + (smoothed_clusters_offset(er,s)*(c-1))),...
                        'color',cmap_cluster{s}(c,:),'linewidth',5);
                else
                    legend_lines(1,g) = plot(lb_merge{er,1}(time_window{set_token}(1)):...
                        lb_merge{er,1}(time_window{set_token}(2)+1),...
                        (smoothed_clusters{er,s}(c,:,g) + (smoothed_clusters_offset(er,s)*(c-1))),...
                        'color',cmap{set_token}(g,:),'linewidth',5);
                end     
            end
            
            plot([lb_merge{er,1}(time_window{set_token}(1)),lb_merge{er,1}(time_window{set_token}(2)+1)],...
                [(smoothed_clusters_offset(er,s)*(c-1)),(smoothed_clusters_offset(er,s)*(c-1))],'color','k','linewidth',1);
            
        end
        
        % Find the top & bottom
        scrap(1,1) = smoothed_clusters_offset(er,s)*size(smoothed_clusters{er,s},1);
        scrap(2,1) = 0;
        
        % Night patches
        a = 1;
        for n = 1:size(nights{set_token},2) % For each night
            r(a) = rectangle('Position',[lb_merge{er,1}(nights_crop{set_token}(nights{set_token}(n))) ...
                min(scrap(2,:)) - (min(scrap(2,:))*0.05)...
                (lb_merge{er,1}(nights_crop{set_token}(nights{set_token}(n))+1)-1) - ...
                lb_merge{er,1}(nights_crop{set_token}(nights{set_token}(n)))...
                max(scrap(1,:)) + (max(scrap(1,:))*0.05)],...
                'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
            uistack(r(a),'bottom'); % Send to back
            a = a + 1;
        end
        
        % Axis etc
        axis([(lb_merge{er,1}(time_window{set_token}(1)) + (1 + time_bins))...
            (lb_merge{er,1}(time_window{set_token}(2)+1) - time_bins) ...
            scrap(2,1) scrap(1,1)]);
        box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
        xlabel('Time (Days/Nights)','Fontsize',32);
        set(gca, 'XTick', []);
        set(gca,'YTick',[]);
        ylabel('Mean Frequency','Fontsize',32);
        
    end
end

clear er set_token s c g legend_lines scrap a n r 

%% Bout Probabilities Order
% Sort based on WT mean day probability 
er = 1; 
set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings

[~,bp_order(1,:)] = sort(nanmean(nanmean(bout_proportions{1,1}...
    (i_experiment_reps == er,:,days_crop{set_token}(days{set_token})),3)),'descend'); % active clusters
bp_order(2,1:numComp(2)) = 1:numComp(2); % inactive clusters (keep sorted by length)
bp_order(bp_order == 0) = NaN; % remove zeros

clear er set_token 

%% Bout Probabilities Histograms WT - Scatter
% Plots Individual Fish
% Plots Day vs Night on the Same Axis

for er = 1 % for the WT experiments
    figure;
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings
    for s = 1:2 % for active/inactive clusters
        ax = subplot(2,1,s); hold on; set(gca,'FontName','Calibri');
        col = 1; % start colour counter
        sep = 0.75/2; clear spread_cols; 
        for g = 1:max(i_group_tags(i_experiment_reps == er)) % For each group
            
            % Plot Spread
            % Day
            spread_cols(1,:) = plotSpread(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                days_crop{set_token}(days{set_token})),3),...
                'distributionColors',cmap_2{set_token}(col,:),'spreadWidth',sep,'showMM',2,'xValues',1:numComp(s));
            spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change Mean properties
            spread_cols{2}.MarkerSize = 12;
            
            % Night
            spread_cols(2,:) = plotSpread(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                nights_crop{set_token}(nights{set_token})),3),...
                'distributionColors',cmap_2{set_token}(col+1,:),'spreadWidth',sep,'showMM',2,'xValues',(1:numComp(s)) + sep);
            spread_cols{2,2}.LineWidth = 3; spread_cols{2,2}.Color = [1 0.5 0]; % Change Mean properties
            spread_cols{2,2}.MarkerSize = 12;
            
        end
        
        % Axis etc
        set(findall(ax,'type','line'),'markersize',15); % change marker sizes
        axis tight
        box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
        xlabel(horzcat(strings{s},' Clusters'),'Fontsize',32);
        set(gca, 'XTick', (1:numComp(s))+sep/2);
        set(gca,'XTickLabels',bp_order(s,:),'Fontsize',32);
        ylabel('Probability','Fontsize',32);
        
        if er == 1 && s == 1
            scrap = get(gca,'Children');
            [~,icons,plots,~] = legend([scrap(end) scrap(2)],'Day','Night');
            legend('boxoff'); set(plots(1:2,1),'MarkerSize',15); set(icons(1:2),'Fontsize',32) ;
        end
    end
    
end

clear er set_token s ax col sep g spread_cols scrap icons plots

%% Bout Probabilities Histograms - Conditions 
    % Plots Day & Night on seperate axes 
    % Plots a mean + SEM 
    % Sorted based on the WT sorting above 
    
for er = 2:max(experiment_reps) % for each group of experiments
    figure; clear legend_lines
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings
    for s = 1:2 % for active/inactive clusters
        subplot(2,2,s); hold on; set(gca,'FontName','Calibri'); clear l_top; 
        for g = 1:max(i_group_tags(i_experiment_reps == er)) % For each group
            % Means
            legend_lines(1,g) = plot(nanmean(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                days_crop{set_token}(days{set_token})),3))','color',...
                cmap{set_token}(g,:),'linewidth',5);
            % Std 
            clear l_mean l_std; 
            l_mean = nanmean(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                days_crop{set_token}(days{set_token})),3)); 
            l_std = nanstd(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                days_crop{set_token}(days{set_token})),3)); 
            plot([1:numComp(s) ; 1:numComp(s)],...
                [(l_mean - l_std) ; l_mean + l_std],'color',...
                cmap{set_token}(g,:),'linewidth',2.5);
            % Store Top 
            l_top(g) = max(l_mean + l_std) + (max(l_mean + l_std)*0.05); 
        end
        
        % Axis etc
        axis([0.5 numComp(s)+0.5 0 max(l_top)]);
        box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
        set(gca,'XTick', 1:numComp(s)); 
        set(gca,'XTickLabels',bp_order(s,:),'Fontsize',32);
        if s == 1
            ylabel('Probability','Fontsize',32);
        end
        
        % Legend
        if s == 2
            [~,icons,plots,~] = legend(legend_lines,geno_list{set_token}.colheaders,...
                'location','northeast');
            legend('boxoff'); set(icons(1:size(legend_lines,1)),'Fontsize',32) ; set(plots,'LineWidth',5);
        end
        
        subplot(2,2,s+2); hold on; set(gca,'FontName','Calibri'); clear l_top; 
        for g = 1:max(i_group_tags(i_experiment_reps == er)) % For each group
            % Means
            plot(nanmean(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                nights_crop{set_token}(nights{set_token})),3))','color',...
                cmap{set_token}(g,:),'linewidth',5);
            % Std
            clear l_mean l_std;
            l_mean = nanmean(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                nights_crop{set_token}(nights{set_token})),3));
            l_std = nanstd(nanmean(bout_proportions{s,1}...
                (i_experiment_reps == er & i_group_tags == g,bp_order(s,1:numComp(s)),...
                nights_crop{set_token}(nights{set_token})),3));
            plot([1:numComp(s) ; 1:numComp(s)],...
                [(l_mean - l_std) ; l_mean + l_std],'color',...
                cmap{set_token}(g,:),'linewidth',2.5);
            % Store Top
            l_top(g) = max(l_mean + l_std) + (max(l_mean + l_std)*0.05);
        end
        
        % Night patch
        axis([0.5 numComp(s)+0.5 0 max(l_top)]);
        a = 1;
        r(a) = rectangle('Position',[0.5 0 (numComp(s)+0.5) max(l_top)],...
            'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        
        % Axis etc
        box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
        xlabel(horzcat(strings{s},' Clusters'),'Fontsize',32);
        set(gca, 'XTick', 1:numComp(s));
        set(gca,'XTickLabels',bp_order(s,:),'Fontsize',32);
        if s == 1
            ylabel('Probability','Fontsize',32);
        end
    end
end

clear er set_token s g legend_lines l_mean l_std l_top icons plots a r 

%% Bouts in PCA Space 

% Load Active GMModel 
load(horzcat('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\New\Z_Score_Wake\1s_2000000d_1r_',...
    num2str(numComp(1)),'k.mat'),'GMModels'); 

% Sample Equally from each cluster 
sample_size = 1000; % hard coded sample size 
sample = []; sample_tags = []; % strcutures
for c = min(idx_numComp_sorted{1,1}):max(idx_numComp_sorted{1,1}) % for each active cluster 
    sample = [sample ; datasample(score(idx_numComp_sorted{1,1}==c,1:2),sample_size,...
        1,'replace',false,'weights',P{1,1}(idx_numComp_sorted{1,1}==c))]; % sample bouts  
    sample_tags = [sample_tags ; repmat(c,[sample_size,1])]; % store cluster tags 
end 

% Calculate pdf 
density = pdf(GMModels{1},sample); 

% 3D Scatter Plot 
figure; cols = 1; % counts colours 
for c = min(idx_numComp_sorted{1,1}):max(idx_numComp_sorted{1,1}) % for each active cluster 
    scatter3(sample(sample_tags==c,1),sample(sample_tags==c,2),...
        density(sample_tags==c,1),90,...
        'markerfacecolor',cmap_cluster{1,1}(cols,:),...
        'markeredgecolor','k','linewidth',0.01);
    cols = cols + 1; % add to colour counter 
    hold on; 
end

% Figure Settings 
axis tight; box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'ZTickLabels',[]); grid off;  
xlabel('PC 1','Fontsize',32); ylabel('PC 2','Fontsize',32); zlabel('Density','Fontsize',32); 

% 2D Scatter Plot 
figure; cols = 1;
for c = min(idx_numComp_sorted{1,1}):max(idx_numComp_sorted{1,1}) % for each active cluster 
    scatter(sample(sample_tags==c,1),sample(sample_tags==c,2),18,...
        'markerfacecolor',cmap_cluster{1,1}(cols,:),...
        'markeredgecolor','k');
    cols = cols + 1;
    hold on; 
end

clear c cols 

%% Bout Shapes Data 
 
% Cumulative fish tags 
for er = 1:max(experiment_reps) % for each group of experiments 
    fish_tags_cm(er,1) = size(raw_data{er,1},1); % fill in number of fish 
end
fish_tags_cm = cumsum(fish_tags_cm); % cumulative number of fish across experiment groups 

% Allocate experiment reps to every active bout
% experiment_reps_long{1,1} = experiment_tags{1,1}; % temp assign
% for er = 1:max(experiment_reps) % for each repeat
%     found = find(experiment_reps == er); % find experiments
%     
%     for f = found % for each experiment in the repeat
%         experiment_reps_long{1,1}(experiment_reps_long{1,1} == f,1) = er; % tag with grouping variable
%     end
%     
% end

clear er  

%% Sampling from a single fish 
sample_size = 100; % hard coded sample size 

% Sample a fish with enough bouts in each cluster 
a = 0; 
while a == 0
    f = datasample(find(i_experiment_reps == 1),1); % sample a WT fish
    for c = find(ai_states == 1) % for each active cluster
        if sum(fish_tags{1,1} == f & idx_numComp_sorted{1,1} == c) < sample_size
            % if there are too few of this bout type 
           continue % select another fish 
        end 
    end
    a = 1; 
end

counter = 1; % counts active clusters 
for c = find(ai_states == 1) % for each active bout type
    
    clear sample 
    
    % Weighted Sample bouts
    [sample,~] = datasample(...
        wake_cells(fish_tags{1,1}==f & idx_numComp_sorted{1,1}==c,:),...
        sample_size,1,'replace',false,'weights',P{1,1}(fish_tags{1,1}==f & idx_numComp_sorted{1,1}==c,:));
    
    % Allocate Space
    bouts{1,counter} = zeros(sample_size,max(sample(:,3)),'single'); % number of bouts x longest (in frames)
    bouts{2,counter} = bouts{1,counter}; % states (just to check you get the correct raw data)  
   
    % Fill in data
    for b = 1:sample_size % for each bout
        
        % Extract raw data
        bouts{1,counter}(b,1:sample(b,3)) = raw_data{i_experiment_reps(f),1}(f,...
            (sample(b,1)+offset(i_experiment_reps(f),i_experiment_tags(f))):...
            (sample(b,2)+offset(i_experiment_reps(f),i_experiment_tags(f))));
        bouts{2,counter}(b,1:sample(b,3)) = states{i_experiment_reps(f),1}(f,...
            (sample(b,1)+offset(i_experiment_reps(f),i_experiment_tags(f))):...
            (sample(b,2)+offset(i_experiment_reps(f),i_experiment_tags(f))));
    end
    
    counter = counter + 1; % counts active clusters 
end

clear a f c counter sample b 

%% Sampling From Multiple Fish 
% sample_size = 100; % hard coded sample size 
% counter = 1; % counts active clusters 
% for c = find(ai_states == 1) % for each active bout type
%     
%     clear sample sample_tags sample_er sample_fish sample_et 
%     
%     % Weighted Sample bouts
%     [sample,sample_tags] = datasample(...
%         wake_cells(idx_numComp_sorted{1,1}==c,:),...
%         sample_size,1,'replace',false,'weights',P{1,1}(idx_numComp_sorted{1,1}==c,:));
%     sample_tags = sample_tags'; % flip 
%     
%     % Allocate Space
%     bouts{1,counter} = zeros(sample_size,max(sample(:,3)),'single'); % number of bouts x longest (in frames)
%     bouts{2,counter} = bouts{1,counter}; 
%     
%     % Sample tags
%     sample_er = experiment_reps_long{1,1}(idx_numComp_sorted{1,1}==c,:); % experiment groups 
%     sample_fish = fish_tags{1,1}(idx_numComp_sorted{1,1}==c,:); % fish i.d.
%     sample_et = experiment_tags{1,1}(idx_numComp_sorted{1,1}==c,:); % experiment tags 
%     
%     % Fill in data
%     for b = 1:sample_size % for each bout
%         
%         % determine bout indexing numbers
%         er = sample_er(sample_tags(b,1)); % sample experiment groups 
%         if er == 1 % for the first experiment group
%             fish = sample_fish(sample_tags(b,1)); % sample fish tags 
%         else % for other experiment groups 
%             fish = sample_fish(sample_tags(b,1)) - fish_tags_cm(er - 1); 
%             % determine correct fish id
%         end
%         et = find(find(experiment_reps == er) == sample_et(sample_tags(b,1))); % sample experiment tags 
%         
%         % Extract raw data
%         bouts{1,counter}(b,1:sample(b,3)) = raw_data{er,1}(fish,...
%             (sample(b,1)+offset(er,et)):(sample(b,2)+offset(er,et)));
%         bouts{2,counter}(b,1:sample(b,3)) = states{er,1}(fish,...
%             (sample(b,1)+offset(er,et)):(sample(b,2)+offset(er,et)));
%     end
%     
%     counter = counter + 1; % counts active clusters 
% end

% clear er f number counter c sample sample_tags sample_er ...
%     sample_fish sample_et b fish et found 

%% Bout Shapes Figure

% Plotting variables 
for b = 1:size(bouts,2) % for each cluster 
    bs_top(b) = max(bouts{1,b}(:)); % find the max value 
    bs_l(b) = size(bouts{1,b},2); % find the max length 
end

bs_top = max(bs_top) + (max(bs_top)*0.05); % add a bit of space
bs_l = max(bs_l); % max length 

figure; hold on; set(gca,'FontName','Calibri');
for b = 1:size(bouts,2) % for each active bout type
    
    % Padding
    plot([zeros(size(bouts{1,b},1),1) bouts{1,b} zeros(size(bouts{1,b},1),1)]' + bs_top*(b-1),...
        'color',cmap_cluster{1,1}(b,:)+(1-cmap_cluster{1,1}(b,:))*(1-(1/(1)^.5)));
    scrap = nanmean(bouts{1,b}); scrap(scrap < 1) = [];
    plot([0 scrap 0]' + bs_top*(b-1),...
        'k','linewidth',3);
    
    %     % No padding
    %     plot((bouts{1,b}' + bs_top*(b-1)),'color',...
    %         cmap_cluster{1,1}(b,:)+(1-cmap_cluster{1,1}(b,:))*(1-(1/(1)^.5)));
    
    plot([1 bs_l + 2],[bs_top*(b-1) bs_top*(b-1)],'k','linewidth',1);
    
end

box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
axis([1 bs_l 0 bs_top*(size(bouts,2))]); % hard coded
set(gca,'XTick',2:2:bs_l); 
set(gca,'XTickLabels',{(round((1:2:bs_l)/fps{1},2,'decimals'))}); % hard coded
xlabel('Time (seconds)','Fontsize',32);
set(gca,'YTick',[]); 
ylabel('Delta Px','Fontsize',32);

clear b bs_top bs_l b scrap 

%% Hierachical Sequences (Local Version) 
    %https://github.com/aexbrown/Behavioural_Syntax
    %http://www.sequitur.info/

% Constructing the grammar library 

% % Settings 
% sMax = max(idx_numComp_sorted{1,1}) + 1; % Maximum states + 1 (first new symbol) 
% nMax = 10; % Maximum n-grams to consider 
% 
% % Allocate 
% grammar = cell(size(threads,1),2); % fish x test/control 
% compVec = cell(size(threads,1),2); % fish x test/control  
% totSavings = zeros(size(threads,1),2); % fish x test/control 
% gTermCell = cell(size(threads,1),2); % cell for storing n-grams - fish x test/control 
% uniqueSeqs = cell(1,2); % 1 x test/control 
% 
% % Generate a Grammer for each fish + Controls 
% % 16hours for 124 fish - 4days,3nights 
% tic
% parfor f = 1:find(i_experiment_reps == 1,1,'last') % For each fish size(threads,1)
%     
%     for tc = 1:2 % for real vs shuffled data
%         % Compress
%         [grammar{f,tc}, compVec{f,tc}, totSavings(f,tc)] = ...
%             compressSequenceNFast(threads{f,1,tc}',sMax,nMax);
%         
%         % Get terminal symbols
%         gTermCell{f,tc} = getGrammarTerminals(grammar{f,tc});
%     end
%     
%     % Report progress
%     disp(horzcat('Compressed Fish ',num2str(f)));
% end
% disp('Finished Compressing Fish'); 
% toc

%% Load in Legion Data

load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180223.mat'); 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\Compression_Results_Final.mat',...
    'gTermCell','totSavings'); 

%% Sequence Lengths 

% Find Sequence Lengths 
for tc = 1:size(gTermCell,2) % for real & shuffled data
    for f = 1:size(gTermCell,1) % for each fish
        seq_lengths{f,tc} = zeros(size(gTermCell{f,tc},1),1,'single'); % fish x shuffles {sequences x 1} 
        grammar_size(f,tc) = size(gTermCell{f,tc},1); % fish x shuffles 
        
        for s = 1:size(gTermCell{f,tc},1) % for each sequence
            seq_lengths{f,tc}(s,1) = size(gTermCell{f,tc}{s,1},2); % find it's length
        end
        
        sq_l(f,tc,:) = minmax(seq_lengths{f,tc}')'; % fish x shuffles x min/max 
        
    end
end
    
% Fit pdfs  
for tc = 1:size(gTermCell,2) % for real vs shuffled data
    for f = 1:size(gTermCell,1) % for each fish
        pd = fitdist(seq_lengths{f,tc},'kernel','Width',1); % Fit pdf
        seq_lengths_pd(f,:,tc) = pdf(pd,min(sq_l(:)):max(sq_l(:))); % fish x pdf x real & shuffles
    end
end

clear tc f s sq_l pd 

%% Sequence Lengths Figure
    
clear legend_lines; 
figure; box off; set(gca, 'Layer','top'); hold on; 
set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format
for tc = 2:size(seq_lengths_pd,3) % for each shuffle 
    legend_lines(2) = plot(2:(size(seq_lengths_pd,2)+1),...
        nanmean(seq_lengths_pd(:,:,tc)),...
        'color','k','linewidth',5); % plot an average  
end
legend_lines(1) = plot(2:(size(seq_lengths_pd,2)+1),...
    nanmean(seq_lengths_pd(:,:,1)),...
    'color',cmap{1,1},'linewidth',5); % plot the real data

% Settings 
[~,icons,plots,~] = legend(legend_lines,'Real Data','Shuffled Data',...
    'location','northeast');
legend('boxoff'); set(icons(1:size(legend_lines,1)),'Fontsize',32) ; set(plots,'LineWidth',5);
xlabel('Motif Length','Fontsize',32); ylabel('Probability','Fontsize',32); % Y Labels
axis([2 (size(seq_lengths_pd,2)+1) ylim]); 

%% Constructing a common Grammar 
    % Roughly 1 minute per shuffle (629 fish)
    
tic
% Merge all the grammers and keep only unique elements
for tc = 1:size(gTermCell,2) % for real & shuffled data
    uniqueSeqs{1,tc} = []; % structure
    
    for f = 1:size(gTermCell,1) % for each fish
        % append grammar terminals to the total set
        uniqueSeqs{1,tc} = vertcat(uniqueSeqs{1,tc}, gTermCell{f,tc});
    end
    
    % Determine unique sequences from this (larger) set
    [uniqueSeqs{1,tc}, ~] = countUniqueRows(uniqueSeqs{1,tc});
    
    % Report Progress
    disp(horzcat('Merged Grammar ',num2str(tc)));
    
end
disp('Merged all grammers');
toc

clear tc f gTermCell

%% Remove Sequences with NaN (tagged as zeros) 

for tc = 1:size(uniqueSeqs,2) % for real vs shuffled data 
    nan_locs = cellfun(@(s) ismember(0, s), uniqueSeqs{1,tc}); % find sequences with nans 
    uniqueSeqs{1,tc}(nan_locs,:) = []; % remove these 
    clear nan_locs 
end 

clear tc nan_locs 

%% Construct Real Grammar Matrix  

grammar_mat{1,1} = nan(size(uniqueSeqs{1,1},1),size(seq_lengths_pd,2)+1,'single'); % sequences x max length 
    
for s = 1:size(uniqueSeqs{1,1},1) % for each sequence
    grammar_mat{1,1}(s,1:size(uniqueSeqs{1,1}{s,1},2)) = uniqueSeqs{1,1}{s,1}; % fill in sequence
end

clear s 

%% Grammar Matrix Figure 
    % Note: 180226 - Maybe the most interesting figure is actually a
    % grammar with only the sequences above chance? 
    
figure; 
grammar_mat_sorted = flip(sortrows(grammar_mat{1,1})); % sort rows of grammar_mat  
ax = imagesc(grammar_mat_sorted,'AlphaData',isnan(grammar_mat_sorted)==0); % imagesc with nan values in white 
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap  
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(ax,'CDataMapping','direct'); 
c = colorbar; c.Label.String = 'Cluster'; 
xlabel('Position in Motif','Fontsize',32); 
ylabel('Motif','Fontsize',32);

clear ax grammar_mat_sorted c 

%% Compressibility 
% The compressibility of a sequence of uncompressed length l is given by the sum of the savings S 
% at each iteration divided by l (Gomez-Marin et al.,2016) 

% compressibility = zeros(size(totSavings),'single'); % fish x t/c
% compressibility_rel = zeros(size(totSavings,1),1,'single'); % fish x 1 
% 
% for f = 1:size(threads,1) % for each fish 
%     compressibility(f,:) = totSavings(f,:)/size(threads{f,1,1},1); % calculate compressibility 
%     compressibility_rel(f,1) = (compressibility(f,1) - nanmean(compressibility(f,2:end)))/...
%         nanstd(compressibility(f,2:end)); % relative compressibility (Z-Score) 
% end 
% 
% clear f 

%% Compressibility Figure 

% figure;
% for er = 1:max(experiment_reps) % for each group of experiments
%     set_token = find(experiment_reps == er,1,'first'); % settings
%     subplot(2,3,er); counter = 1; % counts groups for plots
%     hold on; set(gca,'FontName','Calibri'); clear scrap;
%     
%     for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
%         clear data; 
%         data = [repmat(compressibility(i_experiment_reps == er & i_group_tags == g,1),(size(compressibility,2)-1),1) ...
%             reshape(compressibility(i_experiment_reps == er & i_group_tags == g,2:end),[],1)]; 
%         plot([counter,counter+1],data,...
%             'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
%         errorbar([counter,counter+1],nanmean(data),nanstd(data),...
%             'color',cmap{set_token}(g,:),'linewidth',3);
%         counter = counter + 2; % add to counter
%         
%         scrap(1,g) = min(min(compressibility(i_experiment_reps == er & i_group_tags == g,:))); 
%         scrap(2,g) = max(max(compressibility(i_experiment_reps == er & i_group_tags == g,:))); 
%     end
%     
%     box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
%     if er == 1 % for the WT Data 
%         set(gca, 'XTick', [1 2]); % set X-ticks
%         set(gca,'XTickLabels',{'Data','Shuffled'}); % X Labels
%     else % for the other experiments 
%         set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
%     end
%     ylabel('Compressibility','Fontsize',32); % Y Labels
%     axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 ...
%         (min(scrap(1,:)) - (min(scrap(1,:))*0.05)) (max(scrap(2,:)) + (max(scrap(2,:))*0.05))]);
% end
% 
% clear er set_token g scrap counter data 

%% Relative Compressibility Figure 
    
% figure;
% for er = 1:max(experiment_reps) % for each group of experiments
%     set_token = find(experiment_reps == er,1,'first'); % settings
%     ax = subplot(2,3,er);
%     hold on; set(gca,'FontName','Calibri'); clear scrap;
%     
%     for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
%         clear data;
%         data = compressibility_rel(i_experiment_reps == er & i_group_tags == g,:);
%         % To use a mean for each fish instead use
%         %data = nanmean(compressibility_rel(i_experiment_reps == er & i_group_tags == g,:),2);
%         
%         spread_cols = plotSpread((data(:)')',...
%             'distributionColors',cmap{set_token}(g,:),'showMM',2,'XValue',g);
%         set(findall(ax,'type','line'),'markersize',15); % change marker sizes
%         spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = 'k'; % Change Mean properties
%         spread_cols{2}.MarkerSize = 12;
%         
%         scrap(1,g) = min(min(data));
%         scrap(2,g) = max(max(data));
%     end
%     
%     box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
%     if er == 1 % for the WT Experiments 
%         set(gca, 'XTick', [1 2]); % set X-ticks
%         set(gca,'XTickLabels',{'Data','Shuffled'}); % X Labels
%     else % for the rest 
%         set(gca,'XTick',1:max(i_group_tags(i_experiment_reps == er)));
%         set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
%     end
%     ylabel({'Relative' ; 'Compressibility' ; '(Z-Score)'},'Fontsize',32); % Y Labels
%     axis([0.5 (max(i_group_tags(i_experiment_reps == er))+.5) ...
%         min(scrap(1,:)) max(scrap(2,:))]);
% end
% 
% clear ax er set_token g scrap data

%% Compressibility N-Way Anova

% for er = 1:max(i_experiment_reps) % for each group of experiments
%     
%     clear anova_group anova_tc anova_experiment data
%     
%     % Grouping Variables
%     anova_group = repmat(i_group_tags(i_experiment_reps==er),...
%         [size(threads,3),1])'; % groups
%     anova_tc = ones(size(anova_group)); % real (0) vs shuffled data (1) 
%     anova_tc(1:size(i_group_tags(i_experiment_reps==er),1)) = 0;
%     anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
%         [size(threads,3),1])'; % experiments
%     
%     % Data to Compare 
%     if er == 1 % WT - compare real vs shuffled data 
%         data = compressibility(i_experiment_reps==er,:); % grab data
%         data = data(:)'; % vectorise
%     else % compare relative compressibility 
%         data = compressibility_rel(i_experiment_reps == er,:); 
%         data = data(:)'; % vectorise 
%         anova_group = anova_group(1:size(data,2)); % trim out real data indicies  
%         anova_experiment = anova_experiment(1:size(data,2)); % trim 
%         anova_tc = anova_tc(1:size(data,2)); % trim 
%     end
%     
%     % Comparison 
%     [twa.compress.p{1,er},~,twa.compress.stats{1,er}] = anovan(data,...
%         {anova_group,anova_tc,anova_experiment},...
%         'display','off','model','full');
%     
% end
% 
% clear er anova_group anova_tc anova_experiment data

%% Grammar in Time 

% Convert real unique seqs to singles before passing to Legion 
for s = 1:size(uniqueSeqs{1,1},1) % for each sequence
    uniqueSeqs{1,1}{s,1} = single(uniqueSeqs{1,1}{s,1}); % convert to single
end

%% Local Version (Parallel - 75mins for 125 fish, 4days,3nights) 

% gCount = cell(size(threads,1),2); % counts - fish x t/c 
% gFreq = cell(size(threads,1),2); % frequency - fish x t/c 
% for tc = 1:2 % for real/control data
%     for f = 1:size(threads,1) % for each fish     
%         gCount{f,tc} = zeros(size(uniqueSeqs{1,tc},1),max(parameter_indicies{1,1}),'single'); % {f,t/c} - uniqueSeqs x time windows 
%         gFreq{f,tc} = zeros(size(uniqueSeqs{1,tc},1),max(parameter_indicies{1,1}),'single'); % {f,t/c} - uniqueSeqs x time windows 
%     end
% end
% t_one = min(parameter_indicies{1,1}); % first time window 
% t_two = max(parameter_indicies{1,1})+1; % 2nd time window + 1 
% 
% tic
% parfor f = 1:size(threads,1) % for each fish
%     for tc = 1:2 % for real vs shuffled data
%         for i = 1:size(uniqueSeqs{1,tc},1) % For each sequence
%             % Find each sequence and count it's time windows
%             % Note that strfind(str,pattern) outputs the starting index of each
%             % occurrence of pattern in str. This is then used to index the
%             % time windows of these occurances, and then fed to histcounts
%             % which counts the occurances in each time window
%             gCount{f,tc}(i,:) = histcounts(threads{f,3,tc}...
%                 (strfind(threads{f,1,tc}',uniqueSeqs{1,tc}{i,1})),...
%                 t_one:t_two);
%             % Calculate frequency
%             if sum(gCount{f,tc}(i,:),2) > 0 % if fish (f) uses pattern (i)
%                 gFreq{f,tc}(i,:) = gCount{f,tc}(i,:)./sum(gCount{f,tc}(i,:));
%                 % calculate it's frequency in each time window
%             end
%         end
%     end
%     disp(horzcat('Finished Grammar in Time for fish ',num2str(f)));
% end
% disp('Overall Time Taken = '); 
% toc 
% 
% clear tc f t_one t_two i c 

%% Load & Reformat Legion Data  
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\Grammar_Results_Final.mat',...
    'gCount_norm');

temp = gCount_norm; clear gCount_norm; 

% Allocate 
for tc = 1:size(temp,2) % for real & shuffled data
    gCount_norm{1,tc} = zeros(size(temp{1,1},1),size(temp{1,1},2),size(temp,1),...
        'single'); % {real/shuffled} sequences x time windows x fish
end

% calculate inf replacements 
scrap = zeros(1,size(temp,2),'single'); % 1 x real & shuffled data  
scrap(1,1) = 1; 
inf_r = (scrap(1,1) - nanmean(scrap(1,2:end)))/nanstd(scrap); 

% Fill data 
for tc = 1:size(temp,2) % for real & shuffled data
    for f = 1:size(temp,1) % for each fish
        gCount_norm{1,tc}(:,:,f) = temp{f,tc}; % fill data
    end
    
    % Replace -Inf & Inf Values 
    gCount_norm{1,tc}(isinf(gCount_norm{1,tc}) & ...
        gCount_norm{1,tc} < 0) = -1*inf_r; % -ve inf 
    gCount_norm{1,tc}(isinf(gCount_norm{1,tc}) & ...
        gCount_norm{1,tc} > 0) = inf_r; % +ve inf 
      
end

clear temp tc scrap inf_r f

%% Real vs Shuffled Z-Scores Pdf (Rounded) 
    % Note: 180313 - may be interesting to compare Kurtosis? 
    % Note: 180406 - Sum from lower values onwards 
    
for tc = 1:size(gCount_norm,2) % for each shuffle
    tb(:,tc) = minmax(round(gCount_norm{1,tc}(:)')); % find it's max & min z-score
end
tb = min(tb(:)):max(tb(:)); % vector from min-max z-score
tb_z = find(tb == 0); % zero location

tc_pdf = NaN(size(gCount_norm{1,1},3),size(tb,2),...
    size(gCount_norm,2),'single'); % fish x z-score range x real/shuffled data
tc_pdf_binned = nan(size(gCount_norm{1,1},3),21,size(gCount_norm,2),'single');
% fish x hard coded bin x real/shuffled data

tic
for f = 1:size(gCount_norm{1,1},3) % for each fish
    for tc = 1:size(gCount_norm,2) % for each shuffle
        clear data pd; 
        data = round(gCount_norm{1,tc}(:,:,f)); data = data(:);
        pd = fitdist(data,'kernel','Width',1); % Fit pdf
        tc_pdf(f,:,tc) = pdf(pd,tb(1):tb(end)); % all data
        
        % Binned data
        tc_pdf_binned(f,2:(end-1),tc) = tc_pdf(f,(tb_z-9):(tb_z+9),tc); % around zero 
        tc_pdf_binned(f,1,tc) = sum(tc_pdf(f,1:(tb_z-10),tc)); % binned 
        tc_pdf_binned(f,end,tc) = sum(tc_pdf(f,(tb_z+10):end,tc)); % binned
    end
    
    % Report Progress
    if mod(f,25) == 0
        disp(horzcat(num2str(f),' of ',num2str(size(gCount_norm{1,1},3))));
    end
    
end
toc

clear tc f data pd

%% Calculate Skewness 
skewness(tc_pdf(f,:,tc)'); 

%% Load Data
gCount_norm(:,2:end) = []; % remove excess shuffled data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180227.mat'); 

%% Real vs Shuffled Z-Scores Figure 
    % Note: 180406 - Sum from lower values onwards 

er = 1; 
set_token =  find(experiment_reps == er,1,'first'); % settings
figure; 
hold on; clear legned_lines
plot(-10:10,tc_pdf_binned(i_experiment_reps == er,:,1),'color',...
    cmap_2{set_token}(1,:)+(1-cmap_2{set_token}(1,:))*(1-(1/(5)^.5)),...
    'linewidth',1.5)
for tc = 2:11
    plot(-10:10,tc_pdf_binned(i_experiment_reps == er,:,tc),'color',...
    ([1 1 1]*(1-(1/(5)^.5))),'linewidth',1.5);
end
legend_lines(2) = plot(-10:10,nanmean(nanmean(tc_pdf_binned(i_experiment_reps == er,:,2:end)),3),...
    'k','linewidth',3); 
legend_lines(1) = plot(-10:10,nanmean(tc_pdf_binned(i_experiment_reps == er,:,1)),...
    'color',cmap_2{set_token}(1,:),'linewidth',3); 

axis tight 
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(gca,'XTick',-10:2:10)
xlabel('Z-Score'); ylabel('Probability'); 
legend(legend_lines,'Real Data','Shuffled Data'); 
legend('boxoff'); 

% Note: Remember to add in greater & less than symols to the ends 

%% "Common-ness" of Grammar
    % Note: 180228 this may be more interesting for only the sequences that
    % occur above chance 
    
% Logical Table 
% common_table = zeros(size(uniqueSeqs{1,1},1),size(gFreq_merge{1,1},3),'single'); % sequences x fish  
% 
% for s = 1:size(common_table,1) % for each sequence 
%     common_table(s,:) = sum(squeeze(gFreq_merge{1,1}(s,:,:))); 
%     % binary 1 (fish uses sequence) or zero (fish doesn't) 
% end 

%% Common-ness of Grammar Figures 

% % Overall Percentage of Fish who use each sequence 
% figure; 
% box off; set(gca, 'Layer','top'); hold on; 
% set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format   
% histogram((sum(common_table,2)/size(common_table,2))*100,'Normalization','probability','binwidth',5,...
%     'EdgeColor','none','FaceColor',cmap{1}(1,:));
% xlabel('% of Fish','Fontsize',32); ylabel('Probability','Fontsize',32); % Axis Labels
% axis tight 
% 
% % Number of Sequences used by each fish  
% figure; 
% for er = 1:max(experiment_reps) % for each group of experiments 
%     set_token = find(experiment_reps == er,1,'first'); % settings
%     ax = subplot(2,3,er); box off; set(gca, 'Layer','top'); hold on;
%     set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format
%     spread_cols = plotSpread(sum(common_table(:,i_experiment_reps==er))',...
%         'distributionIdx',i_group_tags(i_experiment_reps == er),...
%         'distributionColors',cmap{set_token},'showMM',2);
%     set(findall(ax,'type','line'),'markersize',15); % change marker sizes
%     spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = 'k'; % Change Mean properties
%     spread_cols{2}.MarkerSize = 12;
%     ylabel('Number of Sequences','Fontsize',32); % Axis Labels
%     set(gca,'xticklabel',geno_list{set_token}.colheaders,'Fontsize',32); % Name each group
% end 
% 
% clear ax er 

%% Normalise Counts by number of frames per time window
% gCount_norm = gCount_merge(1,1); % sequences x time windows x fish
% tw_length_corrections = []; % allocate 
% 
% for er = 1:max(experiment_reps) % for each group of experiments
%     for e = find(experiment_reps == er) % for each experiment
%         
%         % Check if the water was topped up
%         if sum(sleep_cells_nan_track(experiment_tags{2,1}==e,1)) == 0 % no
%             gCount_norm{1,1}(:,1:(size(lb{e},1)-1),i_experiment_tags==e) = ...
%                 gCount_norm{1,1}(:,1:(size(lb{e},1)-1),i_experiment_tags==e)./diff(lb{e})';
%             % divide each time window by the number of frames
%         
%         else % yes
%             
%             for f = find(i_experiment_tags(i_experiment_reps==er)==e)' % for each fish
%                 
%                 tw_lengths = histcounts(find(states{er,1}(f,:) == 0),lb{e});
%                 % find top-up periods and which time window they fall into 
%                 
%                 gCount_norm{1,1}(:,1:(size(lb{e},1)-1),f+fish_tags_cm(er-1,1)) = ...
%                     gCount_norm{1,1}(:,1:(size(lb{e},1)-1),f+fish_tags_cm(er-1,1))./...
%                     (diff(lb{e})' - tw_lengths); % divide using the corrected number of frames 
%                 
%                 tw_length_corrections = [tw_length_corrections ; tw_lengths]; % store corrections 
%                 clear  tw_lengths   
%             end     
%         end
%     end
% end
% 
% clear er e f 

%% Identifying Interesting Sequences (is)

% 180328 Notes 
% 1. Need to zscore mRMR_data
    % Peng "Each feature variable in the raw data was preprocessed to have 
    % zero mean-value and unit variance (i.e., transformed to their
    % z-scores). 

% 180215 Notes
% 1. https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-mRMR_data.html#d119e2598
% 2. https://uk.mathworks.com/help/stats/fitcdiscr.html

% Settings
comps = 250; % number of sequences to identify

% mRMR Approach
for er = 1:max(experiment_reps) % for each experiment repeat
    set_token =  find(experiment_reps == er,1,'first'); % settings
    
    % Grab Data 
    if er == 1 % for the WT fish
        % Organised to be:
        % rows = all day mRMR_data, then all night mRMR_data 
        % columns = sequences 
        mRMR_data{er,1} = ...
            double([reshape(gCount_norm{1,1}(:,days_crop{set_token}(days{set_token}),...
            i_experiment_reps==er),size(grammar_mat{1,1},1),[])' ; ...
            reshape(gCount_norm{1,1}(:,nights_crop{set_token}(nights{set_token}),...
            i_experiment_reps==er),size(grammar_mat{1,1},1),[])']);
        
        mRMR_tw{er,1} = ones(size(mRMR_data{er,1},1),1)*2; 
        mRMR_tw{er,1}(1:size(mRMR_data{er,1},1)/2) = 1; % day (1) vs night (2) 
        
        %mRMR_data{er,1}(mRMR_data{er,1} < 0) = 0; % Remove negative values for now
        
    else
        % Organised to be:
        % rows = fish 1 (all time windows), fish 2 (all time windows) etc 
        % columns = sequences 
        mRMR_data{er,1} = ...
            double(reshape(gCount_norm{1,1}(:,time_window{set_token}(1):time_window{set_token}(2),...
            i_experiment_reps == er),size(uniqueSeqs{1,1},1),[]))';
        mRMR_tw{er,1} = repmat(i_group_tags(i_experiment_reps == er),1,...
            size(time_window{set_token}(1):time_window{set_token}(2),2))';
        mRMR_tw{er,1} = (mRMR_tw{er,1}(:)')';
    end
        
    % Pairwise Comparisons
    tic
    counter = 1; % start a counter
    for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group
        for g_two = (g_one + 1):max(mRMR_tw{er,1}) % for each comparison
            
            % mRMR
            [comps_v{er,1}(counter,:)] = mrmr_miq_d(...
                zscore(mRMR_data{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,:)),...
                mRMR_tw{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,:),comps);
            
            % Classifiers 
            for s = 1:comps % for each comp sequence
                    % Fit a linear classifier as you add features
                    % Using 10 fold cross validation
                    % Hold 10% of mRMR_data back by default
                    Mdl = fitcdiscr(...
                        zscore(mRMR_data{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,...
                        comps_v{er,1}(counter,1:s))),...
                        mRMR_tw{er,1}(mRMR_tw{er,1} == g_one | mRMR_tw{er,1} == g_two,:),...
                        'DiscrimType','linear','CrossVal','on');
                Mdl_loss{er,1}(counter,s) = kfoldLoss(Mdl);
                Mdl_loss{er,2}(counter,s) = nanstd(kfoldLoss(Mdl,'Mode','individual')); 
                %disp(num2str(s)); 
            end
            
            % Minimal Feature Space
            if er == 1 % for the WT data
                mRMR_ms(er,counter) = find(smooth(Mdl_loss{er,1}(counter,:),3) == ...
                    min(smooth(Mdl_loss{er,1}(counter,:),3)),1,'first');
            else
                try
                    mRMR_ms(er,counter) = find(Mdl_loss{er,1}(counter,:) < 0.05,1,'first');
                catch
                    mRMR_ms(er,counter) = find(smooth(Mdl_loss{er,1}(counter,:),3) == ...
                        min(smooth(Mdl_loss{er,1}(counter,:),3)),1,'first');
                end
            end
            
            counter = counter + 1; % add to counter
        end
    end
    disp(horzcat('Finished mRMR Comparisons ',num2str(er),' of ',...
        num2str(max(experiment_reps)))); % report progress
    toc
       
end

clear er set_token s Mdl

%% Sorted Raw Data Example
    
er = 1; % settings 
set_token =  find(experiment_reps == er,1,'first'); % settings
clear scrap data; 
scrap =  double([reshape(gCount_norm{1,1}(:,days_crop{set_token}(days{set_token}),...
            i_experiment_reps==er),size(grammar_mat{1,1},1),[])' ; ...
            reshape(gCount_norm{1,1}(:,nights_crop{set_token}(nights{set_token}),...
            i_experiment_reps==er),size(grammar_mat{1,1},1),[])']); 
scrap = sortrows(scrap','descend','ComparisonMethod','auto')'; % sort motifs   
data = [sortrows(scrap(1:(size(scrap,1)/2),:),'descend'); ...
    sortrows(scrap(((size(scrap,1)/2)+1):end,:),'descend')]; % sort fish (seperatly day/night)
figure; 
imagesc(zscore(data),[-0.5 0.5]); % imagesc 

% colormap 
n = 15; % number of colours 
CT = [linspace(cmap_2{1,1}(1,1),1,n)'...
    linspace(cmap_2{1,1}(1,2),1,n)'...
    linspace(cmap_2{1,1}(1,3),1,n)']; 
CT = [CT ; [linspace(1,cmap_2{1,1}(2,1),n)'...
    linspace(1,cmap_2{1,1}(2,2),n)'...
    linspace(1,cmap_2{1,1}(2,3),n)']]; 
CT(n,:) = []; % remove repeat 
CT = flip(CT); % flip colormap 
colormap(CT); 
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
c = colorbar; c.Label.String = 'Z-Score';
xlabel('Motif','Fontsize',32); 
set(gca,'YTick',[]); 
set(gca,'YTick',[size(mRMR_data{er,1},1)/4 size(mRMR_data{er,1},1)*(3/4)]); 
set(gca,'YTickLabels',{'Day' ; 'Night'},'Fontsize',32); 
ylabel('Fish ID','Fontsize',32); 

clear er set_token data n CT c 

%% Model Loss Figure 
er = 1; % settings 
set_token =  find(experiment_reps == er,1,'first'); % settings
figure; hold on; 
clear data; 
data = smooth(Mdl_loss{er,1}(1,:),3)'; 

a = plot(Mdl_loss{er,1}(1,:),'color',...
    cmap_2{set_token}(1,:)+(1-cmap_2{set_token}(1,:))*(1-(1/(5)^.5)),'linewidth',3); % raw data 
b = plot(data,'color',cmap{set_token}(1,:),'linewidth',3); % smoothed data
plot([mRMR_ms(er,1) mRMR_ms(er,1)],...
    [0 data(1,mRMR_ms(er,1))],...
    '--k','linewidth',1.5); 
plot([0 mRMR_ms(er,1)],...
    [data(1,mRMR_ms(er,1)) data(1,mRMR_ms(er,1))],...
    '--k','linewidth',1.5); 

set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
xlabel('Motifs','Fontsize',32); 
ylabel('Classification Error (%)','Fontsize',32); 
set(gca,'YTickLabels',{(get(gca,'YTick')*100)},'Fontsize',32); % convert labels to percentages 
legend([a b],'Raw Error','Smoothed Error','location','best'); 
legend('boxoff'); 
axis([1 size(data,2) ylim]); 

clear er set_token data a b 

%% Pairwise Classifier Figure
    % diag(Mdl_loss{5,1}(2:end,mRMR_ms(end,2:end)))

figure;
er = 5; % set interest 
set_token =  find(experiment_reps == er,1,'first'); % settings
clear scrap; counter = 1; % start a counter 
scrap = nan(max(mRMR_tw{er,1}),max(mRMR_tw{er,1}),'single'); 
for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group
    for g_two = (g_one + 1):max(mRMR_tw{er,1}) % for each comparison
        scrap(g_one,g_two) = mRMR_ms(er,counter);
        counter = counter + 1; % add to counter 
    end
end
ax = imagesc(scrap,'AlphaData',isnan(scrap)==0); % imagesc with nan values in white
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
c = colorbar; c.Label.String = 'Number of Motifs';
set(gca,'XTick',1:size(scrap,1)); set(gca,'XTickLabel',geno_list{set_token}.colheaders,'Fontsize',32);
set(gca,'YTick',1:size(scrap,1)); set(gca,'YTickLabel',geno_list{set_token}.colheaders,'Fontsize',32);

clear er set_token scrap g_one g_two counter ax c

%% Calculating tSNE Spaces 
for er = 1:max(experiment_reps) % for each experiment

    if er == 1
        mRMR_tsne{er,1} = tsne(mRMR_data{er,1}(:,comps_v{er,1}(1,1:mRMR_ms(er,1))),...
            'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
            'Perplexity',30,'Standardize',1,'Verbose',0);
    else
        mRMR_tsne_dim{er,1} = [];
        for c = 1:size(comps_v{er,1},1) % for each comparison
            mRMR_tsne_dim{er,1} = [mRMR_tsne_dim{er,1} comps_v{er,1}(c,1:mRMR_ms(er,c))];
        end
        mRMR_tsne_dim{er,1} = unique(mRMR_tsne_dim{er,1},'Stable');
        mRMR_tsne{er,1} = tsne(mRMR_data{er,1}(:,mRMR_tsne_dim{er,1}),...
            'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
            'Perplexity',30,'Standardize',1,'Verbose',0);
    end
end

clear er set_token  

%% Interesting Motifs Figure 
% Settings 
er = 1; % set interest
set_token =  find(experiment_reps == er,1,'first'); % settings

% Motifs
figure;
subplot(1,3,1); % subplot
clear scrap;
scrap = grammar_mat{1,1}(comps_v{er,1}(1,1:mRMR_ms(er,1)),:); % grab motifs
scrap(:,sum(isnan(scrap)) == mRMR_ms(er,1)) = [];
ax = imagesc(scrap,'AlphaData',isnan(scrap)==0); % imagesc with nan values in white
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(ax,'CDataMapping','direct');
%c = colorbar; c.Label.String = 'Cluster'; c.Location = 'westoutside';
xlabel('Position in Motif','Fontsize',32);
ylabel('Motif','Fontsize',32);

% WT Constraint/Enrichment 
subplot(1,3,2); hold on; set(gca,'Ydir','reverse'); 
plot([0 0],[(1-0.5) (mRMR_ms(er,1)+0.5)],'-k','linewidth',1.5); clear scrap;
for s = 1:mRMR_ms(er,1) % for each contextual motif 
    clear data; 
    data = squeeze(gCount_norm{1,1}(comps_v{er,1}(1,s),...
        time_window{set_token}(1):time_window{set_token}(2),...
        i_experiment_reps == er))';
    
    for g = 1:2 % for day/night
        errorbar(nanmean(reshape(data(:,g:2:end),[],1)),s,...
            nanstd(reshape(data(:,g:2:end),[],1)),'horizontal','-o',...
            'markersize',3,'MarkerEdgeColor',cmap_2{set_token}(g,:)...)
            ,'MarkerFaceColor',cmap_2{set_token}(g,:),'linewidth',1.5,'color',cmap_2{set_token}(g,:))
    end
    
end 
axis tight
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
xlabel('Z-Score','Fontsize',32);
set(gca,'YTick',[]); 

% tSNE Plot 
subplot(1,3,3); hold on; 
for g = 1:max(mRMR_tw{er,1}) % for each group
    if er == 1 % for the wt mRMR_data
        scatter(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,1),mRMR_tsne{er,1}(mRMR_tw{er,1} == g,2),...
            'markerfacecolor',cmap_2{set_token}(g,:),...
            'markeredgecolor',cmap_2{set_token}(g,:));
    else
        scatter(mRMR_tsne{er,1}(mRMR_tw{er,1} == g,1),mRMR_tsne{er,1}(mRMR_tw{er,1} == g,2),...
            'markerfacecolor',cmap{set_token}(g,:),...
            'markeredgecolor',cmap{set_token}(g,:));
    end
      
end
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
scrap = get(gca,'Children');
if er == 1 % for the WT data
    [~,icons,plots,~] = legend([scrap(2) scrap(1)],'Day','Night','location','best');
else
    
end
legend('boxoff'); 
set(gca,'XTick',[]); set(gca,'YTick',[]); 
xlabel('tSNE 1','Fontsize',32); 
ylabel('tSNE 2','Fontsize',32); 

clear er set_token scrap ax c s data g icons plots     

%% PlotSpread of an IS. 
% Settings
er = 1; % set experiment of interest  
s = comps_v{er,1}(1,1); % choose sequence of interest
set_token =  find(experiment_reps == er,1,'first'); % settings 
sep = [0 0.75/2]; % day/night spread 

figure; 
hold on; set(gca,'FontName','Calibri'); set(gca,'Fontsize',32);
col = 1; 

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    
    % Days
    spread_cols(1,:) = plotSpread(reshape(squeeze(gCount_norm{1,1}(s,days_crop{set_token}(days{set_token}),...
        i_experiment_reps == er & i_group_tags == g)),[],1),'xValues',g + sep(1),...
        'distributionColors',cmap_2{set_token}(col,:),'spreadWidth',sep(2),'showMM',2);
    col = col + 1;
    
    % Nights
    spread_cols(2,:) = plotSpread(reshape(squeeze(gCount_norm{1,1}(s,nights_crop{set_token}(nights{set_token}),...
        i_experiment_reps == er & i_group_tags == g)),[],1),'xValues',g + sep(2),...
        'distributionColors',cmap_2{set_token}(col,:),'spreadWidth',sep(2),'showMM',2);
    
    % Change Format 
    set(findall(gca,'type','line'),'markersize',30); % change marker sizes
    spread_cols{1,2}.LineWidth = 6; spread_cols{1,2}.Color = 'k'; % Change Mean properties
    spread_cols{1,2}.MarkerSize = 24;
    spread_cols{2,2}.LineWidth = 6; spread_cols{2,2}.Color = 'k'; % Change Mean properties
    spread_cols{2,2}.MarkerSize = 24;
    col = col + 1;

end

x_lims = xlim; 
r(1) = plot([(x_lims(1)+(x_lims(1)*0.05)) (x_lims(2)+(x_lims(2)*0.05))],...
    [0,0],'--','color','k','linewidth',1.5); % plot zero line  
uistack(r(1),'bottom'); % Send to back

axis tight
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
set(gca, 'XTick', (1:g)+(sep(2)/2));
set(gca,'XTickLabels',1:g,'Fontsize',32);
set(gca,'XTickLabels',geno_list{set_token}.colheaders)
ylabel('Z-Score','Fontsize',32);     

clear er s set_token sep col g t spread_cols x_lims r 

%% Grammar Comparison Figure  

% Settings
ts = 5; % number of sequences to show 
cmap_cluster_merge = [cmap_cluster{2,1} ; cmap_cluster{1,1}]; % merged colormap

for er = 1:max(experiment_reps) % for each group of experiments
    figure; clear scrap;
    counter = 1; % counts subplots
    set_token = find(experiment_reps == er,1,'first'); % settings
    
    for s = comps_v{er,1}(1:ts,1)' % for each sequence to show 
        
        % Plot Lines
        subplot(2,ts,counter); hold on; set(gca,'FontName','Calibri');
        for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
            clear data; data = ...
                squeeze(gCount_norm{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g))';
           
            errorbar(nanmean(data),(nanstd(data)/...
                sqrt(sum(i_experiment_reps == er & i_group_tags == g))),...
                'color',cmap{set_token}(g,:),'linewidth',3);
            
            scrap(1,g) = min(nanmean(data) - (nanstd(data)/...
                sqrt(sum(i_experiment_reps == er & i_group_tags == g))));
            scrap(2,g) = max(nanmean(data) + (nanstd(data)/...
                sqrt(sum(i_experiment_reps == er & i_group_tags == g))));        
        end
        
        % Add night patches
        y_lims = [min(scrap(1,:)) max(scrap(2,:))]; % Add a bit of space either side
        
        a = 1; night_start = first_night{set_token}; % Start counters
        for n = 1:size(nights{set_token},2) % For each night
            r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
                1 (y_lims(2)-y_lims(1))],...
                'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
            uistack(r(a),'bottom'); % Send to back
            a = a + 1; night_start = night_start + 2; % Add to counters
        end
        
        % Figure Formatting
        axis([0.5 size([days_crop{set_token}(days{set_token}) nights_crop{set_token}(nights{set_token})],2)+0.5 ...
            y_lims]); % Set axis
        box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
        xlabel('Time (Days/Nights)','Fontsize',12); % X Labels
        set(gca, 'XTick', []); % Turn off X-Ticks
        ylabel({'Z-Score'},'Fontsize',12); % Y Labels
        
        % Legend 
        if counter == ts
            [~,~,~,~] = legend(geno_list{set_token}.colheaders,'Location','best'); % Generate axis
            legend('boxoff'); % Turn legend off
        end
        
        % Plot Sequence - from randomly chosen fish
        subplot(2,ts,counter + ts); hold on; set(gca,'FontName','Calibri');
        if er == 1 % for the WT fish
            exampleFish = datasample(find(sum(squeeze(gCount_norm{1,1}...
                (s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er))) >= 1),...
                1,'replace',false); % find a fish who uses this sequence
        else
            exampleFish = datasample(find(sum(squeeze(gCount_norm{1,1}...
                (s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er))) >= 1),...
                1,'replace',false) + fish_tags_cm(er - 1); % find a fish who uses this sequence
        end
        Tlocs = datasample(strfind(threads{exampleFish,1,1}',uniqueSeqs{1,1}{s,1}),...
            1,'replace',false); % sample an instance of this sequence (position in threads)
        Rlocs = threads{exampleFish,2,1}(Tlocs,1):...
            threads{exampleFish,2,1}(Tlocs+(size(uniqueSeqs{1,1}{s,1},2)-1),2); % start-end (position in frames)
        Rlocs = Rlocs + offset(er,(i_experiment_tags(exampleFish) - ...
            (min(unique(i_experiment_tags(i_experiment_reps == er))) - 1))); % position in offset data
        if er ~= 1 % fish tag for raw_data
            exampleFish = exampleFish - fish_tags_cm(er - 1);
        end
        ax = imagesc([1,size(raw_data{er,1}(exampleFish,Rlocs),2)],...
            [0,max(raw_data{er,1}(exampleFish,Rlocs))],...
            repmat(states{er,1}(exampleFish,Rlocs),[max(raw_data{er,1}(exampleFish,Rlocs)),1]));
        hold on; colormap(cmap_cluster_merge); set(gca,'Ydir','Normal');
        plot(raw_data{er,1}(exampleFish,Rlocs),'k','linewidth',3); hold on;
        
        % Figure Formatting
        axis([1 size(raw_data{er,1}(exampleFish,Rlocs),2) 0 max(raw_data{er,1}(exampleFish,Rlocs))]);
        set(gca,'XTickLabels',xticks/fps{1});
        box off; set(gca,'FontName','Calibri'); set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
        xlabel('Time (Seconds)','FontSize',12);
        ylabel('Delta Px','FontSize',12);
        set(ax,'CDataMapping','direct');
        
        counter = counter + 1; % add to counter
        
    end
    
end

clear er counter set_token s g data scrap y_lims a n r ...
    exampleFish Tlocs Rlocs ax  

%% Localising IS in Time 

for er = 1:max(experiment_reps) % for each group of experiments
    for t = 1:ts
        timeSeqs{er,t} = sparse(sum(i_experiment_reps == er),...
            size(raw_data{er,1},2));
        % experiment reps x is {fish x frames)
    end
end

clear scrap;
for er = 1:max(experiment_reps) % for each group of experiments
    counter_2 = 1; % start a counter (counts sequences)
    
    for s = comps_v{er,1}(1:ts,1)' % for each sequence to show 
        counter = 1; % start a counter (counts fish)
        
        for f = find(i_experiment_reps == er)' % for each fish
            
            Tlocs = strfind(threads{f,1,1}',uniqueSeqs{1,1}{s,1}); % find pattern starts (positions in threads)
            Rlocs = threads{f,2,1}(Tlocs,1); % find pattern starts in real time (frames)
            Rlocs = Rlocs + offset(er,(i_experiment_tags(f) - ...
                (min(unique(i_experiment_tags(i_experiment_reps == er))) - 1))); % correct for offset (frames)
            
            timeSeqs{er,counter_2}(counter,Rlocs) = 1; % fill in ones
            counter = counter + 1; % add to counter (counts fish)
            
        end
        
        counter_2 = counter_2 + 1; % add to counter (counts sequences)
        
    end
end

%% IS in Time Figure 

figure; 
er = 1; 
for s = 1:ts
    subplot(1,ts,s); hold on;
    clear scrap; scrap = full(timeSeqs{er,s}); 
    for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
        plot(smooth(sum(scrap(i_group_tags(i_experiment_reps==er)==g,:))/...
            sum(sum(timeSeqs{er,s}(i_group_tags(i_experiment_reps==er)==g,:))),25*60*15),...
            'color',cmap{set_token}(g,:)); 
        
    end
end

%% IS Transition Probabilities 
steps = 1; % Markov Steps to Take 

% Pre-allocate
transitions = zeros(all_states,all_states,size(threads,1),...
    7,'single');
% states x states x time windows
    
for er = [1 4 5] % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    
    for f = find(i_experiment_reps==er)' % for each fish 
        
        for t = time_window{set_token}(1):time_window{set_token}(2)
            transitions(:,:,f,t) = accumarray(...
                [threads{f,1,1}(find(threads{f,3,1} == t,1,'first'):...
                find(threads{f,3,1} == t,1,'last')-steps,1),...
                threads{f,1,1}(find(threads{f,3,1} == t,1,'first')+steps:...
                find(threads{f,3,1} == t,1,'last'),1)],...
                1,[all_states,all_states]); % Accumulate transitions
            transitions(:,:,f,t) = bsxfun(@rdivide,transitions(:,:,f,t),...
                sum(transitions(:,:,f,t),2)); % Convert to probabilities
        end
    end
end

%% Statistical Comparisons of Counts 

% % Allocate 
% twa.gs = cell(1,max(experiment_reps)); 
% 
% for er = 1:max(experiment_reps) % for each group of experiments
%     set_token =  find(experiment_reps == er,1,'first'); % settings 
%     
%     % Grouping variables
%     anova_group = repmat(i_group_tags(i_experiment_reps==er),...
%         [size([days{set_token} nights{set_token}],2),1])'; % groups
%     anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
%         [size([days{set_token} nights{set_token}],2),1])'; % experiments
%     anova_time = [];
%     for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
%         anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
%         % Allocate alternating zeros and ones to each time window
%     end
%     anova_time = anova_time';
%     
%     % Statistical Comparisons 
%     tic
%     for s = 1:size(uniqueSeqs{1,1},1) % for each sequence
%         
%         data = squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
%             i_experiment_reps == er))';
%         data = data(:)'; % vectorise
%         
%         [twa.gs{1,er}(s,:)] = single(anovan(data,...
%             {anova_group,anova_time,anova_experiment},...
%             'display','off','model','full'));
%         
%         % Allocate 
%         if s == 1 % for the first sequence 
%            twa.gs{1,er}(2:size(uniqueSeqs{1,1},1),1:size(twa.gs{1,er},2)) = single(NaN); 
%            %twa.gs{2,er} = [twa.gs{2,er} ; cell(size(uniqueSeqs{1,1},1)-1,1)]; 
%         end 
%         
%         % Report Progress 
%         if mod(s,1000) == 0
%            disp(horzcat('Completed ',num2str(s),' Comparisons of ',...
%                num2str(size(uniqueSeqs{1,1},1)),' er = ',num2str(er)));
%         end 
%     end
%     toc 
%     clear set_token anova_group anova_experiment anova_time data 
% end

%% Correct P-values for multiple comparisons 

%% Correlations between fish 
% scrap = corrcoef(squeeze(sum(gCount_merge{1,1},2)));
% 
% er = 5; 
% 
% imAlpha=ones(size(scrap(i_experiment_reps == er,i_experiment_reps == er)),'single'); 
% imAlpha(tril(scrap(i_experiment_reps == er,i_experiment_reps == er)) > 0) = 0; % find nan values 
% imagesc(scrap(i_experiment_reps == er,i_experiment_reps == er),'AlphaData',imAlpha);
% colormap(cmap_cluster_merge); colorbar; 
