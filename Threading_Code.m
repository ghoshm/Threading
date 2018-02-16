% Threading Code  

% Aims to use data from State_Space_4 (clustered active and inactive bouts)
    % To study transition dynamics 
    
%% Required scripts 

%% Notes 

% The code for calculating the transition probabilities is adapted from 
 %https://uk.mathworks.com/matlabcentral/answers/57877-seeking-help-creating-a-transition-probability-matrix-for-a-markov-chain

 %% Load data (Post State_Space_4)
    % Note that running State_Space_4 first is necessary to re-order the
        % Clusters by a suitable metric 
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 

% Load  
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat')

% Load Delta_px_sq data 
for f = 1:size(filename,2) %For each file
    delta_px_sq{1,f} = load(strcat(pathname,filename{f}),'delta_px_sq'); % load delta_px_sq data  
end 

%% Threading 

% Pre-allocation
threads = cell(max(fish_tags{1,1}),3,2); % fish x (clusters,times (start & stop),time windows) x (samples vs controls) 
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
    
    % Control 
    threads{f,1,2} = nan(size(threads{f,1,1}),'single'); % clusters    
    threads{f,3,2} = nan(size(threads{f,1,1}),'single'); % time windows
    
    % Deterime starting state (a = wake, b = sleep) 
    if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active 
        a = 1; b = 2; 
    else % If the fish starts inactive 
        a = 2; b = 1; 
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
   
    % Generate Control data 
    % Clusters
    clear scrap; scrap = idx_numComp_sorted{1,1}...
        (fish_tags{1,1} == f,1); % Take active bout data 
    threads{f,1,2}(a:2:end,1) = scrap(randperm(length(scrap))); % Mix  
    clear scrap; scrap = idx_numComp_sorted{2,1}...
        (fish_tags{2,1} == f,1); % Take inactive bout data 
    threads{f,1,2}(b:2:end,1) = scrap(randperm(length(scrap))); % Mix 
    % Time windows 
    threads{f,3,2} = threads{f,3,1}; % Break up using the same time windows  
    
    % Report progress 
    if mod(f,100) == 0 % every 100 fish 
        disp(horzcat('Threaded Fish ',num2str(f),' of ',...
            num2str(max(fish_tags{1,1})))); % Report progress
    end
end 
toc 

clear f a b scrap  

%% Filling in Data 
tic
% Calculate time off-sets between repeats of experiments 
c = nan(max(experiment_reps),1,'single'); % experiment groups x 1 
for er = 1:max(experiment_reps) % for each group of experiments
    lb_merge{er,1} = []; % structure 
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

%% Distribution in time Figure (Run this one figure @ a time) 

all_states = double(max(idx_numComp_sorted{1,1})); % possible states
    % Note that subplot can't take single values
ai_states(1:all_states) = 2; % specify if states are active (1) or inactive (2)
ai_states(min(idx_numComp_sorted{1,1}):end) = 1; % specify active states (1)
time_bins = fps{1}*60*5; % set smoothing window (note assumes a constant frame rate across experiments)

for er = 1:max(experiment_reps) % for each group of experiments
    clear scrap legend_lines legend_cell
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings
    figure; % make a figure for this experiment 
    counter = 1; % start counter (counts subplots)
    for s = 1:2 % for active/inactive
        counter_2 = 1; % start counter_2 (counts clusters) 
        for c = find(ai_states == s) % for each cluster 
            subplot(4,4,counter); hold on; clear scrap; set(gca,'FontName','Calibri');
            title(horzcat(strings{s},' - ',num2str(counter_2)));
            
            for g = 1:max(i_group_tags(i_experiment_reps == er)) % For each group
                if er == 1 % for the WT Experiments
                    legend_lines(1,g) = plot(lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1),...
                        smooth(sum(states{er,1}(i_group_tags(i_experiment_reps == er) == g,...
                        lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1))==c)./...
                        sum(isnan(states{er,1}(i_group_tags(i_experiment_reps == er) == g,...
                        lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1)))==0),...
                        time_bins),'color',cmap_cluster{s}(counter_2,:),'linewidth',1.5);
                else
                    legend_lines(1,g) = plot(lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1),...
                        smooth(sum(states{er,1}(i_group_tags(i_experiment_reps == er) == g,...
                        lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1))==c)./...
                        sum(isnan(states{er,1}(i_group_tags(i_experiment_reps == er) == g,...
                        lb_merge{er,1}(time_window{set_token}(1)):lb_merge{er,1}(time_window{set_token}(2)+1)))==0),...
                        time_bins),'color',cmap{set_token}(g,:),'linewidth',1.5);
                end
                
                % Find the top & bottom
                % Excluding values that smooth outside of the data
                scrap(1,g) = max(legend_lines(1,g).YData((1+time_bins):...
                    (size(legend_lines(1,g).YData,2) - time_bins)));
                scrap(2,g) = min(legend_lines(1,g).YData((1+time_bins):...
                    (size(legend_lines(1,g).YData,2) - time_bins)));
                
                % Legend
                if s == 1 && counter_2 == 4 
                    legend_cell{g} = horzcat(geno_list{set_token}.colheaders{g},', n = ',...
                        num2str(sum(i_group_tags(i_experiment_reps == er) == g)));
                    % Append the group size to each group name
                end
            end
            
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
            
            % Legend
            if s == 1 && counter_2 == 4 % For the 4th state - add a legend
                [~,~,~,~] = legend(legend_cell,'Location','northeast');
                legend('boxoff');
            end
            
            % Axis etc
            axis([(lb_merge{er,1}(time_window{set_token}(1)) + (1 + time_bins))...
                (lb_merge{er,1}(time_window{set_token}(2)+1) - time_bins) ...
                min(scrap(2,:)) - (min(scrap(2,:))*0.05) ...
                max(scrap(1,:)) + (max(scrap(1,:))*0.05)]);
            box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12);
            xlabel('Time (Days/Nights)','Fontsize',12);
            set(gca, 'XTick', []);
            yticks(round(max(scrap(1,:)),2,'Significant')); % just label the max frequency 
            yticklabels(num2str(round(max(scrap(1,:)),2,'Significant')));
            ylabel('Frequency','Fontsize',12);
            
            counter = counter + 1; % (counts subplots)
            counter_2 = counter_2 + 1; % (counts clusters)
        end
    end
end

clear er scrap legend_lines legend_cell set_token counter counter_2 s c g scrap a n r  

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
figure; cols = 1;
for c = min(idx_numComp_sorted{1,1}):max(idx_numComp_sorted{1,1}) % for each active cluster 
    scatter3(sample(sample_tags==c,1),sample(sample_tags==c,2),...
        density(sample_tags==c,1),90,...
        'markerfacecolor',cmap_cluster{1,1}(cols,:),...
        'markeredgecolor','k','linewidth',0.01);
    cols = cols + 1;
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
experiment_reps_long{1,1} = experiment_tags{1,1}; % temp assign 
for er = 1:max(experiment_reps) % for each repeat 
    found = find(experiment_reps == er); % find experiments
    
    for f = found % for each experiment in the repeat 
        experiment_reps_long{1,1}(experiment_reps_long{1,1} == f,1) = er; % tag with grouping variable 
    end 
    
end 

number = 1000; % number of bouts to sample from each active cluster 
counter = 1; % counts active clusters 
for c = find(ai_states == 1) % for each active bout type
    
    clear sample sample_tags sample_er sample_fish sample_et 
    
    % Weighted Sample bouts
    [sample,sample_tags] = datasample(...
        wake_cells(idx_numComp_sorted{1,1}==c,:),...
        number,1,'replace',false,'weights',P{1,1}(idx_numComp_sorted{1,1}==c,:));
    sample_tags = sample_tags'; % flip 
    
    % Allocate Space
    bouts{1,counter} = zeros(number,max(sample(:,3)),'single');
    
    % Sample tags
    sample_er = experiment_reps_long{1,1}(idx_numComp_sorted{1,1}==c,:); % experiment groups 
    sample_fish = fish_tags{1,1}(idx_numComp_sorted{1,1}==c,:); % fish i.d.
    sample_et = experiment_tags{1,1}(idx_numComp_sorted{1,1}==c,:); % experiment tags 
    
    % Fill in data
    for b = 1:number % for each bout
        
        % determine bout indexing numbers
        er = sample_er(sample_tags(b,1)); % sample experiment groups 
        if er == 1 % for the first experiment group
            fish = sample_fish(sample_tags(b,1)); % sample fish tags 
        else % for other experiment groups 
            fish = sample_fish(sample_tags(b,1)) - fish_tags_cm(er - 1); 
            % determine correct fish id
        end
        et = find(find(experiment_reps == er) == sample_et(sample_tags(b,1))); % sample experiment tags 
        
        % Extract raw data
        bouts{1,counter}(b,1:sample(b,3)) = raw_data{er,1}(fish,...
            (sample(b,1)+offset(er,et)):(sample(b,2)+offset(er,et)));
    end
    
    counter = counter + 1; % counts active clusters 
end

clear er f number counter c sample sample_tags sample_er ...
    sample_fish sample_et b fish et found 

%% Bout Shapes Figure 

figure; hold on; set(gca,'FontName','Calibri');
for b = 1:size(bouts,2) % for each active bout type
    
    legend_lines(b) = shadedErrorBar((0:(size(bouts{1,b},2)+1))/fps{1},[0 nanmean(bouts{1,b}) 0],...
        [0 nanstd(bouts{1,b})./sqrt(size(bouts{1,b},1)) 0],'lineprops',...
        {'color',cmap_cluster{1,1}(b,:)}); 
            
    legend_lines(b).mainLine.LineWidth = 1.5;
    legend_cols(b) = legend_lines(b).mainLine; % Store color
    legend_cell{b} = horzcat('Cluster ',num2str(b)); 
    
end
box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
axis([0 8/fps{1} ylim]); % hard coded 
set(gca,'XTick',(1/fps{1}):(2/fps{1}):(8/fps{1})); % hard coded 
set(gca,'XTickLabels',{(round((1/fps{1}):(2/fps{1}):(8/fps{1}),2,'decimals'))}); % hard coded 
xlabel('Time (seconds)','Fontsize',32);
ylabel('Delta Px','Fontsize',32);
[~,icons,plots,~] = legend(legend_cols,legend_cell,'Location','northeast');
legend('boxoff'); set(icons(1:size(bouts,2)),'Fontsize',16) ; set(plots,'LineWidth',2);

clear b icons plots legend_lines legend_cols legend_cell  

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

%% Load in Legion Data (START HERE)
% Ease memory 
clear raw_data_all cells 

clear compVec gTermCell grammar totSavings
load('D:\Behaviour\SleepWake\Re_Runs\Threading\Compression_Results_Final.mat')

%% Constructing a common Grammar 

tic % 19s for 629 fish
% Merge all the grammers and keep only unique elements 
for tc = 1:2 % for real vs shuffled data 
    for f = 1:size(gTermCell,1) % for each fish
        % append grammar terminals to the total set
        uniqueSeqs{1,tc} = vertcat(uniqueSeqs{1,tc}, gTermCell{f,tc});
    end
end
disp('Merged all grammers'); 
toc 

tic % 106s for 629 fish 
% Determine unique sequences from larger set
[uniqueSeqs{1,1}, ~] = countUniqueRows(uniqueSeqs{1,1});
[uniqueSeqs{1,2}, ~] = countUniqueRows(uniqueSeqs{1,2});
disp('Determined unique sequences'); 
toc 

%% Remove Sequences with NaN (tagged as zeros) 

for tc = 1:2 % for real vs shuffled data 
    nan_locs = cellfun(@(s) ismember(0, s), uniqueSeqs{1,tc}); % find sequences with nans 
    uniqueSeqs{1,tc}(nan_locs,:) = []; % remove these 
end 

clear tc nan_locs 

%% Sequence Lengths 

seq_lengths{1,1} = zeros(size(uniqueSeqs{1,1},1),1,'single'); % tc - seqs x 1 
seq_lengths{1,2} = zeros(size(uniqueSeqs{1,2},1),1,'single'); % tc - seqs x 1 

for tc = 1:2 % for real vs shuffled data
    
    for s = 1:size(uniqueSeqs{1,tc},1) % for each sequence
        seq_lengths{1,tc}(s,1) = size(uniqueSeqs{1,tc}{s,1},2); % find it's length
    end
    
end

clear tc s 

%% Sequence Lengths Figure 

figure; subplot(1,2,1); box off; set(gca, 'Layer','top'); 
set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format
title('Data','Fontsize',32); hold on; 
histogram(seq_lengths{1,1},'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap{1}(1,:));
xlabel('Terminal Length','Fontsize',32); ylabel('Probability','Fontsize',32); % Y Labels
axis([(min([seq_lengths{1,1} ; seq_lengths{1,2}])-0.5) (max([seq_lengths{1,1} ; seq_lengths{1,2}])+0.5) 0 0.65]) % hard coded 
subplot(1,2,2); box off; set(gca, 'Layer','top'); 
set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format 
title('Shuffled Data','Fontsize',32); hold on; 
histogram(seq_lengths{1,2},'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap{1}(1,:));
xlabel('Terminal Length','Fontsize',32); ylabel('Probability','Fontsize',32); % Y Labels
axis([(min([seq_lengths{1,1} ; seq_lengths{1,2}])-0.5) (max([seq_lengths{1,1} ; seq_lengths{1,2}])+0.5) 0 0.65])

%% Construct Grammar Matrix  

grammar_mat{1,1} = zeros(size(uniqueSeqs{1,1},1),max(seq_lengths{1,1}),'single'); % sequences x max length 
grammar_mat{1,2} = zeros(size(uniqueSeqs{1,2},1),max(seq_lengths{1,1}),'single'); % sequences x max length 
    % Note that the shuffled matrix (grammar_mat{1,2}, must have the same
    % number of columns as the real data) also ismember won't take nan
    % values, therefore these are zero for now 
    
for tc = 1:2
    for s = 1:size(uniqueSeqs{1,tc},1) % for each sequence
        grammar_mat{1,tc}(s,1:size(uniqueSeqs{1,tc}{s,1},2)) = uniqueSeqs{1,tc}{s,1}; % fill in sequence
    end
end

%% "Translating" between the two grammars 
[~,grammar_mat_trans] = ismember(grammar_mat{1,1},grammar_mat{1,2},'rows'); % real sequences x 1 

grammar_mat{1,1}(grammar_mat{1,1}==0) = NaN; % replace zeros with nan 
grammar_mat{1,2}(grammar_mat{1,2}==0) = NaN; % replace zeros with nan 

%% Grammar Matrix Figure 

figure; 
grammar_mat_sorted = flip(sortrows(grammar_mat{1,2})); % sort rows of grammar_mat  
ax = imagesc(grammar_mat_sorted,'AlphaData',isnan(grammar_mat_sorted)==0); % imagesc with nan values in white 
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap  
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(ax,'CDataMapping','direct'); 
c = colorbar; c.Label.String = 'Cluster'; 
xlabel('Position in Sequence','Fontsize',32); 
ylabel('Sequence','Fontsize',32);
clear ax grammar_mat_sorted c 

%% Compressibility 
% The compressibility of a sequence of uncompressed length l is given by the sum of the savings S 
% at each iteration divided by l (Gomez-Marin et al.,2016) 

compressibility = zeros(size(threads,1),2,'single'); % fish x t/c
for f = 1:size(threads,1) % for each fish 
    compressibility(f,:) = totSavings(f,:)./size(threads{f,1,1},1);
end 

clear f 

%% Compressibility N-Way Anova  

for er = 1:max(i_experiment_reps) % for each group of experiments
    
    clear anova_group anova_tc anova_experiment data
    
    % Grouping Variables
    anova_group = repmat(i_group_tags(i_experiment_reps==er),[2,1])'; % groups
    anova_tc = ones(size(anova_group)); % real vs shuffled data
    anova_tc(1:size(i_group_tags(i_experiment_reps==er),1)) = 0;
    anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),[2,1])'; % experiments
    
    % Data to Compare 
    if er == 1 % WT - compare real vs shuffled data 
        data = compressibility(i_experiment_reps==er,:); % grab data
        data = data(:)'; % vectorise
    else % compare difference in compression between groups 
        data = compressibility(i_experiment_reps==er,1) - ... 
            compressibility(i_experiment_reps==er,2); % difference in compressibility 
        data = data(:)'; % vectorise 
        anova_group = anova_group(1:size(data,2)); % trim  
        anova_tc = anova_tc(1:size(data,2)); % trim 
        anova_experiment = anova_experiment(1:size(data,2)); % trim 
    end
    
    % Comparison 
    [twa.compress.p{1,er},~,twa.compress.stats{1,er}] = anovan(data,...
        {anova_group,anova_tc,anova_experiment},...
        'display','off','model','full');
    
end

clear er anova_group anova_tc anova_experiment data

%% Compressibility Figure 

figure;
for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % settings
    subplot(2,3,er); counter = 1; clear scrap; 
    hold on; set(gca,'FontName','Calibri');
    
    for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
        plot([counter,counter+1],compressibility(i_experiment_reps == er & i_group_tags == g,:)',...
            'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
        errorbar([counter,counter+1],nanmean(compressibility(i_experiment_reps == er & i_group_tags == g,:)),...
            nanstd(compressibility(i_experiment_reps == er & i_group_tags == g,:)),...
            'color',cmap{set_token}(g,:),'linewidth',3);
        counter = counter + 2; % add to counter
        
        scrap(1,g) = min(min(compressibility(i_experiment_reps == er & i_group_tags == g,:))); 
        scrap(2,g) = max(max(compressibility(i_experiment_reps == er & i_group_tags == g,:))); 
    end
    
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    if er == 1
        set(gca, 'XTick', [1 2]); % set X-ticks
        set(gca,'XTickLabels',{'Data','Shuffled'}); % X Labels
    else
        set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
    end
    ylabel('Compressibility','Fontsize',32); % Y Labels
    axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 ...
        (min(scrap(1,:)) - (min(scrap(1,:))*0.05)) (max(scrap(2,:)) + (max(scrap(2,:))*0.05))]);
end

clear er set_token g scrap counter 

%% Grammar in Time 

% Convert data to singles 
for tc = 1:2 % for real/control data
    for i = 1:size(uniqueSeqs{1,tc},1) % for each sequence
        uniqueSeqs{1,tc}{i,1} = single(uniqueSeqs{1,tc}{i,1}); % convert to single
    end
    for f = 1:size(threads,1) % for each fish
        for c = 1:size(threads,2) % for each column
            threads{f,c,tc} = single(threads{f,c,tc}); % convert to single 
        end
    end
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

%% Load Grammar Freqs From Legion 

load('D:\Behaviour\SleepWake\Re_Runs\Threading\Grammar_Results_Final.mat',...
    'gCount','gFreq')

%% Merge Legion Data 

gFreq_merge = cell(1,2); % allocate
gFreq_merge{1,1} = zeros(size(gFreq{1,1},1),size(gFreq{1,1},2),...
    size(gFreq,1),'single'); % {t,c} - seqs x time windows x fish
gFreq_merge{1,2} = zeros(size(gFreq{1,2},1),size(gFreq{1,2},2),...
    size(gFreq,1),'single'); % {t,c} - seqs x time windows x fish
gCount_merge = gFreq_merge; % {t,c} - seqs x time windows x fish 

for tc = 1:2 % for data/control
    for f = 1:size(gFreq,1) % for each fish
        gFreq_merge{1,tc}(:,:,f) = gFreq{f,tc}; % fill data
        gCount_merge{1,tc}(:,:,f) = gCount{f,tc}; % fill data 
    end
end

clear tc f gFreq gCount

 %% Grammar in Time Figure (VERY SLOW...) 
% 
% % Generate dn_colormap
% colors_p = [linspace(cmap_2{1}(1,1),cmap_2{1}(2,1),20)',...
%     linspace(cmap_2{1}(1,2),cmap_2{1}(2,2),20)',...
%     linspace(cmap_2{1}(1,3),cmap_2{1}(2,3),20)'];
% 
% figure;
% for tc = 1:2 % for real vs shuffled data
%     subplot(1,2,tc); hold on; set(gca,'FontName','Calibri');
%     if tc == 1 
%         title('Data','Fontsize',32); 
%     else 
%         title('Shuffled Data','Fontsize',32); 
%     end 
%     
%     scrap = nanmean(gFreq_merge{1,tc},3); % data 
%     counter = 1; % start a counter
%     for thrs = 1:-0.1:-1 % for each threshold (day to night)
%         try
%             plot(scrap(dn_score{1,tc} > (thrs-0.1) & dn_score{1,tc} <= thrs,:)',...
%                 'color',colors_p(counter,:),'linewidth',0.1) % Note will skip any == -1.0
%             % only really a problem if these exist 
%         catch
%         end
%         counter = counter + 1;
%     end
%     
%     if tc == 1 % use the same top for both plots 
%         y_lims = ylim; % find top
%     end 
%     a = 1; night_start = first_night{1}; % Start counters
%     for n = 1:size(nights_crop{1},2) % For each night
%         r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
%             1 (y_lims(2)-y_lims(1))],...
%             'FaceColor',night_color{1},'Edgecolor',[1 1 1]);
%         uistack(r(a),'bottom'); % Send to back
%         a = a + 1; night_start = night_start + 2; % Add to counters
%     end
%     
%     axis([0.5 7.5 0 y_lims(2)]); % hard coded 
%     box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
%     xlabel('Time (Days/Nights)','Fontsize',32); % X Labels
%     set(gca, 'XTick', []); % Turn off X-Ticks
%     ylabel('Mean Frequency','Fontsize',32); % Y Labels
%     
%     clear a night_start n r scrap counter
% end
%  
% clear y_lims 

%% "Common-ness" of Grammar

% Logical Table 
common_table = zeros(size(uniqueSeqs{1,1},1),size(gFreq_merge{1,1},3),'single'); % sequences x fish  

for s = 1:size(common_table,1) % for each sequence 
    common_table(s,:) = sum(squeeze(gFreq_merge{1,1}(s,:,:))); 
    % binary 1 (fish uses sequence) or zero (fish doesn't) 
end 

%% Common-ness of Grammar Figures 

% Overall Percentage of Fish who use each sequence 
figure; 
box off; set(gca, 'Layer','top'); hold on; 
set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format   
histogram((sum(common_table,2)/size(common_table,2))*100,'Normalization','probability','binwidth',5,...
    'EdgeColor','none','FaceColor',cmap{1}(1,:));
xlabel('% of Fish','Fontsize',32); ylabel('Probability','Fontsize',32); % Axis Labels
axis tight 

% Number of Sequences used by each fish  
figure; 
for er = 1:max(experiment_reps) % for each group of experiments 
    set_token = find(experiment_reps == er,1,'first'); % settings
    ax = subplot(2,3,er); box off; set(gca, 'Layer','top'); hold on;
    set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format
    spread_cols = plotSpread(sum(common_table(:,i_experiment_reps==er))',...
        'distributionIdx',i_group_tags(i_experiment_reps == er),...
        'distributionColors',cmap{set_token},'showMM',2);
    set(findall(ax,'type','line'),'markersize',15); % change marker sizes
    spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = 'k'; % Change Mean properties
    spread_cols{2}.MarkerSize = 12;
    ylabel('Number of Sequences','Fontsize',32); % Axis Labels
    set(gca,'xticklabel',geno_list{set_token}.colheaders,'Fontsize',32); % Name each group
end 

clear ax er 

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

%% Normalise Counts by shuffled data 
gCount_norm = gCount_merge(1,1); % sequences x time windows x fish

for s = 1:size(grammar_mat{1,1},1) % for each real sequence
    if grammar_mat_trans(s) > 0 % if the sequence is found in the shuffled data
       gCount_norm{1,1}(s,:,:) = gCount_merge{1,1}(s,:,:) - ...
           gCount_merge{1,2}(grammar_mat_trans(s),:,:); 
    end 
end 

%% Identifying Interesting Sequences (is) 

% 180215 Notes 
    % 1. Need to normalise for chance 
    % 2. https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-data.html#d119e2598
    % 3. https://uk.mathworks.com/help/stats/fitcdiscr.html
    
% Settings 
comps = 100; % number of sequences to identify 

% mRMR Approach 
for er = 1:max(experiment_reps) % for each experiment repeat
    set_token =  find(experiment_reps == er,1,'first'); % settings
    
    if er == 1 % for the WT fish
        clear data tw; data = ...
            double([squeeze(nansum(gCount_norm{1,1}(:,days_crop{set_token}(days{set_token}),i_experiment_reps==er),2))' ; ...
            squeeze(nansum(gCount_norm{1,1}(:,nights_crop{set_token}(nights{set_token}),i_experiment_reps==er),2))']); 
        tw = ones(size(data,1),1); tw(1:size(data,1)/2) = 0; % day vs night data 
    else
%         clear data tw; data = ...
%             double(squeeze(nansum(gCount_merge{1,1}(:,days_crop{set_token}(days{set_token}),i_experiment_reps==er),2) - ...
%             nansum(gCount_merge{1,1}(:,nights_crop{set_token}(nights{set_token}),i_experiment_reps==er),2))'); 
%         tw = i_group_tags(i_experiment_reps == er); 
    end
    
    %data = zscore(data); % z-score data 
    
    tic
    [comps_v{er,1}(:,1)] = mrmr_miq_d(data,tw, comps); 
    toc

    data = data(:,comps_v{er,1}); % reduce to samples x comp sequences 

    tic
    for s = 1:size(comps_v{er,1},1) % for each comp sequence
        Mdl = fitcdiscr(data(:,1:s),tw,'DiscrimType','quadratic','CrossVal','on');
        L(s) = kfoldLoss(Mdl);
    end
    toc

end

clear er set_token data tw 

%% WT Day Night Score 
    % Currently an alternative to the PCA approach for the WT data 
    
% dn_score = cell(1,2); % allocate - high dn score = day, low dn_score = night (1 -> -1)
% dn_score{1,1} = zeros(size(uniqueSeqs{1,1},1),1,'single'); % {t/c} - seqs x 1
% dn_score{1,2} = zeros(size(uniqueSeqs{1,2},1),1,'single'); % {t/c} - seqs x 1
% 
% for tc = 1:2 % for data/control
%     for i = 1:size(dn_score{1,tc},1) % for each sequence
%         dn_score{1,tc}(i,1) = ...
%             sum(nanmean(gFreq_merge{1,tc}(i,days_crop{1}(days{1}),i_experiment_reps==1),3)) - ...
%             sum(nanmean(gFreq_merge{1,tc}(i,nights_crop{1}(nights{1}),i_experiment_reps==1),3));
%     end
% end
% 
% dn_count(:,1) = nanmean(nanmean(gCount_merge{1,1}(:,days_crop{1}(days{1}),i_experiment_reps==1),2),3); 
% dn_count(:,2) = nanmean(nanmean(gCount_merge{1,1}(:,nights_crop{1}(nights{1}),i_experiment_reps==1),2),3); 
% 
% clear exampleSeqs scrap; % number of comparisons to show
% [~,exampleSeqs(1:comps/2)] = maxk(dn_count(:,1).*dn_score{1,1},comps/2); % find most day like sequence 
% [~,scrap] = mink(dn_count(:,2).*dn_score{1,1},comps/2); % find most night like sequences 
% exampleSeqs = [exampleSeqs flip(scrap)']; clear scrap; 
% 
% scatter(dn_count(:,1),dn_score{1,1},36,'markerfacecolor',cmap_2{1}(1,:),...
%     'markeredgecolor',cmap_2{1}(1,:)); 
% hold on; scatter(dn_count(:,2),dn_score{1,1},36,'markerfacecolor',cmap_2{1}(2,:),...
%     'markeredgecolor',cmap_2{1}(2,:)); 
% scatter(dn_count(exampleSeqs(1:comps/2),1),dn_score{1,1}(exampleSeqs(1:comps/2),1),36,'markerfacecolor',cmap_2{1}(1,:),...
%     'markeredgecolor',[1 0.5 0])
% scatter(dn_count(exampleSeqs((comps/2)+1:end),2),dn_score{1,1}(exampleSeqs((comps/2)+1:end),1),36,'markerfacecolor',cmap_2{1}(2,:),...
%     'markeredgecolor',[1 0.5 0])

%% Set WT Interesting Sequences to be based on dn_score/probability 

% comps_v{1,1}(:,1) = exampleSeqs; 

%% Interesting Sequences Figure 

figure; 
for er = 1:max(experiment_reps) % for each group of experiments 
    subplot(2,3,er); % subplot 
    clear scrap;
    scrap = grammar_mat{1,1}(comps_v{er,1}(:,1),:); % grab sequences 
    scrap = scrap(:,1:find(sum(isnan(scrap))~=comps,1,'last')); % trim to longest length 
    ax = imagesc(scrap,'AlphaData',isnan(scrap)==0); % imagesc with nan values in white
    colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap
    set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
    set(ax,'CDataMapping','direct');
    %c = colorbar; c.Label.String = 'Cluster';
    xlabel('Position in Sequence','Fontsize',32);
    ylabel('Sequence','Fontsize',32);
    clear ax grammar_mat_sorted c
end

%% Grammar Comparison Figure  

% Settings
cmap_cluster_merge = [cmap_cluster{2,1} ; cmap_cluster{1,1}]; % merged colormap

for er = 1:max(experiment_reps) % for each group of experiments
    figure; clear scrap;
    counter = 1; % counts subplots
    set_token = find(experiment_reps == er,1,'first'); % settings
    
    for s = comps_v{er,1}' % for each sequence
        
        % Plot Lines
        subplot(2,comps,counter); hold on; set(gca,'FontName','Calibri');
        for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
%                         if er == 1
%                         plot(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
%                             i_experiment_reps == er & i_group_tags == g)),...
%                             'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
%                         end
            errorbar(nanmean(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g)),2),...
                (nanstd(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g))')/sqrt(sum(i_group_tags(i_experiment_reps==er) == g))),...
                'color',cmap{set_token}(g,:),'linewidth',3);
            
            scrap(1,g) = min(nanmean(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g)),2)' - ...
                (nanstd(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g))')/sqrt(sum(i_group_tags(i_experiment_reps==er) == g))));
            scrap(2,g) = max(nanmean(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g)),2)' + ...
                (nanstd(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                i_experiment_reps == er & i_group_tags == g))')/sqrt(sum(i_group_tags(i_experiment_reps==er) == g))));
        end
        
        % Add night patches
        y_lims = [(min(scrap(1,:)) - min(scrap(1,:))*0.05) ...
            (max(scrap(2,:)) + max(scrap(2,:))*0.05)]; % Add a bit of space either side
        
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
        ylabel({'Mean' ; 'Sequence Counts'},'Fontsize',12); % Y Labels
        
        if counter == comps
            [~,~,~,~] = legend(geno_list{set_token}.colheaders,'Location','best'); % Generate axis
            legend('boxoff'); % Turn legend off
        end
        
        % Plot Sequence - from randomly chosen fish
        subplot(2,comps,counter + comps); hold on; set(gca,'FontName','Calibri');
        if er == 1 % for the WT fish
            exampleFish = datasample(find(nansum(squeeze(gFreq_merge{1,1}(s,:,i_experiment_reps == er))) == 1),...
                1,'replace',false); % find a fish who uses this sequence
        else
            exampleFish = datasample(find(nansum(squeeze(gFreq_merge{1,1}(s,:,i_experiment_reps == er))) == 1),...
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

clear ax er set_token g scrap counter seqs

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

%% Localising IS in Time 

for er = 1:max(experiment_reps) % for each group of experiments
    for s = 1:comps
        timeSeqs{er,s} = sparse(size(raw_data{er,1},1),size(raw_data{er,1},2));
        % experiment reps x is {fish x frames)
    end
end

clear scrap;
for er = 1:max(experiment_reps) % for each group of experiments
    counter_2 = 1; % start a counter (counts sequences)
    
    for s = comps_v{er,1}'  % for each interesting sequence
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
for s = 1:comps
    subplot(1,comps,s); hold on;
    clear scrap; scrap = full(timeSeqs{er,s}); 
    for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
        plot(smooth(sum(scrap(i_group_tags(i_experiment_reps==er)==g,:))/...
            sum(sum(timeSeqs{er,s}(i_group_tags(i_experiment_reps==er)==g,:))),25*60*15),...
            'color',cmap{set_token}(g,:)); 
        
    end
end

%% PlotSpread of an IS. 

figure; clear scrap;
set_token = find(experiment_reps == er,1,'first'); % settings
hold on; set(gca,'FontName','Calibri'); set(gca,'Fontsize',32);
col = 1;

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    
    for t = 1:2
        spread_cols = plotSpread(squeeze(gCount_merge{1,1}(s,time_window{set_token}(t),...
            i_experiment_reps == er & i_group_tags == g)),'xValues',g,...
            'distributionColors',cmap_2{set_token}(col,:),'showMM',2);
        set(findall(gca,'type','line'),'markersize',30); % change marker sizes
        spread_cols{2}.LineWidth = 6; spread_cols{2}.Color = 'k'; % Change Mean properties
        spread_cols{2}.MarkerSize = 24;
        col = col + 1;
    end
    
end
set(gca,'XTick',1:g);
set(gca,'XTickLabels',geno_list{set_token}.colheaders)
ylabel({'Mean' ; 'Sequence Counts'},'Fontsize',32); % Y Labels
xlabel('Melatonin','Fontsize',32); 

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

%% Old Working 


%% Grammar in time - WT Example Sequences  

% comps_p = 3; % number of comparisons to show
% 
% % randomly choose a fish who all both sequences 
% exampleFish = datasample(find(sum(common_table(exampleSeqs,i_experiment_reps == 1))==...
%     size(exampleSeqs,2)),1,'replace',false); 
% 
% % Figure 
% %figure; cmap_cluster_merge = [cmap_cluster{2,1} ; cmap_cluster{1,1}]; 
% for e = 1:size(exampleSeqs,2) % for each example sequence 
%     subplot(2,comps,e);     
%     Tlocs = datasample(strfind(threads{exampleFish,1,1}',uniqueSeqs{1,1}{exampleSeqs(e),1}),...
%         1,'replace',false); 
%     Rlocs = threads{exampleFish,2,1}(Tlocs,1):...
%         threads{exampleFish,2,1}(Tlocs+(size(uniqueSeqs{1,1}{exampleSeqs(e),1},2)-1),2); 
%     Rlocs = Rlocs + offset(1,i_experiment_tags(exampleFish)); 
%     
%     ax = imagesc([1,size(Rlocs,2)],[0,max(raw_data{1,1}(exampleFish,Rlocs))],...
%         repmat(states{1,1}(exampleFish,Rlocs),[max(raw_data{1,1}(exampleFish,Rlocs)),1]));
%     hold on; colormap(cmap_cluster_merge); set(gca,'Ydir','Normal'); set(ax,'CDataMapping','direct'); 
%     plot(raw_data{1,1}(exampleFish,Rlocs),'k','linewidth',3); hold on; 
%     axis([1 size(Rlocs,2) 0 (max(raw_data{1,1}(exampleFish,Rlocs)) + 0.5)]);
%     set(gca,'XTickLabels',xticks/fps{1})
%     box off; set(gca,'FontName','Calibri'); set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
%     xlabel('Time (Seconds)'); 
%     ylabel('Delta Px'); 
%     title(horzcat('dn Score = ',num2str(dn_score{1,1}(exampleSeqs(e),1)))); 
% end 
% 
% clear ax exampleFish 