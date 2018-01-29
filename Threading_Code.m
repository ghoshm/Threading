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
imagesc(states{1,1}(:,find(sum(isnan(states{1,1}))==0,1,'first'):...
    find(sum(isnan(states{1,1}))==0,1,'last'))); % plot times where all experiments have data 
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap  
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
c = colorbar; c.Label.String = 'Cluster'; 
x = 0:(fps{1}*60*60*12):size((find(sum(isnan(states{1,1}))==0,1,'first'):...
find(sum(isnan(states{1,1}))==0,1,'last')),2); % set x ticks every 12 hours 
set(gca,'XTick',x); 
set(gca,'XTickLabel',{(0:size(x,2)-1)*12}); 
xlabel('Time (Hours)','Fontsize',32); 
ylabel('Fish ID','Fontsize',32);
 
clear c x 

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

% Scatter Plot 
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

clear er fish_tags_cm f number counter c sample sample_tags sample_er ...
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

%% Load in Legion Data 
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

grammar_mat = nan(size(uniqueSeqs{1,1},1),max(seq_lengths{1,1}),'single'); % sequences x max length 

for s = 1:size(uniqueSeqs{1,1},1) % for each sequence 
    grammar_mat(s,1:size(uniqueSeqs{1,1}{s,1},2)) = uniqueSeqs{1,1}{s,1}; % fill in sequence  
end 

%% Grammar Matrix Figure 

figure; 
grammar_mat_sorted = flip(sortrows(grammar_mat)); % sort rows of grammar_mat  
imAlpha=ones(size(grammar_mat_sorted),'single'); imAlpha(isnan(grammar_mat_sorted))=0; % find nan values 
imagesc(grammar_mat_sorted,'AlphaData',imAlpha); % imagesc with nan values in white 
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap  
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
c = colorbar; c.Label.String = 'Cluster'; 
xlabel('Position in Sequence','Fontsize',32); 
ylabel('Sequence','Fontsize',32);

clear grammar_mat_sorted c 

%% Compressibility 
% The compressibility of a sequence of uncompressed length l is given by the sum of the savings S 
% at each iteration divided by l (Gomez-Marin et al.,2016) 

compressibility = zeros(size(threads,1),2,'single'); % fish x t/c
for f = 1:size(threads,1) % for each fish 
    compressibility(f,:) = totSavings(f,:)./size(threads{f,1,1},1);
end 

clear f 

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

%% Grammar in Time - dn score (WT)

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

dn_score = cell(1,2); % allocate - high dn score = day, low dn_score = night (1 -> -1)
dn_score{1,1} = zeros(size(uniqueSeqs{1,1},1),1,'single'); % {t/c} - seqs x 1
dn_score{1,2} = zeros(size(uniqueSeqs{1,2},1),1,'single'); % {t/c} - seqs x 1

for tc = 1:2 % for data/control
    for i = 1:size(dn_score{1,tc},1) % for each sequence
        dn_score{1,tc}(i,1) = ...
            sum(nanmean(gFreq_merge{1,tc}(i,days_crop{1}(days{1}),i_experiment_reps == 1),3)) - ...
            sum(nanmean(gFreq_merge{1,tc}(i,nights_crop{1}(nights{1}),i_experiment_reps == 1),3));
    end
end

clear tc f i gFreq gCount

 %% Grammar in Time Figure
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

%% Grammar in time - Example Sequences  

exampleSeqs(1) = find(dn_score{1,1} == max(dn_score{1,1})); % find most day like sequence 
exampleSeqs(2) = find(dn_score{1,1} == min(dn_score{1,1})); % find most night like sequence  

% randomly choose a fish who uses both sequences 
exampleFish = datasample(find(sum(squeeze(gFreq_merge{1,1}(exampleSeqs(1),:,i_experiment_reps == 1))) == 1 & ...
    sum(squeeze(gFreq_merge{1,1}(exampleSeqs(2),:,i_experiment_reps == 1))) == 1),1,'replace',false); 

% Figure 
figure; cmap_cluster_merge = [cmap_cluster{2,1} ; cmap_cluster{1,1}]; 
for e = 1:size(exampleSeqs,2) % for each example sequence 
    subplot(2,1,e);     
    Tlocs = datasample(strfind(threads{exampleFish,1,1}',uniqueSeqs{1,1}{exampleSeqs(e),1}),...
        1,'replace',false); 
    Rlocs = threads{exampleFish,2,1}(Tlocs,1):...
        threads{exampleFish,2,1}(Tlocs+(size(uniqueSeqs{1,1}{exampleSeqs(e),1},2)-1),2); 
    Rlocs = Rlocs + offset(1,i_experiment_tags(exampleFish)); 
    
    imagesc(repmat(states{1,1}(exampleFish,Rlocs),[max(raw_data{1,1}(exampleFish,Rlocs)),1]));
    hold on; colormap(cmap_cluster_merge); set(gca,'Ydir','Normal');  
    plot(raw_data{1,1}(exampleFish,Rlocs),'k','linewidth',3); hold on; 
    set(gca,'XTickLabels',xticks/fps{1})
    box off; set(gca,'FontName','Calibri'); set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    if e == size(exampleSeqs,2)
    xlabel('Time (Seconds)'); 
    end 
    ylabel('Delta Px'); 
    title(horzcat('dn Score = ',num2str(dn_score{1,1}(exampleSeqs(e),1)))); 
end 

clear exampleFish 

%% "Common-ness" of Grammar

% Logical Table 
common_table = zeros(size(uniqueSeqs{1,1},1),size(gFreq_merge{1,1},3),'single'); 

for s = 1:size(common_table,1) % for each sequence 
    common_table(s,:) = sum(squeeze(gFreq_merge{1,1}(s,:,:))); 
end 

%% Common-ness of Grammar Figure 

figure; 
box off; set(gca, 'Layer','top'); hold on; 
set(gca,'FontName','Calibri'); set(gca,'Fontsize',32); % Format   
histogram((sum(common_table,2)/size(common_table,2))*100,'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap{1}(1,:));
xlabel('% of Fish','Fontsize',32); ylabel('Probability','Fontsize',32); % Axis Labels
axis([0 100 0 0.05]); 

%% Statistical Comparisons of Counts 

% Allocate 
twa.gs = cell(1,max(experiment_reps)); 

for er = 1:max(experiment_reps) % for each group of experiments
    set_token =  find(experiment_reps == er,1,'first'); % settings 
    
    % Grouping variables
    anova_group = repmat(i_group_tags(i_experiment_reps==er),...
        [size([days{set_token} nights{set_token}],2),1])'; % groups
    anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
        [size([days{set_token} nights{set_token}],2),1])'; % experiments
    anova_time = [];
    for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
        anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
        % Allocate alternating zeros and ones to each time window
    end
    anova_time = anova_time';
    
    % Statistical Comparisons 
    tic
    for s = 1:size(uniqueSeqs{1,1},1) % for each sequence
        
        data = squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
            i_experiment_reps == er))';
        data = data(:)'; % vectorise
        
        [twa.gs{1,er}(s,:)] = single(anovan(data,...
            {anova_group,anova_time,anova_experiment},...
            'display','off','model','full'));
        
        % Allocate 
        if s == 1 % for the first sequence 
           twa.gs{1,er}(2:size(uniqueSeqs{1,1},1),1:size(twa.gs{1,er},2)) = single(NaN); 
           %twa.gs{2,er} = [twa.gs{2,er} ; cell(size(uniqueSeqs{1,1},1)-1,1)]; 
        end 
        
        % Report Progress 
        if mod(s,1000) == 0
           disp(horzcat('Completed ',num2str(s),' Comparisons of ',...
               num2str(size(uniqueSeqs{1,1},1)),' er = ',num2str(er)));
        end 
    end
    toc 
    clear set_token anova_group anova_experiment anova_time data 
end

%% Correct P-values for multiple comparisons 

%% Grammar Comparison Figure 

% Settings 
comps = 5; % Set the number of sequences to show 

for er = 1:max(experiment_reps) % for each group of experiments
    figure; clear scrap;
    counter = 1; % counts subplots  
    set_token = find(experiment_reps == er,1,'first'); % settings
    
    % Find the top comp significant sequences 
    [~,comps_v(:,1)] = mink(twa.gs{1,er}(:,1),comps);
    %[~,comps_v(:,2)] = mink(twa.gs{1,er}(:,4),comps); 
    
    for c = 1:size(comps_v,2) % for each variable
        for s = comps_v(:,c)' % for each sequence
            subplot(2,5,counter); hold on; set(gca,'FontName','Calibri');
            counter_2 = 1; % counts groups 
            
            for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
                plot([counter_2,counter_2+1],squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                    i_experiment_reps == er & i_group_tags == g)),...
                    'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
                errorbar([counter_2,counter_2+1],nanmean(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                    i_experiment_reps == er & i_group_tags == g)),2),...
                    nanstd(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                    i_experiment_reps == er & i_group_tags == g))'),...
                    'color',cmap{set_token}(g,:),'linewidth',3);
                
                counter_2 = counter_2 + 2; % add to counter
                
                scrap(1,g) = min(min(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                    i_experiment_reps == er & i_group_tags == g))));
                scrap(2,g) = max(max(squeeze(gCount_merge{1,1}(s,time_window{set_token}(1):time_window{set_token}(2),...
                    i_experiment_reps == er & i_group_tags == g))));
            end
            
            box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
            set(gca,'XTickLabels',geno_list{set_token}.colheaders); % X Labels
            ylabel('Counts','Fontsize',12); % Y Labels
            axis([0.5 (max(i_group_tags(i_experiment_reps == er))*2)+.5 ...
                (min(scrap(1,:)) - (min(scrap(1,:))*0.05)) (max(scrap(2,:)) + (max(scrap(2,:))*0.05))]);
            counter = counter + 1; % add to counter
            
        end
    end

end

clear er set_token g scrap counter 

%% Old Working 

% %% Bout Features Figure 
% figure; 
% 
% % Active 
% subplot(1,2,1); title('Active Bout Features'); 
% a = find(states(2,:) == 13,1,'last'); % Hard chosen bout 
% hold on; clear scrap; 
% scrap = raw_data(2,a-7:a+4); 
% plot(scrap,'color',cmap(1,:),'linewidth',3); 
% scatter(7,max(scrap),72,'k','filled'); 
% text(7.5,double(max(scrap)),'Maximum','Fontsize',16); 
% scatter(8,min(scrap(scrap>0)),72,'k','filled'); 
% text(8.5,double(min(scrap(scrap>0))),'Minimum','Fontsize',16); 
% plot([6 8], [mean(scrap(scrap > 0)) mean(scrap(scrap > 0))],'-k',...
%     'linewidth',3); 
% plot([6 6],[mean(scrap(scrap > 0))-0.5 mean(scrap(scrap > 0))+0.5],'k',...
%     'linewidth',3); 
% plot([8 8],[mean(scrap(scrap > 0))-0.5 mean(scrap(scrap > 0))+0.5],'k',...
%     'linewidth',3); 
% text(8.5,10.5,{'Length','Mean','Variance','Total'},'Fontsize',16); 
% set(gca,'Fontsize',16); 
% axis([1 12 0 15.5]); 
% xticks(1:2:12); 
% set(gca,'XTickLabel',{(1:2:12)/fps}); 
% xlabel('Time (seconds)'); 
% ylabel('Delta Px (a.u)'); 
% 
% % Inactive 
% subplot(1,2,2); hold on; title('Inactive Bout Features'); set(gca,'Fontsize',16);  
% a = find(states(2,:) == 2,1,'last'); % Hard chosen bout 
% hold on; clear scrap;
% scrap = raw_data(2,a-18:a+1);
% plot(scrap,'color',cmap(1,:),'linewidth',3); 
% axis([1 20 ylim]); 
% plot([2 19],[7 7],'-k','linewidth',3); 
% plot([2 2],[6.5 7.5],'k','linewidth',3); 
% plot([19 19],[6.5 7.5],'k','linewidth',3); 
% text(10,7.5,'Length','Fontsize',16); 
% xticks(1:2:20); 
% set(gca,'XTickLabel',{(1:2:20)/fps}); 
% xlabel('Time (seconds)'); 
% ylabel('Delta Px (a.u)'); 
% %% Tsne  
% 
% % Downsample data 
% sample = []; k = 1000; 
% for b = 1:size(bouts,2) % for each active bout type 
%     sample = [sample ; datasample(wake_cells(idx_numComp_sorted{1,1}==...
%         (b+max(idx_numComp_sorted{2,1})),3:end),k)];
% end
% 
% % Settings 
%     % Reduce to 2 dimensions (for visulaisation) 
%     % Start with PCA to your previously determined number of dimensions
%     % Try a couple of perplexities
%     knee_dim = 2; % define dimensisons to reduce to 
%     perplex = [30 300 3000]; % perplexities to try 
%     mappedX = cell(size(perplex,2),1);
%     counter = 1; % start a counter 
%     for per = perplex % For a range of possible perplexities
%         mappedX{counter} = tsne(sample,[],2,knee_dim,per);
%         counter = counter + 1; % add to counter 
%     end
%     
% %% T-Sne Figure
% 
% % Choice of Perplexity
% knee_perplexity = 3; 
% 
% figure; hold on; 
% a = 1; 
% for c = 1:size(bouts,2) % For each cluster 
%     scatter(mappedX{knee_perplexity,1}(a:a+k-1,1),...
%         mappedX{knee_perplexity,1}(a:a+k-1,2),...
%         'markerfacecolor',cmap_cluster{1,1}(c,:),...
%         'markeredgecolor',cmap_cluster{1,1}(c,:)); hold on; 
%     a = a + k;
% end 
% 
% %% Transitions 
% 
% % Pre-allocation
% % Hard
% steps = [1 2]; % markov steps to take
% 
% % Soft
% transitions = cell(max(fish_tags{1,1}),size(steps,2),size(threads,3)); % fish x steps
% all_states = max(idx_numComp_sorted{1,1}); % possible all_states
% 
% for sc = 1:size(threads,3) % For data/controls
%     for f = 1:max(fish_tags{1,1}) % For each fish
%         for s = 1:size(steps,2) % For each step size
%             % Pre-allocate
%             transitions{f,s,sc} = zeros(all_states,all_states,time_window(2),'single');
%             % states x states x time windows
%             
%             for t = time_window(1):time_window(2) % For each time window
%                 transitions{f,s,sc}(:,:,t) = accumarray(...
%                     [threads{f,1,sc}(find(threads{f,3,sc} == t,1,'first'):...
%                     find(threads{f,3,sc} == t,1,'last')-steps(s),1),...
%                     threads{f,1,sc}(find(threads{f,3,sc} == t,1,'first')+steps(s):...
%                     find(threads{f,3,sc} == t,1,'last'),1)],...
%                     1,[all_states,all_states]); % Accumulate transitions
%                 transitions{f,s,sc}(:,:,t) = bsxfun(@rdivide,transitions{f,s,sc}(:,:,t),...
%                     sum(transitions{f,s,sc}(:,:,t),2)); % Convert to probabilities
%             end
%             transitions{f,s,sc}(isnan(transitions{f,s,sc})) = 0; % Convert nans to zeros
%         end
%     end
% end
% clear f s t sc
% 
% %% Sorting transitions 
%     % For ease of averaging/stats 
%     
% % Pre-allocate
% sorted_transitions = cell(2,size(steps,2),2); % day/night x steps x data/control
%     % Then states x states x fish
% 
%     for sc = 1:size(threads,3) % for data/control
%         for f = 1:max(fish_tags{1,1}) % For each fish
%             for s = 1:size(steps,2) % For each step size
%                 sorted_transitions{1,s,sc}(:,:,f) = nanmean(transitions{f,s,sc}(:,:,days),3);
%                 sorted_transitions{2,s,sc}(:,:,f) = nanmean(transitions{f,s,sc}(:,:,nights),3);
%             end
%         end
%     end
% 
% %% New Working 
% figure; 
% subplot(2,2,1); title('Day Mean','Fontsize',18); hold on; 
% imagesc(nanmean(sorted_transitions{1,1,1},3) + nanmean(sorted_transitions{1,2,1},3),[0 0.5]); 
% 
% subplot(2,2,2); title('Day Shuffled Mean','Fontsize',18); hold on; 
% imagesc(nanmean(sorted_transitions{1,1,2},3) + nanmean(sorted_transitions{1,2,2},3),[0 0.5]); 
% 
% subplot(2,2,3); title('Night Mean','Fontsize',18); hold on; 
% imagesc(nanmean(sorted_transitions{2,1,1},3) + nanmean(sorted_transitions{2,2,1},3),[0 0.5]); 
% 
% subplot(2,2,4); title('Night Shuffled Mean','Fontsize',18); hold on; 
% imagesc(nanmean(sorted_transitions{2,1,2},3) + nanmean(sorted_transitions{2,2,2},3),[0 0.5]); 
% 
% set(gca,'YDir','reverse')
% plot([0.5,all_states+0.5],[max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],...
%     'k','linewidth',3);
% plot([max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],[0.5,all_states+0.5],...
%     'k','linewidth',3);
% set(gca, 'XTick', [median([min(idx_numComp_sorted{2,1})-0.5...
%     max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
%     max(idx_numComp_sorted{1,1})+0.5])]);
% set(gca,'XTickLabels',{'Inactive','Active'},'Fontsize',12)
% xlabel('T_{2}','Fontsize',14); % X Labels
% set(gca, 'YTick', [median([min(idx_numComp_sorted{2,1})-0.5...
%     max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
%     max(idx_numComp_sorted{1,1})+0.5])]);
% set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12);
% ylabel('T_{1}','Fontsize',14); % X Labels
% axis([0.5 max(idx_numComp_sorted{1,1})+0.5 0.5 max(idx_numComp_sorted{1,1})+0.5]) 
% c = colorbar; c.Label.String = 'Probability'; c.Label.FontSize = 14;
%     
% subplot(2,2,3); 
% hold on; 
% a = plot(squeeze(sorted_transitions{1,1,1}(1,7:end,:)),'color',cmap(1,:)); 
% b = plot(squeeze(sorted_transitions{1,1,2}(1,7:end,:)),'color',[255 128 0]/255); 
% box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16);
% xlabel('T_{2} Active Clusters','Fontsize',16);
% ylabel('T_{1} = Inactive 1','Fontsize',18);
% lgd = legend([a(1) b(1)],'Data','Shuffled Data','Location','Northeast');
% lgd.LineWidth = 300; 
% legend('boxoff');
% 
% subplot(2,2,4); 
% hold on; 
% a = plot(squeeze(sorted_transitions{1,1,1}(2,7:end,:)),'color',cmap(1,:)); 
% b = plot(squeeze(sorted_transitions{1,1,2}(2,7:end,:)),'color',[255 128 0]/255); 
% box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16);
% xlabel('T_{2} Active Clusters','Fontsize',16);
% ylabel('T_{1} = Inactive 2','Fontsize',18);
% lgd = legend([a(1) b(1)],'Data','Shuffled Data','Location','Northeast'); 
% lgd.LineWidth = 3; 
% legend('boxoff');
% 
% %% Statistics - chance Mask 
%     
% chance_mask = zeros(all_states,all_states,'single');
% scrap = []; d = 0;
% for e = 1:max(i_experiment_tags) % for each experiment
%     for g = 1:max(i_group_tags) % for each group
%         % for each fish
%         for t = time_window(1):time_window(2) % for each time window
%             for s = 1:size(transitions,2) % for each step size
%                 for trans = 1:max(idx_numComp_sorted{1,1}) % for each transition
%                     for cs = 1:size(transitions,3) % for data vs control
%                         for f = find(i_experiment_tags == e & i_group_tags == g)'
%                             if s == 1
%                                 if trans <= max(idx_numComp_sorted{2,1})
%                                     scrap = [scrap ;...
%                                         transitions{f,s,cs}(trans,min(idx_numComp_sorted{1,1}):end,t)];
%                                 else
%                                     scrap = [scrap ;...
%                                         transitions{f,s,cs}(trans,1:max(idx_numComp_sorted{2,1}),t)];
%                                 end
%                             else
%                                 if trans <= max(idx_numComp_sorted{2,1})
%                                     scrap = [scrap ;...
%                                         transitions{f,s,cs}(trans,1:max(idx_numComp_sorted{2,1}),t)];
%                                 else
%                                     scrap = [scrap ;...
%                                         transitions{f,s,cs}(trans,min(idx_numComp_sorted{1,1}):end,t)];
%                                 end
%                             end
%                         end
%                     end
%                     
%                     anova_group(1:size(scrap,1)) = 1;
%                     anova_group(1:size(scrap,1)/2) = 0;
%                     
%                     try
%                         d = manova1(scrap,anova_group);
%                     catch
%                     end
%                     
%                     if d >= 1
%                         if s == 1
%                             chance_mask(trans,min(idx_numComp_sorted{1,1}):end)...
%                                 = 1;
%                         else
%                             chance_mask(trans,1:max(idx_numComp_sorted{2,1}))...
%                                 = 1;
%                         end
%                     end
%                     scrap = []; d = 0;
%                     clear anova_group;
%                 end
%             end
%         end
%     end
% end
% 
% %% Statistics Working 
% scrap = [squeeze(sorted_transitions{1,2,1}(16,7:18,:))...
%     squeeze(sorted_transitions{1,2,2}(16,7:18,:))]'; 
% 
% anova_group(1:size(scrap,1)) = 1; 
% anova_group(1:size(scrap,1)/2) = 0; 
% 
% d = manova1(scrap,anova_group)
% 
% subplot(1,2,1)
% imagesc([nanmean(sorted_transitions{1,1,1},3) + nanmean(sorted_transitions{1,2,1},3)])
% colorbar
% subplot(1,2,2)
% imagesc([nanmean(sorted_transitions{1,1,2},3) + nanmean(sorted_transitions{1,2,2},3)])
% 
% %% Statistics 
% % Data
% data_vec = [];
% for f = 1:size(transitions,1) % For each fish
%     data_vec = [data_vec reshape(transitions{f,1} + transitions{f,2},...
%         [1,(all_states^2)*size([days nights],2)])];
% end
% 
% % Anova parameters 
% anova_group = [];
% anova_experiment = []; 
% anova_time = []; 
% at = ones(1,size(transitions{1,1},3));
% at(2:2:size(at,2)) = 2; 
% for f = 1:size(transitions,1) % For each fish 
%     anova_group = [anova_group...
%         repmat(i_group_tags(f),[1,size(transitions{f,1},3)])]; 
%     anova_experiment = [anova_experiment...
%         repmat(i_experiment_tags(f),[1,size(transitions{f,1},3)])]; 
%     anova_time = [anova_time at];
% end 
% clear at; 
% 
% % Calculation 
% for t = 1:all_states^2 % For each transition
%     clear anova_vec;
%     anova_vec = data_vec(t:all_states^2:end); % Grab data for this transition 
%     [twa.st.p(:,t),~,twa.st.stats{t}] = anovan(anova_vec,...
%                 {anova_group,anova_time,anova_experiment},...
%                 'display','off','model','full');
% end
% 
% %% Image SC Figure - With joined two all_states  
% 
% % Chance Filter 
% chance_filter = cell(2,1); 
% for t = 1:2 % For day / night 
%     clear scrap 
%     chance_filter{t,1} = nanmean(sorted_transitions{t,1},3) + ...
%         nanmean(sorted_transitions{t,2},3); 
%     chance_filter{t,1}(chance_filter{t,1}(:,1:max(idx_numComp_sorted{2,1}))...
%         < 1/max(idx_numComp_sorted{2,1})) = NaN;
%     scrap = chance_filter{t,1}(:,max(idx_numComp_sorted{2,1})+1:end); 
%     scrap(scrap < 1/(max(idx_numComp_sorted{1,1}) - max(idx_numComp_sorted{2,1})))...
%         = NaN; 
%     chance_filter{t,1}(:,max(idx_numComp_sorted{2,1})+1:end) = scrap; 
% end 
% chance_filter_logical = isnan(chance_filter{1,1}) + isnan(chance_filter{2,1}); 
% chance_filter_logical(chance_filter_logical == 2) = NaN;  
% chance_filter_logical(isnan(chance_filter_logical) == 0) = 1; 
% 
% % Figure 
% top = 0.4; % Hard coded 
% 
% figure;
% for t = 1:2
%     ax(t) = subplot(1,3,t); hold on; set(gca,'Layer','top'); box on;
%     imAlpha=ones(size(chance_filter{t,1})); imAlpha(isnan(chance_filter{t,1}))=0;
%     imagesc(chance_filter{t,1},'AlphaData',imAlpha,...
%         [1/(max(idx_numComp_sorted{1,1}) - max(idx_numComp_sorted{2,1})) top]); % remove NaN values
%     axis([0.5 all_states+0.5 0.5 all_states+0.5])
%     set(gca,'YDir','reverse')
%     plot([0.5,all_states+0.5],[max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],...
%         'k','linewidth',3);
%     plot([max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],[0.5,all_states+0.5],...
%         'k','linewidth',3);
%     
%     set(gca, 'XTick', [median([min(idx_numComp_sorted{2,1})-0.5...
%         max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
%         max(idx_numComp_sorted{1,1})+0.5])]);
%     set(gca,'XTickLabels',{'Inactive','Active'},'Fontsize',12)
%     xlabel('T_{2}','Fontsize',14); % X Labels
%     
%     set(gca, 'YTick', [median([min(idx_numComp_sorted{2,1})-0.5...
%         max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
%         max(idx_numComp_sorted{1,1})+0.5])]);
%     set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12)
%     ylabel('T_{1}','Fontsize',14); % X Labels
%     c = colorbar; c.Label.String = 'Probability'; c.Label.FontSize = 14; 
% 
% end
% 
% ax(1).Title.String = 'Day'; ax(1).Title.FontSize = 18; 
% ax(2).Title.String = 'Night'; ax(2).Title.FontSize = 18;
% 
% %% Filtered ImageSC
% 
% difference = (nanmean(sorted_transitions{1,1},3) + nanmean(sorted_transitions{1,2},3))...
%     - (nanmean(sorted_transitions{2,1},3) + nanmean(sorted_transitions{2,2},3)); 
% difference = difference(:)'; 
% found = find(twa.st.p(2,:) < 0.05 & twa.st.p(6,:) > 0.05); % filter 
% clear scrap; scrap = nan(1,size(twa.st.p,2)); % nan vector 
% scrap(found) = difference(found); 
% scrap = reshape(scrap,[all_states all_states]); % reshape to match transitions matrix 
% 
% % Figure
% subplot(1,3,3); hold on; set(gca,'Layer','top'); box on; set(gca,'Fontsize',18); 
% title('Day - Night Transitions','Fontsize',18); 
% imAlpha=ones(size(scrap)); imAlpha(isnan(scrap))=0; % Filter for sig
% imAlpha(isnan(chance_filter_logical)) = 0; % Filter for chance 
% imagesc(scrap,'AlphaData',imAlpha,[-0.02 0.02]); % remove NaN values 
% c = colorbar; c.Label.String = 'Probability (difference)'; c.Label.FontSize = 14;
% % Details 
% axis([0.5 all_states+0.5 0.5 all_states+0.5])
% set(gca,'YDir','reverse')
% plot([0.5,all_states+0.5],[max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],...
% 'k','linewidth',3);
% plot([max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],[0.5,all_states+0.5],...
% 'k','linewidth',3);
% set(gca, 'XTick', [median([min(idx_numComp_sorted{2,1})-0.5...
% max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
% max(idx_numComp_sorted{1,1})+0.5])]);
% set(gca,'XTickLabels',{'Inactive','Active'},'Fontsize',12)
% xlabel('T_{2}','Fontsize',14); % X Labels
% set(gca, 'YTick', [median([min(idx_numComp_sorted{2,1})-0.5...
% max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
% max(idx_numComp_sorted{1,1})+0.5])]);
% set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12)
% ylabel('T_{1}','Fontsize',14); % X Labels
% 
% %% % of Sig Transitions 
% size(find(isnan(chance_filter{1,1}) == 0),1)/all_states^2;
% size(find(isnan(chance_filter{2,1}) == 0),1)/all_states^2;
% sum(imAlpha(:)')/all_states^2; 
% 
% 
% %% Sequitur 
% 
% % Info 
%     %https://github.com/aexbrown/Behavioural_Syntax
%     %http://www.sequitur.info/
% 
%     % Run Sequitur 
%     tic
%     [grammar, compVec, totSavings] = ...
%         compressSequenceNFast(threads{1,1,2}', [], 5);
%     toc
%     
%     % Expand Grammer 
%     tic 
%     [gExp, ~] = expandGrammar(grammar, compVec);
%     toc
%            
%     % Rule Lengths 
%     for i = 1:size(gExp,1) % for each rule 
%         rule_lengths(i) = size(gExp{i},2);  
%     end
%     
%     rule_counts = histcounts(rule_lengths,1:max(rule_lengths)+1);
%     
%     figure; 
%     plot(rule_counts/nansum(rule_counts,2),'color',cmap_2(1,:),'linewidth',3); 
%     axis([2 max(rule_lengths) ylim]); 
%     
%     % Compression 
%     compressability = totSavings/size(threads{1,1,1}',2);
%     % gExp{1}{1,end} - most compressive 
%     % gExp{find(rule_lengths == 8,1,'last')}{1,end} - most nested 
% 
%     % Grammar in Time
% 
% % Variables
% gTime = zeros(size(gExp,1),max(times),'single'); % time winodws
% gI = zeros(size(gExp,1),1,'single'); % instances 
% 
% % Loop
%     times = threads{1,3,2}'; % time windows
% tic 
% for g = 1:size(gExp,1) % for each rule
%     clear terminalInds ruleLength
%     
%     % Locate
%     [terminalInds, ruleLength] = ...
%         getTerminalInds(compVec, gExp, gExp{g}{1});
%     
%     % Find time windows/Transitions
%     gTime(g,:) = ...
%         histcounts(times(terminalInds),min(days):max(days)+1)/ruleLength;
%     
%     % Normalise 
%     gTime(g,:) = gTime(g,:)./nansum(gTime(g,:),2); 
%     
%     disp(horzcat('Located symbol ',num2str(g),' of ',num2str(size(gExp,1))))
% end
% toc 
%     % Tag Patterns as Day or Night 
%     dn_score = nansum(gTime(:,days),2) - nansum(gTime(:,nights),2); 
%     
%     dn_tag = zeros(size(gExp,1),1); 
%     dn_tag(dn_score > 0.25) = 1; 
%     dn_tag(dn_score < - 0.25) = 2; 
%     
%     figure; axis([1 7 0 1]); hold on; 
%     a = plot(gTime(dn_tag == 1,:)','color',cmap_2(1,:)); 
%     b = plot(gTime(dn_tag == 2,:)','color',cmap_2(2,:)); 
%     c = plot(gTime(dn_tag == 0,:)','color',[1 0.5 0]); 
% 
%     legend([a(1),b(1),c(1)],...
%         horzcat('Day - ',num2str(round(size(a,1)/size(gExp,1),2,'Significant'))),...
%         horzcat('Night - ',num2str(round(size(b,1)/size(gExp,1),2,'Significant'))),...
%         horzcat('Non - ',num2str(round(size(c,1)/size(gExp,1),2,'Significant'))),...
%         'location','best');
%     
%     % Counts Figure 
%     figure; 
%     plot(gI,'color',cmap_2(1,:),'linewidth',3)
%     set(gca,'XScale','log'); % set log axis
%     
