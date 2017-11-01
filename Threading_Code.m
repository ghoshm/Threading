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

% Old 
% load('F:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\WT_Data.mat');
% clear Sleep Wake

% Temp 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\WT_Data.mat', 'fish_tags')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\WT_Data.mat', 'idx_numComp_sorted')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\WT_Data.mat', 'wake_cells')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\WT_Data.mat', 'sleep_cells')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\WT_Data.mat', 'parameter_indicies')

cmap(1,:) = [135 206 250]/255; % light sky blue
cmap_2(1,:) = cmap;
cmap_2(2,:) = [25 25 112]/255; % midnight blue 

% Day vs Night Score
days = [1 3 5 7];
nights = [2 4 6];
%% Threading 

% Pre-allocation
threads = cell(max(fish_tags{1,1}),3,2); % fish x (clusters,times,time windows) x samples vs controls 
idx_numComp_sorted{1,1} = idx_numComp_sorted{1,1} + max(idx_numComp_sorted{2,1}); 
    % Assign higher numbers to the wake bouts 

for f = 1:max(fish_tags{1,1}) % For each fish 
    
    % Pre-allocate
    % Data
    threads{f,1,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % clusters
    threads{f,2,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),2,'single'); % times (start and stop)
    threads{f,3,1} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % time windows
    
    % Control 
    threads{f,1,2} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % clusters    
    threads{f,3,2} = nan(size(find(fish_tags{1,1} == f),1) + ...
        size(find(fish_tags{2,1} == f),1),1,'single'); % time windows
    
    % Deterime starting state 
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
    disp(horzcat('Threading Fish ',num2str(f),' of ',...
        num2str(max(fish_tags{1,1})))); % Report progress 
end 

clear f a b 

%% Filling in states x time 

% Time off-set 
[~,c] = find(lb == max(lb(2,:))); % Find the longest start window 

for e = 1:max(i_experiment_tags) % For each experiment 
    offset(e) = lb(2,c) - lb(2,e); % Calculate the offset  
end 

% Check lb alignment 
figure; hold on; cmap_scrap = hsv(max(i_experiment_tags));
for e = 1:max(i_experiment_tags) % for each experiment 
    plot([lb(:,e)+offset(e),lb(:,e)+offset(e)],...
        [0 1],'color',cmap_scrap(e,:),'linewidth',3)
end 
clear cmap_scrap

% Pre-allocate 
states = nan(max(fish_tags{1,1}),max(lb(end,:)),'single'); % fish x time 
raw_data = nan(max(fish_tags{1,1}),max(lb(end,:)),'single'); % fish x time 

for f = 1:max(fish_tags{1,1}) % for each fish 
    
    for b = 1:size(threads{f,1,1},1) % for each bout 
        states(f,threads{f,2,1}(b,1)+offset(i_experiment_tags(f)):...
            threads{f,2,1}(b,2)+offset(i_experiment_tags(f))) = ...
                threads{f,1,1}(b,1); % Fill in cluster number  
    end
    
    disp(horzcat('Filled States for fish ',num2str(f),' of ',...
        num2str(max(fish_tags{1,1})))); % Report progress 
end 

counter = 1; % start counter
for e = 1:max(i_experiment_tags) % for each experiment
    for f = 1:size(delta_px_sq{e},2) % for each fish
        raw_data(counter,offset(e)+1:size(delta_px_sq{e},1)+offset(e)) =...
            delta_px_sq{1,e}(:,f);
        
        disp(horzcat('Filled raw data for fish ',num2str(counter),' of ',...
            num2str(max(fish_tags{1,1})))); % Report progress
        
        counter = counter + 1; % add to counter
        
    end
end

lb = lb(:,c); 
clear c e f delta_px_sq 

%% States Matrix 
    % Note (170930 - would be cool to split into sleep/wake subplots) 
figure; 
imAlpha = ones(size(states),'single'); imAlpha(isnan(states))=0;
imagesc(states,'AlphaData',imAlpha); 
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); 
c = colorbar; c.Label.String = 'Cluster'; 
set(gca,'XTick',[]); 
set(gca,'Fontsize',18);
xlabel('Time (frames)','Fontsize',32); 
ylabel('Fish ID','Fontsize',32);
 
%% Distribution in time Figure 
all_states = double(max(idx_numComp_sorted{1,1})); % possible states  
    % Note that subplot can't take single values  
ai_states(1:all_states) = 2; % specify if states are active or inactive 
ai_states(min(idx_numComp_sorted{1,1}):end) = 1; % specify active states  
ai_string = {'Active' ; 'Inactive'}; % define states 
time_bins = fps*60; % set smoothing window 

figure; 
ai = 1; % start counters 
for s = 1:all_states % for each state
    subplot(5,4,s); hold on; clear scrap; 
    title(horzcat(ai_string{ai_states(s)},' - ',num2str(ai))); 
    
    for g = 1:max(i_group_tags) % For each group
        legend_lines(1,g) = plot(smooth(sum(states...
            (i_group_tags == g,:) == s)./...
            sum(isnan(states(i_group_tags == g,:))==0),time_bins),...
            'color',cmap(g,:));
        
        % Find the top & bottom
        % Excluding values that smooth outside of the data
        scrap(1,g) = max(legend_lines(1,g).YData((1+time_bins):...
            (size(legend_lines(1,g).YData,2) - time_bins)));
        scrap(2,g) = min(legend_lines(1,g).YData((1+time_bins):...
            (size(legend_lines(1,g).YData,2) - time_bins)));
        
        % Legend 
        if s == 4
            legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
                num2str(size(find(i_group_tags == g),1)),')');
            % Append the group size to each group name
        end
    end
    
    % Night patches     
    a = 1;
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[lb(nights_crop(nights(n))) ...
            min(scrap(2,:)) - (min(scrap(2,:))*0.05)...
            (lb(nights_crop(nights(n))+1)-1) - lb(nights_crop(nights(n)))...
            max(scrap(1,:)) + (max(scrap(1,:))*0.05)],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1;
    end
    
    % Legend 
    if s == 4 % For the 4th state - add a legend
        [~,~,~,~] = legend(legend_cell,'Location','northwest');
        legend('boxoff');
    end
    
    % Axis etc 
    axis([(1 + time_bins) (size(states,2) - time_bins) ...
        min(scrap(2,:)) - (min(scrap(2,:))*0.05) ...
        max(scrap(1,:)) + (max(scrap(1,:))*0.05)]);
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12);
    xlabel('Time (Days/Nights)','Fontsize',12);
    set(gca, 'XTick', []);
    yticks(round(max(scrap(1,:)),2,'Significant')); 
    yticklabels(num2str(round(max(scrap(1,:)),2,'Significant'))); 
    ylabel('Frequency','Fontsize',12);
    
    % State Counter 
    if s == find(ai_states == 2,1,'last')
        ai = 1;
    else
        ai = ai + 1;
    end
    
end 

%% Bout Shapes 

counter = 1; 
for s = find(ai_states == 1,1,'first'):find(ai_states == 1,1,'last') % for each active bout type 
    clear scrap_wc scrap_ft; 
    
    % Pre-alloacte 
    bouts{1,counter} = zeros(size(find(idx_numComp_sorted{1,1}==s),1),...
        max(wake_cells(idx_numComp_sorted{1,1} == s,3)),'single');
    
    scrap_wc = wake_cells(idx_numComp_sorted{1,1} == s,1:3); % take these wake cells
    scrap_ft = fish_tags{1,1}(idx_numComp_sorted{1,1} == s,1); % take their fish tags
    
    for b = 1:size(bouts{1,counter},1) % for each bout 
        bouts{1,counter}(b,1:scrap_wc(b,3)) = raw_data(scrap_ft(b,1),...
            (scrap_wc(b,1)+offset(i_experiment_tags(scrap_ft(b,1)))):...
            (scrap_wc(b,2)+offset(i_experiment_tags(scrap_ft(b,1))))); 
        % Extract raw data 
    end 
    
    disp(horzcat('Calculated Bout Shapes for Cluster ',num2str(counter),...
        ' of ',num2str(size(find(ai_states == 1),2)))); % Report progress 
    counter = counter + 1; % Add to counter
end 

clear s b counter

%% Bout Shapes Figure 
figure; 
subplot(1,2,1); title('Active Bout Shapes'); hold on;
for b = 1:size(bouts,2) % for each active bout type
    plot((0:size([nanmean(bouts{1,b}) 0],2))/fps,...
        [0 nanmean(bouts{1,b}) 0],'color',cmap_cluster{1,1}(b,:),...
        'linewidth',3)
    
    if b == size(bouts,2) % For the last bout type 
    set(gca,'Fontsize',32); 
    axis([0 10/fps ylim]); 
    set(gca,'XTick',(1/fps):(2/fps):(10/fps)); 
    xticklabels(num2str(((1/fps):(2/fps):(10/fps))')); 
    xlabel('Time (seconds)','Fontsize',32); 
    ylabel('Delta Px (a.u)','Fontsize',32);
    end 
    
end  

subplot(1,2,2); title('Inactive Bout Length'); hold on;
for b = find(ai_states == 2,1,'first'):find(ai_states == 2,1,'last') % for each inactive bout  
    plot([(min(sleep_cells(idx_numComp_sorted{2,1} == b,3))/fps)...
        (max(sleep_cells(idx_numComp_sorted{2,1} == b,3))/fps)],...
        [b b],'color',cmap_cluster{2,1}(b,:),...
        'linewidth',3) % Plot average inactive bout length, note that this is 
        % Already converted to seconds 
        
        if b == max(idx_numComp_sorted{2,1})
            set(gca,'Fontsize',32)
            axis([0 100 ylim]);
            xlabel('Time (seconds)','Fontsize',32);
            ylabel('Cluster Number','Fontsize',32);
        end
end 

%% Bout Features Figure 
figure; 

% Active 
subplot(1,2,1); title('Active Bout Features'); 
a = find(states(2,:) == 13,1,'last'); % Hard chosen bout 
hold on; clear scrap; 
scrap = raw_data(2,a-7:a+4); 
plot(scrap,'color',cmap(1,:),'linewidth',3); 
scatter(7,max(scrap),72,'k','filled'); 
text(7.5,double(max(scrap)),'Maximum','Fontsize',16); 
scatter(8,min(scrap(scrap>0)),72,'k','filled'); 
text(8.5,double(min(scrap(scrap>0))),'Minimum','Fontsize',16); 
plot([6 8], [mean(scrap(scrap > 0)) mean(scrap(scrap > 0))],'-k',...
    'linewidth',3); 
plot([6 6],[mean(scrap(scrap > 0))-0.5 mean(scrap(scrap > 0))+0.5],'k',...
    'linewidth',3); 
plot([8 8],[mean(scrap(scrap > 0))-0.5 mean(scrap(scrap > 0))+0.5],'k',...
    'linewidth',3); 
text(8.5,10.5,{'Length','Mean','Variance','Total'},'Fontsize',16); 
set(gca,'Fontsize',16); 
axis([1 12 0 15.5]); 
xticks(1:2:12); 
set(gca,'XTickLabel',{(1:2:12)/fps}); 
xlabel('Time (seconds)'); 
ylabel('Delta Px (a.u)'); 

% Inactive 
subplot(1,2,2); hold on; title('Inactive Bout Features'); set(gca,'Fontsize',16);  
a = find(states(2,:) == 2,1,'last'); % Hard chosen bout 
hold on; clear scrap;
scrap = raw_data(2,a-18:a+1);
plot(scrap,'color',cmap(1,:),'linewidth',3); 
axis([1 20 ylim]); 
plot([2 19],[7 7],'-k','linewidth',3); 
plot([2 2],[6.5 7.5],'k','linewidth',3); 
plot([19 19],[6.5 7.5],'k','linewidth',3); 
text(10,7.5,'Length','Fontsize',16); 
xticks(1:2:20); 
set(gca,'XTickLabel',{(1:2:20)/fps}); 
xlabel('Time (seconds)'); 
ylabel('Delta Px (a.u)'); 
%% Tsne  

% Downsample data 
sample = []; k = 1000; 
for b = 1:size(bouts,2) % for each active bout type 
    sample = [sample ; datasample(wake_cells(idx_numComp_sorted{1,1}==...
        (b+max(idx_numComp_sorted{2,1})),3:end),k)];
end

% Settings 
    % Reduce to 2 dimensions (for visulaisation) 
    % Start with PCA to your previously determined number of dimensions
    % Try a couple of perplexities
    knee_dim = 2; % define dimensisons to reduce to 
    perplex = [30 300 3000]; % perplexities to try 
    mappedX = cell(size(perplex,2),1);
    counter = 1; % start a counter 
    for per = perplex % For a range of possible perplexities
        mappedX{counter} = tsne(sample,[],2,knee_dim,per);
        counter = counter + 1; % add to counter 
    end
    
%% T-Sne Figure

% Choice of Perplexity
knee_perplexity = 3; 

figure; hold on; 
a = 1; 
for c = 1:size(bouts,2) % For each cluster 
    scatter(mappedX{knee_perplexity,1}(a:a+k-1,1),...
        mappedX{knee_perplexity,1}(a:a+k-1,2),...
        'markerfacecolor',cmap_cluster{1,1}(c,:),...
        'markeredgecolor',cmap_cluster{1,1}(c,:)); hold on; 
    a = a + k;
end 

%% Transitions 

% Pre-allocation
% Hard
steps = [1 2]; % markov steps to take

% Soft
transitions = cell(max(fish_tags{1,1}),size(steps,2),size(threads,3)); % fish x steps
all_states = max(idx_numComp_sorted{1,1}); % possible all_states

for sc = 1:size(threads,3) % For data/controls
    for f = 1:max(fish_tags{1,1}) % For each fish
        for s = 1:size(steps,2) % For each step size
            % Pre-allocate
            transitions{f,s,sc} = zeros(all_states,all_states,time_window(2),'single');
            % states x states x time windows
            
            for t = time_window(1):time_window(2) % For each time window
                transitions{f,s,sc}(:,:,t) = accumarray(...
                    [threads{f,1,sc}(find(threads{f,3,sc} == t,1,'first'):...
                    find(threads{f,3,sc} == t,1,'last')-steps(s),1),...
                    threads{f,1,sc}(find(threads{f,3,sc} == t,1,'first')+steps(s):...
                    find(threads{f,3,sc} == t,1,'last'),1)],...
                    1,[all_states,all_states]); % Accumulate transitions
                transitions{f,s,sc}(:,:,t) = bsxfun(@rdivide,transitions{f,s,sc}(:,:,t),...
                    sum(transitions{f,s,sc}(:,:,t),2)); % Convert to probabilities
            end
            transitions{f,s,sc}(isnan(transitions{f,s,sc})) = 0; % Convert nans to zeros
        end
    end
end
clear f s t sc

%% Sorting transitions 
    % For ease of averaging/stats 
    
% Pre-allocate
sorted_transitions = cell(2,size(steps,2),2); % day/night x steps x data/control
    % Then states x states x fish

    for sc = 1:size(threads,3) % for data/control
        for f = 1:max(fish_tags{1,1}) % For each fish
            for s = 1:size(steps,2) % For each step size
                sorted_transitions{1,s,sc}(:,:,f) = nanmean(transitions{f,s,sc}(:,:,days),3);
                sorted_transitions{2,s,sc}(:,:,f) = nanmean(transitions{f,s,sc}(:,:,nights),3);
            end
        end
    end

%% New Working 
figure; 
subplot(2,2,1); title('Day Mean','Fontsize',18); hold on; 
imagesc(nanmean(sorted_transitions{1,1,1},3) + nanmean(sorted_transitions{1,2,1},3),[0 0.5]); 

subplot(2,2,2); title('Day Shuffled Mean','Fontsize',18); hold on; 
imagesc(nanmean(sorted_transitions{1,1,2},3) + nanmean(sorted_transitions{1,2,2},3),[0 0.5]); 

subplot(2,2,3); title('Night Mean','Fontsize',18); hold on; 
imagesc(nanmean(sorted_transitions{2,1,1},3) + nanmean(sorted_transitions{2,2,1},3),[0 0.5]); 

subplot(2,2,4); title('Night Shuffled Mean','Fontsize',18); hold on; 
imagesc(nanmean(sorted_transitions{2,1,2},3) + nanmean(sorted_transitions{2,2,2},3),[0 0.5]); 

set(gca,'YDir','reverse')
plot([0.5,all_states+0.5],[max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],...
    'k','linewidth',3);
plot([max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],[0.5,all_states+0.5],...
    'k','linewidth',3);
set(gca, 'XTick', [median([min(idx_numComp_sorted{2,1})-0.5...
    max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
    max(idx_numComp_sorted{1,1})+0.5])]);
set(gca,'XTickLabels',{'Inactive','Active'},'Fontsize',12)
xlabel('T_{2}','Fontsize',14); % X Labels
set(gca, 'YTick', [median([min(idx_numComp_sorted{2,1})-0.5...
    max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
    max(idx_numComp_sorted{1,1})+0.5])]);
set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12);
ylabel('T_{1}','Fontsize',14); % X Labels
axis([0.5 max(idx_numComp_sorted{1,1})+0.5 0.5 max(idx_numComp_sorted{1,1})+0.5]) 
c = colorbar; c.Label.String = 'Probability'; c.Label.FontSize = 14;
    
subplot(2,2,3); 
hold on; 
a = plot(squeeze(sorted_transitions{1,1,1}(1,7:end,:)),'color',cmap(1,:)); 
b = plot(squeeze(sorted_transitions{1,1,2}(1,7:end,:)),'color',[255 128 0]/255); 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16);
xlabel('T_{2} Active Clusters','Fontsize',16);
ylabel('T_{1} = Inactive 1','Fontsize',18);
lgd = legend([a(1) b(1)],'Data','Shuffled Data','Location','Northeast');
lgd.LineWidth = 300; 
legend('boxoff');

subplot(2,2,4); 
hold on; 
a = plot(squeeze(sorted_transitions{1,1,1}(2,7:end,:)),'color',cmap(1,:)); 
b = plot(squeeze(sorted_transitions{1,1,2}(2,7:end,:)),'color',[255 128 0]/255); 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',16);
xlabel('T_{2} Active Clusters','Fontsize',16);
ylabel('T_{1} = Inactive 2','Fontsize',18);
lgd = legend([a(1) b(1)],'Data','Shuffled Data','Location','Northeast'); 
lgd.LineWidth = 3; 
legend('boxoff');

%% Statistics - chance Mask 
    
chance_mask = zeros(all_states,all_states,'single');
scrap = []; d = 0;
for e = 1:max(i_experiment_tags) % for each experiment
    for g = 1:max(i_group_tags) % for each group
        % for each fish
        for t = time_window(1):time_window(2) % for each time window
            for s = 1:size(transitions,2) % for each step size
                for trans = 1:max(idx_numComp_sorted{1,1}) % for each transition
                    for cs = 1:size(transitions,3) % for data vs control
                        for f = find(i_experiment_tags == e & i_group_tags == g)'
                            if s == 1
                                if trans <= max(idx_numComp_sorted{2,1})
                                    scrap = [scrap ;...
                                        transitions{f,s,cs}(trans,min(idx_numComp_sorted{1,1}):end,t)];
                                else
                                    scrap = [scrap ;...
                                        transitions{f,s,cs}(trans,1:max(idx_numComp_sorted{2,1}),t)];
                                end
                            else
                                if trans <= max(idx_numComp_sorted{2,1})
                                    scrap = [scrap ;...
                                        transitions{f,s,cs}(trans,1:max(idx_numComp_sorted{2,1}),t)];
                                else
                                    scrap = [scrap ;...
                                        transitions{f,s,cs}(trans,min(idx_numComp_sorted{1,1}):end,t)];
                                end
                            end
                        end
                    end
                    
                    anova_group(1:size(scrap,1)) = 1;
                    anova_group(1:size(scrap,1)/2) = 0;
                    
                    try
                        d = manova1(scrap,anova_group);
                    catch
                    end
                    
                    if d >= 1
                        if s == 1
                            chance_mask(trans,min(idx_numComp_sorted{1,1}):end)...
                                = 1;
                        else
                            chance_mask(trans,1:max(idx_numComp_sorted{2,1}))...
                                = 1;
                        end
                    end
                    scrap = []; d = 0;
                    clear anova_group;
                end
            end
        end
    end
end

%% Statistics Working 
scrap = [squeeze(sorted_transitions{1,2,1}(16,7:18,:))...
    squeeze(sorted_transitions{1,2,2}(16,7:18,:))]'; 

anova_group(1:size(scrap,1)) = 1; 
anova_group(1:size(scrap,1)/2) = 0; 

d = manova1(scrap,anova_group)

subplot(1,2,1)
imagesc([nanmean(sorted_transitions{1,1,1},3) + nanmean(sorted_transitions{1,2,1},3)])
colorbar
subplot(1,2,2)
imagesc([nanmean(sorted_transitions{1,1,2},3) + nanmean(sorted_transitions{1,2,2},3)])

%% Statistics 
% Data
data_vec = [];
for f = 1:size(transitions,1) % For each fish
    data_vec = [data_vec reshape(transitions{f,1} + transitions{f,2},...
        [1,(all_states^2)*size([days nights],2)])];
end

% Anova parameters 
anova_group = [];
anova_experiment = []; 
anova_time = []; 
at = ones(1,size(transitions{1,1},3));
at(2:2:size(at,2)) = 2; 
for f = 1:size(transitions,1) % For each fish 
    anova_group = [anova_group...
        repmat(i_group_tags(f),[1,size(transitions{f,1},3)])]; 
    anova_experiment = [anova_experiment...
        repmat(i_experiment_tags(f),[1,size(transitions{f,1},3)])]; 
    anova_time = [anova_time at];
end 
clear at; 

% Calculation 
for t = 1:all_states^2 % For each transition
    clear anova_vec;
    anova_vec = data_vec(t:all_states^2:end); % Grab data for this transition 
    [twa.st.p(:,t),~,twa.st.stats{t}] = anovan(anova_vec,...
                {anova_group,anova_time,anova_experiment},...
                'display','off','model','full');
end

%% Image SC Figure - With joined two all_states  

% Chance Filter 
chance_filter = cell(2,1); 
for t = 1:2 % For day / night 
    clear scrap 
    chance_filter{t,1} = nanmean(sorted_transitions{t,1},3) + ...
        nanmean(sorted_transitions{t,2},3); 
    chance_filter{t,1}(chance_filter{t,1}(:,1:max(idx_numComp_sorted{2,1}))...
        < 1/max(idx_numComp_sorted{2,1})) = NaN;
    scrap = chance_filter{t,1}(:,max(idx_numComp_sorted{2,1})+1:end); 
    scrap(scrap < 1/(max(idx_numComp_sorted{1,1}) - max(idx_numComp_sorted{2,1})))...
        = NaN; 
    chance_filter{t,1}(:,max(idx_numComp_sorted{2,1})+1:end) = scrap; 
end 
chance_filter_logical = isnan(chance_filter{1,1}) + isnan(chance_filter{2,1}); 
chance_filter_logical(chance_filter_logical == 2) = NaN;  
chance_filter_logical(isnan(chance_filter_logical) == 0) = 1; 

% Figure 
top = 0.4; % Hard coded 

figure;
for t = 1:2
    ax(t) = subplot(1,3,t); hold on; set(gca,'Layer','top'); box on;
    imAlpha=ones(size(chance_filter{t,1})); imAlpha(isnan(chance_filter{t,1}))=0;
    imagesc(chance_filter{t,1},'AlphaData',imAlpha,...
        [1/(max(idx_numComp_sorted{1,1}) - max(idx_numComp_sorted{2,1})) top]); % remove NaN values
    axis([0.5 all_states+0.5 0.5 all_states+0.5])
    set(gca,'YDir','reverse')
    plot([0.5,all_states+0.5],[max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],...
        'k','linewidth',3);
    plot([max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],[0.5,all_states+0.5],...
        'k','linewidth',3);
    
    set(gca, 'XTick', [median([min(idx_numComp_sorted{2,1})-0.5...
        max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
        max(idx_numComp_sorted{1,1})+0.5])]);
    set(gca,'XTickLabels',{'Inactive','Active'},'Fontsize',12)
    xlabel('T_{2}','Fontsize',14); % X Labels
    
    set(gca, 'YTick', [median([min(idx_numComp_sorted{2,1})-0.5...
        max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
        max(idx_numComp_sorted{1,1})+0.5])]);
    set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12)
    ylabel('T_{1}','Fontsize',14); % X Labels
    c = colorbar; c.Label.String = 'Probability'; c.Label.FontSize = 14; 

end

ax(1).Title.String = 'Day'; ax(1).Title.FontSize = 18; 
ax(2).Title.String = 'Night'; ax(2).Title.FontSize = 18;

%% Filtered ImageSC

difference = (nanmean(sorted_transitions{1,1},3) + nanmean(sorted_transitions{1,2},3))...
    - (nanmean(sorted_transitions{2,1},3) + nanmean(sorted_transitions{2,2},3)); 
difference = difference(:)'; 
found = find(twa.st.p(2,:) < 0.05 & twa.st.p(6,:) > 0.05); % filter 
clear scrap; scrap = nan(1,size(twa.st.p,2)); % nan vector 
scrap(found) = difference(found); 
scrap = reshape(scrap,[all_states all_states]); % reshape to match transitions matrix 

% Figure
subplot(1,3,3); hold on; set(gca,'Layer','top'); box on; set(gca,'Fontsize',18); 
title('Day - Night Transitions','Fontsize',18); 
imAlpha=ones(size(scrap)); imAlpha(isnan(scrap))=0; % Filter for sig
imAlpha(isnan(chance_filter_logical)) = 0; % Filter for chance 
imagesc(scrap,'AlphaData',imAlpha,[-0.02 0.02]); % remove NaN values 
c = colorbar; c.Label.String = 'Probability (difference)'; c.Label.FontSize = 14;
% Details 
axis([0.5 all_states+0.5 0.5 all_states+0.5])
set(gca,'YDir','reverse')
plot([0.5,all_states+0.5],[max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],...
'k','linewidth',3);
plot([max(idx_numComp_sorted{2,1})+0.5,max(idx_numComp_sorted{2,1})+0.5],[0.5,all_states+0.5],...
'k','linewidth',3);
set(gca, 'XTick', [median([min(idx_numComp_sorted{2,1})-0.5...
max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
max(idx_numComp_sorted{1,1})+0.5])]);
set(gca,'XTickLabels',{'Inactive','Active'},'Fontsize',12)
xlabel('T_{2}','Fontsize',14); % X Labels
set(gca, 'YTick', [median([min(idx_numComp_sorted{2,1})-0.5...
max(idx_numComp_sorted{2,1})+0.5]) median([min(idx_numComp_sorted{1,1})-0.5...
max(idx_numComp_sorted{1,1})+0.5])]);
set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12)
ylabel('T_{1}','Fontsize',14); % X Labels

%% % of Sig Transitions 
size(find(isnan(chance_filter{1,1}) == 0),1)/all_states^2;
size(find(isnan(chance_filter{2,1}) == 0),1)/all_states^2;
sum(imAlpha(:)')/all_states^2; 


%% Sequitur 

% Info 
    %https://github.com/aexbrown/Behavioural_Syntax
    %http://www.sequitur.info/

    % Run Sequitur 
    tic
    [grammar, compVec, totSavings] = ...
        compressSequenceNFast(threads{1,1,2}', [], 5);
    toc
    
    % Expand Grammer 
    tic 
    [gExp, ~] = expandGrammar(grammar, compVec);
    toc
           
    % Rule Lengths 
    for i = 1:size(gExp,1) % for each rule 
        rule_lengths(i) = size(gExp{i},2);  
    end
    
    rule_counts = histcounts(rule_lengths,1:max(rule_lengths)+1);
    
    figure; 
    plot(rule_counts/nansum(rule_counts,2),'color',cmap_2(1,:),'linewidth',3); 
    axis([2 max(rule_lengths) ylim]); 
    
    % Compression 
    compressability = totSavings/size(threads{1,1,1}',2);
    % gExp{1}{1,end} - most compressive 
    % gExp{find(rule_lengths == 8,1,'last')}{1,end} - most nested 

    % Grammar in Time

% Variables
gTime = zeros(size(gExp,1),max(times),'single'); % time winodws
gI = zeros(size(gExp,1),1,'single'); % instances 

% Loop
    times = threads{1,3,2}'; % time windows
tic 
for g = 1:size(gExp,1) % for each rule
    clear terminalInds ruleLength
    
    % Locate
    [terminalInds, ruleLength] = ...
        getTerminalInds(compVec, gExp, gExp{g}{1});
    
    % Find time windows/Transitions
    gTime(g,:) = ...
        histcounts(times(terminalInds),min(days):max(days)+1)/ruleLength;
    
    % Normalise 
    gTime(g,:) = gTime(g,:)./nansum(gTime(g,:),2); 
    
    disp(horzcat('Located symbol ',num2str(g),' of ',num2str(size(gExp,1))))
end
toc 
    % Tag Patterns as Day or Night 
    dn_score = nansum(gTime(:,days),2) - nansum(gTime(:,nights),2); 
    
    dn_tag = zeros(size(gExp,1),1); 
    dn_tag(dn_score > 0.25) = 1; 
    dn_tag(dn_score < - 0.25) = 2; 
    
    figure; axis([1 7 0 1]); hold on; 
    a = plot(gTime(dn_tag == 1,:)','color',cmap_2(1,:)); 
    b = plot(gTime(dn_tag == 2,:)','color',cmap_2(2,:)); 
    c = plot(gTime(dn_tag == 0,:)','color',[1 0.5 0]); 

    legend([a(1),b(1),c(1)],...
        horzcat('Day - ',num2str(round(size(a,1)/size(gExp,1),2,'Significant'))),...
        horzcat('Night - ',num2str(round(size(b,1)/size(gExp,1),2,'Significant'))),...
        horzcat('Non - ',num2str(round(size(c,1)/size(gExp,1),2,'Significant'))),...
        'location','best');
    
    % Counts Figure 
    figure; 
    plot(gI,'color',cmap_2(1,:),'linewidth',3)
    set(gca,'XScale','log'); % set log axis
    
%% For multiple fish + Shuffled Data  

% Settings 
nMax = 10; % Maximum n-grams to consider 
sMax = 21; % Maximum states (max + 1)  
night_color = [0.9608 0.9608 0.9608]; % For background  
first_night = 2; 

% Allocate 
grammar = cell(size(threads,1),2); % fish x test/control 
compVec = cell(size(threads,1),2); % fish x test/control  
totSavings = zeros(size(threads,1),2); % fish x test/control 
gTermCell = cell(size(threads,1),2); % cell for storing n-grams
uniqueSeqs = cell(1,2); % 1 x test/control 

% Generate a Grammer for each fish + Controls
tic
parfor f = 1:3 % For each fish (size(threads,1))
    for tc = 1:2 % for real vs shuffled data
        % Compress
        [grammar{f,tc}, compVec{f,tc}, totSavings(f,tc)] = ...
            compressSequenceNFast(threads{f,1,tc}',sMax,nMax);
        
        % Get terminal symbols
        gTermCell{f,tc} = getGrammarTerminals(grammar{f,tc});
        disp(horzcat('Compressed Fish ',num2str(f),' of ',...
            '3, tc = ',num2str(tc))); 
    end
end
toc

% Merge all the grammers and keep only unique elements 
for tc = 1:2 % for real vs shuffled data 
    for f = 1:3 %size(gTermCell,1) % for each fish
        % append grammar terminals to the total set
        uniqueSeqs{1,tc} = vertcat(uniqueSeqs{1,tc}, gTermCell{f,tc});
    end
end

% Determine unique sequences from larger set
[uniqueSeqs{1,1}, ~] = countUniqueRows(uniqueSeqs{1,1});
[uniqueSeqs{1,2}, ~] = countUniqueRows(uniqueSeqs{1,2});

% Grammar in Time
% Allocate 
gCount{1,1} = zeros(size(uniqueSeqs{1,1},1),max(days),size(threads,1),'single'); % t/c - uniqueSeqs x time windows x fish 
gCount{1,2} = zeros(size(uniqueSeqs{1,2},1),max(days),size(threads,1),'single'); % t/c - uniqueSeqs x time windows x fish
gFreq{1,1} = zeros(size(uniqueSeqs{1,1},1),max(days),size(threads,1),'single'); % t/c - uniqueSeqs x time windows x fish
gFreq{1,2} = zeros(size(uniqueSeqs{1,2},1),max(days),size(threads,1),'single'); % t/c - uniqueSeqs x time windows x fish

for tc = 1:2 % for real vs shuffled data
    for f = 1:3 % for each fish (size(threads,1))
        for i = 1:size(uniqueSeqs{1,tc},1) % For each sequence
            % Find each sequence and count it's time windows
            gCount{1,tc}(i,:,f) = histcounts(threads{f,3,tc}...
                (strfind(threads{f,1,tc}',uniqueSeqs{1,tc}{i,1})),min(days):max(days)+1);
            % Calculate frequency
            if sum(gCount{1,tc}(i,:,f),2) > 0
                gFreq{1,tc}(i,:,f) = gCount{1,tc}(i,:,f)./sum(gCount{1,tc}(i,:,f));
            end
        end
    end
end

% Tag Patterns as Day or Night
dn_score{1,1} = nan(size(gFreq{1,1},1),size(gFreq{1,1},3),'single'); % t/c - uniqueSeqs x fish 
dn_score{1,2} = nan(size(gFreq{1,2},1),size(gFreq{1,2},3),'single'); % t/c - uniqueSeqs x fish 
dn_tag{1,1} = zeros(size(dn_score{1,1},1),1); % t/c - uniqueSeqs x 1  
dn_tag{1,2} = zeros(size(dn_score{1,2},1),1); % t/c - uniqueSeqs x 1 

for tc = 1:2 % for real vs shuffled data
    for f = 1:3 % for each fish
        dn_score{1,tc}(:,f) = nansum(gFreq{1,tc}(:,days,f),2) - ...
            nansum(gFreq{1,tc}(:,nights,f),2);
    end
    dn_tag{1,tc}(nanmean(dn_score{1,tc},2) > 0.25) = 1;
    dn_tag{1,tc}(nanmean(dn_score{1,tc},2) < - 0.25) = 2;
end

% Figure 
figure;
for tc = 1:2 % for real vs shuffled data
    subplot(1,2,tc); axis([0.5 7.5 0 0.025]); hold on;
    lines{1} = plot(nanmean(gFreq{1,tc}(dn_tag{1,tc} == 1,:,:),3)','color',cmap_2(1,:));
    lines{2} = plot(nanmean(gFreq{1,tc}(dn_tag{1,tc} == 2,:,:),3)','color',cmap_2(2,:));
    lines{3} = plot(nanmean(gFreq{1,tc}(dn_tag{1,tc} == 0,:,:),3)','color',[1 0.5 0]);

    y_lims = ylim; 
    a = 1; night_start = first_night; % Start counters
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
            1 (y_lims(2)-y_lims(1))],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1; night_start = night_start + 2; % Add to counters
    end
     
    if size(lines{2},1) > 0
        legend([lines{1}(1),lines{2}(1),lines{3}(1)],...
            horzcat('Day - ',num2str(round(size(lines{1},1)/size(gFreq{1,tc},1),2,'Significant'))),...
            horzcat('Night - ',num2str(round(size(lines{2},1)/size(gFreq{1,tc},1),2,'Significant'))),...
            horzcat('Non - ',num2str(round(size(lines{3},1)/size(gFreq{1,tc},1),2,'Significant'))),...
            'location','northeast');
    else
        legend([lines{1}(1),lines{3}(1)],...
            horzcat('Day - ',num2str(round(size(lines{1},1)/size(gFreq{1,tc},1),2,'Significant'))),...
            horzcat('Non - ',num2str(round(size(lines{3},1)/size(gFreq{1,tc},1),2,'Significant'))),...
            'location','northeast');
    end
    
    box off; legend boxoff; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format 
    xlabel('Time (Days/Nights)','Fontsize',12); % X Labels 
    set(gca, 'XTick', []); % Turn off X-Ticks 
    ylabel('Mean Frequency','Fontsize',12); % Y Labels 
    
    clear lines y_lims a night_start n r  
end
 
% Sequence lengths 
seq_lengths{1,1} = zeros(size(uniqueSeqs{1,1},1),1,'single'); % tc - seqs x 1 
seq_lengths{1,2} = zeros(size(uniqueSeqs{1,2},1),1,'single'); % tc - seqs x 1 

for tc = 1:2 % for real vs shuffled data
    
    for s = 1:size(uniqueSeqs{1,tc},1) % for each sequence
        seq_lengths{1,tc}(s,1) = size(uniqueSeqs{1,tc}{s,1},2);
    end
end

% Figure 
figure; subplot(1,2,1); box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
title('Data','Fontsize',18); hold on; 
histogram(seq_lengths{1},'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap(1,:));
xlabel('Terminal Length','Fontsize',12); ylabel('Probability','Fontsize',12); % Y Labels
axis([1 nMax+1 0 1]); subplot(1,2,2); box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format 
title('Shuffled Data','Fontsize',18); hold on; 
histogram(seq_lengths{2},'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap(1,:));
xlabel('Terminal Length','Fontsize',12); ylabel('Probability','Fontsize',12); % Y Labels
axis([1 nMax+1 0 1]);

% Compressibility 
compressibility = zeros(size(threads,1),2,'single'); % fish x t/c
for f = 1:3 % for each fish (1:size(threads,1)) 
    compressibility(f,:) = totSavings(f,:)./size(threads{f,1,1},1);  
end 

% figure 
figure; hold on; 
plot(compressibility(1:3,:)',...
    'color',cmap(1,:)+(1-cmap(1,:))*(1-(1/(5)^.5)),'linewidth',1.5);
errorbar(nanmean(compressibility(1:3,:)),...
    nanstd(compressibility(1:3,:))./sqrt(3),'color',cmap(1,:),'linewidth',3); 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
set(gca, 'XTick', [1 2]); % Turn off X-Ticks
set(gca,'XTickLabels',{'Data','Shuffled'}); % X Labels 
ylabel('Compressibility','Fontsize',12); % Y Labels
axis([0.5 2.5 ylim]); 

% "Common-ness" of Grammar 
    % How many fish utilise each grammar sequence 
uniqueSeqs_common{1,1} = nan(size(uniqueSeqs{1,1},1),1,'single'); 
uniqueSeqs_common{1,2} = nan(size(uniqueSeqs{1,2},1),1,'single'); 

for tc = 1:2 % for real vs shuffled data
    for c = 1:size(uniqueSeqs{1,tc},1) % for each fish (1:size(threads,1))
        uniqueSeqs_common{1,tc}(c,1) = size(find(sum(squeeze(gCount{1,tc}(c,:,:)))>0),2)/3;
    end
end

% Figure 
figure; subplot(1,2,1); 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
title('Data','Fontsize',18); hold on;    
histogram(uniqueSeqs_common{1,1},'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap(1,:));
xlabel('Fish (%)','Fontsize',12); ylabel('Probability','Fontsize',12); % Axis Labels
axis([0 1 0 1]); subplot(1,2,2); 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
title('Shuffled Data','Fontsize',18); hold on; 
histogram(uniqueSeqs_common{1,2},'Normalization','probability','EdgeColor','none',...
    'FaceColor',cmap(1,:));
xlabel('Fish (%)','Fontsize',12); ylabel('Probability','Fontsize',12); % Axis Labels
axis([0 1 0 1]);