%% Hour by Hour Threading Code  

% An adaptation of Threading_Code to 
    % Seperatly compress each hour 
    % Count motifs every hour 
    
%% Required scripts 

%% Notes 

%% Settings 
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 
set(0,'defaultfigurecolor',[1 1 1]); % white background

%% Load Data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180324_Hours.mat'); 

%% Chunk Compression 
% Compress cluster chunks

% settings
step = 500;
sMax = max(idx_numComp_sorted{1,1}) + 1; % Maximum states + 1 (first new symbol)
nMax = 10; % Maximum n-grams to consider

% allocate
chunks = cell(size(threads,1),1); % fish x 1
totSavings_cells = cell(size(threads,1),1); % fish x 1

tic
for f = 1:size(threads,1) % for each fish
    chunks{f,1} = find(threads{f,3,1} == 1,1,'first'):...
        step:find(isnan(threads{f,3,1}) == 0,1,'last');
    
    for t = 1:(length(chunks{f,1})-1) % for each chunk
        [~,~, totSavings_cells{f,1}(1,t)] = ...
            compressSequenceNFast(threads{f,1,1}(chunks{f,1}(1,t):(chunks{f,1}(t+1)-1))',...
            sMax,nMax); % compress returning totSavings
    end
    
    chunks{f,1}(end) = []; % remove the spare edge
    
    disp(num2str(f));
    
end
toc 

clear sMax nMax f t

%% Reshape Chunks & Calculate Compression  
    % The compressibility of a sequence of uncompressed length l is given by the sum of the savings S
    % at each iteration divided by l (Gomez-Marin et al.,2016)
    
% settings
tw = 48; % hard coded maximum number of time windows
dn_hour = ones(1,tw)*2; % day and night hours
dn_hour([1:14 25:38]) = 1; % day (1) and night (2) hours

% allocate
compressibility = nan(size(threads,1),tw,'single'); % fish x max time hour bins

% fill data
for f = 1:124 %size(threads,1) % for each fish
    for h = 1:max(threads{1,3,1}) % for each hour
        try
            compressibility(f,h) = nanmean(totSavings_cells{f,1}(...
                threads{f,3,1}(chunks{f,1}) == h)/step); 
        catch
        end
    end
    
end

clear f h 

%% Relative Compressibility Figure
er = 1; 
set_token = find(experiment_reps == er,1,'first'); % settings
figure;
hold on; set(gca,'FontName','Calibri'); clear scrap;

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    clear data;
    data = compressibility(i_experiment_reps == er & i_group_tags == g,:);
    
    if er == 1 % for the WT experiments 
    plot(data',...
        'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
    end 
    
    errorbar(nanmean(data),nanstd(data),...
        'color',cmap{set_token}(g,:),'linewidth',3);
    
    scrap(1,g) = min(data(:));
    scrap(2,g) = max(data(:));
end

y_lims = [(min(scrap(1,:)) - min(scrap(1,:))*0.05) ...
    (max(scrap(2,:)) + max(scrap(2,:))*0.05)]; % Add a bit of space either side
   
% Night Patches 
a = 1; night_start = 15; % hard coded counter
for n = 1:size(nights{set_token},2) % For each night
    r(a) = rectangle('Position',[(night_start) y_lims(1)...
        9 (y_lims(2)-y_lims(1))],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r(a),'bottom'); % Send to back
    a = a + 1; night_start = night_start + 24; % Add to counters
end
    
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
xlabel('Time (Hours)','Fontsize',32); % X Labels 
ylabel('Compressibility','Fontsize',32); % Y Labels
axis([1 size(data,2)+0.5 y_lims]); 

clear er set_token g data scrap y_lims a night_start n r 

%% WT Relative Compressibility Day vs Night 

figure; hold on; 
er = 1; g = 1; % settings
set_token = find(experiment_reps == er,1,'first'); % settings

scrap = compressibility(i_experiment_reps == er & i_group_tags == g,:);

data(:,1) = nanmean(scrap(:,dn_hour == 1),2); 
data(:,2) = nanmean(scrap(:,dn_hour == 2),2); 

plot(data',...
    'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
errorbar(nanmean(data),nanstd(data),...
    'color',cmap{set_token}(g,:),'linewidth',3);

y_lims = [(min(data(:)) - min(data(:))*0.05) ...
    (max(data(:)) + max(data(:))*0.05)]; % Add a bit of space either side

% Night Patch 
a = 1;
r(a) = rectangle('Position',[1.5 y_lims(1)...
    1 (y_lims(2)-y_lims(1))],...
    'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
uistack(r(a),'bottom'); % Send to back

box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
set(gca,'XTick',1:2); 
set(gca,'XTickLabels',{'Day','Night'}); 
ylabel({'Relative' ; 'Compressibility' ;'(Z-Score)'},'Fontsize',32); % Y Labels
axis([0.5 size(data,2)+0.5 y_lims]);

clear er g scrap dn_hour y_lims a r 

%% WT Relative Compressibility Two Way ANOVA
er = 1; % settings

% Grouping Variables
anova_group = repmat(i_group_tags(i_experiment_reps==er),...
    [size(data,2),1])'; % groups

anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
    [size(data,2),1])'; % experiments

anova_time = [];
for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
    anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
    % Allocate alternating zeros and ones to each time window
end
anova_time = anova_time';
%anova_time = anova_time(1:size(anova_time,2)/2); 

% Development Grouping Variable
if size(days_crop{set_token}(days{set_token}),2) == ...
        size(nights_crop{set_token}(nights{set_token}),2) ...
        && size(days_crop{set_token}(days{set_token}),2) > 1 % If there are an equal number of windows (>1)
    
    anova_development = []; % development
    anova_development = zeros(1,size(anova_group,2)); % Pre-allocate
    d = 1:size(anova_development,2)/(size(time_window{set_token}(1):...
        time_window{set_token}(2),2)/2):...
        size(anova_development,2); % divide into "24h" windows
    for t = 1:size(d,2)-1
        anova_development(d(t):d(t+1)-1) = t;
    end
else
    anova_development = ones(size(anova_experiment)); % use all ones
end
anova_development = ones(size(anova_experiment)); 

% Comparison
data = data(:)'; 

[twa.rc.p{er}(:,1),~,twa.rc.stats{er}] = anovan(data,...
    {anova_group,anova_time,anova_development,anova_experiment},...
    'display','off','model','full');

clear er anova_group anova_experiment anova_time t anova_development d data

%% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat');

%% Generate Hour Indexing Variable
% Note that as this is it will work far better for experiments where
% you cut out a middle section

parameter_indicies_hours{1,1} = []; % allocate active
parameter_indicies_hours{2,1} = []; % allocate inactive

for e = 1:size(experiment_reps,2) % for each experiment
    
    gap = round((lb{e}(time_window{e}(2)+1) - lb{e}(time_window{e}(1)))/...
        (size(days{e},2)*24)); % frames per hour 
    
    edges = [1 lb{e}(time_window{e}(1)):gap:...
        lb{e}(time_window{e}(2)+1) lb{e}(end,1)]; % hour times 
    
    if size(edges,2) ~= (size(days{e},2)*24) + 3 % if there are not the correct number of edges 
        edges = [1 lb{e}(time_window{e}(1)):gap:...
        lb{e}(time_window{e}(2)+1) lb{e}(time_window{e}(2)+1) lb{e}(end,1)]; 
        % add a bin with slighlty different number of frames
    end
    
    % Active
    parameter_indicies_hours{1,1} = [parameter_indicies_hours{1,1} ; ...
        single(discretize(wake_cells(experiment_tags{1,1} == e,1),...
        edges))]; % bin active data into hours 
    
    % Inactive
    parameter_indicies_hours{2,1} = [parameter_indicies_hours{2,1} ; ...
        single(discretize(sleep_cells(experiment_tags{2,1} == e,1),...
        edges))]; % bin inactive data into hours 
    
end

clear e gap edges

%% Threading Every Hour With Multiple Shuffles of Data 
% Very Slow, Roughly 16 hours for 629 fish, each with 10 shuffles 
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
    
    % Fill in Real Data
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
    threads{f,3,1}(a:2:end,1) = parameter_indicies_hours{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active time windows
    threads{f,3,1}(b:2:end,1) = parameter_indicies_hours{2,1}...
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
        for tw = min(parameter_indicies_hours{1,1}(fish_tags{1,1} == f)):...
                max(parameter_indicies_hours{1,1}(fish_tags{1,1} == f)) % for each time window
            data = []; scrap = []; % empty structures
            
            % take real clusters from this time window
            scrap{1,1} = idx_numComp_sorted{1,1}(fish_tags{1,1} == f & parameter_indicies_hours{1,1} == tw,1); % active
            scrap{2,1} = idx_numComp_sorted{2,1}(fish_tags{2,1} == f & parameter_indicies_hours{2,1} == tw,1); % inactive
            
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
    
    % Crop out extra time windows
    threads{f,3,1}(threads{f,3,1} == 1) = NaN; % remove first window 
    threads{f,3,1}(threads{f,3,1} == max(threads{f,3,1})) = NaN; % remove last window 
    threads{f,3,1} = threads{f,3,1} - 1; % shift windows back one 
        
    % Report progress
    if mod(f,50) == 0 % every 50 fish
        disp(horzcat('Threaded Fish ',num2str(f),' of ',...
            num2str(max(fish_tags{1,1})))); % Report progress
    end
    
end
toc

clear f a b tc tw data scrap

save('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180324_Hours','-v7.3');

%% Load Grammar Data From Legion 

% load grammar freq data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\Grammar_Hours_Results_Final.mat','gCount_norm')

% variables 
dn_hour = ones(1,size(gCount_norm{1,1},2))*2; % day and night hours 
dn_hour([1:14 25:38]) = 1; % day (1) and night (2) hours 
shuffles = 11; 

%% Reformat Grammar Data 
temp = gCount_norm; clear gCount_norm; 

% Allocate 
for tc = 1:size(temp,2) % for real & shuffled data
    gCount_norm{1,tc} = zeros(size(temp{1,1},1),size(temp{1,1},2),size(temp,1),...
        'single'); % {real/shuffled} sequences x time windows x fish
end

% calculate inf replacements 
scrap = zeros(1,shuffles,'single'); % 1 x real & shuffled data  
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

clear temp tc scrap inf_r f shuffles

%% Load threaded data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180227.mat'); 
clear states raw_data

%% Identifying Interesting Sequences (is)

% 180215 Notes
% 1. https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-mRMR_data.html#d119e2598
% 2. https://uk.mathworks.com/help/stats/fitcdiscr.html

% Settings
comps = 100; % number of sequences to identify

% mRMR Approach
for er = 1:max(experiment_reps) % for each experiment repeat
    
    % Grab Data
    if er == 1 % for the WT fish
        % Organised to be:
        % rows = all day mRMR_data, then all night mRMR_data
        % columns = sequences
        mRMR_data{er,1} = ...
            double(reshape(gCount_norm{1,1}(:,:,i_experiment_reps == er),...
            size(gCount_norm{1,1},1),[])');
        
        mRMR_tw{er,1} = reshape((ones(sum(i_experiment_reps == er),size(gCount_norm{1,1},2))) .* ...
            [1:size(gCount_norm{1,1},2)/2 1:size(gCount_norm{1,1},2)/2],[],1);
        
        %mRMR_data{er,1}(mRMR_data{er,1} < 0) = 0; % Remove negative values for now
        
        mRMR_data{er,1}(mRMR_tw{er,1} > 14,:) = [];
        mRMR_tw{er,1}(mRMR_tw{er,1} > 14,:) = [];
        
        
        
    else
        % Organised to be:
        % rows = fish 1 (all time windows), fish 2 (all time windows) etc
        % columns = sequences
%         mRMR_data{er,1} = ...
%             double(reshape(gCount_norm{1,1}(:,time_window{set_token}(1):time_window{set_token}(2),...
%             i_experiment_reps == er),size(uniqueSeqs{1,1},1),[]))';
%         mRMR_tw{er,1} = repmat(i_group_tags(i_experiment_reps == er),1,...
%             size(time_window{set_token}(1):time_window{set_token}(2),2))';
%         mRMR_tw{er,1} = (mRMR_tw{er,1}(:)')';
    end
    
    % Pairwise Comparisons
    tic
    counter = 1; 
    for g_one = min(mRMR_tw{er,1}):max(mRMR_tw{er,1}) % for each group
        
        tags = ones(size(mRMR_tw{er,1}))*2; % all data 
        tags(mRMR_tw{er,1} == g_one) = 1; % group of interest 
        
        % mRMR
        [comps_v{er,1}(counter,:)] = mrmr_miq_d(...
            zscore(mRMR_data{er,1}),tags,comps);

        % Classifiers
        for s = 1:comps % for each comp sequence
            % Fit a linear classifier as you add features
            % Using 10 fold cross validation
            % Hold 10% of mRMR_data back by default
            Mdl = fitcdiscr(...
                zscore(mRMR_data{er,1}(:,...
                comps_v{er,1}(counter,1:s))),...
                tags,...
                'DiscrimType','linear','CrossVal','on');
            Mdl_loss{er,1}(counter,s) = kfoldLoss(Mdl);
            Mdl_loss{er,2}(counter,s) = nanstd(kfoldLoss(Mdl,'Mode','individual')); 
        end
        
        % Minimal Feature Space
        if er == 1 % for the WT data
            mRMR_ms(er,counter) = find(islocalmin(smooth(...
                Mdl_loss{er,1}(counter,:),3)) == 1,1,'first');
        else
            mRMR_ms(er,counter) = find(Mdl_loss{er,1}(counter,:) < 0.05,1,'first');
        end
        
        disp(num2str(counter)); 
        counter = counter + 1; 
    end
    disp(horzcat('Finished mRMR Comparisons ',num2str(er),' of ',...
        num2str(max(experiment_reps)))); % report progress
    toc
    
end

clear er set_token s Mdl

%% Interesting Motifs Figure
% Left hand side = classification error + std
% Right hand side = best motifs
% Settings
er = 1; % set interest
set_token =  find(experiment_reps == er,1,'first'); % settings

figure; hold on;

% Motifs
clear scrap;
scrap = grammar_mat{1,1}(comps_v{er,1}(:,1),:); % grab motifs
scrap(:,sum(isnan(scrap)) == size(comps_v{er,1},1)) = [];
ax = imagesc(scrap,'AlphaData',isnan(scrap)==0); % imagesc with nan values in white
colormap([cmap_cluster{2,1} ; cmap_cluster{1,1}]); % merged colormap
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(ax,'CDataMapping','direct');
c = colorbar; c.Label.String = 'Cluster';
xlabel('Position in Motif','Fontsize',32);
ylabel('Motif','Fontsize',32);

clear er set_token scrap ax c

%% Minimal Feature Space Figure
er = 1; % set interest
set_token = find(experiment_reps == er,1,'first'); % settings

if er == 1
    motifs = []; 
    for m = 1:size(mRMR_ms,2)
        motifs = [motifs comps_v{er,1}(m,1:mRMR_ms(er,m))];
    end 
    motifs = unique(motifs); 
  
    mRMR_tsne{er,1} = tsne(mRMR_data{er,1}(:,motifs),...
        'Algorithm','barneshut','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
        'Perplexity',30,'Standardize',1,'Verbose',0);
else
%     scrap = [];
%     for c = 1:size(comps_v{er,1},1) % for each comparison
%         scrap = [scrap comps_v{er,1}(c,1:mRMR_ms(er,c))]; 
%     end
%     scrap = unique(scrap); 
%     mRMR_tsne{er,1} = tsne(mRMR_data{er,1}(:,scrap),...
%         'Algorithm','exact','Exaggeration',4,'NumDimensions',2,'NumPCAComponents',0,...
%         'Perplexity',30,'Standardize',1,'Verbose',0);
end

figure; hold on;
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

