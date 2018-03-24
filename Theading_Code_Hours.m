%% Hour by Hour Threading Code  

% An adaptation of Threading_Code to 
    % Seperatly compress each hour 
    % Count motifs every hour 
    
%% Required scripts 

%% Notes 

%% Settings 
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 
set(0,'defaultfigurecolor',[1 1 1]); % white background

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

%% Load Threaded Data & Seperately Compressed Data From Legion 

% load threaded data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\Compression_Hours_Results_Final.mat'); % load compressed data 

% 

%% Calculate Relative Compressibility Each Day/Night
% The compressibility of a sequence of uncompressed length l is given by the sum of the savings S 
% at each iteration divided by l (Gomez-Marin et al.,2016) 

% Uncompressed lengths 
tw = 48; % hard coded maximum number of time windows 
seq_lengths = nan(size(threads,1),tw,'single'); % fish x time windows  

for f = 1:size(threads,1) % for each fish 
    seq_lengths(f,1:max(threads{f,3,1})) = histcounts(threads{f,3,1}); 
    % count sequence length per time window 
end 

% Compressibility 
compressibility = zeros(size(totSavings),'single'); % fish x time windows x t/c
compressibility_rel = zeros(size(totSavings,1),tw,'single'); % fish x time windows  

for f = 1:size(threads,1) % for each fish 
    compressibility(f,:,:) = totSavings(f,:,:)./seq_lengths(f,:); % calculate compressibility 
    compressibility_rel(f,:) = (compressibility(f,:,1) - nanmean(compressibility(f,:,2:end),3))./...
        nanstd(compressibility(f,:,2:end),0,3); % relative compressibility (Z-Score) 
end 

clear f 

%% WT Relative Compressibility Figure

set_token = find(experiment_reps == er,1,'first'); % settings
figure;
hold on; set(gca,'FontName','Calibri'); clear scrap;

for g = 1:max(i_group_tags(i_experiment_reps == er)) % for each group
    clear data;
    data = compressibility_rel(i_experiment_reps == er,time_window{set_token}(1):time_window{set_token}(2));
    plot(data',...
        'color',cmap{set_token}(g,:)+(1-cmap{set_token}(g,:))*(1-(1/(5)^.5)),'linewidth',1.5);
    errorbar(nanmean(data),nanstd(data),...
        'color',cmap{set_token}(g,:),'linewidth',3);
    
    scrap(1,g) = min(data(:));
    scrap(2,g) = max(data(:));
end

y_lims = [(min(scrap(1,:)) + min(scrap(1,:))*0.05) ...
    (max(scrap(2,:)) + max(scrap(2,:))*0.05)]; % Add a bit of space either side
   
% Night Patches 
a = 1; night_start = first_night{set_token}; % Start counters
for n = 1:size(nights{set_token},2) % For each night
    r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
        1 (y_lims(2)-y_lims(1))],...
        'FaceColor',night_color{set_token},'Edgecolor',[1 1 1]);
    uistack(r(a),'bottom'); % Send to back
    a = a + 1; night_start = night_start + 2; % Add to counters
end
    
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
set(gca,'XTick',[]); 
xlabel('Time (Days/Nights)','Fontsize',32); % X Labels 
ylabel({'Relative' ; 'Compressibility' ;'(Z-Score)'},'Fontsize',32); % Y Labels
axis([0.5 size(data,2)+0.5 y_lims]); 

clear er set_token g data scrap y_lims a n r 

%% WT Relative Compressibility Two Way ANOVA
er = 1;
set_token = find(experiment_reps == er,1,'first'); % settings

% Grouping Variables
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

% Comparison
scrap = compressibility_rel(i_experiment_reps == er,time_window{set_token}(1):time_window{set_token}(2));
scrap = scrap(:)'; % Vectorise

[twa.rc.p{er}(:,1),~,twa.rc.stats{er}] = anovan(scrap,...
    {anova_group,anova_time,anova_development,anova_experiment},...
    'display','off','model','full');

clear er set_token anova_group anova_experiment anova_time anova_development scrap

%% Load & Reformat Hour Grammar Data From Legion 

%% Identifying Interesting Sequences (is)
    % To seperate each hour 
    
%% Interesting Motifs Figure 

%% Minimal Feature Space Figure
