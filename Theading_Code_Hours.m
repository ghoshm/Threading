%% Hour by Hour Threading Code  

% An adaptation of Threading_Code to 
    % Seperatly compress each day/night 
    % Count motifs every hour 
    
%% Required scripts 

%% Notes 

%% Settings 
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 
set(0,'defaultfigurecolor',[1 1 1]); % white background

%% Load data 
load('D:\Behaviour\SleepWake\Re_Runs\Threading\New\180223.mat')

%% Load Seperately Compressed Data From Legion 

%% Calculate Relative Compressibility Each Day/Night
% The compressibility of a sequence of uncompressed length l is given by the sum of the savings S 
% at each iteration divided by l (Gomez-Marin et al.,2016) 

%% Relative  Compressibility Figure 

%% Generate Hour Indexing Variable

parameter_indicies_hours{1,1} = []; % allocate active
parameter_indicies_hours{2,1} = []; % allocate inactive

% New
for e = 1:size(experiment_reps,2) % for each experiment
    
    % Active
    parameter_indicies_hours{1,1} = [parameter_indicies_hours{1,1} ; ...
        single(discretize(wake_cells(experiment_tags{1,1} == e,1),...
        lb{e}(time_window{e}(1)):(fps{e}*60*60):...
        (lb{e}(time_window{e}(2)+1) + (fps{e}*60*60))))];
    
    % Inactive
    parameter_indicies_hours{2,1} = [parameter_indicies_hours{2,1} ; ...
        single(discretize(sleep_cells(experiment_tags{2,1} == e,1),...
        lb{e}(time_window{e}(1)):(fps{e}*60*60):...
        (lb{e}(time_window{e}(2)+1) + (fps{e}*60*60))))];
    
end

clear e

%% Threads_Hours 
    % Replace the day/night tags in threads with hour tags 
threads_hours = threads; 

for f = 1:max(fish_tags{1,1}) % For each fish
    
    % Deterime starting state (a = active, b = inactive)
    if wake_cells(find(fish_tags{1,1} == f,1,'first'),1) == 1 % If the fish starts active
        a = 1; b = 2; % start active
    else % If the fish starts inactive
        a = 2; b = 1; % start inactive
    end
    
    % Time Windows
    threads_hours{f,3,1}(a:2:end,1) = parameter_indicies_hours{1,1}...
        (fish_tags{1,1} == f,1); % Fill in active time windows
    threads_hours{f,3,1}(b:2:end,1) = parameter_indicies_hours{2,1}...
        (fish_tags{2,1} == f,1); % Fill in inactive time windows
    
end 

clear f a b 

%% Load & Reformat Grammar Data From Legion 

%% Identifying Interesting Sequences (is)
    % To seperate each hour 
    
%% Interesting Motifs Figure 

%% Minimal Feature Space Figure
