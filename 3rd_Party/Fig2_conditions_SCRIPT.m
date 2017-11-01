% script to use compressive algorithm to identify potentially informative
% sequences from input worm videos.  These are then combined into a total
% dictionary which is used to compare wild isolates.
%
% Reproduces part of figure 2

% should repeats be included?
repeats = false;
nMax = 15; % the maximum length sequence to consider when compressing
postureNum = 30;

% set the root directory
directory = '/Users/abrown/Andre/wormVideos/robynTracking/';

% load the representative postures
load(['/Users/abrown/Andre/wormVideos/results-12-05-10/Laura Grundy/'...
    'gene_NA/allele_NA/N2/on_food/XX/30m_wait/'...
    'postures_' num2str(postureNum) ...
    '-centers_20-files_5000-framesPerFile.mat'])
postures = postures';

% get the file names
[fileList, ~] = ...
    dirSearch(directory, ...
    ['stateSequence_' num2str(postureNum) 'means_N2-1000_ds-5.mat']);
% remove optogenetics experiments
dropInds = zeros(numel(fileList), 1);
for ii = 1:numel(fileList)
    if ~isempty(strfind(fileList{ii}, '/optoG/'))
        dropInds(ii) = 1;
    end
end
fileList(logical(dropInds)) = [];

% hack the worm names
wormNames = cell(numel(fileList), 1);
for ii = 1:numel(fileList)
    % get the positions of forward slashes in file name
    slashPositions = strfind(fileList{ii}, '/');
    
    % get the strain name
    wormNames{ii} = fileList{ii}(slashPositions(6)+1:slashPositions(7)-1);
end

% get the unique worm names
[uniqueNames, ~, classInts] = unique(wormNames);

% initialisation and model set up
seqCell = cell(numel(uniqueNames), 1); % cell for storing state sequences
gTermCell = cell(numel(uniqueNames), 1); % cell for storing n-grams
freqCell = cell(numel(uniqueNames), 1); % frequency distributions for each strain
uniqueSeqs = cell(0);

% get the grammars for each worm
for jj = 1:numel(uniqueNames)
    
    % get the file names of the current strain
    currentInds = find(strcmp(wormNames, uniqueNames{jj}));
    % loop through files to populate seqCell and get unique n-grams
    for ii = 1:numel(currentInds)
        disp([jj/numel(uniqueNames), ii/numel(currentInds)])
        
        % import the data
        stateData = cell2mat(struct2cell(load(fileList{currentInds(ii)})));
        % ds is 1 for expansion because although data are downsampled, the
        % data used to train the language model were also downsampled
        if repeats
            expandedSeq = expandSequence(stateData, 1);
        else
            expandedSeq = stateData;
        end
        seqCell{jj, ii} = expandedSeq(:, 1)';
        
        % get the grammar for the current iteration
        [grammar, ~] = compressSequenceNFast(expandedSeq(:, 1)', [], nMax);
        
        % get the terminal symbols
        if ~isempty(grammar)
            gTerm = getGrammarTerminals(grammar);
            gTermCell{jj, ii} = gTerm;
            
            % append grammar terminals to the total set
            uniqueSeqs = vertcat(uniqueSeqs, gTerm);
            
            % re-calculate unique sequences from larger set
            [uniqueSeqs, ~] = countUniqueRows(uniqueSeqs);
        end
    end
end


% make a matrix of counts for each of the unique sequences in each
% individual worm's postural sequence
countMat = NaN(numel(wormNames), size(uniqueSeqs, 1));
freqMat = NaN(numel(wormNames), size(uniqueSeqs, 1));
nn = 1;
% loop over strains
for ii = 1:numel(uniqueNames)
    disp(ii/numel(uniqueNames))
    % loop over non-empty state sequences for current strain
    for jj = 1:sum(~cellfun(@isempty, seqCell(ii, :)))
        % loop over unique n-grams
        for kk = 1:size(uniqueSeqs, 1)
            % count the number of occurences of the current unique
            % sequences in the current state sequence
            counts = ...
                numel( strfind(seqCell{ii, jj}, uniqueSeqs{kk, :}) );
            countMat(nn, kk) = counts;
            freqMat(nn, kk) = counts / numel(seqCell{ii, jj});
        end
        nn = nn + 1;
    end
end

% convert unique sequences to cell array of strings for use with ismember
% below
uniqueSeqsString = cellfun(@num2str, uniqueSeqs, 'UniformOutput', false);

% do pair-wise comparisons of strains
pValCell = cell(numel(uniqueNames), numel(uniqueNames));

% numHits = 100; % how many sequences should be kept from each comparison
for ii = 1:numel(uniqueNames)-1
    % get the n-grams of the ii'th strain
    gTerms1 = vertcat(gTermCell{ii, :});
    
    % get name matches of ii'th strain
    nameMatches1 = find(strcmp(wormNames, uniqueNames{ii}));
    
    for jj = ii+1:numel(uniqueNames)
        disp([ii/numel(uniqueNames), jj/numel(uniqueNames)])
        % get the grammar terminals of the jj'th strain
        gTerms2 = vertcat(gTermCell{jj, :});
        
        % get the unique grammar terminals of the current pair
        [uniqueSeqsPair, ~] = countUniqueRows([gTerms1; gTerms2]);
        
        % find the indices of the current unique sequences in the total set
        % of unique sequences.  To use ismember, must first convert to cell
        % of stings
        uniqueSeqsPairString = ...
            cellfun(@num2str, uniqueSeqsPair, 'UniformOutput', false);
        [~, matchInds] = ...
            ismember(uniqueSeqsPairString, uniqueSeqsString);
        
        % get the submatrix from count mat for the current pair
        nameMatches2 = find(strcmp(wormNames, uniqueNames{jj}));
        rowInds = [nameMatches1; nameMatches2];
        countMatPair = countMat(rowInds, matchInds);
        freqMatPair = freqMat(rowInds, matchInds);
        
        %         % drop columns corresponding to very rare behaviours
        %         dropInds = sum(countMatPair) < 5;
        %         countMatPair(:, dropInds) = [];
        %         freqMatPair(:, dropInds) = [];
        %         matchInds(dropInds) = [];
        
        % calculate the F-stats for each behaviour
        classIntsPair = [ones(numel(nameMatches1), 1); ...
            ones(numel(nameMatches2), 1) + 1];
        
        pVals = NaN(size(freqMatPair, 2), 1);
        for kk = 1:size(freqMatPair, 2)
            % also do rank sum tests
            p = ranksum(freqMatPair(classIntsPair == 1, kk), ...
                freqMatPair(classIntsPair == 2, kk));
            
            % get the p-value
            pVals(kk) = p;
        end
        
        % add to cell
        pValCell{ii, jj} = pVals;
    end
end


% get p-values and correct for multiple comparisons
pValsTotal = [];
for ii = 1:size(pValCell, 1)-1
    for jj = ii+1:size(pValCell, 2)
        pValsTotal = [pValsTotal; pValCell{ii, jj}];
    end
end
[~, pCrit, bhFDR] = fdr_bh(pValsTotal, 0.05 , 'dep');

% for the significantly modulated behaviours, record the comparison
% statistics
comparisonStats = cell(numel(uniqueNames), numel(uniqueNames));
informativeBehaviours = cell(numel(uniqueNames), numel(uniqueNames));
for ii = 1:numel(uniqueNames)-1
    % get the n-grams of the ii'th strain
    gTerms1 = vertcat(gTermCell{ii, :});
    
    % get name matches of ii'th strain
    nameMatches1 = find(strcmp(wormNames, uniqueNames{ii}));
    
    for jj = ii+1:numel(uniqueNames)
        disp([ii/numel(uniqueNames), jj/numel(uniqueNames)])
        % get the grammar terminals of the jj'th strain
        gTerms2 = vertcat(gTermCell{jj, :});
        
        % get the unique grammar terminals of the current pair
        [uniqueSeqsPair, ~] = countUniqueRows([gTerms1; gTerms2]);
        
        % find the indices of the current unique sequences in the total set
        % of unique sequences.  To use ismember, must first convert to cell
        % of stings
        uniqueSeqsPairString = ...
            cellfun(@num2str, uniqueSeqsPair, 'UniformOutput', false);
        [~, matchInds] = ...
            ismember(uniqueSeqsPairString, uniqueSeqsString);
        
        % get the submatrix from count mat for the current pair
        nameMatches2 = find(strcmp(wormNames, uniqueNames{jj}));
        rowInds = [nameMatches1; nameMatches2];
        countMatPair = countMat(rowInds, matchInds);
        freqMatPair = freqMat(rowInds, matchInds);
        
        %         % drop columns corresponding to very rare behaviours
        %         dropInds = sum(countMatPair) < 5;
        %         countMatPair(:, dropInds) = [];
        %         freqMatPair(:, dropInds) = [];
        %         matchInds(dropInds) = [];
        
        % calculate the F-stats for each behaviour
        classIntsPair = [ones(numel(nameMatches1), 1); ...
            ones(numel(nameMatches2), 1) + 1];
        
        
        % control false discovery rate
        hitInds = find(pValCell{ii, jj} < pCrit);
        
        % sort hitInds by bhFDR
        [~, hitSortInds] = sort(pValCell{ii, jj}(hitInds));
        hitInds = hitInds(hitSortInds);
        
        % get the frequencies for each strain
        freq1 = mean(freqMatPair(classIntsPair == 1, :));
        [~, freqSortInds1] = sort(freq1, 'descend');
        freq2 = mean(freqMatPair(classIntsPair == 2, :));
        [~, freqSortInds2] = sort(freq2, 'descend');
        
        % for the hits, record the following:
        % frequency in strain 1, frequency in strain 2, rank in
        % strain 1, and rank in strain 2, p-values for comparisons
        comparisonMat = NaN(length(hitInds), 5);
        
        % add data
        comparisonMat(:, 1) = freq1(hitInds)'; % frequencies in strain 1
        comparisonMat(:, 2) = freq2(hitInds)'; % frequencies in strain 2
        [~, freqRank1] = ismember(hitInds, freqSortInds1);
        comparisonMat(:, 3) = freqRank1; % rank in strain 1
        [~, freqRank2] = ismember(hitInds, freqSortInds2);
        comparisonMat(:, 4) = freqRank2; % rank in strain 2
        comparisonMat(:, 5) = pValCell{ii, jj}(hitInds); % p-values for each comparison
        
        
        % add to cell
        comparisonStats{ii, jj} = comparisonMat;
        
        % also record informative n-grams
        informativeBehaviours{ii, jj} = ...
            uniqueSeqs(matchInds(hitInds), :);
    end
end

% get all of the hits and the corresponding comparison stats
allHits = vertcat(informativeBehaviours{:});
comparisonStatMat = vertcat(comparisonStats{:});


% count the hits
[uniqueHits, counts] = countUniqueRows(allHits);



% also make a heat map that shows how the frequencies of the hits vary
% across the strains

% get all the frequencies for the hits
uniqueHitsString = ...
    cellfun(@num2str, uniqueHits, 'UniformOutput', false);
hitInds = find(ismember(uniqueSeqsString, uniqueHitsString));
freqMean = NaN(numel(uniqueNames), numel(hitInds));

for ii = 1:numel(uniqueNames)
    for jj = 1:numel(hitInds)
        freqMean(ii, jj) = mean(freqMat(classInts == ii, hitInds(jj)));
    end
end

% normalise the frequencies
freqMeanMean = mean(freqMean);
freqMeanStd = std(freqMean);
freqMeanNorm = ...
    (freqMean - repmat(freqMeanMean, size(freqMean, 1), 1)) ./ ...
    repmat(freqMeanStd, size(freqMean, 1), 1);

% plot the state densities.  Cluster to aid visualisation.
clf = clustergram(freqMeanNorm, 'Linkage', 'complete');

% it's easier to control the colour using image than clustergram, so take
% row labels and re-plot
rowOrder = str2num(char(clf.RowLabels));
colOrder = str2num(char(clf.ColumnLabels));

figure
imagesc(freqMeanNorm(rowOrder(end:-1:1), colOrder), [-2, 2])
cmap = cbrewer('div', 'RdBu', 50); % RdYlBu
colormap(cmap(end:-1:1, :))


% also make box plots for each of the strains for the different hits
for ii = 1:size(uniqueHits, 1)
    % get the index of the current hit in the total set of unique n-grams
    matchInd = ...
        ismember(uniqueSeqsString, uniqueHitsString(colOrder(ii)));
    
    % use anova1 function to plot box plot of frequencies in all strains
    [~, table, stats] = anova1(freqMat(:, matchInd), wormNames);
    
    
    % recalculate p-values for this n-gram to add bars to the plot
    sigDiffMap = zeros(numel(uniqueNames));
    %     offset = 5;
    for jj = 1:numel(uniqueNames)-1
        % get name matches of jj'th strain
        nameMatches1 = find(strcmp(wormNames, uniqueNames{jj}));
        freq1 = freqMat(nameMatches1, matchInd);
        
        for kk = jj+1:numel(uniqueNames)
            % get the submatrix from count mat for the current pair
            nameMatches2 = find(strcmp(wormNames, uniqueNames{kk}));
            freq2 = freqMat(nameMatches2, matchInd);
            
            % calculate the p-value
            p = ranksum(freq1, freq2);
            
            % for hits, add a bar to the box plot to indicate significant
            % difference
            if p <= pCrit
                % add the significantly different points to the sigDiffMap
                sigDiffMap(jj, kk) = 1;
                sigDiffMap(kk, jj) = 1;
                %                 ylim([-5e-4, (offset+0.5) * stats.s])
                %                 line([jj, kk], [offset * stats.s, offset * stats.s], ...
                %                     'LineWidth', 2, 'Color', [0.1 0.3 0.9])
                %                 offset = offset + 0.5;
            end
        end
    end
    
    
    
    title(['Number of total hits: ' num2str(counts(colOrder(ii))) '.'], ...
        'FontSize', 14)
    
    % save the box plot
    saveas(gcf, ['./boxplots-conditions-nMax_15_k-' num2str(postureNum) ...
        '/boxplot_' num2str(ii) '.pdf'])
    close all
    
    % plot and save the map of significant differences
    imagesc(sigDiffMap)
    colormap([1,1,1; 0,0,0])
    set(gca,'XTick',0.5:numel(uniqueNames)+0.5)
    set(gca,'YTick',0.5:numel(uniqueNames)+0.5)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'GridLineStyle', '-')
    grid on
    
    saveas(gcf, ['./boxplots-conditions-nMax_15_k-' num2str(postureNum) ...
        '/sigDiffMap_' num2str(ii) '.pdf'])
    close all
    
    
    % also plot the behaviour
    figure
    plotSequence(postures, uniqueHits{colOrder(ii)}, 'r', 'b')
    text(0, 0.6, num2str(uniqueHits{colOrder(ii)}), 'FontSize', 16)
    
    % save the behaviour plot
    saveas(gcf, ['./boxplots-conditions-nMax_15_k-' num2str(postureNum) ...
        '/posturePlot_' num2str(ii) '.pdf'])
    close all
end



% for choosing which sample sequences to plot, determine which sequence is
% the most compressive in each of the conditions.
bestCompressibility = [0, 0, 0];
bestSeq = cell(1, 3);
for jj = 1:numel(uniqueNames)
    
    % get the file names of the current strain
    currentInds = find(strcmp(wormNames, uniqueNames{jj}));
    % loop through files to populate seqCell and get unique n-grams
    for ii = 1:numel(currentInds)
        disp([jj/numel(uniqueNames), ii/numel(currentInds)])
        
        % get the most compressive sequence and the savings
        [outSeq, ~, savings] = compressiveNFast(seqCell{jj, ii}, nMax);
        
        % get compression ratio
        compressiblity = savings/length(seqCell{jj, ii});
        
        % update best savings
        if compressiblity > bestCompressibility(jj) && ...
                any(strcmp(uniqueHitsString, num2str(outSeq)))
            bestCompressibility(jj) = compressiblity;
            bestSeq{jj} = outSeq;
        end
    end
end

% convert the most compressive hits to their indices in the boxplots folder
compressMatchInds = NaN(1, numel(uniqueNames));
for ii = 1:numel(uniqueNames)
    % find the best index in uniqueHits
    currentInd = find(strcmp(uniqueHitsString, num2str(bestSeq{ii})));
    
    % convert to plot index based on column order
    compressMatchInds(ii) = find(colOrder == currentInd);
end

