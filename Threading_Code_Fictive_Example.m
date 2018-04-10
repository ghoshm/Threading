%% Fictive Threading Example 

% Data 
% Clusters 
c{1,1} = [0 0 0];
c{2,1} = [1 1];
c{3,1} = [3 2.5 1]; 

vec_dp = [c{1,1} c{2,1} c{1,1} c{2,1} c{1,1} c{2,1} c{1,1} c{2,1} c{1,1} c{3,1}...
    c{1,1} c{2,1} c{1,1} c{2,1} c{1,1} c{2,1} c{1,1} c{2,1} c{1,1} c{3,1} c{1,1}];
vec_c = [1 2 1 2 1 2 1 2 1 3 1 2 1 2 1 2 1 2 1 3 1]; 
vec_c2 = [1 4 4 3 1 4 4 3 1];
vec_c3 = [5 5 1]; 
vec_s = [1 1 1 2 2 1 1 1 2 2 1 1 1 2 2 1 1 1 2 2 1 1 1 3 3 3 ...
    1 1 1 2 2 1 1 1 2 2 1 1 1 2 2 1 1 1 2 2 1 1 1 3 3 3 1 1 1];
m1 = [4 2 1 2 1]; 
m2 = [5 1 4 4 3]; 
motif_1 = {'4', '=', '2', '1', '2', '1'};
motif_2 = {'5', '=', '1','4','4','3'};
[grammar, compVec, totSavings] = compressSequenceNFast(vec_c,(max(vec_c)+1),(length(vec_c)-1)); 

cmap_cluster = flip(lbmap(max(compVec),'RedBlue')); % generate a colormap for the clusters/motifs

%% Fictive Example Figure 
    % Circles & Numbers 
    
% Figure Formatting
figure; hold on; 
axis([0 55 -.5 6]);
ax1 = gca; 
yruler = ax1.YRuler; yruler.Axle.Visible = 'off';
xruler = ax1.XRuler; xruler.Axle.Visible = 'off';
set(gca,'TickLength',[0 0]); set(gca,'XTick',[]);
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(gca,'YTick',0:5);
yticklabels({'Sequence:','Motif 2:','Sequence:','Motif 1:','Sequence:','Delta Px:',''})

% Plotting 
plot((vec_dp/3)+5,'k','linewidth',3); % delta px 
text((1:length(vec_c))-.5,ones(size(vec_c))*4,string(vec_c),'FontName','Calibri','Fontsize',32); 
text((1:length(motif_1))-.5,ones(size(motif_1))*3,motif_1,'FontName','Calibri','Fontsize',32); 
text((1:length(vec_c2))-.5,ones(size(vec_c2))*2,string(vec_c2),'FontName','Calibri','Fontsize',32); 
text((1:length(motif_2))-.5,ones(size(motif_2)),motif_2,'FontName','Calibri','Fontsize',32); 
text((1:length(vec_c3))-.5,ones(size(vec_c3))*0,string(vec_c3),'FontName','Calibri','Fontsize',32); 

% Most Compressive Rectangles
r(1) = rectangle('Position',[1.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(2) = rectangle('Position',[5.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(3) = rectangle('Position',[11.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(4) = rectangle('Position',[15.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(5) = rectangle('Position',[.5 1.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(5,:),'EdgeColor',[0 0 0]);
r(6) = rectangle('Position',[4.5 1.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(5,:),'EdgeColor',[0 0 0]);
uistack(r(:),'bottom');

% Clusters 
gscatter(1:size(vec_s,2),ones(size(vec_s))*4.7,vec_s,cmap_cluster(1:3,:),'.',60); % dp 
gscatter(1:size(vec_c,2),ones(size(vec_c))*3.7,vec_c,cmap_cluster(1:3,:),'.',60); % vec_c
gscatter([1 3 4 5 6],ones(size(m1))*2.7,m1,[cmap_cluster(1:2,:) ; cmap_cluster(4,:)],'.',60); % m1
gscatter(1:size(vec_c2,2),ones(size(vec_c2))*1.7,vec_c2,[cmap_cluster(1,:) ; cmap_cluster(3:4,:)],'.',60); % vec_c2
gscatter([1 3 4 5 6],ones(size(m2))*.7,m2,[cmap_cluster(1,:) ; cmap_cluster(3:5,:)],'.',60); % m1
gscatter(1:size(vec_c3,2),zeros(size(vec_c3))-.3,vec_c3,[cmap_cluster(1,:) ; cmap_cluster(5,:)],'.',60); % vec_c2

% Legend 
legend off

%% Fictive Example Figure 
    % Circles Only  
    
% Figure Formatting
figure; hold on; 
axis([0 55 -.5 6]);
ax1 = gca; 
yruler = ax1.YRuler; yruler.Axle.Visible = 'off';
xruler = ax1.XRuler; xruler.Axle.Visible = 'off';
set(gca,'TickLength',[0 0]); set(gca,'XTick',[]);
set(gca,'FontName','Calibri'); box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
set(gca,'YTick',0:5);
yticklabels({'Sequence:','Motif 2:','Sequence:','Motif 1:','Sequence:','Delta Px:',''})

% Plotting 
plot((vec_dp/3)+5,'k','linewidth',3); % delta px 

% Most Compressive Rectangles
r(1) = rectangle('Position',[1.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(2) = rectangle('Position',[5.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(3) = rectangle('Position',[11.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(4) = rectangle('Position',[15.5 3.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(4,:),'EdgeColor',[0 0 0]);
r(5) = rectangle('Position',[.5 1.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(5,:),'EdgeColor',[0 0 0]);
r(6) = rectangle('Position',[4.5 1.85 4 0.35],'Curvature',0.6,...
    'linewidth',1.5,'FaceColor',cmap_cluster(5,:),'EdgeColor',[0 0 0]);
uistack(r(:),'bottom');

% Clusters 
gscatter(1:size(vec_s,2),ones(size(vec_s))*4.7,vec_s,cmap_cluster(1:3,:),'.',60); % dp 
gscatter(1:size(vec_c,2),ones(size(vec_c))*4,vec_c,cmap_cluster(1:3,:),'.',60); % vec_c
gscatter([1 3 4 5 6],ones(size(m1))*3,m1,[cmap_cluster(1:2,:) ; cmap_cluster(4,:)],'.',60); % m1
gscatter(1:size(vec_c2,2),ones(size(vec_c2))*2,vec_c2,[cmap_cluster(1,:) ; cmap_cluster(3:4,:)],'.',60); % vec_c2
gscatter([1 3 4 5 6],ones(size(m2)),m2,[cmap_cluster(1,:) ; cmap_cluster(3:5,:)],'.',60); % m1
gscatter(1:size(vec_c3,2),zeros(size(vec_c3)),vec_c3,[cmap_cluster(1,:) ; cmap_cluster(5,:)],'.',60); % vec_c2

% Arrows 
text(1.5,3,'\rightarrow','FontName','Calibri','Fontsize',20,'FontWeight','bold'); 
text(1.5,1,'\rightarrow','FontName','Calibri','Fontsize',20,'FontWeight','bold'); 

% Legend 
legend off