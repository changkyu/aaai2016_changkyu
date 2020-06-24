%--------------------------------------------------------------------------
% SfM Synthetic data result plot (Cube)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2012/03/15)
%--------------------------------------------------------------------------
fontsize = 14;

load('result/sfm_cube/result_cube_all.mat');

% Get average over runs w/ different initializations
SAx = mean(SAs, 4);

% Separate out each angles
SAx1 = SAx(:,1,:);
SAx2 = SAx(:,2,:);

% Prepare matrix to plot - each one has (noises) x (runs)
SAx1t = zeros(20, 11);
SAx2t = zeros(20, 11);
for idn = 1:size(SAx,1)
    for idr = 1:size(SAx,3)
        SAx1t(idr, idn) = SAx1(idn,1,idr) * 180 / pi;
        SAx2t(idr, idn) = SAx2(idn,1,idr) * 180 / pi;
    end
end
    
% Prepare vectors to plot - mean of each box (noises) x 1
SAx1c = mean(SAx1t, 1);
SAx2c = mean(SAx2t, 1);

% Prepare figure
figure;
h = axes;
hold on;
hl = xlabel('Noise level (%)'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Subspace angle (degree)'); 
set(hl, 'FontSize',fontsize);
set(h,'FontSize',fontsize);

% Plot Boxes
boxplot(SAx1t,'labels',{'0','1','2','3','4','5','6','7','8','9','10'},...
    'colors','k','symbol','g+','Position',1:1:11,'widths',0.2); 
set(gca,'XTickLabel',{' '});
boxplot(SAx2t,'labels',{'0','1','2','3','4','5','6','7','8','9','10'},...
    'colors','b','symbol','r+','Position',1.3:1:11.3,'widths',0.2); 

drawnow;
