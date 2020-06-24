%--------------------------------------------------------------------------
% SfM Synthetic data result plot (Cube)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2015/04/03)
%--------------------------------------------------------------------------
fontsize = 14;

load(sprintf('results/sfm_%s/%s_all.mat', model, model));

[noises, ~, runs] = size(SAs);

% Separate out each angles
SAx0 = SAs(:,1,:);
SAx1 = SAs(:,2,:);
SAx2 = SAs(:,3,:);
SAx3 = SAs(:,4,:);
SAx4 = SAs(:,5,:);
SAx5 = SAs(:,6,:);
SAx6 = SAs(:,7,:);

% Prepare matrix to plot - each one has (noises) x (runs)
SAx0t = zeros(runs, noises);
SAx1t = zeros(runs, noises);
SAx2t = zeros(runs, noises);
SAx3t = zeros(runs, noises);
SAx4t = zeros(runs, noises);
SAx5t = zeros(runs, noises);
SAx6t = zeros(runs, noises);
for idn = 1 : noises
    for idr = 1 : runs
        SAx0t(idr, idn) = SAx0(idn,1,idr) * 180 / pi;
        SAx1t(idr, idn) = SAx1(idn,1,idr) * 180 / pi;
        SAx2t(idr, idn) = SAx2(idn,1,idr) * 180 / pi;
        SAx3t(idr, idn) = SAx3(idn,1,idr) * 180 / pi;
        SAx4t(idr, idn) = SAx4(idn,1,idr) * 180 / pi;
        SAx5t(idr, idn) = SAx5(idn,1,idr) * 180 / pi;
        SAx6t(idr, idn) = SAx6(idn,1,idr) * 180 / pi;
    end
end
    
% Prepare vectors to plot - mean of each box (noises) x 1
SAx0c = mean(SAx0t, 1);
SAx1c = mean(SAx1t, 1);
SAx2c = mean(SAx2t, 1);
SAx3c = mean(SAx3t, 1);
SAx4c = mean(SAx4t, 1);
SAx5c = mean(SAx5t, 1);
SAx6c = mean(SAx6t, 1);

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
% boxplot(SAx0t, ...
%     'colors','y','symbol','y+','Position',1:1:11,'widths',0.08); 
% set(gca,'XTickLabel',{' '})
% boxplot(SAx1t,...
%     'colors','r','symbol','r+','Position',1.1:1:11.1,'widths',0.08); 
% set(gca,'XTickLabel',{' '})
% boxplot(SAx2t,...
%     'colors','g','symbol','g+','Position',1.2:1:11.2,'widths',0.08); 
% set(gca,'XTickLabel',{' '})
% boxplot(SAx3t,...
%     'colors','m','symbol','m+','Position',1.3:1:11.3,'widths',0.08); 
% set(gca,'XTickLabel',{' '})
% boxplot(SAx4t,...
%     'colors','b','symbol','b+','Position',1.4:1:11.4,'widths',0.08); 
% set(gca,'XTickLabel',{' '})
% boxplot(SAx5t,...
%     'colors','k','symbol','k+','Position',1.5:1:11.5,'widths',0.08); 
% set(gca,'XTickLabel',{' '})
% boxplot(SAx6t,'label',{'0','1','2','3','4','5','6','7','8','9','10'},...
%     'colors','c','symbol','c+','Position',1.6:1:11.6,'widths',0.08); 
%ylim([0 45]);

boxplot(SAx1t,...
    'colors','k','symbol','m+','Position',1.3:1:11.3,'widths',0.1); 
set(gca,'XTickLabel',{' '})
boxplot(SAx3t,...
    'colors','b','symbol','g+','Position',1.5:1:11.5,'widths',0.1); 
set(gca,'XTickLabel',{' '})
boxplot(SAx4t,'label',{'0','1','2','3','4','5','6','7','8','9','10'},...
    'colors','r','symbol','c+','Position',1.7:1:11.7,'widths',0.1); 

title(['PPCA (black), D-PPCA (blue), D-PPCA-ANT (red)']);
% title('SVD (y), PPCA (r), D-PPCA (g), BPCA (m), D-BPCA (b), PPCA (IR, k), VBPCA (IR, c)');

drawnow;
