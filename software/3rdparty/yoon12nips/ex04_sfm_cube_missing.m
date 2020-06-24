%--------------------------------------------------------------------------
% SfM Synthetic data demonstration (MAR and MNAR)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2015/02/08)
%--------------------------------------------------------------------------

clear;

model = 'cube';

missings = 10; % Number of independent runs of missing values (MAR)
repeats = 20; % Number of independent runs of different initial values
miss_rate = 0.2;

% Results are stored here with each row with different missing value set
% SAc will contain subspace angles for centralized PPCA
SAc = zeros(missings+1, 1);
% SAs will contain subspace angles for some independent runs of D-PPCA
SAs = zeros(missings+1, repeats);

%--------------------------------------------------------------
% Data
% Frequency for intermediate object function value output
objfreq_c = 0;
objfreq_d = 0;

% Our data is fixed this time. We have to remove the first point from the
% cube data because it is occluded over all frames, thus never observed.
% This will not harm the validity of our result since this will only
% increase difficulties for MNAR settings. We use all points for MAR since
% all points have same probability of missing.
load(sprintf('data/%s/data_missing.mat', model));
mm_mat = measurement_matrix;
mm_mat_MAR = mm_mat;
GT_MAR = GT;
mm_mat_MNAR = mm_mat(:,2:end);
GT_MNAR = GT(2:end,:);

% get dimension of measurement matrix
[D, N_MAR] = size(mm_mat_MAR);
N_MNAR = size(mm_mat_MNAR, 2);

% Measurement matrix should be in the specific form
if mod(D, 2) ~= 0
    error('Measurement matrix should be (2 x #frames)x(# points) form!');
end

% Prepare missing value matrices
% First MissIDX = MNAR, other missing value sets (indices) of it are MAR
MissIDXs = cell(missings+1,1);
MissIDXs{1} = MissIDX(:,2:end);
for idx = 2:(missings+1)
    MissIDXs{idx} = gen_sfm_missing_mar(mm_mat, miss_rate);
end

% Create result output directory if needed
if ~exist(sprintf('result/sfm_%s_m', model), 'dir');
    mkdir(sprintf('result/sfm_%s_m', model));
end

for idm = 1:length(MissIDXs)
    % Pick current missing value set
    MissIDX = MissIDXs{idm};
    
    if idm == 1
        mm_mat = mm_mat_MNAR;
        GT = GT_MNAR;
    else
        mm_mat = mm_mat_MAR;
        GT = GT_MAR;
    end

    % translate origin of measurement matrix to zero
    % NOTE: In the real setting, the mean over columns are actually
    %  the mean over rows as we transposed. Thus, there's no difference
    %  in actual implementation as we have access to all points in each
    %  node.
    centroid = mean(mm_mat, 2);
    if idm == 1
        mm_trans = mm_mat - repmat(centroid, [1, N_MNAR]);
    else
        mm_trans = mm_mat - repmat(centroid, [1, N_MAR]);
    end
    
    mm_trans_t = mm_trans';
    MissIDX_t = double(MissIDX');
    
    %% Centralized setting
    M = 3; % latent space dimension
    THRESHc = 1e-3; % convergence precision

    % Get initial values
    m_init_c = get_init_value(mm_trans_t', M, ...
        'ModelType', 'sfm_c', 'VarianceFactor', 10, ...
        'SfMInitPerturb', 0.5, 'MissingIndex', MissIDX_t');

    % Note: we cannot use SVD-base PPCA for data with missing values.
    cm_c = cppca_em(mm_trans_t, M, 'Threshold', THRESHc, ...
            'InitModel', m_init_c, 'ShowObjPer', objfreq_c, ...
            'MissIDX', MissIDX_t);
    
    % Print Result
    fprintf('*** Result (Missing IDX: %d) ***\n', idm);

    % We consider centralized SVD solution as ground truth
    fprintf('### Subspace angle (S/A)             :1e-123456789    15\n');

    % Centralized settting result
    fprintf('(1) GT and (2) PPCA SfM              : %.15f\n', subspace(GT, cm_c.W));
    
    SAc(idm) = subspace(GT, cm_c.W);

    %% Distributed setting
    ETA = 10; % learning parameter
    THRESHd = 1e-3; % convergence precision

    % Run with 5 nodes
    J = 5; % number of nodes

    parfor idr2 = 1:repeats % independent runs (initial value)
        m_init = get_init_value(mm_trans_t', M, ...
                'ModelType', 'sfm_d', 'VarianceFactor', 10, ...
                'SfMInitPerturb', 0.5, 'MissingIndex', MissIDX_t', ...
                'NumberOfNodes', J);

        % D-PPCA
        V = get_sample_assign(J, D, 'IsSfM', true); % node assignment of samples(=frames)
        E = get_adj_graph(J); E = E{2}; % graph topology (1 = compete, 2 = ring)
        cm_dppca = dppca(mm_trans_t, M, V, E, 'InitModel', m_init, ...
               'Eta', ETA, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
               'ZeroMean', true);

        % Distributed settting result
        max_angle = subspace(GT, cm_dppca.W(:,:,1));
        for idj = 2:J
            if max_angle < subspace(GT, cm_dppca.W(:,:,idj))
                max_angle = subspace(GT, cm_dppca.W(:,:,idj));
            end
        end

        fprintf('(1) GT and (2) D-PPCA    (Nodes:%d)   : %.15f (Init: %02d / %f)\n', ...
            J, max_angle, idr2, m_init.VAR);

        save_var_parfor1(sprintf('result/sfm_%s_m/angle_%03d_%03d.mat', model, idm, idr2), max_angle);
        save_var_parfor2(sprintf('result/sfm_%s_m/result_%03d_%03d.mat', model, idm, idr2), cm_dppca, m_init);
    end % parfor: independent runs (different initializations)

    %% Store Result
    for idr2 = 1:repeats
        load(sprintf('result/sfm_%s_m/angle_%03d_%03d.mat', model, idm, idr2));
        SAs(idm, idr2) = cm;
        delete(sprintf('result/sfm_%s_m/angle_%03d_%03d.mat', model, idm, idr2));
    end
end

%% Output overall result
% Note that if we compute 
save(sprintf('result/sfm_%s_m/result_%s_all.mat', model, model), 'SAc', 'SAs');
mean_dist = mean(SAs, 2);

mar_cent = mean(SAc(2:end)) * 180 / pi;
mar_dist = mean(mean_dist(2:end)) * 180 / pi;

mnar_cent = SAc(1) * 180 / pi;
mnar_dist = mean_dist(1) * 180 / pi;

fprintf('### Result MAR  (subspace angle b/w ground truth and centralized PPCA): %f (degrees)\n', mar_cent);
fprintf('### Result MNAR (subspace angle b/w ground truth and centralized PPCA): %f (degrees)\n', mnar_cent);
fprintf('### Result MAR  (subspace angle b/w ground truth and distributed PPCA): %f (degrees)\n', mar_dist);
fprintf('### Result MNAR (subspace angle b/w ground truth and distributed PPCA): %f (degrees)\n', mnar_dist);
