%--------------------------------------------------------------------------
% SfM Synthetic data demonstration (Cube)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2015/04/02)
%--------------------------------------------------------------------------
clear; close all;

% Noise levels
noises = [1e-5, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
% Number of independent runs of perturbed data
reps = 20;
model = 'cube';

if ~exist('data', 'dir'), mkdir('data'); end
if ~exist(sprintf('data/%s', model), 'dir')
    error('Please run sfm_gen_cube_data.m first!');
end
if ~exist('results', 'dir'), mkdir('data'); end
if ~exist(sprintf('results/sfm_%s', model), 'dir')
    mkdir(sprintf('results/sfm_%s', model));
end

% store subspace angles (i.e. FINAL RESULT)
SAs = zeros(length(noises), 7, reps);

% for each noise level
for idn = 1:length(noises)
    % Load data
    measurement_matrices = cell(reps, 1);
    for idr = 1 : reps
        datastr = sprintf('data/%s/%03d_%03d.mat', model, idn, idr);
        load(datastr);
        measurement_matrices{idr} = measurement_matrix;
    end

    parfor idr = 1 : reps % independent runs (per noise level)
        datastr = sprintf('data/%s/%03d_%03d.mat', model, idn, idr);
        disp(['[[ Current Model: ' datastr ' ]]']);
        mm_mat = measurement_matrices{idr};

        % get dimension of measurement matrix
        [D, N] = size(mm_mat);

        % Measurement matrix should be stack of 2D point frames
        if mod(D, 2) ~= 0
            error('Measurement matrix should be (2 x #frames)x(# points) form!');
        end

        % translate origin of measurement matrix to zero
        centroid = mean(mm_mat, 2);
        mm_trans = mm_mat - repmat(centroid, [1, N]);
        X = mm_trans';

        %--------------------------------------------------------------
        % Reference run
        [sfm_M, sfm_S, sfm_b, U3, W3, V3] = sfm_affine(mm_mat);

        % Print Result
        fprintf('*** Result (Run: %d) ***\n', idr);

        % Reference centralized SVD solution
        fprintf('### Subspace angle (S/A)             :1e-123456789    15\n');

        % Centralized settting result
        angle0 = subspace(GT, sfm_S');
        fprintf('(1) GT and (2) SVD SfM               : %.15f\n', angle0);

        %------------------------------------------------------------------
        % Options
        M = 3; % latent space dimension is x-y-z coordinate of structure
        NV = 5; % number of cameras
        ETA = 10; % learning parameter
        
        THRESHc = 1e-3; % convergence precision
        THRESHd = 1e-3; % convergence precision
        
        objfreq_c = 0; % intermediate object function output (centralized)
        objfreq_d = 0; % intermediate object function output (distributed)
        
        %------------------------------------------------------------------
        % PPCA, BPCA, D-PPCA, D-BPCA
        
        % node assignment of samples(=frames)
        V = get_sample_assign(NV, D, 'IsSfM', true);
        % graph topology (1 = compete, 2 = ring)
        E = get_adj_graph(NV); E = E{2}; 
        
        m_init_c = get_init_value(X, M, 'ModelType', 'cppca');
        m_init_d = get_init_value(X, M, 'ModelType', 'sfm_d', ...
            'NumberOfNodes', NV, 'SampleAssignVec', V);  

        [cm1, cm2] = expr_run_cppca(X, M, THRESHc, objfreq_c, m_init_c, true);
        [cm3, cm4] = expr_run_dppca(X, M, V, E, ETA, THRESHd, objfreq_d, m_init_d, true);
        [cm5, cm6] = expr_run_vbpca(X, M);

        % Centralized settting result
        angle1 = subspace(GT, cm1.W);
        %angle2 = subspace(GT, cm2.mW);

        % VBPCA (Ilin & Raiko)
        %angle5 = subspace(GT, cm5.W);
        %angle6 = subspace(GT, cm6.mW);
        
        % Distributed settting result
        angle3 = subspace(GT, cm3.W{1});
        angle4 = subspace(GT, cm4.W{1});
        for idj = 2:NV
            if angle3 < subspace(GT, cm3.W{idj})
                angle3 = subspace(GT, cm3.W{idj});
            end
            if angle4 < subspace(GT, cm4.W{idj})
               angle4 = subspace(GT, cm4.W{idj});
            end
        end

        fprintf('(1) GT and (2) PPCA                   : %.15f (Init: %02d)\n', ...
            angle1, idr);
%         fprintf('(1) GT and (2) BPCA                   : %.15f (Init: %02d)\n', ...
%             angle2, idr);
        
        fprintf('(1) GT and (2) D-PPCA    (Nodes:%d)   : %.15f (Init: %02d)\n', ...
            NV, angle3, idr);
        fprintf('(1) GT and (2) D-PPCA-ANT(Nodes:%d)   : %.15f (Init: %02d)\n', ...
            NV, angle4, idr);
        
%         fprintf('(1) GT and (2) PPCA (IR)              : %.15f (Init: %02d)\n', ...
%             angle5, idr);
%         fprintf('(1) GT and (2) BPCA (IR)              : %.15f (Init: %02d)\n', ...
%             angle6, idr);        
        
%         angles = [angle0 angle1 angle2 angle3 angle4 angle5 angle6];
        angles = [angle0 angle1 0 angle3 angle4 0 0];

        parfor_var_save9( ...
            sprintf('results/sfm_%s/%03d_%03d.mat', model, idn, idr), ...
            cm1, cm2, cm3, cm4, cm5, cm6, m_init_c, m_init_d, angles);
    end % for: independent runs (noise)
    
    %----------------------------------------------------------------------
    % Store Result
    for idr = 1 : reps
        tmp = load(sprintf('results/sfm_%s/%03d_%03d.mat', model, idn, idr));
        SAs(idn, :, idr) = tmp.angles;
    end
end % noise level

save(sprintf('results/sfm_%s/%s_all.mat', model, model), 'SAs');

%% ------------------------------------------------------------------------
% Show result
test6_sfm_cube_show;
