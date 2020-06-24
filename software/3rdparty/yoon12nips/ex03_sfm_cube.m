%--------------------------------------------------------------------------
% SfM Synthetic data demonstration (Cube)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2012/03/15)
%--------------------------------------------------------------------------

model = 'cube';
noises = [1e-5, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];

if ~exist(sprintf('data/%s', model), 'dir')
    error('Please run cube data generation script first!');
elseif ~exist(sprintf('result/sfm_%s', model), 'dir')
    mkdir(sprintf('result/sfm_%s', model));
end

datagens = 20; % Number of independent runs of perturbed data
repeats = 5; % Number of independent runs of different initial values

SAs = zeros(length(noises), 2, datagens, repeats);

for idn = 1:length(noises) % noise level
    %--------------------------------------------------------------
    % Data
    % Frequency for intermediate object function value output
    objfreq_c = 0;
    objfreq_d = 0;

    measurement_matrices = cell(20, 1);
    for idr = 1:datagens
        datastr = sprintf('data/%s/%03d_%03d.mat', model, idn, idr);
        load(datastr);
        measurement_matrices{idr} = measurement_matrix;
    end

    for idr = 1:datagens % independent runs (noise level)
        datastr = sprintf('data/%s/%03d_%03d.mat', model, idn, idr);
        disp(['[[ Current Model: ' datastr ' ]]']);
        mm_mat = measurement_matrices{idr};
        mm_mat_t = mm_mat';

        % get dimension of measurement matrix
        [D, N] = size(mm_mat);

        % Measurement matrix should be in the specific form
        if mod(D, 2) ~= 0
            error('Measurement matrix should be (2 x #frames)x(# points) form!');
        end

        % translate origin of measurement matrix to zero
        % NOTE: In the real setting, the mean over columns are actually
        %  the mean over rows as we transposed. Thus, there's no difference
        %  in actual implementation as we have access to all points.
        centroid = mean(mm_mat, 2);
        mm_trans = mm_mat - repmat(centroid, [1, N]);
        mm_trans_t = mm_trans';

        %--------------------------------------------------------------
        % Centralized setting
        [sfm_M, sfm_S, sfm_b, U3, W3, V3] = affine_sfm(mm_mat);

        M = 3; % latent space dimension
        THRESHc = 1e-3; % convergence precision

        % Print Result
        fprintf('*** Result (Run: %d) ***\n', idr);

        % We consider centralized SVD solution as ground truth
        fprintf('### Subspace angle (S/A)             :1e-123456789    15\n');

        % Centralized settting result
        fprintf('(1) GT and (2) SVD SfM               : %.15f\n', ...
            subspace(GT, sfm_S'));

        %--------------------------------------------------------------
        % Distributed setting
        ETA = 10; % learning parameter
        THRESHd = 1e-3; % convergence precision

        % Run with 5 nodes
        J = 5; % number of nodes

        parfor idr2 = 1:repeats % independent runs (initial value)
            m_init = get_init_value(mm_trans_t', M, ...
                'ModelType', 'sfm_d', 'VarianceFactor', 10, ...
                'SfMInitPerturb', 0.5, 'NumberOfNodes', J);

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

            save_var_parfor1(sprintf('result/sfm_%s/angle_%03d_%03d_%03d.mat', ...
                model, idn, idr, idr2), max_angle);
            save_var_parfor2(sprintf('result/sfm_%s/result_%03d_%03d_%03d.mat', ...
                model, idn, idr, idr2), cm_dppca, m_init);
        end % parfor: independent runs (initial value)

        %----------------------------------------------------------
        % Store Result
        for idr2 = 1:repeats
            SAs(idn, 1, idr, idr2) = subspace(GT, sfm_S');
            load(sprintf('result/sfm_%s/angle_%03d_%03d_%03d.mat', ...
                model, idn, idr, idr2));
            SAs(idn, 2, idr, idr2) = cm;
            delete(sprintf('result/sfm_%s/angle_%03d_%03d_%03d.mat', ...
                model, idn, idr, idr2));
        end
    end % for: independent runs (noise)
end % noise level

save(sprintf('result/sfm_%s/result_%s_all.mat', model, model), 'SAs');

show_result_3_sfm_cube;
