function [X, n_X, dim_X, dim_Z, GT] = generate_dataset_synthetic( varargin )
% Generate Synthetic Data for ex01
% 
% Description
%  Generate Sytheic Data for experiment 01
%
% Input
%  params = structure(N, M, D, VAR, ...);
%  N     : Number of Data
%  M     : Dimension of subspace
%  D     : Dimension of original data
%  VAR   : Variance of noise
%
% Output
%  data = structure(X, W, Z, E, ...);
%  X     : Observation
%  W     : Projection
%  Z     : Latent variable
%  E     : Noise

N   = varargin{1};
D   = varargin{2};
M   = varargin{3};
VAR = varargin{4};

env_setting;
dir_data = fullfile(dir_dataset,'synthetic/');

filename_data = sprintf('dataset_ex01_synthetic_N%d_D%d_M%d_VAR%1.2f.mat',N,D,M,VAR);
filepath_data = fullfile(dir_data, filename_data);
if exist( filepath_data, 'file' )
    load(filepath_data);
else    
    % random seed
    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s);
    reset(s,0);

    % Z (dim=M) comes from N(0,I)
    Z = randn(s,N,M);

    % W (dim=M) comes from N(0,I)
    % NOTE: W can be an arbitrary matrix, i.e. W = rand(s,D,M);
    W = randn(s,D,M);

    % E (dim=D) comes from N(0,VAR*I)
    E = (sqrt(VAR) .* randn(s,N,D))';

    % Our PPCA model is X = W * Z + E
    X = W * Z' + E;
    GT    = W;
    n_X   = N;
    dim_Z = M;
    dim_X = D;

    if( ~exist(dir_data, 'dir') )
        mkdir(dir_data);
    end
    save(filepath_data, 'X', 'n_X', 'dim_X', 'dim_Z', 'GT');
end

end

