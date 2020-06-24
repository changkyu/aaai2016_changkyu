function model_init = get_init_value( X, M, varargin )
% GET_INIT_VALUE    Get initial value model for PPCA/D-PPCA
%
% Input
%  X          : Input data (measurement matrix)
%  M          : Scalar value for latent space dimension
%  [Optional Parameters]
%  ModelType      : One of model types below
%   cppca         - General centralized PPCA parameters (Default)
%   dppca         - General D-PPCA parameters
%   dppca_dup     - Simple duplicate for given set of parameters
%   sfm_c         - Centralized Structure from Motion 
%                   (initialize by x and y coordinates of the first frame)
%   sfm_d         - Distributed Structure from Motion
%                   (same as sfm_c but use the first frame of each camera)
%   sfm_dsvd      - Distributed Structure from Motion
%                   (do svd locally to initialize; FrameAssignVec required)
%   sfm_drnd      - Distributed Structure from Motion
%                   (use the same random initialization for all cameras)
%  VarianceFactor : Variance factor for random variance (Default 1)
%  NumberOfNodes  : Number of nodes in network
%                   (Only for : dppca, dppca_dup)
%  Init           : Initial model to be duplicated
%                   (Only for : dppca_dup)
%  MissingIndex   : Binary matrix containing missing value information.
%                   0 means available feature, 1 means missing feature.
%                   (Default all ones)
%  SfMInitPerturb : Special parameter for SfM initialization (decides how
%                   much we perturbe initial values regardless of noise)
%  FrameAssignVec : A vector containing frame assignments to nodes.
%                   Required if you use sfm_dsvd as ModelType.
%
% Output
%  model_init : Model structure with appropriate initial parameters
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.01.30 (last modified on 2014/05/13)

% Parse parameters
p = inputParser;
p.StructExpand = false;

defaultType = 'cppca';
expectedTypes = {'cppca', 'dppca', 'dppca_dup', ...
    'sfm_c', 'sfm_d', 'sfm_dsvd', 'sfm_drnd'};

D   = size(X, 1);
W   = rand(D, M)./100;
MU  = zeros(D, 1);
VAR = rand(1)/10;
defaultMODEL = structure(W, MU, VAR);
defaultVARI = 1;
defaultJ = 1;
defaultMissingIndex = zeros(size(X));
defaultSfMInitPerturb = 1;
defaultFrameAssignVec = zeros(D,1);

addParamValue(p,'ModelType',defaultType, ...
    @(x) any(validatestring(x,expectedTypes)));
addParamValue(p,'Init',defaultMODEL);
addParamValue(p,'VarianceFactor',defaultVARI,@isnumeric);
addParamValue(p,'NumberOfNodes',defaultJ,@isnumeric);
addParamValue(p,'MissingIndex',defaultMissingIndex,@isnumeric);
addParamValue(p,'SfMInitPerturb',defaultSfMInitPerturb,@isnumeric);
addParamValue(p,'FrameAssignVec',defaultFrameAssignVec);

parse(p,varargin{:});

% Create initial model
if strcmp(p.Results.ModelType, 'cppca')
    %----------------------------------------------------------------------
    % Centralized PPCA
    %----------------------------------------------------------------------
    model_init.W   = rand(D, M)./100;
    model_init.MU  = zeros(D, 1);
    model_init.VAR = rand(1) / p.Results.VarianceFactor;
elseif strcmp(p.Results.ModelType, 'dppca')
    %----------------------------------------------------------------------
    % Distributed PPCA
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    model_init.W   = repmat(rand(D, M)/100, [1 1 J]);
    model_init.MU  = zeros(D, 1);
    model_init.VAR = 100;
elseif strcmp(p.Results.ModelType, 'dppca_dup')
    %----------------------------------------------------------------------
    % Simply duplicate a model over all nodes
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    model_init.W   = repmat(p.Results.Init.W, [1 1 J]);
    model_init.MU  = p.Results.Init.MU;
    model_init.VAR = p.Results.Init.VAR;
elseif strcmp(p.Results.ModelType, 'sfm_c')
    %----------------------------------------------------------------------
    % Centralized PPCA for SfM
    %----------------------------------------------------------------------
    % get number of points
    N = size(X, 2);
    vari = p.Results.Init.VAR;
    SfMInitPerterb = p.Results.SfMInitPerturb;
    % W is initialized x-y coordinates of the first frame
    % with small variance from uniform distribution scaled by some value
    model_init.W = zeros(N, 3);
    for idx = 1:N
        if p.Results.MissingIndex(1, idx) ~= 1
            model_init.W(idx,:) = [ X(1:2,idx) + ones(2,1)/SfMInitPerterb; 1]';
        else
            model_init.W(idx,:) = [ ones(2,1)/SfMInitPerterb; 1 ]';
        end
    end
    % MU should be initialized as zero
    model_init.MU  = zeros(N, 1);
    % VAR is just a small number
    model_init.VAR = rand(1)/vari;    
elseif strcmp(p.Results.ModelType, 'sfm_d')
    %----------------------------------------------------------------------
    % Distributed PPCA for SfM
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    vari = p.Results.Init.VAR;
    SfMInitPerturb = p.Results.SfMInitPerturb;
    % calculate frames per node
    [FF, N] = size(X);
    fpnode = floor(FF / (2 * J)); 
    % W is initialized x-y coordinates of the first frame AT EACH NODE
    % with small variance from uniform distribution scaled by some value
    model_init.W = zeros(N, 3, J); 
    for idj = 1:J
        past_frames = (fpnode * 2) * (idj - 1);
        model_init.W(:,:,idj) = ...
            [ X( past_frames + 1 : past_frames + 2, : ) + ...
            ones(2,N)/SfMInitPerturb; ones(1,N)]';
    end
    % MU should be initialized as zero
    model_init.MU = zeros(N, 1);
    % VAR is just a small number
    model_init.VAR = rand(1)/vari;  
elseif strcmp(p.Results.ModelType, 'sfm_dsvd')
    %----------------------------------------------------------------------
    % Distributed PPCA for SfM
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    vari = p.Results.Init.VAR;
    % check correct frame assignment
    V = p.Results.FrameAssignVec;
    assert(sum(V) ~= 0, 'Wrong frame assignment vector!');
    % W is initialized by doing local SfM in each node
    [~, N] = size(X);
    model_init.W = zeros(N, 3, J); 
    for idj = 1:J
        [~, tmp, ~, ~, ~, ~] = affine_sfm(X(V==idj,:));
        model_init.W(:,:,idj) = tmp';
    end
    % MU should be initialized as zero
    model_init.MU = zeros(N, 1);
    % VAR is just a small number
    model_init.VAR = rand(1)/vari;
elseif strcmp(p.Results.ModelType, 'sfm_drnd')
    %----------------------------------------------------------------------
    % Distributed PPCA for SfM
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    vari = p.Results.Init.VAR;
    % W is initialized by doing local SfM in each node
    [~, D] = size(X);
    model_init.W = repmat(rand(D, M)/100, [1 1 J]);
    % MU should be initialized as zero
    model_init.MU = zeros(D, 1);
    % VAR is just a small number
    model_init.VAR = rand(1)/vari;
else
    error('invalid param');
end

end
