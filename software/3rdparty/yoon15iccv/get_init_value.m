function model_init = get_init_value( X, M, varargin )
% GET_INIT_VALUE    Get initial value model for PPCA/D-PPCA
%
% Input
%  X          : Input data (measurement matrix)
%  M          : Scalar value for latent space dimension
%  [Optional Parameters]
%  ModelType      : One of model types below
%   cppca         - General centralized PPCA initialization (Default)
%   dppca         - General distributed PPCA initialization
%                   (use PPCA on each local node; FrameAssignVec required)
%   sfm_c         - Use first frame's points to initialize W matrix
%   sfm_d         - Use first frame's points to initialize W matrix in each
%                   node. The first frame is the first one in each node.
%  NumberOfNodes  : Number of nodes in network
%                   (Required for : dppca)
%  SampleAssignVec: A vector containing frame assignments to nodes.
%                   (Required for : dppca)
%
% Output
%  model_init : Model structure with appropriate initial parameters
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.01.30 (last modified on 2015/03/23)

% Parse parameters
p = inputParser;
p.StructExpand = false;

defaultType = 'cppca';
expectedTypes = {'cppca', 'dppca', 'dppca_r' ,'sfm_c', 'sfm_d'};

[D, N] = size(X);
defaultJ = 1;
defaultFrameAssignVec = ones(N, 1);

addParameter(p,'ModelType',defaultType, ...
    @(x) any(validatestring(x,expectedTypes)));
addParameter(p,'NumberOfNodes',defaultJ,@isnumeric);
addParameter(p,'SampleAssignVec',defaultFrameAssignVec);

parse(p,varargin{:});

% Create initial model
if strcmp(p.Results.ModelType, 'cppca')
    %----------------------------------------------------------------------
    % Centralized PPCA
    %----------------------------------------------------------------------
    model_init.W = orth(randn(D, M));
    model_init.MU = zeros(D, 1);
    model_init.VAR = 1;
elseif strcmp(p.Results.ModelType, 'dppca')
    %----------------------------------------------------------------------
    % Distributed PPCA using local PPCA for initialization
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    V = p.Results.SampleAssignVec;
    model_init.W = cell(J, 1);
    model_init.MU = cell(J, 1);
    model_init.VAR = cell(J, 1);
    for j = 1 : J
        % When there's missing data, use EM-PPCA. Otherwise, use SVD-PPCA.
        if sum(sum(isnan(X))) > 0
            cm = cppca_em(X(:,V == j), M, 'ShowObjPer', 0);
        else
            cm = cppca(X(:,V == j), M);
        end
        model_init.W{j} = cm.W;
        model_init.MU{j} = cm.MU;
        model_init.VAR{j} = cm.VAR;
    end
elseif strcmp(p.Results.ModelType, 'dppca_r')
    %----------------------------------------------------------------------
    % Distributed PPCA using fully random initialization
    %----------------------------------------------------------------------
    J = p.Results.NumberOfNodes;
    model_init.W = cell(J, 1);
    model_init.MU = cell(J, 1);
    model_init.VAR = cell(J, 1);
    for j = 1 : J
        model_init.W{j} = orth(randn(D, M));
        model_init.MU{j} = zeros(D, 1);
        model_init.VAR{j} = 1;
    end    
elseif strcmp(p.Results.ModelType, 'sfm_c')
    %----------------------------------------------------------------------
    % Centralized PPCA (Transposed) for SfM
    %----------------------------------------------------------------------
    X = X';
    % Parameters
    VARI = 10;
    % get number of points
    N = size(X, 2);
    % W is initialized x-y coordinates of the first frame
    % with small variance from uniform distribution scaled by some value
    model_init.W = zeros(N, 3);
    for n = 1 : N
        if ~isnan(X(1, n))
            model_init.W(n,:) = [ X(1:2,n) + ones(2,1)/VARI; 1]';
        else
            model_init.W(n,:) = [ ones(2,1)/VARI; 1 ]';
        end
    end
    % MU should be initialized as zero
    model_init.MU  = zeros(N, 1);
    % VAR is just a small number
    model_init.VAR = rand(1);
elseif strcmp(p.Results.ModelType, 'sfm_d')
    %----------------------------------------------------------------------
    % Distributed PPCA for SfM
    %----------------------------------------------------------------------
    X = X';
    % Parameters
    VARI1 = 10;
	VARI2 = 10;
	VARI3 = 2;
    D_PARAM = 0.5;
    % Number of nodes
    J = p.Results.NumberOfNodes;
    % calculate frames per node
    [FF, N] = size(X);
    fpnode = floor(FF / (2 * J)); 
    % W is initialized x-y coordinates of the first frame AT EACH NODE
    % with small variance from uniform distribution scaled by some value
    model_init.W = cell(J, 1);
    model_init.MU = cell(J, 1);
    model_init.VAR = cell(J, 1);    
    for j = 1 : J
        model_init.W{j} = zeros(N, 3); 
        past_frames = (fpnode * 2) * (j - 1);
        model_init.W{j} = ...
            [ X( past_frames + 1 : past_frames + 2, : ) + ones(2,N)/VARI1; ones(1,N)]';
        % cancel out if missing
        NcI = isnan(X(past_frames + 1, :));
        model_init.W{j}(NcI,:) = [ ones(2,sum(NcI))/VARI2; ones(1,sum(NcI)) ]';
        % MU should be initialized as zero
        model_init.MU{j} = zeros(N, 1);
        % VAR is just a small number
        model_init.VAR{j} = rand(1)/VARI3 + D_PARAM;    
    end
end

end
