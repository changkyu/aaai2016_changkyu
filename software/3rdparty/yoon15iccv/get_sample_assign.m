function V = get_sample_assign( J, N, varargin )
% GET_SAMPLE_ASSIGN   Get node assignment vector of samples for the
%                     distributed models
%
% INPUT
%  J       Number of nodes
%  N       Number of frames
%  [Optional Parameters]
%  IsSfM   false for general case, true for SfM case (N is a multiple of 2)
%          (Default: false)
%
% OUTPUT
%  V       N x 1 dimensional vector containing node assignment of samples
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.01.30 (last modified on 2014/02/14)

% Parse parameters
p = inputParser;
p.StructExpand = false;

defaultIsSfM = false;

addParamValue(p,'IsSfM',defaultIsSfM,@islogical);

parse(p,varargin{:});

% Node assignment
if nargin > 2 && p.Results.IsSfM
    % ... to each frame (x,y)
    frames = N / 2;
    if frames < J
        error('cannot process because %d frames < %d nodes!', frames, J);
    else
        frame_per_node = floor(frames / J);
        for idy = 1:(J - 1)
            V_tmp(frame_per_node*(idy-1)+1:frame_per_node*idy,1) = idy;
        end
        V_tmp(frame_per_node*(J-1)+1:frames,1) = J;
    end
    V = zeros(frames * 2, 1);
    V(1:2:frames*2-1,1) = V_tmp;
    V(2:2:frames*2,1) = V_tmp;
else
    % ... to each sample
    samples_per_node = floor(N / J);
    for idx = 1:(J - 1)
        V(samples_per_node*(idx-1)+1:samples_per_node*idx,1) = idx;
    end
    V(samples_per_node*(J-1)+1:N,1) = J;
end

end
