function Networks = get_adj_graph_ex( J )
% GET_ADJ_GRAPH   Get Adjacent Graphs
%
% Description
%  Returns 4-5 adjacent matrices for different network topology structures 
%  as a cell array given J nodes. Network topologies are (1) complete (2) 
%  ring (3) star (4) chain (5) clusters. Cluster is only available when 
%  the condition J > 4 && mod(J,4) == 0 holds.
%
% Input
%  J    : Number of nodes in network
%
% Output
%  Networks : 4 or 5 element cell array each with an adjacent matrix
%
% Implemented/Modified
%  by     Changkyu Song (changkyu.song@rutgers.edu)
%
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2011.10.07 (last modified on 2014/02/14)

% exception handling
if J < 2
    error('J < 2');    
end

if J > 4 && mod(J, 4) == 0
    n_Networks = 6;
else    
    n_Networks = 4;
end
Networks = cell(n_Networks,1);

% complete
E1 = ones(J) - eye(J);
Networks{1} = struct('name', 'complete', 'adj', E1);

% ring
E2 = diag(ones(1,J-1),1); E2(J,1) = 1; E2 = E2+E2'; 
Networks{2} = struct('name', 'ring', 'adj', E2);

% star
E3 = [ones(1,J); zeros(J-1,J)]; E3 = E3 + E3';
E3(1,1) = 0;
Networks{3} = struct('name', 'star', 'adj', E3);

% chain
E4 = diag(ones(1,J-1),1); E4 = E4+E4'; 
Networks{4} = struct('name', 'chain', 'adj', E4);

% cluster
if J > 4 && mod(J, 4) == 0
    E5 = [ones(1,J/2) zeros(1,J/2); zeros(J/2-1,J); ...
          zeros(1,J/2) ones(1,J/2) ; zeros(J/2-1,J)]; 
    E5(1,J/2+1) = 1; E5 = E5 + E5';
    Networks{5} = struct('name', 'cluster1', 'adj', E5);
    
    E6 = [ones(J/2, J/2) zeros(J/2, J/2); ...
          zeros(J/2, J/2) ones(J/2, J/2);]; 
    E6(J/2,J/2+1) = 1;
    E6(J/2+1,J/2) = 1;
    Networks{6} = struct('name', 'cluster2', 'adj', E6);    
end

end
