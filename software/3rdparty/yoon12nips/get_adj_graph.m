function Earr = get_adj_graph( J )
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
%  Earr : 4 or 5 element cell array each with an adjacent matrix
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2011.10.07 (last modified on 2014/02/14)

% exception handling
if J == 1
    Earr = {[], [], [], [], []};
    return;
end

% complete
E1 = ones(J) - eye(J);

% ring
E2 = diag(ones(1,J-1),1); E2(J,1) = 1; E2 = E2+E2'; 

% star
E3 = [ones(1,J); zeros(J-1,J)]; E3 = E3 + E3';
E3(1,1) = 0;

% chain
E4 = diag(ones(1,J-1),1); E4 = E4+E4'; 

% cluster
if J > 4 && mod(J, 4) == 0
    E5 = [ones(1,J/2) zeros(1,J/2); zeros(J/2-1,J); ...
          zeros(1,J/2) ones(1,J/2) ; zeros(J/2-1,J)]; 
    E5(1,J/2+1) = 1; E5 = E5 + E5';
    Earr = {E1, E2, E3, E4, E5};
else
    Earr = {E1, E2, E3, E4};
end

end
