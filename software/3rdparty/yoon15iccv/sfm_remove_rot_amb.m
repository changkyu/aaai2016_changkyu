%% remove_rot_amb
function [M, S] = sfm_remove_rot_amb(M, S)

% Remove ambiguity by finding nullspace (= solving the least square)
QQt = getQQt(M);
[Uq, Dq, ~] = svd(QQt);
Q = Uq * sqrt(Dq);

% find final SfM solution
M = M * Q;
S = Q \ S;

end

%% getQQt
function QQt = getQQt(sfm_M, F)

if nargin < 2
    F = size(sfm_M, 1);
end

A = [];
b = repmat([1; 1; 0], [F/2, 1]);

for idx = 1:ceil(F/2)
    r11 = sfm_M(idx*2-1,1:3)'*sfm_M(idx*2-1,1:3);  % r1*r1'
    r22 = sfm_M(idx*2,1:3)'*sfm_M(idx*2,1:3);      % r2*r2'
    r12 = sfm_M(idx*2-1,1:3)'*sfm_M(idx*2,1:3);    % r1*r2'

    A = [ A; r11(:)'; r22(:)'; r12(:)' ];
end

% Since QQt is symmetric, only 6 variables are independent
idx = [1 2 3; 2 4 5; 3 5 6];
idx = idx(:);

% Form the matrix of independent coefficients
B = zeros(3*F/2,6);
for i=1:6,
  B(:,i) = sum(A(:,idx==i),2);
end

% find QQt
Y = B \ b;
QQt = reshape(Y(idx), [3, 3])';

end
