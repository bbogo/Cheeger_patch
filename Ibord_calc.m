function Ibord = Ibord_calc(p,T)

% simple boundary detection
% idea: find points which lie on edges belonging to only one triangle

ap = size(p,1);

A = min(sparse(T(:,1),T(:,2),1,ap,ap)+sparse(T(:,2),T(:,3),1,ap,ap)+sparse(T(:,3),T(:,1),1,ap,ap),1);
A = min(A+A',1);

B = A^2.*A==1;
Ibord = find(sum(B,2)>0);
