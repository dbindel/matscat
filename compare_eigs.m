% [ll,diff] = compare_eigs(ll1, ll2, tol)
%
% Compare the results of two eigenvalue vectors (ll1 and ll2).
% Inputs:
%   ll1, ll2 - eigenvalue lists
%   tol      - matching tolerance
% Outputs:
%   ll   - all elements of ll1 within tol of some element of ll2
%   diff - diff(i) is min(abs(ll1(i) - ll2));

function [ll,diff] = compare_eigs(ll1, ll2, tol)

if size(ll1,1) == 1, ll1 = ll1.'; end
if size(ll2,1) ~= 1, ll2 = ll2.'; end
e1 = 0*ll1 + 1;
e2 = 0*ll2 + 1;
diff = min(abs(ll1*e2 - e1*ll2), [], 2);
ll   = ll1(find(diff < tol));
