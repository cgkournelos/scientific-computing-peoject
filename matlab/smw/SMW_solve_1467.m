% Author: Christos Gkournelos
% Date: 16/01/2019 
% 
% SMW_solve function for solving linear systems 
% using the Sherman Morrison formula 
%
function x = SMW_solve_1467(A, b, M, P, Q, sdir)
  % 
  %  Input error checking 
  % 
  if (nargin < 5 || nargin   > 6)
    error('Wrong input. \nThe function requires 5 or 6 inputs', -1) 
  end %if
  %
  % Input params M, P and Q 
  % according sdir input 
  %
  n = length(A);
  if(exist('sdir','var'))
    if (isequal(sdir,'colwise'))
      M = diag(diag(A));
      P = (A - M)';
      Q = eye(n);
    elseif (isequal(sdir,'rowwise'))
      M = diag(diag(A));
      Q = (A - M)';
      P = eye(n);
    end %if
  end %if
  %
  % Initialization
  %
  x = M\b;
  
  P = M\P;
  %
  % Recursion until n-1
  %
  for i = 1 : (n-1)
    x = x - P(:,i) * ( Q(:,i)' * x ) / ( 1 + Q(:,i)' * P(:,i) );
    for j = (i + 1) : n
      P(:,j) = P(:,j) - ( ( Q(:,i)' * P(:,j) ) / ( 1 + Q(:,i)' * P(:,i) ) ) * P(:,i);
    end %  for j
  end % for i
  %
  % Final step for computing x
  %
  x = x - P(:,n) * (Q(:,n)' * x) / (1 + Q(:,n)' * P(:,n));
end %SMW_solve_1467