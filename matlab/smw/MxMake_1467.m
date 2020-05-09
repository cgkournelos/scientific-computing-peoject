% Author: Christos Gkournelos
% Date: 15/01/2019 
% 
% Mx_Make function for providing a known matrix
%
function A = MxMake_1467(mx_id,n)
  % 
  %  Input error checking 
  % 
  if nargin ~= 2
    error('Wrong input. \nThe function requires 2 inputs', -1) 
  end
  %
  % Hadamard 
  %
  if(isequal('had', mx_id) == 1)
    A = hadamard(n);
  end %if
  %
  % Upper Triagonal Hadamard 
  %
  if(isequal('trihad', mx_id) == 1)
    A = triu(hadamard(n));
  end %if
  %
  % Triagonal Toeplitz 
  %
  if(isequal('toep', mx_id) == 1)
    if( n < 3)
      error('Wrong input. \nFor tridiagonal Toeplitz must nsize must be higher to 3', -1) 
    end %if
    A = toeplitz([4,-1,zeros(1,n-2)]);
  end %if
  %
  % Dense matrix
  %
  if(isequal('mc', mx_id) == 1)
    for i = 1 : n
        for j = 1 : n
          if i==j
            A(i,i)=1+i;
            break;
          end %if
          A(i,j)=1/(abs(i+j)^2);
        end %for j
     end %for i 
  end %if
  %
  % Sparse matrix wathen
  %
  if(isequal('wathen', mx_id) == 1)   
    A = gallery('wathen', n, 11);
  end %if
  %
  % Sparse matrix CollegeMsg
  %
  if(isequal('CollegeMsg', mx_id) == 1)   
    load(websave("CollegeMsg.mat",'https://sparse.tamu.edu/mat/SNAP/CollegeMsg.mat'));
    A = Problem.A;
  end %if
end %for
