% Author: Christos Gkournelos
% Date: 15/01/2019 
% 
% A script for testing Google's PageRank
%
clc; clear;
%
% Load the SNAP/CollegeMsg structure 
%
load(websave("CollegeMsg.mat",'https://sparse.tamu.edu/mat/SNAP/CollegeMsg.mat'));
A = Problem.A;
% 
%  Create a plot for visualizing A
%
figure(1)
spy(A,'r')
% set(gcf, 'Position', get(0, 'Screensize'));
ac = gca;
ac.FontSize = 18;

% PageRank impl
%
n = length(A);
[from_index,to_index] = find(A);
G = sparse(from_index,to_index, 1, n, n);
% fix the dangling nodes values
for j=1:n
  if(sum(G(:,j)) == 0)
    G(:,j) = 1/n;
  end % if
end % for
c = sum(G);
D = spdiags(1./c',0,n,n);
p =0.85;
delta = (1-p) / n;
e =  ones(n,1);
I = speye(n,n);
x = (I - p*G*D)\(delta*e); 
%
% Find top 20 nodes
%
[x_sort,x_index] = sort(x,'descend');
top20_table = table(x_index(1:20),x_sort(1:20),'VariableNames',{'node','value_PR'});
%
% Evalutate the condition value (k) of 
% (I - p*G*D) matrix for multiple p 
%
p = [0.25;0.45;0.65;0.85;0.95;0.99];
for i = 1:6
  k(i) = cond((I - p(i)*full(G)*full(D)),Inf);
end % for
