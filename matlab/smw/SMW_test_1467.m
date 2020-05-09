% Author: Christos Gkournelos
% Date: 17/01/2019 
% 
% A script for testing SMW_solve function
%
clc; clear;
%
% Compute the initial x_sol
%
for i = 1: (1899/2)
  x_sol(2*i-1)=1;
  x_sol(2*i)=((-1)^(i+1)*1/(2*i));
end %for
x_sol(1899)=1;
%
% For A as a Hadamard Matrix
%
A.had = MxMake_1467('had',64);
b.had = A.had * x_sol(1:length(A.had))';
x.had = SMW_solve_1467(A.had, b.had, 1, 1, 1, 'colwise');
k.had = cond(A.had,inf);
k_est.had = condest(A.had);
a_post_error.had = norm(b.had-A.had * x.had) / ((norm(A.had) * norm(x.had)) + norm(b.had));
implicit_fw_error.had = 2 * k.had * a_post_error.had;
explicit_error.had = norm(x.had - x_sol, inf) / norm(x_sol, inf);
% Matlab impl
mat_x.had = A.had \ b.had;
mat_a_post_error.had = norm(b.had-A.had * mat_x.had) / ((norm(A.had) * norm(mat_x.had)) + norm(b.had));
mat_implicit_fw_error.had = 2 * k.had * mat_a_post_error.had;
mat_explicit_error.had = norm(mat_x.had - x_sol, inf) / norm(x_sol, inf);
%
% For A as an Upper Triagonal Hadamard 
%
A.trihad = MxMake_1467('trihad',64);
b.trihad = A.trihad * x_sol(1:length(A.trihad))';
x.trihad = SMW_solve_1467(A.trihad, b.trihad, 1, 1, 1, 'colwise');
k.trihad = cond(A.trihad,inf);
k_est.trihad = condest(A.trihad);
a_post_error.trihad = norm(b.trihad-A.trihad * x.trihad) / ((norm(A.trihad) * norm(x.trihad)) + norm(b.trihad));
implicit_fw_error.trihad = 2 * k.trihad * a_post_error.trihad;
explicit_error.trihad = norm(x.trihad - x_sol, inf) / norm(x_sol, inf);
% Matlab impl
mat_x.trihad = A.trihad \ b.trihad;
mat_a_post_error.trihad = norm(b.trihad-A.trihad * mat_x.trihad) / ((norm(A.trihad) * norm(mat_x.trihad)) + norm(b.trihad));
mat_implicit_fw_error.trihad = 2 * k.trihad * mat_a_post_error.trihad;
mat_explicit_error.trihad = norm(mat_x.trihad - x_sol, inf) / norm(x_sol, inf);
%
% For A as an Triagonal Toeplitz 
%
A.toep = MxMake_1467('toep',64);
b.toep = A.toep * x_sol(1:length(A.toep))';
x.toep = SMW_solve_1467(A.toep, b.toep,1, 1, 1, 'colwise');
k.toep = cond(A.toep,inf);
k_est.toep = condest(A.toep);
a_post_error.toep = norm(b.toep-A.toep * x.toep) / ((norm(A.toep) * norm(x.toep)) + norm(b.toep));
implicit_fw_error.toep = 2 * k.toep * a_post_error.toep;
explicit_error.toep = norm(x.toep - x_sol, inf) / norm(x_sol, inf);
% Matlab impl
mat_x.toep = A.toep \ b.toep;
mat_a_post_error.toep = norm(b.toep-A.toep * mat_x.toep) / ((norm(A.toep) * norm(mat_x.toep)) + norm(b.toep));
mat_implicit_fw_error.toep = 2 * k.toep * mat_a_post_error.toep;
mat_explicit_error.toep = norm(mat_x.toep - x_sol, inf) / norm(x_sol, inf);
%
% For A as an Dense matrix
%
A.mc = MxMake_1467('mc',400);
b.mc = A.mc * x_sol(1:length(A.mc))';
x.mc = SMW_solve_1467(A.mc, b.mc, 1, 1, 1, 'colwise');
k.mc = cond(A.mc,inf);
k_est.mc = condest(A.mc);
a_post_error.mc = norm(b.mc-A.mc * x.mc) / ((norm(A.mc) * norm(x.mc)) + norm(b.mc));
implicit_fw_error.mc = 2 * k.mc * a_post_error.mc;
explicit_error.mc = norm(x.mc - x_sol, inf) / norm(x_sol, inf);
% Matlab impl
mat_x.mc = A.mc \ b.mc;
mat_a_post_error.mc = norm(b.mc-A.mc * mat_x.mc) / ((norm(A.mc) * norm(mat_x.mc)) + norm(b.mc));
mat_implicit_fw_error.mc = 2 * k.mc * mat_a_post_error.mc;
mat_explicit_error.mc = norm(mat_x.mc - x_sol, inf) / norm(x_sol, inf);
%
%
% For A as an sparse matrix wathen 
% matrix size n = 443
%
A.wathen = MxMake_1467('wathen',12);
b.wathen = A.wathen * x_sol(1:length(A.wathen))';
x.wathen = SMW_solve_1467(A.wathen, b.wathen, 1, 1, 1, 'colwise');
k.wathen = cond(A.wathen,inf);
k_est.wathen = condest(A.wathen);
a_post_error.wathen = normest(b.wathen-A.wathen * x.wathen) / ((normest(A.wathen) * normest(x.wathen)) + normest(b.wathen));
implicit_fw_error.wathen = 2 * k.wathen * a_post_error.wathen;
explicit_error.wathen = normest(x.wathen - x_sol, inf) / normest(x_sol, inf);
% Matlab impl
mat_x.wathen = A.wathen \ b.wathen;
mat_a_post_error.wathen = normest(b.wathen-A.wathen * mat_x.wathen) / ((normest(A.wathen) * normest(mat_x.wathen)) + normest(b.wathen));
mat_implicit_fw_error.wathen = 2 * k.wathen * mat_a_post_error.wathen;
mat_explicit_error.wathen = normest(mat_x.wathen - x_sol, inf) / normest(x_sol, inf);
%
% For A as an sparse matrix CollegeMsg 
%
A.coll = MxMake_1467('CollegeMsg', 1899);
A.coll = eye(1899) - 0.85*A.coll;
b.coll = A.coll * x_sol';
x.coll = SMW_solve_1467(A.coll, b.coll, 1, 1, 1, 'colwise');
k.coll = cond(A.coll,inf);
k_est.coll = condest(A.coll);
a_post_error.coll = norm(b.coll-A.coll * x.coll) / ((norm(A.coll) * norm(x.coll)) + norm(b.coll));
implicit_fw_error.coll = 2 * k.coll * a_post_error.coll;
explicit_error.coll = norm(x.coll - x_sol, inf) / norm(x_sol, inf);
% Matlab impl
mat_x.coll = A.coll \ b.coll;
mat_a_post_error.coll = norm(b.coll-A.coll * mat_x.coll) / ((norm(A.coll) * norm(mat_x.coll)) + norm(b.coll));
mat_implicit_fw_error.coll = 2 * k.coll * mat_a_post_error.coll;
mat_explicit_error.coll = norm(mat_x.coll - x_sol, inf) / norm(x_sol, inf);