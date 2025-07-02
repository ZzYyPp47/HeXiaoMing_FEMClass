% Do Newton's iteration
function uh_new = Newton_iteration(iteration_num,tol,A,b,u0,B_m_edges,B_f_edges_u,Nf_u,Nf_p,Nm,P,T,Pb_u,Pb_p,Nb_u,basis_u,Tb_u,func_boundary1,func_boundary2,func_real_p)
uh_new = u0;
for ii = 1:1:iteration_num
    %% Modify to meet your needs
    uh_old = uh_new;
    [uh_old_1,uh_old_2,~] = deal(uh_old(1:Nf_u,:),uh_old(Nf_u + 1:2 * Nf_u,:),uh_old(2 * Nf_u + 1:end,:));
    AN1 = assembler_A_Newton(uh_old_1, 2, 1, 0, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, basis_u, basis_u, Tb_u, Tb_u, 0, 0, 0, 0);
    AN2 = assembler_A_Newton(uh_old_1, 2, 0, 0, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, basis_u, basis_u, Tb_u, Tb_u, 1, 0, 0, 0);
    AN3 = assembler_A_Newton(uh_old_2, 2, 0, 0, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, basis_u, basis_u, Tb_u, Tb_u, 0, 1, 0, 0);
    AN4 = assembler_A_Newton(uh_old_1, 2, 0, 1, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, basis_u, basis_u, Tb_u, Tb_u, 0, 0, 0, 0);
    AN5 = assembler_A_Newton(uh_old_2, 2, 1, 0, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, basis_u, basis_u, Tb_u, Tb_u, 0, 0, 0, 0);
    AN6 = assembler_A_Newton(uh_old_2, 2, 0, 1, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, basis_u, basis_u, Tb_u, Tb_u, 0, 0, 0, 0);
    O1 = sparse(Nf_p,Nf_p);
    O2 = sparse(Nf_u,Nf_p);
    O3 = sparse(Nf_u,Nf_p);
    AN = [AN1 + AN2 + AN3,AN4,O2;AN5,AN6 + AN2 + AN3,O3;O2',O3',O1];
    bN1 = assembler_b_Newton(uh_old_1, uh_old_1, 0, 0, 1, 0, 2, Nf_u, Nm, P, T, Tb_u, Nb_u, Nb_u, basis_u, basis_u, Tb_u, 0, 0);
    bN2 = assembler_b_Newton(uh_old_2, uh_old_1, 0, 0, 0, 1, 2, Nf_u, Nm, P, T, Tb_u, Nb_u, Nb_u, basis_u, basis_u, Tb_u, 0, 0);
    bN3 = assembler_b_Newton(uh_old_1, uh_old_2, 0, 0, 1, 0, 2, Nf_u, Nm, P, T, Tb_u, Nb_u, Nb_u, basis_u, basis_u, Tb_u, 0, 0);
    bN4 = assembler_b_Newton(uh_old_2, uh_old_2, 0, 0, 0, 1, 2, Nf_u, Nm, P, T, Tb_u, Nb_u, Nb_u, basis_u, basis_u, Tb_u, 0, 0);
    O = sparse(Nf_p,1);
    bN = [bN1 + bN2;bN3 + bN4;O];

    [A_sum,b_sum] = deal(A + AN,b + bN);
    [A_sum,b_sum] = boundary_adjust(A_sum,b_sum,P,T,Pb_u,Tb_u,B_f_edges_u,B_m_edges,Nf_u,Nb_u,Nb_u,func_boundary1,func_boundary2,[],[],[],[],basis_u(1),basis_u(2));% 调用调整器
    [A_sum,b_sum] = fixed_presure(A_sum,b_sum,func_real_p,Nf_u + Nf_u,1,Pb_p);

    uh_new = A_sum \ b_sum;
    test_error = norm(uh_old - uh_new);
    fprintf("Iteration %d: test_error = %.5e\n",ii,test_error);
    if test_error < tol
        fprintf("In iteration %d, we get expected test_error = %.5e, which is lower than tol = %.5e,\n",ii,test_error,tol);
        fprintf("So that we terminate the iteration.\n");
        break;
    end
end