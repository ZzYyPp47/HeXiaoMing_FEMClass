%% THIS PROCEDURE IS ONLY USED TO CALCULATE NON-TIME-DEPENDENT A.
% If option is set to "all", the solution will be recorded over the entire time domain.
% Otherwise, only the solution at the final time will be recorded.
% Output: U is the solution determined by option, uh_end is the solution at the final time (column vector).
% Input: u0 = [u1;u2;p] at t = 0.
function [U,uh_end] = t_fdm(iteration_num, tol, M, A, b_all, Time, u0, option, P, T, Pb_u, Pb_p, Tb_u, Nf_u, Nf_p, Nm, B_f_edges_u, B_m_edges, Nb_u, func_boundary1, func_boundary2, func_real_p, func_c, func_p, func_q, func_r, basis_u)
%% Modify to meet your needs
if option == "all"
    U = zeros([size(u0,1),size(Time,2)]); % 全局答案存储
else
    U = zeros(size(u0)); % 只存储最终时间的解
end
U(:,1) = u0;
ht = Time(2) - Time(1);

for iteration = 1:length(Time) - 1
    fprintf("Now we calculate time = %.5f,the initial time is %.5f,the end time is %.5f.\n",Time(iteration + 1),Time(1),Time(end))
    if iteration == 1
        u_current = U(:,1);
    else
        if option == "all"
            u_current = U(:,iteration);
        else
            u_current = U;
        end
    end

    u_next = Newton_iteration(iteration_num,tol,ht,A,M,b_all{iteration + 1},u_current,B_m_edges,B_f_edges_u,Nf_u,Nf_p,Nm,P,T,Pb_u,Pb_p,Nb_u,basis_u,Tb_u,func_boundary1,func_boundary2,func_real_p,Time(iteration + 1));

    if option == "all"
        U(:,iteration + 1) = u_next;
    else
        U(:) = u_next;
    end
end
if option == "all"
    uh_end = U(:,end);
else
    uh_end = U;
end
end