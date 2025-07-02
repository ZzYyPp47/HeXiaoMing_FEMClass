%% THIS PROCEDURE IS ONLY USED TO CALCULATE NON-TIME-DEPENDENT A.
% If option is set to "all", the solution will be recorded over the entire time domain.
% Otherwise, only the solution at the final time will be recorded.
% Output: U is the solution determined by option, uh_end is the solution at the final time (column vector).
% Input: u0 = [u1;u2;p] at t = 0.
function [U,uh_end] = t_fdm(M, A, b_all, theta, Time, u0, option, P, T, Pb_u, Pb_p, Tb_u, Nf_u, B_f_edges_u, B_m_edges, Nb_u, func_boundary1, func_boundary2, func_real_p, func_c, func_p, func_q, func_r, dim, type)
%% Modify to meet your needs
    if option == "all"
        U = zeros([size(u0,1),size(Time,2)]); % 全局答案存储
    else
        U = zeros(size(u0)); % 只存储最终时间的解
    end
    U(:,1) = u0;
    ht = Time(2) - Time(1);
    M_ht = M ./ ht;
    theta_A = theta * A;
    one_minus_theta_A = (1 - theta) * A;    
    left_hand = M_ht + theta_A;% 差分方程左手侧
    for iteration = 1:length(Time) - 1
        if iteration == 1
            u_current = U(:,1);
        else
            if option == "all"
                u_current = U(:,iteration);
            else
                u_current = U;
            end
        end       
        right_hand = (M_ht - one_minus_theta_A) * u_current + ...
                     theta * b_all{iteration + 1} + ...
                     (1 - theta) * b_all{iteration};% 差分方程右手侧
        [left_hand, right_hand] = boundary_adjust(left_hand, right_hand, P, T, Pb_u, Tb_u, B_f_edges_u, B_m_edges, Nf_u, Nb_u, Nb_u, func_boundary1, func_boundary2, func_c, func_p, func_q, func_r, dim, type, Time(iteration + 1));% 调用边界调整器
        [left_hand, right_hand] = fixed_presure(left_hand,right_hand,func_real_p,Nf_u + Nf_u,1,Pb_p,Time(iteration + 1));
        u_next = left_hand \ right_hand;
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