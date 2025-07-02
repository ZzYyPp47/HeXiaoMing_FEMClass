% 边界条件调整器
function [A,b] = boundary_adjust(A,b,P,T,Pb_trial,Tb,B_f_edges,B_m_edges,Nf,Nb_test,Nb_trial,func_boundary1,func_boundary2,func_c,func_p,func_q,func_r,dim,type,t0)
% [A,b] = check_robin(A,b,P,T,Tb,B_m_edges,Nb_test,Nb_trial,func_c,func_q,func_r,dim,type);
% b = check_neu(b,P,T,Tb,B_m_edges,Nb_test,func_c,func_p,dim,type);
[A,b] = check_dirich(A,b,Nf,Pb_trial,B_f_edges,func_boundary1,func_boundary2,t0);
end

function [A,b] = check_dirich(A,b,Nf,Pb_trial,B_f_edges,func_boundary1,func_boundary2,t0)
idx = (B_f_edges(1,:) == -1);
row = B_f_edges(2,idx);
% eqa 1
A(row,:) = 0;
diag_A = diag(A);% 提取对角线
diag_A(row) = 1; % 将对应对角元置换为1
A = spdiags(diag_A,0,A);% 替换主对角元
b(row) = func_boundary1(Pb_trial(1,row),Pb_trial(2,row),t0);
% eqa 2
A(row + Nf,:) = 0;
diag_A = diag(A);% 提取对角线
diag_A(row + Nf) = 1; % 将对应对角元置换为1
A = spdiags(diag_A,0,A);% 替换主对角元
b(row + Nf) = func_boundary2(Pb_trial(1,row),Pb_trial(2,row),t0);
end

function b = check_neu(b,P,T,Tb,B_m_edges,Nb_test,func_c,func_p,dim,type)
for ii = 1:1:size(B_m_edges,2)
    if B_m_edges(1,ii) == -2
        jj = B_m_edges(2,ii);
        for beta = 1:1:Nb_test
            p1 = P(:,B_m_edges(3,ii));
            p2 = P(:,B_m_edges(4,ii));
            int_func = @(x,y) func_c(x,y) .* func_p(x,y) .* basis_local_function(x, y, dim, type, beta, 0, 0, P(:,T(:,jj)));
            if p1(1) == p2(1)%  x = c0型
                result = gauss_2d_line(int_func, p1(2), p2(2), p1, p2,'v');
            else
                result = gauss_2d_line(int_func, p1(1), p2(1), p1, p2,'other');
            end
            b(Tb(beta,jj)) = b(Tb(beta,jj)) + result;
        end
    end
end
end

function [A,b] = check_robin(A,b,P,T,Tb,B_m_edges,Nb_test,Nb_trial,func_c,func_q,func_r,dim,type)
for ii = 1:1:size(B_m_edges,2)
    if B_m_edges(1,ii) == -3
        jj = B_m_edges(2,ii);
        for beta = 1:1:Nb_test
            p1 = P(:,B_m_edges(3,ii));
            p2 = P(:,B_m_edges(4,ii));
            int_func = @(x,y) func_c(x,y) .* func_q(x,y) .* basis_local_function(x, y, dim, type, beta, 0, 0, P(:,T(:,jj)));
            if p1(1) == p2(1)%  x = c0型
                result = gauss_2d_line(int_func, p1(2), p2(2), p1, p2,'v');
            else
                result = gauss_2d_line(int_func, p1(1), p2(1), p1, p2,'other');
            end
            b(Tb(beta,jj)) = b(Tb(beta,jj)) + result;
            for alpha = 1:1:Nb_trial
                int_func = @(x,y) func_c(x,y) .* func_r(x,y)...
                    .* basis_local_function(x, y, dim, type, alpha, 0, 0, P(:,T(:,jj)))...
                    .* basis_local_function(x, y, dim, type, beta, 0, 0, P(:,T(:,jj)));
                if p1(1) == p2(1)%  x = c0型
                    result = gauss_2d_line(int_func, p1(2), p2(2), p1, p2,'v');
                else
                    result = gauss_2d_line(int_func, p1(1), p2(1), p1, p2,'other');
                end
                A(Tb(beta,jj),Tb(alpha,jj)) = A(Tb(beta,jj),Tb(alpha,jj)) + result;
            end
        end
    end
end
end
