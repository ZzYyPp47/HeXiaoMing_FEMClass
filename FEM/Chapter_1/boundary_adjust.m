% 1D boundary adjuster
function [A,b] = boundary_adjust(A, b, Pb, B_info, func_boundary, func_c)
for ii = 1:1:2
    if B_info(1,ii) == -1
        jj = B_info(2,ii);
        A(jj,:) = 0;
        A(jj,jj) = 1;
        b(jj) = func_boundary(Pb(jj),[]);
    end
    if B_info(1,ii) == -2
        jj = B_info(2,ii);
        b(jj) = b(jj) + B_info(3,ii) * B_info(5,ii) * func_c(Pb(jj),[]);
    end
    if B_info(1,ii) == -3
        jj = B_info(2,ii);
        A(jj,jj) = A(jj,jj) + B_info(3,ii) * B_info(4,ii) * func_c(Pb(jj),[]);
        b(jj) = b(jj) + B_info(3,ii) * B_info(5,ii) * func_c(Pb(jj),[]);
    end
end
end