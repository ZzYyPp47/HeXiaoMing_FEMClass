% make sure apply Dirichlet condition in intersection.
function B_f_edges = check_legal(B_m_edges,B_f_edges,Nx,Ny,Nfx,Nfy)
% Neumann
if (B_m_edges(1,1) == -2 && B_m_edges(1,end) == -1) || (B_m_edges(1,1) == -1 && B_m_edges(1,end) == -2)
    B_f_edges(1,1) = -1;
end
if (B_m_edges(1,Nx) == -2 && B_m_edges(1,Nx + 1) == -1) || (B_m_edges(1,Nx) == -1 && B_m_edges(1,Nx + 1) == -2)
    B_f_edges(1,Nfx) = -1;
end
if (B_m_edges(1,Nx + Ny) == -2 && B_m_edges(1,Nx + Ny + 1) == -1) || (B_m_edges(1,Nx + Ny) == -1 && B_m_edges(1,Nx + Ny + 1) == -2)
    B_f_edges(1,Nfx + Nfy - 1) = -1;
end
if (B_m_edges(1,2 * Nx + Ny) == -2 && B_m_edges(1,2 * Nx + Ny + 1) == -1) || (B_m_edges(1,2 * Nx + Ny) == -1 && B_m_edges(1,2 * Nx + Ny + 1) == -2)
    B_f_edges(1,2 * Nfx + Nfy - 2) = -1;
end
% Robin
if (B_m_edges(1,1) == -3 && B_m_edges(1,end) == -1) || (B_m_edges(1,1) == -1 && B_m_edges(1,end) == -3)
    B_f_edges(1,1) = -1;
end
if (B_m_edges(1,Nx) == -3 && B_m_edges(1,Nx + 1) == -1) || (B_m_edges(1,Nx) == -1 && B_m_edges(1,Nx + 1) == -3)
    B_f_edges(1,Nfx) = -1;
end
if (B_m_edges(1,Nx + Ny) == -3 && B_m_edges(1,Nx + Ny + 1) == -1) || (B_m_edges(1,Nx + Ny) == -1 && B_m_edges(1,Nx + Ny + 1) == -3)
    B_f_edges(1,Nfx + Nfy - 1) = -1;
end
if (B_m_edges(1,2 * Nx + Ny) == -3 && B_m_edges(1,2 * Nx + Ny + 1) == -1) || (B_m_edges(1,2 * Nx + Ny) == -1 && B_m_edges(1,2 * Nx + Ny + 1) == -3)
    B_f_edges(1,2 * Nfx + Nfy - 2) = -1;
end
end
