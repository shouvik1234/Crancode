function [Rsum, Sigma_opt] = find_uplink_NCF_sumrate(G, N_r, N_u, sigmasquared)
numusers = size(G, 2)/N_u;
numrelays = size(G, 1)/N_r;
current_interference_matrix = zeros(size(G, 1));
Sigma_opt = cell(1, numusers);
numiter = 3;
for kk = 1:numusers
    Sigma_opt{kk} = zeros(N_u);
end
for jj = numiter
    optimal_matrix = zeros(size(G, 1));
    for kk = 1:numusers
        Gtmp = G(:, ((kk-1)*N_u + 1):(kk*N_u));
        Q = current_interference_matrix + (sigmasquared + 1)*eye(size(G, 1));
        [~, Sigma_opt{kk}] = find_broadcast_capacity(Gtmp.', Q, 1);
        optimal_matrix = optimal_matrix + (Gtmp*Sigma_opt{kk}*Gtmp.');
        Gtmpnext = G(:, ((mod(kk, numusers))*N_u + 1):((mod(kk, numusers) + 1)*N_u));
        current_interference_matrix = current_interference_matrix + ...
            (Gtmp*Sigma_opt{kk}*Gtmp.') - (Gtmpnext*Sigma_opt{mod(kk, numusers) + 1}*Gtmpnext.');
    end
end
Rsum = 0.5*log2(det((optimal_matrix/(sigmasquared + 1)) + eye(size(G, 1))));
end