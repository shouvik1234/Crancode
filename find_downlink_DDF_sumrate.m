function [Rsum, Sigma_opt] = find_downlink_DDF_sumrate(H, N_r, N_u, sigmasquared)
numusers = size(H, 1)/N_u;
numrelays = size(H, 2)/N_r;
current_interference_matrix = zeros(size(H, 1));
Sigma_opt = cell(1, numrelays);
numiter = 3;
for kk = 1:numrelays
    Sigma_opt{kk} = zeros(N_r);
end
for jj = numiter
    optimal_matrix = zeros(size(H, 1));
    for kk = 1:numrelays
        Htmp = H(:, ((kk-1)*N_r + 1):(kk*N_r));
        Q = current_interference_matrix + (sigmasquared)*eye(size(H, 1));
        [~, Sigma_opt{kk}] = find_broadcast_capacity(Htmp.', Q, 1);
        optimal_matrix = optimal_matrix + (Htmp*Sigma_opt{kk}*Htmp.');
        Htmpnext = H(:, ((mod(kk, numrelays))*N_r + 1):((mod(kk, numrelays) + 1)*N_r));
        current_interference_matrix = current_interference_matrix + ...
            (Htmp*Sigma_opt{kk}*Htmp.') - (Htmpnext*Sigma_opt{mod(kk, numrelays) + 1}*Htmpnext.');
    end
end
Rsum = 0.5*log2(det((optimal_matrix/sigmasquared) + eye(size(H, 1))));
end