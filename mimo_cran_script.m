clearvars
close all
clc
path_loss = [2.5 3.5];
N_u = 4;
N_r = 4;
K = [4 10];
L = [6 15];
scaling = {'linear', 'exp', 'fixedk'};
Csum = [10 20 30 40];
for ii = path_loss
    for kk = 1:length(scaling)
        mimo_uplink_stochgeom_sim_wf(ii, N_u, N_r, scaling{kk});
        mimo_downlink_stochgeom_sim_wf(ii, N_u, N_r, scaling{kk});
    end
    for jj = 1:length(K)
        for kk = Csum
            stochgeom_sim_lin_fixedfronthaul_wf(kk, ii, K(jj), L(jj));
        end
    end
end