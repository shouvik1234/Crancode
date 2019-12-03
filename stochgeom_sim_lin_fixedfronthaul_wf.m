function [] = stochgeom_sim_lin_fixedfronthaul_wf(Csum, path_loss, numusers, numrelays)
% lambdaL = 2lambdaK, i.e., the relay density is twice the user density
side_length = 1; % we consider a square of side length 1
N_u = 2:2:10;
N_r = N_u;
search_set = [.1:.1:5 6:(numrelays*150)];
numsims = 100;
for jj = 1:length(N_u)
    tic
    maxsumratetmp_uplink = zeros(1, numsims);
    maxsumratetmp2_uplink = maxsumratetmp_uplink;
    maxsumratetmp_downlink = zeros(1, numsims);
    maxsumratetmp2_downlink = maxsumratetmp_uplink;
    sigmavec_uplink = zeros(1, numsims);
    sigmavec2_uplink = zeros(1, numsims);
    sigmavec_downlink = zeros(1, numsims);
    sigmavec2_downlink = zeros(1, numsims);
    for ii = 1:numsims
        %% modification: poisson_user_relay_dropping is replaced with ...
        [G_ori, G_normRef] = gen_stochgeom_channel_matrix_fixedKL(path_loss, ...
            side_length, N_u(jj), N_r(jj), numusers, numrelays);
        [G_ori1, G_normRef1] = gen_stochgeom_channel_matrix_fixedKL(path_loss, ...
            side_length, N_u(jj), N_r(jj), numusers, numrelays);
        H_ori = G_ori1.'; H_normRef = G_normRef1.';
        % calculate centralized MIMO sum rate
        % Modified code to change choice of sigma squared -- Shouvik
        Rsumtmp_up_vec = arrayfun(@(x) find_uplink_NCF_sumrate(G_ori, ...
            N_r(jj), N_u(jj), search_set(x)), 1:length(search_set));
        Rsumtmp_down_vec = arrayfun(@(x) find_downlink_DDF_sumrate(H_ori, ...
            N_r(jj), N_u(jj), search_set(x)), 1:length(search_set));
        Rsumtmp2_up_vec = arrayfun(@(x) find_uplink_NCF_sumrate(G_normRef, ...
            N_r(jj), N_u(jj), search_set(x)), 1:length(search_set));
        Rsumtmp2_down_vec = arrayfun(@(x) find_downlink_DDF_sumrate(H_normRef, ...
            N_r(jj), N_u(jj), search_set(x)), 1:length(search_set));
        
        [maxsumratetmp_uplink(ii), sigmavec_uplink(ii)] = max(min(Csum - ...
            (0.5*N_r(jj)*numrelays*log2(1+(1./search_set))), Rsumtmp_up_vec));
        [maxsumratetmp2_uplink(ii), sigmavec2_uplink(ii)]= max(min(Csum - ...
            (0.5*N_r(jj)*numrelays*log2(1+(1./search_set))), Rsumtmp2_up_vec));
        
        [maxsumratetmp_downlink(ii), sigmavec_downlink(ii)] = max(min(Csum, ...
            Rsumtmp_down_vec) - (0.5*N_u(jj)*numusers*log2(1+(1./search_set))));
        [maxsumratetmp2_downlink(ii), sigmavec2_downlink(ii)] = max(min(Csum, ...
            Rsumtmp2_down_vec) - (0.5*N_u(jj)*numusers*log2(1+(1./search_set))));
    end
    % compute the median values over all simulation runs
    Rsumuplink(jj) = median(maxsumratetmp_uplink);
    Rsumuplink2(jj) = median(maxsumratetmp2_uplink);
    Rsumdownlink(jj) = median(maxsumratetmp_downlink);
    Rsumdownlink2(jj) = median(maxsumratetmp2_downlink);
    sigmaopt_uplink(jj) = median(search_set(sigmavec_uplink));
    sigmaopt2_uplink(jj) = median(search_set(sigmavec2_uplink));
    sigmaopt_downlink(jj) = median(search_set(sigmavec_downlink));
    sigmaopt2_downlink(jj) = median(search_set(sigmavec2_downlink));
    toc
end

figure % Original Figure by UCSD
plot(N_u, Rsumuplink, 'r-o', 'LineWidth', 2); hold on;
plot(N_u, Csum*ones(size(N_u)), 'k-', 'LineWidth', 3);

grid on;
xl = xlabel('N_u = N_r');
if path_loss == 2.5
    yl = ylabel('Median rates');
    set(yl, 'FontSize', 35);
end
%s1 = sprintf('Centralized MIMO\nsum rate');
% s1 = sprintf('R_{sum}^{MIMO}');
% s2 = sprintf('R_{sum}^{NCF}');
% lg = legend({'R_{sum}^{MIMO}', 'R_{sum}^{NCF}', 'C^*', '\lambda_{u}log\lambda_{u}/2'});
lg = legend({'R_{sum}^{max}', 'C_{\Sigma}'});
set(lg, 'Location', 'northwest');
set(xl, 'FontSize', 35); set(lg, 'FontSize', 17); set(gca,'FontSize',20);
xlim([N_u(1) N_u(end)]);

eval(['print -depsc beta' num2str(10*path_loss) '_KL' num2str(numusers) '_' ...
    num2str(numrelays) '_Csum_' num2str(Csum) '_linear_uplink_fixedfronthaul.eps'])
eval(['print -djpeg beta' num2str(10*path_loss) '_KL' num2str(numusers) '_' ...
    num2str(numrelays) '_Csum_' num2str(Csum) '_linear_uplink_fixedfronthaul.jpg'])
eval(['savefig(''beta' num2str(10*path_loss) '_KL' num2str(numusers) '_' ...
    num2str(numrelays) '_Csum_' num2str(Csum) '_linear_uplink_fixedfronthaul.fig'');'])
hold off

figure % Original Figure by UCSD
plot(N_u, Rsumdownlink, 'r-o', 'LineWidth', 2); hold on;
plot(N_u, Csum*ones(size(N_u)), 'k-', 'LineWidth', 3);

grid on;
xl = xlabel('N_u = N_r');
if path_loss == 2.5
    yl = ylabel('Median rates');
    set(yl, 'FontSize', 35);
end
%s1 = sprintf('Centralized MIMO\nsum rate');
% s1 = sprintf('R_{sum}^{MIMO}');
% s2 = sprintf('R_{sum}^{NCF}');
% lg = legend({'R_{sum}^{MIMO}', 'R_{sum}^{NCF}', 'C^*', '\lambda_{u}log\lambda_{u}/2'});
lg = legend({'R_{sum}^{max}', 'C_{\Sigma}'});
set(lg, 'Location', 'northwest');
set(xl, 'FontSize', 35); set(lg, 'FontSize', 17); set(gca,'FontSize',20);
xlim([N_u(1) N_u(end)]);

eval(['print -depsc beta' num2str(10*path_loss) '_KL' num2str(numusers) '_' ...
    num2str(numrelays) '_Csum_' num2str(Csum) '_linear_downlink_fixedfronthaul.eps'])
eval(['print -djpeg beta' num2str(10*path_loss) '_KL' num2str(numusers) '_' ...
    num2str(numrelays) '_Csum_' num2str(Csum) '_linear_downlink_fixedfronthaul.jpg'])
eval(['savefig(''beta' num2str(10*path_loss) '_KL' num2str(numusers) '_' ...
    num2str(numrelays) '_Csum_' num2str(Csum) '_linear_downlink_fixedfronthaul.fig'');'])
hold off

if path_loss == 2.5
    figure % Figure with some randomness (Blockage, Shadowing, Fast Fading)
    plot(N_u, Rsumuplink2, 'r-o', 'LineWidth', 2); hold on;
    plot(N_u, Csum*ones(size(N_u)), 'k-', 'LineWidth', 3);
    
    grid on;
    xl = xlabel('N_u = N_r');
    %s1 = sprintf('Centralized MIMO\nsum rate');
    % s1 = sprintf('R_{sum}^{MIMO}');
    % s2 = sprintf('R_{sum}^{NCF}');
    % lg = legend({'R_{sum}^{MIMO}', 'R_{sum}^{NCF}', 'C^*', '\lambda_{u}log\lambda_{u}/2'});
    lg = legend({'R_{sum}^{max}', 'C_{\Sigma}'});
    set(lg, 'Location', 'northwest');
    set(xl, 'FontSize', 35); set(lg, 'FontSize', 17); set(gca,'FontSize',20);
    xlim([N_u(1) N_u(end)]);
    eval(['print -depsc SG_KL' num2str(numusers) '_' ...
        num2str(numrelays) '_Csum_' num2str(Csum) '_linear_uplink_fixedfronthaul.eps'])
    eval(['print -djpeg SG_KL' num2str(numusers) '_' ...
        num2str(numrelays) '_Csum_' num2str(Csum) '_linear_uplink_fixedfronthaul.jpg'])
    eval(['savefig(''SG_KL' num2str(numusers) '_' ...
        num2str(numrelays) '_Csum_' num2str(Csum) '_linear_uplink_fixedfronthaul.fig'');'])
    hold off
    
    figure % Figure with some randomness (Blockage, Shadowing, Fast Fading)
    plot(N_u, Rsumdownlink2, 'r-o', 'LineWidth', 2); hold on;
    plot(N_u, Csum*ones(size(N_u)), 'k-', 'LineWidth', 3);
    
    grid on;
    xl = xlabel('N_u = N_r');
    %s1 = sprintf('Centralized MIMO\nsum rate');
    % s1 = sprintf('R_{sum}^{MIMO}');
    % s2 = sprintf('R_{sum}^{NCF}');
    % lg = legend({'R_{sum}^{MIMO}', 'R_{sum}^{NCF}', 'C^*', '\lambda_{u}log\lambda_{u}/2'});
    lg = legend({'R_{sum}^{max}', 'C_{\Sigma}'});
    set(lg, 'Location', 'northwest');
    set(xl, 'FontSize', 35); set(lg, 'FontSize', 17); set(gca,'FontSize',20);
    xlim([N_u(1) N_u(end)]);
    eval(['print -depsc SG_KL' num2str(numusers) '_' ...
        num2str(numrelays) '_Csum_' num2str(Csum) '_linear_downlink_fixedfronthaul.eps'])
    eval(['print -djpeg SG_KL' num2str(numusers) '_' ...
        num2str(numrelays) '_Csum_' num2str(Csum) '_linear_downlink_fixedfronthaul.jpg'])
    eval(['savefig(''SG_KL' num2str(numusers) '_' ...
        num2str(numrelays) '_Csum_' num2str(Csum) '_linear_downlink_fixedfronthaul.fig'');'])
    hold off
end

matfilename1 = sprintf('Csum%dbeta%dKL%d_%d_linear_fixedfronthaul.mat', Csum, 10*path_loss, numusers, numrelays);
save(matfilename1, 'Csum', 'numusers', 'numrelays', 'N_u' , 'N_r', 'Rsumuplink', ...
    'Rsumuplink2', 'Rsumdownlink', 'Rsumdownlink2', 'sigmaopt_uplink', 'sigmaopt2_uplink', ...
    'sigmaopt_downlink', 'sigmaopt2_downlink');
close all

end