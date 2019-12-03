function [] = mimo_downlink_stochgeom_sim_wf(path_loss, N_u, N_r, scaling)
% lambdaL = 2lambdaK, i.e., the relay density is twice the user density
side_length = 1; % we consider a square of side length 1
lambdaL = 4:15;
switch lower(scaling)
    case 'linear'
        fprintf('Lambda_K is half of lambda_L.\n\n');
        lambdaK = 0.5*lambdaL;
        filename_identifier_string = 'linear';
    case 'exp'
        fprintf('Lambda_K is the square root of lambda_L.\n\n');
        lambdaK = sqrt(lambdaL);
        filename_identifier_string = 'exp';
    case 'fixedk'
        fprintf('Lambda_K is kept fixed at 10.\n\n');
        lambdaK = 10*ones(size(lambdaL));
        filename_identifier_string = 'Kfixed';
    otherwise
        error('Unknown method.')
end
search_set = [.1:.1:5 6:(lambdaL(end)*150)];
num_local_trials = 20;
numsims = 100;
for jj = 1:length(lambdaK)
    tic
    Cstartmp = zeros(1, numsims);
    Rsumtmp = Cstartmp;
    Rsumcentraltmp = inf(size(Rsumtmp));
    
    % initializations -- Shouvik
    Cstartmp2 = zeros(1, numsims);
    Rsumtmp2 = Cstartmp2;
    Rsumcentraltmp2 = inf(size(Rsumtmp2));
    sigmavec = zeros(1, numsims);
    sigmavec2 = zeros(1, numsims);
    for ii = 1:numsims
        %% modification: poisson_user_relay_dropping is replaced with ...
        [ G_ori, G_normRef ] = (gen_stochgeom_channel_matrix_MIMO(lambdaK(jj), lambdaL(jj),...
            path_loss, side_length, N_u, N_r));
        % generate channel matrix H by dropping users and relays according
        % to independent Poisson point processes
        H_ori = G_ori.'; H_normRef = G_normRef.';
        numrelays = size(H_ori, 2)/N_r; numusers = size(H_ori, 1)/N_u;
        for kk = 1:num_local_trials
            if numrelays > 1
                A = sort(rand(1, numrelays-1));
                q = [A(1) A(2:end)-A(1:(end-1)) 1-A(end)];
            else
                q = 1;
            end
            qnew = kron(diag(q), eye(N_r));
            [C, ~] = find_broadcast_capacity(H_ori, qnew, 1);
            Rsumcentraltmp(ii) = min(Rsumcentraltmp(ii), C);
            [C, ~] = find_broadcast_capacity(H_normRef, qnew, 1);
            Rsumcentraltmp2(ii) = min(Rsumcentraltmp2(ii), C);
        end
        % calculate centralized MIMO sum rate
        % Modified code to change choice of sigma squared -- Shouvik
        Cstartmp_vec = arrayfun(@(x) find_downlink_DDF_sumrate(H_ori, ...
            N_r, N_u, search_set(x)), 1:length(search_set));
        Cstartmp2_vec = arrayfun(@(x) find_downlink_DDF_sumrate(H_normRef, ...
            N_r, N_u, search_set(x)), 1:length(search_set));
        Rsumtmp_vec = Cstartmp_vec - 0.5*size(H_ori, 1)*log2(1 + (1./search_set));
        Rsumtmp2_vec = Cstartmp2_vec - 0.5*size(H_normRef, 1)*log2(1 + (1./search_set));
        
        [~, idx1] = max(Rsumtmp_vec);
        [~, idx2] = max(Rsumtmp2_vec);
        Cstartmp(ii) = Cstartmp_vec(idx1);
        % calculate C* according to Theorem 2; choose sigma^2 to maximize
        % R_sum^DDF
        Rsumtmp(ii) = Rsumtmp_vec(idx1);
        % calculate sum rate achievable by DDF
        sigmavec(ii) = search_set(idx1);
        
        %        sigmasquared = search_set(idx1);
        Cstartmp2(ii) = Cstartmp2_vec(idx2);
        % calculate C* according to Theorem 2; choose sigma^2 to maximize
        % R_sum^DDF
        Rsumtmp2(ii) = Rsumtmp2_vec(idx2);
        % calculate sum rate achievable by DDF
        sigmavec2(ii) = search_set(idx2);
    end
    % compute the median values over all simulation runs
    Cstar(jj) = median(Cstartmp);
    Cstar2(jj) = median(Cstartmp2);
    Rsum(jj) = median(Rsumtmp);
    Rsum2(jj) = median(Rsumtmp2);
    Rsumcentral(jj) = median(Rsumcentraltmp);
    Rsumcentral2(jj) = median(Rsumcentraltmp2);
    sigma_opt(jj) = median(sigmavec);
    sigma_opt2(jj) = median(sigmavec2);
    toc
end

figure % Original Figure by UCSD
plot(lambdaL, Rsumcentral, 'r-o', 'LineWidth', 2); hold on;
plot(lambdaL, Rsum, 'b-^', 'LineWidth', 2);
plot(lambdaL, Cstar, 'm-s', 'LineWidth', 2);

grid on;
xl = xlabel('\lambda_{r}');
if path_loss == 2.5
    yl = ylabel('Median rates');
    set(yl, 'FontSize', 35);
end
%s1 = sprintf('Centralized MIMO\nsum rate');
% s1 = sprintf('R_{sum}^{MIMO}');
% s2 = sprintf('R_{sum}^{NCF}');
% lg = legend({'R_{sum}^{MIMO}', 'R_{sum}^{NCF}', 'C^*', '\lambda_{u}log\lambda_{u}/2'});
lg = legend({'R_{sum}^{un--cl}', 'R_{sum}^{DDF}', 'C^*'});
set(lg, 'Location', 'northwest');
%eval(['tl = title(''\beta = ' num2str(path_loss) ', \lambda_{r} = 2\lambda_{u}'');'])
set(xl, 'FontSize', 35); set(lg, 'FontSize', 17); set(gca,'FontSize',20);
xlim([lambdaL(1) lambdaL(end)]);

eval(['print -depsc beta' num2str(10*path_loss) '_NuNr' num2str(N_u) '_' num2str(N_r) '_L' ...
    filename_identifier_string '_downlink.eps'])
eval(['print -djpeg beta' num2str(10*path_loss) '_NuNr' num2str(N_u) '_' num2str(N_r) '_L' ...
    filename_identifier_string '_downlink.jpg'])
eval(['savefig(''beta' num2str(10*path_loss) '_NuNr' num2str(N_u) '_' num2str(N_r) '_L' ...
    filename_identifier_string '_downlink.fig'');'])
hold off

if path_loss == 2.5
    figure % Figure with some randomness (Blockage, Shadowing, Fast Fading)
    plot(lambdaL, Rsumcentral2, 'r-o', 'LineWidth', 3); hold on;
    plot(lambdaL, Rsum2, 'b-^', 'LineWidth', 2);
    plot(lambdaL, Cstar2, 'm-s', 'LineWidth', 2);
    
    grid on;
    xl = xlabel('\lambda_{r}');
    
    lg = legend({'R_{sum}^{un--cl}', 'R_{sum}^{DDF}', 'C^*'});
    set(lg, 'Location', 'northwest');
    %eval(['tl = title(''\beta = ' num2str(path_loss) ', \lambda_{r} = 2\lambda_{u}'');'])
    set(xl, 'FontSize', 35); set(lg, 'FontSize', 17); set(gca,'FontSize',20);
    xlim([lambdaL(1) lambdaL(end)]);
    
    eval(['print -depsc SG_NuNr' num2str(N_u) '_' num2str(N_r) '_L' ...
        filename_identifier_string '_downlink.eps'])
    eval(['print -djpeg SG_NuNr' num2str(N_u) '_' num2str(N_r) '_L' ...
        filename_identifier_string '_downlink.jpg'])
    eval(['savefig(''SG_NuNr' num2str(N_u) '_' num2str(N_r) '_L' ...
        filename_identifier_string '_downlink.fig'');'])
    hold off
end
matfilename1 = sprintf('beta%dNuNr%d_%d_L%s_downlink.mat', 10*path_loss, N_u, ...
    N_r, filename_identifier_string);
save(matfilename1, 'N_u' , 'N_r', 'lambdaL', 'lambdaK', 'Cstar', 'Cstar2', ...
    'Rsum', 'Rsum2', 'Rsumcentral', 'Rsumcentral2', 'sigma_opt', 'sigma_opt2');
close all
