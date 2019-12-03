function [C, optimal_Sigma] = find_broadcast_capacity(H, Q, P1)
% H is KxL channel matrix
% Q is a LxL symmetric PSD matrix.
% constraint: trace(optimal_Sigma) <= 1 
Gtilde = sqrtm(Q)\(H.');
[~,S,V] = svd(Gtilde, 0);
l = length(S);
noise_levels = 1./(diag(S).^2);
noise_levels = noise_levels(~isinf(noise_levels));
[sorted_levels, idx] = sort(noise_levels);
cumulated_levels = cumsum((1:(length(sorted_levels)-1)).*(sorted_levels(2:end) ...
    - sorted_levels(1:(end-1))));
if sum(cumulated_levels >= P1)
    last_filled_level = find(cumulated_levels >= P1, 1);
    water_level = (P1 + sum(sorted_levels(1:last_filled_level)))/last_filled_level;
else
    water_level = (P1 + sum(sorted_levels))/length(sorted_levels);
end
power_alloc_sorted = max(water_level - sorted_levels, 0);
C = 0.5*sum(log2(1 + (power_alloc_sorted./sorted_levels)));
power_alloc_original = zeros(size(power_alloc_sorted));
for ii = 1:length(power_alloc_sorted)
    power_alloc_original(idx(ii)) = power_alloc_sorted(ii);
end
optimal_Sigma = (V*diag([power_alloc_original.' zeros(1, l-length(power_alloc_original))].'))*V.'; 
end