function noise_cutoff = calculate_iterative_noise(t_index)
% Author: Sethuraman Sankaran, as described in Rajasethupathy et al Nature
% 2015.

cutoff_i = 0.3;
noise_std = std(t_index(find(t_index < cutoff_i)));
threshold = 0.02;
n_iter = 1;
while ( (abs(3 * noise_std - cutoff_i) > threshold) && (n_iter <= 200))
    if ( (3 * noise_std) > cutoff_i)
        cutoff_i = 1.1 * cutoff_i;
    else
        cutoff_i = 0.9 * cutoff_i;
    end
    noise_std = std(t_index(find(t_index < cutoff_i)));
    n_iter = n_iter + 1;
end
noise_cutoff = cutoff_i;
    