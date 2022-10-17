format long
clear;
close all;

%%% parameters of frequency hopping signal  
%%% only support this input format 
mod_para = struct("mem0", struct("mod", "msk", "symbol_rate", 5e6), ...
                  "mem1", struct("mod", "msk", "symbol_rate", 5e6), ...
                  "mem2", struct("mod", "msk", "symbol_rate", 5e6));

fs = 610e6;              %%% sample rate
hop_period = 76923;      %%% period of frequency hopping signal (hop/s)
hop_length = round(1 / hop_period * fs);  %%% time -> samples
hop_num = 7;             %%% number of hop 
mem_num = size(fieldnames(mod_para), 1);     %%%  get number of fh signal
net_interval = 30;       %%% minimum frequency between two adjacent signal (in MHz) 

%%%  return link16 class "l" according to the input parameters
l = link16(mem_num, hop_num, net_interval, fs);
freq_pattern = l.freq_pattern;   %%% real frequency pattern 

%%%  real doa
union_doa = 10;
doa_pattern = repmat((1:1:size(fieldnames(mod_para), 1))' .* union_doa, 1, hop_num);

%%%  return fh class "fh_ss" according to the input parameters
%%%  fh_ss contains source frequency hopping signal
fh_ss = fh(fs, mem_num, hop_num, hop_length, net_interval, freq_pattern, doa_pattern, mod_para);

sig_src = [];  %%% save source frequency hopping signal  (vecotr)

if mem_num == 1
    sig_src = f_ss.src_signal;
else
    sig_src = fh_ss.src_signal(1, :);
    for i = 2:1:mem_num
        sig_src = sig_src + fh_ss.src_signal(i, :);
    end
end

%%% draw source signal  (time domain)
subplot(2, 1, 1)
plot(1:1:length(sig_src), sig_src);
axis([0 length(sig_src) + 1000 -4 4]);
xlabel("samples");
ylabel("Amplitude");
title("\fontsize{13}source signal wave");

%%% stft parameters  
win_length = 256;
dft_length = win_length * 2;
win = hann(win_length);
overlap_length = round(0.75 * win_length);

sig_src_tfspec = stft(sig_src, fs, 'FFTLength', dft_length, ...
'Window', win, 'Centered', false, 'OverlapLength', overlap_length);
%%% draw source signal   (time-frequency domain)
subplot(2, 1, 2)
contour(abs(sig_src_tfspec))
xlabel("column length");
ylabel("DFT length");
title("\fontsize{13}source signal time-frequnecy spectrum");


ant_num = 2;     %%% number of receive antenna
snr = 12e0;      

%%%  return rx_signal class "rx"
%%%  rx contains receive signal
rx = rx_signal(ant_num, 0.1, snr, fh_ss);

th = 0.3;
%%%  return tfdec class "tf"
%%% tf contains estimation
tf = tfdec(rx, win, overlap_length, dft_length, fs, th, 0);

%%% draw receive signal  (time domain)
figure;
subplot(2, 1, 1)
plot(1:1:length(rx.receive_signal(1, :)), real(rx.receive_signal(1, :)));
axis([0 length(rx.receive_signal(1, :)) + 1000 -4 4]);
xlabel("samples");
ylabel("Amplitude");
title("\fontsize{13}received signal wave")
%%% draw receive signal  (time-frequency domain)
subplot(2, 1, 2)
contour(abs(tf.stft_tensor(:, :, 1)))
xlabel("column length");
ylabel("DFT length");
title("\fontsize{13}received signal time-frequnecy spectrum");

%%%  get result form tf 
doa_est = tf.doa_est_;        %%% estimation of doa (matrix)
freq_est = tf.freq_est;       %%% estimation of frequency pattern  (matirx)
hop_num_est = round(length(tf.hop_vec) / 2);
hop_vec = tf.hop_vec;
hop_vec_diff = diff(hop_vec);
hop_vec_mod = hop_vec_diff(find(hop_vec_diff > 1));
hoptime_est = mean(hop_vec_mod(2:end - 1)) / fs;     %%% estimation of hopping time
num_est = tf.num_est;         %%% estimation of number of signal


doa_est_vec = [];           
for i = 1:1:size(doa_est, 1)
    doa_est_vec = [doa_est_vec, doa_est(i, :)];
end 
doa_est_vec = sort(doa_est_vec);
doa_delta = union_doa / 2;
label = find(diff(doa_est_vec) > doa_delta);
seg = [1, label, length(doa_est_vec)];
doa_est_vec_ = [];            %%% estimation of doa (vector)
for i = 1:1:length(seg) - 1
    doa_est_vec_(i) = mean(doa_est_vec(seg(i):seg(i+1))); 
end
fh_node = zeros(num_est, hop_num_est);    %%% sorting frequency hopping signal
for i = 1:1:length(doa_est_vec_)
    [r, c] = find(abs(doa_est - doa_est_vec_(i)) < doa_delta);
    for j = 1:1:length(r)
        fh_node(i, j) = freq_est(r(j), c(j));    %%% row: node   column: frequency
    end
end

%%% draw frequency pattern
figure;
if size(freq_pattern) == size(freq_est)
    for i = 1:1:size(freq_pattern, 1)
        scatter(1:1:size(freq_pattern, 2), freq_pattern(i, :), 'b');
        hold on;
        scatter(1:1:size(freq_est, 2), freq_est(i, :), 'r*');
        hold on;
    end
end
legend('real frequency', 'estimated frequency')
xlabel("hop number");
ylabel("frequency(MHz)");
title("\fontsize{13}frequency pattern");

figure;

%%% draw sorted result
for i = 1:1:size(fh_node, 1)
    subplot(length(doa_est_vec_), 1, i);
    scatter(1:1:size(freq_est, 2), fh_node(i, :), 'r*');
    axis([1 hop_num_est 969 - 10 1206 + 10]);
    xlabel("hop number");
    ylabel("frequency (MHz)");
    if i == 1
        title({"\fontsize{13}hopping frequency sorting";...
                "\fontsize{10}node:" + string(i) + "    doa:" + string(doa_est_vec_(i)) + " °"});
    else
        title("\fontsize{10}node:" + string(i) + "    doa:" + string(doa_est_vec_(i)) + " °");

    end
end

%%%  print result
fprintf("     real hopping time: %e (s)   \n", hop_length / fs);
fprintf("estimated hopping time: %e (s)   \n", hoptime_est);
fprintf("\n");
fprintf("     real hopping period: %d (hop/s)\n", hop_period);
fprintf("estimated hopping period: %d (hop/s)\n", round(1 / hoptime_est));
fprintf("\n");
fprintf("     real number of fh nodes: %d  \n", mem_num);
fprintf("estimated number of fh nodes: %d  \n", num_est);
fprintf("\n");
fprintf("     real doa: %f°\n", doa_pattern(:,1)');
fprintf("estimated doa: %f°\n", doa_est_vec_);

% close all
