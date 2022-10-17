format long
clear;
close all;

%%% parameters of frequency hopping signal  
%%% only support this input format 
mod_para = struct("mem0", struct("mod", "msk", "symbol_rate", 5e6), ...
                  "mem1", struct("mod", "msk", "symbol_rate", 5e6), ...
                  "mem2", struct("mod", "msk", "symbol_rate", 5e6));

fs = 610e6;                   %%% sample rate
hop_period = 76923;           %%% period of frequency hopping signal (hop/s)
hop_length = round(1 / hop_period * fs);   %%% time -> samples
hop_num = 14;                              
mem_num = size(fieldnames(mod_para), 1);   %%%  get number of fh signal
net_interval = 30;    %%% minimum frequency between two adjacent signal (in MHz) 

%%%  return link16 class "l" according to the input parameters
l = link16(mem_num, hop_num, net_interval, fs);
freq_pattern = l.freq_pattern;    

%%%  real doa
doa_pattern = repmat((1:1:size(fieldnames(mod_para), 1))' .* 10, 1, hop_num);


%%%  return fh class "fh_ss" according to the input parameters
%%%  fh_ss contains source frequency hopping signal
fh_ss = fh(fs, mem_num, hop_num, hop_length, net_interval, freq_pattern, doa_pattern, mod_para);
fh_ss.src_signal(:, 4*hop_length:5 * hop_length) = 0;    %%%  safe interval
fh_ss.src_signal(:, 9*hop_length:10 * hop_length) = 0;   %%%  safe interval

sig_src = [];   %%% save source frequency hopping signal  (vecotr)

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
title("\fontsize{13}source signal wave")

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


ant_num = 2;
snr = 12e0;
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

hop_vec = tf.hop_vec;
hop_vec_diff =  diff(hop_vec);
hop_vec_mod = hop_vec_diff(find(hop_vec_diff > 1)); 

delata = 1000;
jmp_label = [];
for i = 1:1:length(hop_vec_mod) - 1
    if (abs((hop_vec_mod(i)) - hop_vec_mod(i + 1)) > delata) && (hop_vec_mod(i) > delata)
        jmp_label = [jmp_label, i];
    end
end

interval = diff(jmp_label(2:2:end));
ll = jmp_label(2:2:end);
if length(jmp_label(2:2:end)) == 2
    hop_length_est = mean(hop_vec_mod(ll(1) + 1:ll(2) - 1));
else 
    throw('fail to estimating');
end

time_slot_real= (hop_length * 4) / fs;
time_slot_est= (hop_length_est * interval) / fs;

fprintf("     real time slot length: %e (s)\n", time_slot_real);
fprintf("estimated time slot length: %e (s)\n", time_slot_est);
