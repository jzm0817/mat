%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% class : 
%%%        tfdec
%%% fea: 
%%%        using stft to estimate parameter according to receive signal 
%%% parameter: 
%%%   signal_input: receive signal
%%%            win: stft window
%%% overlap_length: stft overlap length between two adjacent window
%%%     dft_length: stft fft length
%%%             fs: sample rate
%%%             th: stft denoise threshold 
%%%            net: 0: parameter estimation   1: network estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef tfdec

    properties
        dft_length;
        win_length;
        overlap_length;
        fs;
        th;
        win;
        
        num_est;    %% save the number of fh signal (estimation)
        freq_est;                      %% save frequency pattern (estimation)
        doa_est;                       %% save doa (estimation, music)
        doa_est_;                      %% save doa (estimation, bss)
        hop_vec;                       %% save signal segment (estimation)
        hop; 

        f_resolution;                  %% resolution according to length of window
        dft_resolution;                %% resolution according to length of dft
        rx;                            %% receive signal

        label = containers.Map();
        label_d = containers.Map();
        index_delta_min;               %% index_delta_min according to band_width_min
        index_delta_max;               %% index_delta_max according to band_width
        band_width = 7.5e6 * 2;
        band_width_min = 2.5e6 * 2;
        stft_tensor = [];              %% result of stft
        seg = containers.Map();
        loss = containers.Map();
        % p_mu = containers.Map();
        % p_mu;
        network;
    end


    methods
        %% constructor
        function obj = tfdec(signal_input, win, overlap_length, dft_length, fs, th, net)
            obj.fs = fs;
            obj.th = th;
            obj.win = win;
            obj.rx = signal_input;
            obj.win_length = length(win);
            obj.overlap_length = overlap_length;
            obj.dft_length = dft_length;
            obj.f_resolution = obj.fs / obj.win_length * 1e-6;      %% resolution according to length of window (MHz)
            obj.dft_resolution = obj.fs / obj.dft_length * 1e-6;    %% resolution according to length of dft (MHz)

            obj.index_delta_min = round(1 * (obj.band_width_min / (obj.fs / obj.dft_length)));  %% max interval
            obj.index_delta_max = round(1 * (obj.band_width / (obj.fs / obj.dft_length)));      %% min interval

            %% using bss method to estimate doa, at least two antenna are needed
            for i = 1:1:2
                %% stft
                [s, f, t] = stft(signal_input.receive_signal(i, :), obj.fs, 'FFTLength', obj.dft_length, 'Window', obj.win, 'Centered', false, 'OverlapLength', obj.overlap_length);
                %% save stft result
                obj.stft_tensor(:, :, i) = s;
                [obj.label(string(i)), obj.label_d(string(i))] = obj.deduct(s);
            end

            if net == 0
                [obj.freq_est, obj.doa_est_, obj.hop_vec, obj.num_est] = obj.get_res();
            else
                obj.network = obj.net_dec();
            end

        end

        %% find peaks according to stft result
        function [s_label, s_label_d] = deduct(obj, tf_spec)
           
            %% find the maximum in column
            s_col_max = max(abs(tf_spec));
              
            %% denoise according to 'obj.th'
            s_ovth = abs(tf_spec) > obj.th .* (s_col_max);
            s_mod = tf_spec .* s_ovth;

            %% save frequency peak result
            s_label = containers.Map();
            s_label_d = containers.Map();
           
            %% traversing the column
            for i = 1:size(s_mod, 2)

                if sum(s_mod(:, i)' * s_mod(:, i) > 100)

                    %% half search
                    s_mod_diff = diff(abs(s_mod(1:round(size(s_mod, 1) / 2), i)));
                    %% find stagnation point and return row label
                    s_mod_label = find(diff((s_mod_diff > 0) - (s_mod_diff < 0 )) == -2) + 1;

                    if (length(s_mod_label) > 1)
                        s_mod_label_cluster = find(diff(s_mod_label) > obj.index_delta_max);
                    else
                        s_mod_label_cluster = s_mod_label;
                    end

                    s_freq_label = [];

                    %% divide stagnation point segment
                    if length(s_mod_label_cluster) >= 1

                        for kk = 1:1:length(s_mod_label_cluster) + 1
                            if kk == 1
                                l = s_mod_label(1:s_mod_label_cluster(kk));
                            elseif kk == length(s_mod_label_cluster) + 1
                                l = s_mod_label(s_mod_label_cluster(kk-1) + 1:length(s_mod_label));
                            else
                                l = s_mod_label(s_mod_label_cluster(kk-1) + 1:s_mod_label_cluster(kk));
                            end
                       
                        ll = [];
                            %% find real stagnation point 
                            if (0.35 * max(abs(tf_spec(l, i))) > min(abs(tf_spec(l, i))))

                                l = l(find(abs(tf_spec(l, i)) > 0.35 * max(abs(tf_spec(l, i)))));
                            end
                            
                            if isempty(l)
                                disp "zero length warning"
                            end
                            
                            %% find two adjacent real stagnation point 
                            %% convert it to carrier stagnation point 
                            if (~isempty(find(diff(l) <= obj.index_delta_min))) && (~mod(length(l), 2) && (length(l) > 2))

                                for ii = 1:round(length(l) / 2)
                                    ll(ii) = round(mean([l(ii), l(ii+2)]));
                                end
                                
                            else

                                if (~mod(length(l), 2))

                                    for ii = 1:round(length(l) / 2)
                                        ll(ii) = round(mean([l(ii), l(ii+1)]));
                                    end
                                else

                                    if length(l) > 2
                                        ll(1) = round(mean([l(1), l(2)]));
                                        ll(2) = round(mean([l(2), l(end)]));
                                    else
                                        ll = l;
                                    end
                                end
                                
                            end
                        %% label of stagnation point 
                        s_freq_label = [s_freq_label; ll'];

                        end
                    else
                        s_freq_label = round(mean(s_mod_label));
                    end

                    s_label("col" + string(i-1)) = s_freq_label;
                    s_label_d("col" + string(i-1)) = s_mod_label;
                else
                    s_label("col" + string(i-1)) = 0;
                    s_label_d("col" + string(i-1)) = 0;

                end
            end
            
        end


        %% parameter estimation
        function [f_est, doa_est_, hop_vec, num] = get_res(obj)
           
            %% convert map container to matrix
            label = obj.label("1");

            %% find maximum length in map container
            mx_len = length(label("col" + string(0)));
            for i = 1:1:label.length - 1
                if mx_len < length(label("col" + string(i)))
                    mx_len = length(label("col" + string(i)));
                end
            end
            
            %% label matirx
            label_m = zeros(mx_len, length(label));
            
            for i = 1:1:label.length
                label_m(:, i) = [label("col" + string(i-1)); zeros(mx_len - length(label("col" + string(i-1))), 1)];
            end

            %% label_vec contain the number of frequency
            label_vec = [];
            for i = 1:1:size(label_m, 2)
                label_vec(i) = length(find(label_m(:, i) > 0));
            end

            %% the number of frequency
            num_est = round(mean(label_vec(find(label_vec > 0))));
            num = num_est;
            %% empty zone in stft matirx
            label_zero = find(label_vec == 0);
            %% maybe frequency change time
            label_greate = find(label_vec > num_est);

            %% find empty zone column
            if ~isempty(label_zero)

                %% preserve interval
                delta_zero = 5;
                label_zero_temp = [];
                for i = 1:1:length(label_zero) - 1
                    if (label_zero(i + 1) - label_zero(i) > delta_zero)
                        label_zero_temp = [label_zero_temp, i];
                    end
                end

                label_zero_seg = [];

                for i = 1:1:length(label_zero_temp)
                    label_zero_seg = [label_zero_seg, label_zero(label_zero_temp(i)), label_zero(label_zero_temp(i) + 1)];
                end
                label_zero_seg = [label_zero(1), label_zero_seg, label_zero(end)];
            end

            %% find frequency change time
            label_jmp = [];

            %% find frequency change time
            if ~isempty(label_greate)

                delta_greate = 10;
                label_delta = 2;
                %% search range
                span = 10;
                
                for i = 1:1:length(label_greate)

                    if (label_greate(i) < delta_greate) || (label_greate(i) > length(label_vec) - delta_greate)
                        continue;
                    else
                        cnt = 1;
                        while(cnt < span)
                            if (label_vec(label_greate(i) - cnt) == label_vec(label_greate(i) + cnt)) && ...
                                (sum(abs(label_m(1:num_est, label_greate(i) - cnt) - label_m(1:num_est, label_greate(i) + cnt))) > num_est * label_delta);
                                
                                label_jmp = [label_jmp, label_greate(i)];
                                cnt = 20;
                            end
                            cnt = cnt + 1;
                        end
                    end

                end
                %% real frequency change time
                label_jmp = unique([label_jmp, size(label_m, 2)]);
            end

            %% real frequency change time
            if ~isempty(label_jmp)

                jmp_delta = 10;
                label_jmp_diff = diff(label_jmp);
                find(label_jmp_diff > jmp_delta);
                label_jmp = label_jmp(find(label_jmp_diff > jmp_delta));
                label_jmp = [1, label_jmp, size(label_m, 2)];
            end

            label_seg = [];
            %% find start and end within a hop signal
            if (~isempty(label_jmp)) && (~isempty(label_zero))
                label_seg = sort([label_jmp, label_jmp(2:end-1) + 1, label_zero_seg]);
            else
                label_seg = sort([label_jmp, label_jmp(2:end-1) + 1]);
            end


            a_est = [];
            l_temp = [];
            f_est = [];
            %% convert two stagnation point contianer map to matirx
            %% label_d_1 means antenna 1, label_d_2 means antenna 2
            label_d_1 = obj.label_d("1");
            label_d_2 = obj.label_d("2");

            mx_len1 = length(label_d_1("col" + string(0)));
            for i = 1:1:label_d_1.length - 1
                if mx_len1 < length(label_d_1("col" + string(i)))
                    mx_len1 = length(label_d_1("col" + string(i)));
                end
            end

            label_d_1_m = zeros(mx_len1, length(label_d_1));

            for i = 1:1:label_d_1.length
                label_d_1_m(:, i) = [label_d_1("col" + string(i-1)); zeros(mx_len1 - length(label_d_1("col" + string(i-1))), 1)];
            end

            mx_len2 = length(label_d_2("col" + string(0)));
            for i = 1:1:label_d_2.length - 1
                if mx_len2 < length(label_d_2("col" + string(i)))
                    mx_len2 = length(label_d_2("col" + string(i)));
                end
            end

            label_d_2_m = zeros(mx_len2, length(label_d_2));

            for i = 1:1:label_d_2.length
                label_d_2_m(:, i) = [label_d_2("col" + string(i-1)); zeros(mx_len2 - length(label_d_2("col" + string(i-1))), 1)];
            end
            

            for i = 1:2:length(label_seg)

                l_temp = [];
                for ii = label_seg(i):1:label_seg(i+1)

                    if (length(find(label_d_1_m(:, ii) > 0)) == num_est * 2) && (length(find(label_d_2_m(:, ii) > 0)) == num_est * 2)
                        a_est(:, ii) = obj.stft_tensor(label_d_2_m(1:num_est * 2, ii) , ii, 2) ./ obj.stft_tensor(label_d_1_m(1:num_est * 2, ii), ii, 1);
                        l_temp(:, ii) = label_m(1:num_est, ii);
                    end

                end
                                          
                %% delete false column
                del_label = [];

                for j = 1:1:size(l_temp, 2)
                    if l_temp(:, j) == zeros(size(l_temp, 1), 1)
                        del_label(j) = j;
                    end
                    del_label = del_label(find(del_label > 0));
                end
                
                l_temp(:, del_label) = [];
                if ~isempty(l_temp)
                    l_est(:,  round(i / 2)) = round(mean(l_temp, 2));
                else
                    l_est(:,  round(i / 2)) = 0;
                end
                
                a_est(:, del_label) = [];
                a_est_mod = sqrt(a_est(1:2:end, :) .* a_est(2:2:end, :));

                %%% bss method to estimate DOA   
                l = link16(1, 2, 30, obj.fs);
                f_est(:,  round(i / 2)) = l.ifreq_mapping((l_est(:,  round(i / 2)) - 1) * obj.fs / obj.dft_length * 1e-6);
                doa_est_m = asind(atan(imag(a_est_mod) ./ (real(a_est_mod) + 1e-31)) * 3e8 ./ (2 * pi * f_est(:,round(i / 2)) * 1e6 * 0.1));
                

                for k = 1:1:size(doa_est_m, 1)
                    doa_temp = doa_est_m(k, find(doa_est_m(k, :) > 0));
                    doa_temp = doa_temp(find(abs(doa_temp - mean(doa_temp)) < 5));
                    doa_est_(k, round(i / 2)) = mean(doa_temp);
                end
                
                
                l_temp = [];
                
            end

            hop =  diff(label_seg(2:end) * (obj.win_length - obj.overlap_length));

            hop_vec = label_seg;
            
              
            for i = 2:2:length(label_seg) - 1
                hop_vec(i) = hop_vec(i) * (obj.win_length - obj.overlap_length);
                hop_vec(i+1) = hop_vec(i) + 1;
            end
            
            hop_vec(end) = size(obj.rx.receive_signal, 2);
            
            %%% music algorithm to estimate DOA   
            % a_vec = zeros(obj.rx.antenna_num, 1);
            % delta_theta = 2;
            % angle_scan = -90:delta_theta:90 - delta_theta;
            % p_mu = zeros(1, length(angle_scan));

            % doa_est = zeros(size(f_est));

            % eig_val_m = [];
            % p_m = zeros(num_est, length(p_mu), round(length(label_seg) / 2));

            % for ii = 1:2:length(hop_vec)

            %     xx = obj.rx.receive_signal(:, hop_vec(ii) + 2:hop_vec(ii + 1) - 2);
            %     xx_len = size(xx, 2);

            %     rxx = zeros(obj.rx.antenna_num, obj.rx.antenna_num);
            %     rxx = (1 / xx_len) .* (xx * xx');

            %     [V, D] = eig(rxx);
            %     eig_val = diag(V);

            %     eig_val_m(round(ii / 2), :) = eig_val;

            %     Us = V(:, end - num_est + 1:end);
            %     Un = V(:, 1:end - num_est);

                
            %     for i = 1:1:size(f_est, 1)

            %         for theta = 1:length(angle_scan) 
            %             a = exp(1j * 2 * pi * 0.1 / (3e8 / (f_est(i, round(ii / 2)) * 1e6)) * sin(angle_scan(theta) / 180 * pi) * (0 : obj.rx.antenna_num - 1)');
            %             p_mu(theta) = 1 / (a' * Un * Un' * a);
            %             p_m(i, :, round(ii / 2)) = p_mu';
            %         end

            %         doa_est(i, round(ii / 2)) = mean(angle_scan(find(abs(p_mu) == max(abs(p_mu)))));
            %     end

            %     p_scan = p_m;

            % end

        end


        function [net_type] = net_dec(obj)
            label = obj.label("1");

            mx_len = length(label("col" + string(0)));
            for i = 1:1:label.length - 1
                if mx_len < length(label("col" + string(i)))
                    mx_len = length(label("col" + string(i)));
                end
            end
            
            label_m = zeros(mx_len, length(label));
            
            for i = 1:1:label.length
                label_m(:, i) = [label("col" + string(i-1)); zeros(mx_len - length(label("col" + string(i-1))), 1)];
            end

            label_num_m = [];

            for i = 1:1:size(label_m, 2)
                label_num_m(i) = length(find(label_m(:, i) > 0));
            end

            
            label_num_jmp = find(diff(label_num_m) ~= 0);
            label_num_m(label_num_jmp);

            if ~isempty(label_num_jmp)

                delta_greate = 10;
                span = 10;
                
                for i = 1:1:length(label_num_jmp)

                    if (label_num_jmp(i) < delta_greate) || (label_num_jmp(i) > length(label_num_m) - delta_greate)
                        continue;
                    else
                        cnt = 1;
                        while(cnt < span)
                            if (label_num_m(label_num_jmp(i) - cnt) == label_num_m(label_num_jmp(i) + cnt))
                                cnt = 20;
                                label_num_jmp(i) = 0;
                            end
                            cnt = cnt + 1;
                        end
                    end

                end

                if label_num_jmp(end) == length(label_num_m)
                    label_num_jmp(end) = 0;
                end

                label_num_jmp = label_num_jmp(find(label_num_jmp > 0));
            end

            if isempty(label_num_jmp)

                net_type = "syn";
            else
                net_type = "asyn";
            end

        end

    end 
    

end


