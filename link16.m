%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% class : 
%%%        link16
%%% fea: 
%%%        generate frequency pattern
%%% parameter: 
%%%        num: the number of fh signal
%%%    hop_num: the number of hop 
%%% net_interval:frequency interval between two adjacent fh signal (in MHz)
%%%         fs: sample rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef link16

    properties 
        freq_table = [969:3:1008, 1053:3:1065, 1113:3:1206];    %% link 16 fh range
        freq_pattern;
        freq_mapping_base;
        fs;
    end

    methods 
        %% constructor
        function obj = link16(num, hop_num, net_interval, fs)

        obj.fs = fs;
        if num == 1
            obj.freq_pattern = obj.freq_table_samll(hop_num, 30);
        else
            obj.freq_pattern = obj.generate_freq_pattern(num, hop_num, net_interval);

        end
            obj.freq_mapping_base = obj.freq_mapping(obj.freq_table, obj.fs * 1e-6);
        end


        %% generate frequency pattern for single fh signal
        function freq_small = freq_table_samll(obj, num, interval)

            f = randsrc(1, 1, obj.freq_table);

            while length(f) < 2
                f_sel = randsrc(1, 1, obj.freq_table);
                if (abs(f_sel - f(1)) > interval)
                    f = [f, f_sel];
                end
            end

            while length(f) < num
                f_sel = randsrc(1, 1, obj.freq_table);
                f = sort([f, f_sel]);
                label = find(f == f_sel);
                if length(label) > 1
                    label = label(1);
                end


                if label == 1
                    if abs(f(2) - f_sel) < interval
                        f(label) = [];
                    end
                
                elseif label == length(f)
                    if abs(f(label-1) - f_sel) < interval
                        f(label) = [];
                    end

                else
                    if ((abs(f_sel - f(label - 1)) < interval) || (abs(f_sel - f(label + 1)) < interval))
                        f(label) = [];
                    end
                end
            
            end
            
            freq_small = f(randperm((size(f, 2))));
        end

        %% generate frequency pattern for multiple fh signal
        function freq_pattern = generate_freq_pattern(obj, num, N, net_interval)

            f = freq_table_samll(obj, num, 30)';

            if num > 1
                cnt = 1;
                while cnt < N
                    flag = 1;
                    
                    while flag
                        f1 = freq_table_samll(obj, num, 30)';

                        if (sum(abs(f1 - f(:, size(f, 2))) > net_interval) == num)
                            f = [f, f1];
                            flag = 0;
                            cnt = cnt + 1;
                        end
                    end
                end
                
            end

            if num == 1
                freq_pattern = f';
            else
                freq_pattern = f;
            end

        end

        %% frequency mapping : origin freq  -- fs -->  after sampling freq
        function mapping_matrix = freq_mapping(obj, input, sample_rate) 
            fs = sample_rate;
            mapping = zeros(size(input));
            for ii = 1:1:size(input, 1)
                for jj = 1:1:size(input, 2)
                    n = floor(input(ii, jj) / fs);
                    f0 = input(ii, jj)  - (n:1:n+1) * fs;
                    f1 = -input(ii, jj)  + (n:1:n+1) * fs;
                    f0 = f0((f0 < fs / 2) & (f0 > - fs / 2));
                    f1 = f1((f1 < fs / 2) & (f1 > - fs / 2));
                    if ~isempty([abs(f0), abs(f1)])
                        mapping(ii, jj) = unique([abs(f0), abs(f1)]);
                    else
                        mapping(ii, jj) = NaN;
                    end
                end
            end
            mapping_matrix = mapping;
        end

        %% frequency mapping : after sampling freq  -- fs -->  origin freq
        function freq_matrix = ifreq_mapping(obj, input)

            mapping = zeros(size(input));

            for ii = 1:1:size(input, 1)
                for jj = 1:1:size(input, 2)
                    for kk = 1:1:length(obj.freq_table)
                        if abs(input(ii, jj) - obj.freq_mapping_base(kk)) < 1.5
                            mapping(ii, jj) = obj.freq_table(kk);
                        end
                    end
                end
            end
            freq_matrix = mapping;
        end

    end

end




