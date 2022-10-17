%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% class : 
%%%        fh
%%% fea: 
%%%        generate frequency hopping signal
%%% parameter: 
%%%             fs: sample rate
%%%     member_num: the number of fh signal
%%%        hop_num: the number of hop
%%%     hop_length: sample in a hop
%%%   network_interval: frequency interval between two adjacent fh signal (in MHz)
%%%   freq_pattern: freq_pattern 
%%%  theta_pattern: theta_pattern
%%%  modulation_para: fh signal modulation parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef fh
    
    properties 
        freq_pattern;
        theta_pattern
        member_num;
        % modualtion_table;
        hop_num;
        hop_length;
        fs;
        modulation_para;
        modulation;
        network_interval;
        carrier;      %% DUC carrier
        src_signal;

    end

    methods
        %% constructor
        function obj = fh(fs, member_num, hop_num, hop_length, network_interval, freq_pattern, theta_pattern, modulation_para)

        obj.fs = fs;
        obj.member_num = member_num;
        obj.hop_num = hop_num;
        obj.hop_length = hop_length;
        obj.network_interval = network_interval;
        obj.freq_pattern = freq_pattern;
        obj.theta_pattern = theta_pattern;
    
        %% check two numbers are equal
        if ~(obj.member_num == size(fieldnames(modulation_para), 1))
            throw("member_num error!")
        end

        obj.modulation_para = modulation_para;

        %% get member_name
        member_name = string(fieldnames(obj.modulation_para));
        %% generate empty mapping 
        modulation_map = containers.Map();
        
        %% full length carrier matirx
        carrier_all = zeros(obj.member_num, obj.hop_length * obj.hop_num);

        %% save modulated signal
        src_signal = [];

        for i = 1:obj.member_num

            carrier = [];

            %% generate different carrier according to frequency pattern
            for j = 1:obj.hop_num
                %% call complex_exponential_wave class to generate carrier
                %% instantiate complex_exponential_wave class and get member property sample_seq
                carrier = [carrier; complex_exponential_wave(obj.freq_pattern(i, j) * 1e6, obj.fs, obj.hop_length).sample_seq];     
            end

            carrier_all(i, :) = carrier;

            %% save modulation parameter
            kk = obj.modulation_para.(member_name(i));

            %% only support to msk modulation
            if kk.("mod") == "msk"
                %% instantiate msk_modulation class
                modulation_obj = msk_modulation(kk.("symbol_rate"), obj.fs, obj.hop_length * obj.hop_num);
                %% call member method to realize DUC
                modulation_obj = modulation_obj.mixer(carrier, 1);
            end
            
            %% mapping : member_name -> msk_modulation class
            modulation_map(member_name(i)) = modulation_obj;
            %% generate fh signal matirx
            src_signal(i,:) = modulation_obj.modulated_seq_up;


        end
        obj.modulation = modulation_map;
        obj.carrier = carrier_all;
        obj.src_signal = src_signal;

        end

    end


end