%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% class : 
%%%        rx_signal
%%% fea: 
%%%        generate receive signal
%%% parameter: 
%%%    antenna_num: the number of receive antenna
%%%  element_distance: distance between two adjacent antenna
%%%         rx_snr: receive snr
%%%     src_signal: source (multiple) fh signal class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef rx_signal

    properties

        antenna_num;
        element_distance;
        rx_snr; 
        src_signal;
        receive_signal;    %% receive signal
        mix_matrix;        
        theta_pattern;     
        freq_pattern;   
        f_pattern;
    end


    methods
        %% constructor
        function obj = rx_signal(antenna_num, element_distance, rx_snr, src_signal)
            %% source (multiple) fh signal class
            obj.src_signal = src_signal;
            obj.antenna_num = antenna_num;
            obj.element_distance = element_distance;
            obj.rx_snr = rx_snr;
            %% get freq_pattern from source (multiple) fh signal class
            obj.freq_pattern = obj.src_signal.freq_pattern;
            %% get theta_pattern from source (multiple) fh signal class
            obj.theta_pattern = obj.src_signal.theta_pattern;
            obj.mix_matrix = obj.generate_array_structure();
            obj.receive_signal = obj.receive_src_signal();

        end

        %% ULA delay factor
        function ula_row_vec = ula_row_element(obj, element_distance, freq, theta)

            ula_row_vec = exp(1j * 2 * pi  * element_distance * freq .* sind(theta) / (3e8));

        end

        %% generate min_matirx according to freq_pattern and theta_pattern
        %% only for ULA
        function mix_matrix = generate_array_structure(obj)
            
            A = zeros(obj.antenna_num, size(obj.freq_pattern, 1), size(obj.freq_pattern, 2));
            A(1, :, :) = ones(size(obj.freq_pattern, 1), size(obj.freq_pattern, 2));
            A(2, :, :) = ula_row_element(obj, obj.element_distance, obj.freq_pattern * 1e6, obj.theta_pattern);

            if obj.antenna_num > 2
                for iii = 3:obj.antenna_num
                    A(iii, :, :) = power(A(2, :, :), iii-1);
                end
            end

            mix_matrix = A;

        end

        %% generate receive signal according to source signal and mix matrix
        function receive_signal = receive_src_signal(obj)
            %% receive signal matrix
            receive_signal = zeros(obj.antenna_num, obj.src_signal.hop_num * obj.src_signal.hop_length); 
           
            for j = 0:1:obj.src_signal.hop_num - 1
                %% receive signal with agwn
                receive_signal(:, j * obj.src_signal.hop_length + 1: (j+1) * obj.src_signal.hop_length) ...
                    = awgn(obj.mix_matrix(:, :, j+1) * obj.src_signal.src_signal(:, j * obj.src_signal.hop_length + 1: (j+1) * obj.src_signal.hop_length), obj.rx_snr, 'measured');
            
                end

        end

    end

end
