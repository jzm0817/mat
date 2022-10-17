%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% class : 
%%%        msk_modulation
%%% fea: 
%%%        msk modulation
%%% parameter: 
%%%    symbol_rate: symbol rate
%%%    sample_rate: sample rate
%%%  sample_length: length of msk modulated signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef msk_modulation

    properties 
        symbol_rate ;
        sample_rate ;
        sample_per_symbol;   %% equal to sample_rate / symbol_rate
        bit_rate ;           %% for msk, bit_rate = symbol_rate
        modulation_type;     %% 'msk'
        sample_length;
        t;                   %% basic time vector
        initial_bit_length; 
        modulated_seq;       %% base band signal
        modulated_seq_up;    %% frequency band signal
        initial_bit;         %% random symbol
    end

    methods 
        %% constructor
        function obj = msk_modulation(symbol_rate, sample_rate, sample_length)
            obj.symbol_rate = symbol_rate;
            obj.sample_rate = sample_rate;
            obj.bit_rate = symbol_rate;
            obj.sample_per_symbol = sample_rate / symbol_rate;
            obj.sample_length = sample_length;
            obj.initial_bit_length = ceil(obj.sample_length / obj.sample_per_symbol);
            obj.modulation_type = "msk";
            [obj.modulated_seq, obj.initial_bit] = generate_modulated_seq(obj);
            obj.t = (0:1:obj.sample_length-1) * (1 / obj.sample_rate);
            obj.t = obj.t';
        end

        %% generate base band signal
        function [modulated_seq, initial_bit] = generate_modulated_seq(obj)
            %% generate random symbol
            modulate_symbol = randi([0, 1], obj.initial_bit_length, 1);
            initial_bit = modulate_symbol;
            %% call function 'mskmod' to generate base band modulated signal
            modulated_seq = mskmod(modulate_symbol, obj.sample_per_symbol);
            %% cut the base band modulated signal according to parameter 'sample_length'
            modulated_seq = modulated_seq(1:obj.sample_length);
            if obj.initial_bit_length == 1
                modulated_seq = modulated_seq';
            end

            carrier0 = complex_exponential_wave(obj.symbol_rate * 1, obj.sample_rate, obj.sample_length).sample_seq;
            %% base band modulated signal
            modulated_seq = real(modulated_seq .* carrier0);
        end

        %% DUC
        function obj = mixer(obj, carrier, flag)   

            if flag
                obj.modulated_seq_up = real(obj.modulated_seq .* carrier);
            end

        end

    end
    
end