%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% class : 
%%%        complex_exponential_wave
%%% fea: 
%%%        generate complex carrier 
%%% parameter: 
%%%             fc: carrier frequency
%%%             fs: sample rate
%%%  sample_length: length of carrier
%%%       varargin: option, support to set phase,amplitude and output energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef complex_exponential_wave

    properties

        fc;
        fs;
        sample_length;
        amplitude = 1;         %% default :amplitude = 1
        phase = 0;             %% default :phase = 1
        union = 0;             %% default :union = 0, not to normalize output energy
        sample_seq;            %% complex carrier
        t;                     %% basic time vector
        type = "carrier";      %% introduction of class
        
    end
    
    methods 
        %% constructor 
        function obj = complex_exponential_wave(fc, fs, sample_length, varargin)

            obj.fc = fc;
            obj.fs = fs;
            obj.sample_length = sample_length;
            
            %% selsction of optional parameter 
            if ~isempty(varargin)
                optional_para = cell2mat(varargin);
                optional_para_name = string(fieldnames(cell2mat(varargin)));
                option_para_len = size(optional_para_name, 1);
                
                if option_para_len > 3
                    throw("the length of option parameters should less than 3!\n");
                end

                if size(find(optional_para_name == "phase"), 1) > 0
                    obj.phase = optional_para.(optional_para_name(find(optional_para_name == "phase")))
                end

                if size(find(optional_para_name == "union"), 1) > 0
                    obj.union = optional_para.(optional_para_name(find(optional_para_name == "union")))
                end

                if size(find(optional_para_name == "amplitude"), 1) > 0
                    obj.amplitude = optional_para.(optional_para_name(find(optional_para_name == "amplitude")))
                end

            end

            obj.t = (0:1:obj.sample_length-1) * (1 / obj.fs);
            obj.t = obj.t';
            obj.sample_seq = exp(1j * (2*pi * obj.fc * obj.t + obj.phase));    %% generate complex carrier

            %% normalize output energy
            if obj.union == 1
                obj.sample_seq = obj.sample_seq / sqrt(sum((conj(obj.sample_seq)' * obj.sample_length)));
            end
                   
        end
    end

end