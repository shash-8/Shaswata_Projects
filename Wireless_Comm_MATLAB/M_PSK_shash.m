%% Full M-PSK Modulation, Demodulation, BER
clc; clear all; close all;
%%
M_total = [2, 4, 8, 16, 32];          % M-PSK
%M_total = [6,10];
fig_count = 1;
for iter = 1:length(M_total)
        M = M_total(iter);
        phi_M = (0:M-1)*2*pi/M;          % corresponding angles for M symbols
        n_bits = log2(M);               % n_bits = no. of bits in a symbol
        %
        A = 5;
        
        % assuming X axis is sin(2pifct) and Y axis is cos(2pifct)
        % Sk(t) = A*cos(phi_M)*sin(2pifct) + A*sin(phi_M)*cos(2pifct)    [ % modulated signal ]
        coeff_X = A*cos(phi_M);
        coeff_Y = A*sin(phi_M);
        
        % rounding upto 3 decimal points
        coeff_X = round(coeff_X * 1000) / 1000;
        coeff_Y = round(coeff_Y * 1000) / 1000;

        figure(fig_count); fig_count = fig_count+1 ;
   
        scatter(coeff_X, coeff_Y, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 2);
        titleString = sprintf('Constellation Diagram for %d-PSK', M);
        title(titleString);
        grid on
        %
        no_symbols = 400;           % total number of symbols taken into consideration
        no_bits = no_symbols*n_bits;      % n_bits = no. of bits in a symbol which must be multiple of log2(M)
                                    % no_bits = total no of bits taken into consideration
        Rb = 1000;                  % data bit rate
        samples_per_bit = 40;       % samples per bit
        Tb = 1/Rb;                  % bit duration
        fc = 2*Rb;                  % carrier frequency
        Fs = samples_per_bit*Rb;                  % sampling frequency
        Ts = 1/Fs;                  % considered time step
        t = 0:Ts:no_bits*Tb-Ts;           % considered total time interval
            
       
        bits_stream = randi([0,1],1,no_bits);       % 1D row vector
        %
        temp = reshape(bits_stream, [n_bits, no_bits/n_bits])';   % temp dimension = [no of symbols , no of bits in a symbol] 
                                                                  % [400,2] in this case
        symbols_stream = (bi2de(temp, 'left-msb'))';
        symbols_stream_gray = bin2gray(symbols_stream,'psk',M);         % it basically interchanges '2' and '3' [only for qpsk]
        
        stream_phases = symbols_stream_gray*2*pi/M;     % phases for symbols
        % upsampling vector to match timeaxis
        symbol_upsampling = length(t)/length(stream_phases);
        bit_upsampling = length(t)/length(bits_stream);
        %    
        bits_time = reshape(repmat(bits_stream, bit_upsampling, 1), [1, length(t)]);            % bits array matching time dimension
        phi = reshape(repmat(stream_phases, symbol_upsampling, 1), [1, length(t)]);                % phi angle array matching time dimension
        symbols_time = reshape(repmat(symbols_stream, symbol_upsampling, 1), [1, length(t)]);   % symbols array matching time dimension
        
            
        %% Modulation
        modulated_signal = A*sin(2*pi*fc*t + phi);
        
        figure(fig_count); fig_count = fig_count+1 ;
        
        subplot(311)
        plot(t(1:20*samples_per_bit),bits_time(1:20*samples_per_bit), 'LineWidth',2, 'Color', [0.5, 0.25, 0])
        titleString = sprintf('Bit Sequence (%d-PSK)', M);
        title(titleString);
        legend_label = sprintf('Bit Duration %.3f sec', Tb);
        legend(legend_label, 'Location', 'northeast')
        ylim([-0.1,1.1]);
        xlabel('time (seconds)'); ylabel('Labels');
        grid on
        
        subplot(312)
        plot(t(1:20*samples_per_bit),symbols_time(1:20*samples_per_bit),'LineWidth',2, 'Color', [0.5, 0.25, 0])
        %titleString = sprintf('Symbol Sequence (%d-PSK)', M);
        title('Symbol Sequence')
        legend_label = sprintf('Symbol Duration %.3f sec', Tb*n_bits);
        legend(legend_label, 'Location', 'northeast')
        ylim([-1,M]);
        xlabel('time (seconds)'); ylabel('Labels');
        grid on
        
        subplot(313)
        plot(t(1:20*samples_per_bit), modulated_signal(1:20*samples_per_bit),'LineWidth',2, 'Color', [0.5, 0.25, 0])
        title('Modulated Signal')
        % legend_label = sprintf('Bit Duration %.3f sec', Tb);
        legend(legend_label, 'Location', 'northeast')
        ylim([-1-A,1+A]);
        xlabel('time (seconds)'); ylabel('Labels');
        grid on
        
        % PSD
        [pxx,frequency] = pwelch(modulated_signal,250,125,512*2,Fs);  %pwelch(modulated_signal,500,250,512,Fs)
        pxx = pxx/max(pxx);
        figure(fig_count);  fig_count = fig_count+1 ;
        plot(frequency, pxx,'LineWidth', 2, 'Color', [0.5, 0.25, 0]); grid on;
        xlabel('Frequency');
        ylabel('PSD')
        titleString = sprintf('PSD of Modulated Signal (%d-PSK)', M);
        title(titleString);
       
        %% Demodulation and BER vs SNR
        SNR = -25:5:5;         % defining range of SNR (dB) that we will check
        %SNR = [-30];
        BER = [];
        fig_no_LAST = length(M_total)*(length(SNR)+3)+1;   
        symbol_duration_time = 0:Ts:n_bits*Tb - Ts;        % 1 symbol duration over timeaxis, here the upper limit is significant (subtracting Ts)
                                                        
        carrier = sin(2*pi*fc*symbol_duration_time);    % carrier sig for demodulation (taken over 1 symbol)
        carrier_quad = cos(2*pi*fc*symbol_duration_time);  % quadrature carrier sig for demodulation
            
        carrier_len = length(carrier_quad);             % no of samples over 1 symbol
        %    
            for i = 1:length(SNR)
                signal_with_noise = awgn(modulated_signal, SNR(i), 'measured');    % contains the noisy signal over full interval
        
                % suppose, no of symbols=400, and no of samples in 1 symbol is 80.
                % so noisy signal length is (400*80). now we will divide it into
                % 400x80 2D matrix, and multiply each symbol (length of 80) with
                % carrier signal. to match dimension of matrix multiplication,
                % necessaty reshaping is done. Finally we need to divide by
                % 'carrier_len' for averaging.
                
                % Multiply by carrier and do average
                a_i = carrier*reshape(signal_with_noise, [carrier_len, length(symbols_stream)])/carrier_len;      % a_i & a_q will be of 400 length 1D matrix
                a_q = carrier_quad*reshape(signal_with_noise, [carrier_len, length(symbols_stream)])/carrier_len;
                
                % searching distance from constellation points and find the minimum one
                constellation_distance = (A*cos(phi_M')-2*a_i).^2+(A*sin(phi_M')-2*a_q).^2;
                [~, decoded_gray] = min(constellation_distance); % just to get the index
                
                decoded_gray = decoded_gray - 1; % graycode starts from 0 but index starts from 1
                decoded_symbols = bin2gray(decoded_gray,'psk',M);    % basically interchanging '2' and '3' for qpsk
                %decoded_bits = fliplr(de2bi(decoded_symbols, n_bits))';
                decoded_bits = (de2bi(decoded_symbols, n_bits,'left-msb'))';
                decoded_msg = decoded_bits(:)';         % size(decoded_msg) = [1,total_no_bits]
                %size(decoded_msg)
                
                figure(fig_count); fig_count = fig_count+1 ;
                
                subplot(211)
                stairs(bits_stream(1:20),'LineWidth', 2, 'Color',[0.5,0.25,0]);     % plotting first 20 bits
                titleString = sprintf('Transmitted Bits ; SNR=%d ; (%d-PSK)',SNR(i), M);
                title(titleString);
                xlim([1,20]);
                ylim([-0.2,1.2])
                xlabel("Bit Sequence number"); grid on;
                
                subplot(212)
                stairs(decoded_msg(1:20),'LineWidth', 2, 'Color',[0.5,0.25,0]);     % plotting first 20 bits
                title('Decoded Bit Stream');
                xlim([1,20]);
                ylim([-0.2,1.2])
                xlabel("Bit Sequence number"); grid on;
                BER = [BER, mean(abs(bits_stream-decoded_msg))];
            end

            %BER
            figure(fig_no_LAST); % fig_count = fig_count+1 ;
            plot(SNR, BER, 'LineWidth', 2);
            xlabel('SNR in dB');
            ylabel('Bit Error Rate');
            grid on;
            hold on;
                
end   
figure(fig_no_LAST);
legend('2 PSK', '4 PSK', '8 PSK', '16 PSK', '32 PSK');
title('BER vs SNR [M-PSK]');

%% here bin2gray() and gray2bin() has been used
% there is an alternate way to implement both functions

% implement bin2gray() alternate:
% psk_gray = pskmod(psk_symbols,M);
% psk_phases=angle(psk_gray);

% implement gray2bin() alternate:
% constellation_distance = (A*cos(phiM')-2*ai).^2+(A*sin(phiM')-2*aq).^2;
% [~, demod_graycode] = min(constellation_distance);
% demod_graycode = demod_graycode - 1;
% demod_symbols = pskmod(demod_graycode,M);
% demod_symbols=angle(demod_symbols);
% demod_symbols(demod_symbols<0)=demod_symbols(demod_symbols<0)+2*pi;
% demod_symbols=(demod_symbols*M)/(2*pi);
% demod_symbols=int32(demod_symbols);
% demod_bits = fliplr(de2bi(demod_symbols, n_bits))';
% demod_msg = demod_bits(:)';
