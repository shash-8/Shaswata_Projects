
%% just the plotting part of constellation diagram
clc; clear all; close all;


M_total = [4, 8, 16, 32];
%M_total = [8];
for iter = 1:length(M_total)

    M = M_total(iter);              
    x = (0:M-1)';       
    sym = qammod(x,M);  % constellation diagram
    % my addition
    real_sym=real(sym);
    imag_sym=imag(sym);
    
    figure();
    scatter(real_sym,imag_sym,"filled"); grid on;
    titleString = sprintf('Constellation Diagram for %d-QAM', M);
    title(titleString);
end
%% Modulation part + demod
M_total = [4, 8, 16, 32];
fig_count = length(M_total)+1;
for iter = 1:length(M_total)
    M = M_total(iter);          
    %M = 4;
    n_bits = log2(M);           % n_bits = no. of bits in a symbol
    
    no_symbols = 2000; 
    no_bits = no_symbols*n_bits;            % n_bits = no. of bits in a symbol 
                                    % no_bits = total no of bits taken into consideration
    Rb = 1000;                  % data bit rate
    samples_per_bit = 20;                     % samples per bit
    Tb = 1/Rb;                  % 1 bit duration
    fc = 2*Rb;                  % carrier frequency
    Fs = samples_per_bit*Rb;                  % simulation sampling frequency
    Ts = 1/Fs;                  % considered time step (important)
    t = 0:Ts:no_bits*Tb-Ts;           % considered total time interval
    
    bits_stream = randi([0,1],1,no_bits);      

    % Binary to symbol
    temp = reshape(bits_stream, [n_bits, no_bits/n_bits])';
    symbols_stream = (bi2de(temp, 'left-msb'))';
    
    modulated_stream = qammod(symbols_stream, M);           % stream_modulated is complex
    
    % Time sampling the data sequences to match time axis
    bit_upsampling = length(t)/length(bits_stream);
    symbol_upsampling = length(t)/length(symbols_stream);
    
    
    bits_time = reshape(repmat(bits_stream, bit_upsampling, 1), [1, length(t)]);
    symbols_time = reshape(repmat(symbols_stream, symbol_upsampling, 1), [1, length(t)]);
    mod = reshape(repmat(modulated_stream, symbol_upsampling, 1), [1, length(t)]);
   
    % Modulation
    modulated_sig = real(mod).*cos(2*pi*fc*t)+imag(mod).*sin(2*pi*fc*t);
    
    % Plot
    figure(fig_count); fig_count = fig_count++1;
    subplot(311) 
    plot(t(1:20*samples_per_bit),bits_time(1:20*samples_per_bit), 'LineWidth',2, 'Color', [0.5, 0.25, 0])
    ylim([-0.2, 1.2]);
    titleString = sprintf('Bit Sequence (%d-QAM)', M);
    title(titleString);
    
    subplot(312) 
    plot(t(1:20*samples_per_bit),symbols_time(1:20*samples_per_bit), 'LineWidth',2, 'Color', [0.5, 0.25, 0])
    ylim([min(symbols_time)-1, max(symbols_time)+1]);
    titleString = sprintf('Symbol Sequence (%d-QAM)', M);
    title(titleString);
   
    subplot(313) 
    plot(t(1:20*samples_per_bit), modulated_sig(1:20*samples_per_bit),'LineWidth',2, 'Color', [0.5, 0.25, 0])
    titleString = sprintf('Modulated Signal (%d-QAM)', M);
    title(titleString);
     
    %% Power Spectrum
    
    [pxx, frequency] = pwelch(modulated_sig,500,250,512,Fs);
    pxx = pxx/max(pxx);
    
    figure(fig_count); fig_count = fig_count +1;
    
    %subplot(510+M_idx)
    plot(frequency, pxx,'LineWidth', 2, 'Color', [0.5, 0.25, 0]); grid on;
    xlabel('Frequency');
    ylabel('PSD')
    titleString = sprintf('PSD of Modulated Signal (%d-QAM)', M);
    title(titleString);
    
    %% Demodulation
    
    SNR = -25:5:15;
    BER = [];
    
    fig_no_LAST = length(M_total)*(2)+1+length(M_total)+1;

    symbol_duration_time = 0:Ts:n_bits*Tb-Ts;
    carrier = cos(2*pi*fc*symbol_duration_time);    % main carrier
    carrier_quad = sin(2*pi*fc*symbol_duration_time);  % quadrature carrier
    
    carrier_len = length(carrier);
    
    for i = 1:length(SNR)
        signal_with_noise = awgn(modulated_sig, SNR(i), 'measured');
        
        % Multiplication by carrier and averaging
        a_i = carrier*reshape(signal_with_noise, [carrier_len, length(symbols_stream)])/carrier_len;
        a_q  = carrier_quad*reshape(signal_with_noise, [carrier_len, length(symbols_stream)])/carrier_len;
        
        decoded_symbols = qamdemod(a_i + 1i*a_q, M);
        
        %decoded_bits = fliplr(de2bi(decoded_symbols, n_bits))';
        decoded_bits = (de2bi(decoded_symbols, n_bits,'left-msb'))';
        decoded_msg = decoded_bits(:)';
        
        BER = [BER, mean(abs(bits_stream-decoded_msg))];
    end
    
    figure(fig_no_LAST);
    plot(SNR, BER, 'LineWidth', 2);
    xlabel('SNR in dB');
    ylabel('Bit Error Rate')
    titleString = sprintf('BER vs SNR (%d-QAM)', M);
    title(titleString);
    
    grid on;
    hold on;
    
end

figure(fig_no_LAST);
legend('4-QAM', '8-QAM', '16-QAM', '32-QAM');
title('BER vs SNR [M-QAM]');