clear all
%% ------Part A------

Fs = 40000; %Sampling rate
t = (0:2^19)*(1/Fs); %Time-sample points 32767
omega = 2*pi*(5400/Fs); %Convert 5400 Hz into units of omega
n = (0:2^19);

figure(1)

x = (sqrt(2))*sin(omega*n); %The input sinusoid

subplot(2,1,1)
plot(t(1:100),x(1:100),'-o') %Plot of first 100 points of input sinusoid vs time
title('Value of a Sinusoid vs. Time: 100 Points')
xlabel('Time (s)')

X_4096 = fft(x,4096); %DFT of the first 4096 points
X_4096_dB = 20*log10(abs(fftshift(X_4096))); %shifted magnitude of X_4096 in dB
f = (-2048:2047)*(Fs/4096); %Frequency axis values

subplot(2,1,2)
plot(f,X_4096_dB)
title('Magnitude of DFT vs. Frequency')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%% ------Part B------
figure(2)

v = randn(size(x)); %Creates random noise 

subplot(2,1,1)
plot(t(1:1000),v(1:1000)) %plots first 1000 noise points
title('Random Noise vs. Time')
xlabel('Time (s)')

V_4096 = fft(v,4096); %Get DFT of the random noise
V_4096_dB = 20*log10(abs(fftshift(V_4096))); %Shifted magnitude of V_4096
Pnf = 4096; %Noise floor power
Pnf_dB = 10*log10(Pnf); %Convert noise floor power to dB

subplot(2,1,2)
plot(f,V_4096_dB,f,Pnf_dB,'r.')
title('Magnitude of DFT of Random Noise vs. Frequency')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%% ------Part C------

r1 = x+v;
SNR = 10*log10(sum(x.^2)/sum((r1-x).^2)); % SNR value from derivations
SNRline = SNR.*ones(100,1); % creates constant line to show on graph

figure(3)

plot(t(1:100),r1(1:100),t(1:100),x(1:100),t(1:100),SNRline)
legend('sinusoidal signal + noise','sinusoidal signal','SNR')
title('Signal and Noise compared to SNR')
xlabel('Time (s)')
grid

figure(4)

%For SNR=-20dB
SNR_in = -20; % finds r values for SNR_in = -20
b=(sum(x.^2)/sum(v.^2))*10^-(SNR_in/10);
r=x+(v*sqrt(b));

subplot (3,1,1)
plot(t(1:100),r(1:100),t(1:100),x(1:100))
legend('sinusoidal signal + noise','sinusoidal signal')
title('Signal and Noise compared to SNR (-20dB)')
xlabel('Time (s)')
grid

%For SNR=10dB
SNR_in = 10; % finds r values for SNR_in = 10
b=(sum(x.^2)/sum(v.^2))*10^-(SNR_in/10);
r=x+(v*sqrt(b));

subplot (3,1,2)
plot(t(1:100),r(1:100),t(1:100),x(1:100))
legend('sinusoidal signal + noise','sinusoidal signal')
title('Signal and Noise compared to SNR (10dB)')
xlabel('Time (s)')
grid

%For SNR=30dB
SNR_in = 30; % finds r values for SNR_in = 30
b=(sum(x.^2)/sum(v.^2))*10^-(SNR_in/10);
r=x+(v*sqrt(b));

subplot (3,1,3)
plot(t(1:100),r(1:100),t(1:100),x(1:100))
legend('sinusoidal signal + noise','sinusoidal signal')
title('Signal and Noise compared to SNR (30dB)')
xlabel('Time (s)')
grid

%% ------Part D------

SNR_in = 20; % input SNR value
b=(sum(x.^2)/sum(v.^2))*10^-(SNR_in/10); % Equation for b from previous part
r=x+(v*sqrt(b)); 

SNR_out=zeros(6,1);ps=SNR_out;SNR_out_zp=ps; % pre-setting these values before the loop
for i=1:6
    
    figure(5) % This creates the DFT at varying sample sizes
        
        N=2^(9+i); % Increase N value from 1024 by power of two each iteration
        freq=(-N/2:N/2-1)*Fs/N; %creates bounds with N values
        R = fft(r(1:N),N); %calculates DFT of r
        R_dB = 20*log10(abs(fftshift(R))); % converts R to dB
        vi=v(1:N); % limits the v value to this iteration's N
        Noise_floor = 10*log10(sum(b*vi.^2)); % Set value of Noise_floor f
        
        DFT_peak = max(R_dB); % calculates peak of graph
        Nf_line = Noise_floor*ones(N,1); % creates flat line to graph Noise floor
        SNR_out(i) = DFT_peak - Noise_floor; % sets each iteration of SNR_out according to given formula
        
        ps(i)=log2(N); % ps represents power of number of samples on each iteration
        
        %plots without zero-padding
        
        subplot(3,2,i) %Creates next subplot on every iteration
        plot(freq,R_dB,freq,Nf_line,'r')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        axis([-2e4 2e4 -20 120]) % These axes adequately encompass all values we will need
        grid
        title(['DFT with ',num2str(N),' Samples'])
    
    figure(7) %This creates the same graphs but with zero padding 
        
        N_zp=4*N; % N_zp is zero-padded number we want to use
        freq = ((-N_zp/2):((N_zp/2)-1))*Fs/N_zp; %Creates N_zp values of frequency from -2kHz to 2kHz
        R_zp = fft(r(1:N),N_zp); % DFT of r, but will have N_zp values
        R_dB_zp = 20*log10(abs(fftshift(R_zp))); % zero-padded magnitude of the DFT of output signal, r
        Noise_floor_zp = 10*log10(sum(b*vi.^2)/4); % calculates noise floor for zp
        
        DFT_peak_zp = max(R_dB_zp); % These all have same function as before without zero-padding, but with all zp values
        Nf_line = Noise_floor_zp*ones(N_zp,1);
        SNR_out_zp(i) = DFT_peak_zp - Noise_floor_zp; 
        
        subplot(3,2,i) % same process as before without zero padding 
        plot(freq,R_dB_zp,freq,Nf_line,'r')
        xlabel('Frequency (Hz)')
        ylabel('|R| (dB)')
        axis([-2e4 2e4 -20 120])
        title(['Zero-Padded DFT with ',num2str(N),' Samples'])
        grid
end

figure(6) % displays the values of SNR_out vs the power of two without zero padding

plot(ps,SNR_out,'-o')
title('Plot of SNR_{out} before Zero-Padding')
xlabel('Power of 2')
ylabel('SNR_{out} (dB)')


figure(8) % displays the values of SNR_out vs the power of two with zero padding

plot(ps,SNR_out_zp,'-o')
xlabel('Power of 2')
ylabel('SNR_{out} (dB)')
title('Plot of SNR_{out} with Zero-Padding')