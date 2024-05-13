%% compare plots
% STFT related parameters
clear;
addpath("FromThomasDietzen\")
N_STFT = 2048;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));

r = 0.057; % coordinate unit length
c = 343.3; % speed of sound
% microphone system parameters definition
m_pos = [r 0 0; r/2 -sqrt(3)/2*r 0; -r/2 -sqrt(3)/2*r 0; -r 0 0; -r/2 sqrt(3)/2*r 0; r/2 sqrt(3)/2*r 0]';
% sound source parameters definition
azimuth_deg = 0;
elevation_deg = -90;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = 50*r*[cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

for setNr = 1:2
    if setNr == 1
        totalFileNr = 16;
    else
        totalFileNr = 8;
    end
    for fileNr = 1:totalFileNr
        % load data
        % mvdr coefficients calculation data time domain
        [N, fs_noiseOnly] = audioread(sprintf("Temporary/toBeTested/set%d_Noise (%d).flac", setNr, fileNr));
        % mvdr coefficients calculation data frequency domain
        N_stft = calc_STFT(N, fs_noiseOnly, win, N_STFT, R_STFT, 'onesided');
        
        % input data time domain
        [input_t, fs_input] = audioread(sprintf("Temporary/toBeTested/set%d_Recording (%d).flac", setNr, fileNr));
        % input data frequency domain
        input_stft = calc_STFT(input_t, fs_input, win, N_STFT, R_STFT, 'onesided');
        xTickProp = [0, R_STFT/fs_input, 0]; % R_STFT/fs_input, fs_input/R_STFT
        yTickProp = [0, fs_input/(2000*R_STFT), R_STFT/2];
        cRange    = [-45 15];
        % plot input
        fig_in = figure;
        iterations = numel(input_stft(1,1,:));
        for i = 1:iterations
            subplot (iterations, 1, i);
            plotSpec(input_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('Freq (kHz)');
        end
        xlabel("Time");
        
        % initialize the output
        output_stft = zeros(numel(input_stft(:,1,1)),numel(input_stft(1,:,1)),1);
        
        % pipeline
        % calculate psd
        [P_NN_smth, P_NN_mean] = estim_corrmat(N_stft, 1);
        P_NN = squeeze(P_NN_mean);
        n_freq_bins = numel(input_stft(:,1,1));
        
        % delay and sum algorithm
        dasb_delay = s_pos*m_pos/norm(s_pos)/c; % to compensate the delay aka alignment: times "-" to a "-"
        
        % calculate coefficients then apply
        for i = 1:n_freq_bins
            % delay and sum algorithm
            d_dasb = exp(-1j*2*pi*((fs_input/2)/(N_STFT/2+1)*(i-1))*dasb_delay)/numel(m_pos(1,:));
            % mvdr algorithm
            d_mvdr = d_dasb.';
            Phi_NN = squeeze(P_NN(i,:,:));
            w_mvdr = ((Phi_NN\d_mvdr)/((d_mvdr.'/Phi_NN)*d_mvdr));
            w_mvdr = w_mvdr/sum(abs(w_mvdr));
            output_stft(i,:,:) = squeeze(input_stft(i,:,:))*w_mvdr;
        end
        
        % plot output
        fig_out = figure;
        plotSpec(output_stft(:,:,1),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('Freq (kHz)');
        xlabel("Time");
        
        % output in time domain
        output_t = calc_ISTFT(output_stft, win, N_STFT, R_STFT, 'onesided');
        
        % write output signal
        audiowrite(sprintf("Temporary/toBeTested/out_MVDR/set%d_Recording (%d).flac", setNr, fileNr), output_t, fs_input);
    
        % save plots
        saveas(fig_in, sprintf("Temporary/figures/in_MVDR/set%d_Recording (%d).png", setNr, fileNr));
        saveas(fig_out,sprintf("Temporary/figures/out_MVDR/set%d_Recording (%d).png", setNr, fileNr));
        disp("Job done!");
    end
end

%% calculate SNR
clc;
setNr = 1;
fileNr = 1;
input_t = audioread(sprintf("Temporary/toBeTested/set%d_Recording (%d).flac", setNr, fileNr));
N_in = audioread(sprintf("Temporary/toBeTested/set%d_Noise (%d).flac", setNr, fileNr));

output_t = audioread(sprintf("Temporary/toBeTested/out_MVDR/set%d_Recording (%d).flac", setNr, fileNr));
N_out = audioread(sprintf("Temporary/toBeTested/out_MVDR/set%d_Recording (%d).flac", setNr, fileNr));
noise_in_psd = pwelch(N_in);
noise_out_psd = pwelch(N_out);
y_psd = pwelch(input_t);
x_hat_psd = pwelch(output_t);

signal_in_psd = mean(abs(y_psd)) - mean(abs(noise_in_psd));
input_SNR = 10*log10(sum(signal_in_psd)/sum(mean(abs(noise_in_psd))));
disp(["input SNR: ",input_SNR]);

signal_psd_out = mean(abs(x_hat_psd)) - mean(abs(noise_out_psd));
output_SNR = 10*log10(sum(signal_psd_out)/sum(mean(abs(noise_out_psd))));
disp(["output SNR: ",output_SNR]);
snr_improvement = output_SNR - input_SNR;
disp(["SNR improvement: ",snr_improvement]);
