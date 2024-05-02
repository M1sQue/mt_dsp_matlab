%% compare plots
% STFT related parameters
clear;
addpath("FromThomasDietzen\")
N_STFT = 2048;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));

% load data
% mwf coefficients calculation data time domain
[N, fs_noiseOnly] = audioread("Temporary/00N05_MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
[Y, fs_targetPlusNoise]= audioread("Temporary/00Y04_MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
% mwf coefficients calculation data frequency domain
N_stft = calc_STFT(N, fs_noiseOnly, win, N_STFT, R_STFT, 'onesided');
Y_stft = calc_STFT(Y, fs_targetPlusNoise, win, N_STFT, R_STFT, 'onesided');

% input data time domain
[input_t, fs_input] = audioread("Temporary/toBeTested/MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
% input data frequency domain
input_stft = calc_STFT(input_t, fs_input, win, N_STFT, R_STFT, 'onesided');
xTickProp = [0, R_STFT/fs_input, 0]; % R_STFT/fs_input, fs_input/R_STFT
yTickProp = [0, fs_input/(2000*R_STFT), R_STFT/2];
cRange    = [-45 15];
% plot input
figure;
iterations = numel(input_stft(1,1,:));
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(input_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('Frequency (kHz)');
end
xlabel("Time");

% initialize the output
output_stft = zeros(size(input_stft));

% pipeline
% calculate psd
[P_NN_smth, P_NN_mean] = estim_corrmat(N_stft, 1);
[P_YY_smth, P_YY_mean] = estim_corrmat(Y_stft, 1);
P_NN = squeeze(P_NN_mean);
P_YY = squeeze(P_YY_mean);
n_freq_bins = numel(input_stft(:,1,1));

% calculate coefficients then apply
for i = 1:n_freq_bins
    w_mwf = (squeeze(P_YY(i,:,:))-squeeze(P_NN(i,:,:)))/squeeze(P_YY(i,:,:));
%     for j = 1:numel(w_mwf(:,1))
%         w_mwf(j, :) = w_mwf(j, :)/sum(abs(w_mwf(j, :)));
%     end
    output_stft(i,:,:) = squeeze(input_stft(i,:,:))*w_mwf;
end

% plot output
figure;
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(output_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('Frequency (kHz)');
end
xlabel("Time");

% output in time domain
output_t = calc_ISTFT(output_stft, win, N_STFT, R_STFT, 'onesided');

% write output signal
audiowrite("Temporary/00_current_test_result.wav", output_t, fs_input);

%% calculate SNR
N_in = audioread("Temporary/zz_current_noise_in.flac");
N_out = audioread("Temporary/zz_current_noise_out.wav");
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
