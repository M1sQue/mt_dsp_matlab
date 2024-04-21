clear;
load("MatData/w_mwf.mat");

% input time domain
[input_t, fs_input] = audioread("Temporary\SNR_-15.wav");
% STFT parameters
N_STFT = 512;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));
% spectogram figure settings
xTickProp = [0, R_STFT/fs_input, 0];
yTickProp = [0, fs_input/(2000*R_STFT), R_STFT/2]; % fs_input/R_STFT
cRange    = [-45 15];
% plot input
input_stft   = calc_STFT(input_t, fs_input, win, N_STFT, R_STFT, 'onesided');
figure;
iterations = numel(input_stft(1,1,:));
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(input_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
end

% input frequency domain
input_f = fft(input_t);
% output frequency domain
output_f = input_f*(w_mwf.');

% output time domain
output_t = abs(ifft(output_f));
% plot input
output_stft   = calc_STFT(output_t, fs_input, win, N_STFT, R_STFT, 'onesided');
figure;
iterations = numel(output_stft(1,1,:));
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(output_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
end

audiowrite("Temporary\00_current_test_result.wav", output_t, fs_input);
