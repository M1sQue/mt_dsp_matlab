clear;
load("MatData/w_mwf.mat");

% input time domain
[input_t, fs_input] = audioread("Temporary/MO202501-WQB4BVWP-20211117-023500-MULTICHANNEL.flac");
% STFT parameters
N_STFT = 512;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));
% spectogram figure settings
xTickProp = [0, R_STFT/fs_input, 0]; % R_STFT/fs_input, fs_input/R_STFT
yTickProp = [0, fs_input/(2000*R_STFT), R_STFT/2];
cRange    = [-45 15];
% plot input
input_stft = calc_STFT(input_t, fs_input, win, N_STFT, R_STFT, 'onesided');
figure;
iterations = numel(input_stft(1,1,:));
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(input_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
end
output_stft = zeros(numel(input_stft(:,1,1)), numel(input_stft(1,:,1)), flag);

% filter coefficients
f_start = 0;
f_end = 0;
frame_len = numel(input_stft(:,1,1));
iterations = numel(input_stft(1,:,1));
for i=1:numel(f_sound_group)
    f_start = f_end;
    if(i ~= numel(f_sound_group))
        f_end = (f_sound_group(1,i)+f_sound_group(1,i+1))/2;
    else
        f_end = fs_input/2;
    end
    i_start = floor(f_start*frame_len/(fs_input/2)) +1;
    i_end = floor(f_end*frame_len/(fs_input/2));
    w_filter = W_MWF{i};
    for j = 1:iterations
        output_stft(i_start:i_end,j,:) = squeeze(input_stft(i_start:i_end,j,:)) * w_filter.';
    end
end

output_t = calc_ISTFT(output_stft, win, N_STFT, R_STFT, 'onesided');

% plot output
iterations = flag;
figure;
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(output_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
end

audiowrite("Temporary\00_current_test_result.wav", output_t, fs_input);

%% calculate SNR
N_in = audioread("Temporary/SNR_-15_pure_noise.wav");
N_out = audioread("Temporary/result_noise_SNR_-15_2.wav");
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
