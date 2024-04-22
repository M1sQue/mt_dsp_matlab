clear;
load("MatData/w_mwf.mat");

% input time domain
[input_t, fs_input] = audioread("Temporary/SNR_-15.wav");
% STFT parameters
N_STFT = 512;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));
% spectogram figure settings
xTickProp = [0, R_STFT/fs_input, 0];
yTickProp = [0, fs_input/(2000*R_STFT), R_STFT/2]; % fs_input/R_STFT
cRange    = [-45 15];
% plot input
input_stft = calc_STFT(input_t, fs_input, win, N_STFT, R_STFT, 'onesided');
figure;
iterations = numel(input_stft(1,1,:));
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(input_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
end

% % input frequency domain
% input_f = fft(input_t);
% input_len = numel(input_f(:,1));
% f_start = 0;
% f_end = 0;
% output_f = zeros(size(input_f));
% for i=1:numel(f_sound_group)
%     f_start = f_end;
%     if(i ~= numel(f_sound_group))
%         f_end = (f_sound_group(1,i)+f_sound_group(1,i+1))/2;
%     else
%         f_end = fs_input/2;
%     end
%     i_start = round(f_start*input_len/fs_input) +1;
%     i_end = round(f_end*input_len/fs_input);
%     w_mwf = W_MWF{i};
%     output_f(i_start:i_end,:) = input_f(i_start:i_end,:)*w_mwf.';
% end
% 
% output_f = output_f(1:input_len/2,:);
% 
% % output frequency domain
% % output_f = input_f*(w_mwf.');
% 
% % output time domain
% % output_t = ifft(output_f);
% output_t = abs(ifft(output_f));
% % plot input
% output_stft = calc_STFT(output_t, fs_input, win, N_STFT, R_STFT, 'onesided');

output_stft = zeros(size(input_stft));
iterations = numel(input_stft(1,:,1));
for i = 1:iterations
    output_stft(:,i,:) = squeeze(input_stft(:,i,:)) * W_MWF{6}.';
end
output_t = calc_ISTFT(output_stft, win, N_STFT, R_STFT, 'onesided');

% plot output
iterations = numel(input_stft(1,1,:));
figure;
for i = 1:iterations
    subplot (iterations, 1, i);
    plotSpec(output_stft(:,:,i),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
end

audiowrite("Temporary\00_current_test_result.wav", output_t, fs_input);
