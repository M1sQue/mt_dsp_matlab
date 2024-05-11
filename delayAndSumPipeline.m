%% Pipeline for Delay and Sum beamformer
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

% mean directivity pattern for DASB weights compensation

step = 5;
angles = 0:step:355;
num_angles = length(angles);
ir_data = cell(1, num_angles);
fft_data = cell(1, num_angles);

num_polar_freq = N_STFT/2 +1;
polars = zeros(num_polar_freq, num_angles,6);

monName = ["monCQGLL74L" "monG7SMCCBW" "monKD8N255G" "monGHHX" "mon99PJ"];
for channel = 1:6
    for k = 1:5
        for i = 1:num_angles
            %filename = sprintf("Microphone_Impulse_Responses/ShureSM58_125cm_Normalised_IRs/IRs/ShureSM58_125cm_%dDeg.wav", angles(i));
            filename = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_%s_%d_Channel_%d.wav", channel, monName(k), angles(i), channel);
            [ir_data{i}, fs] = audioread(filename);
            N = length(ir_data{i});
            N_FFT = 2^nextpow2(N);
            fft_data{i} = fft(ir_data{i},N_FFT);
            Z = fft_data{i};
            magnitude = abs(Z);
            index_deg = 1+mod(i-1+270/step,360/step);
            index_freq = N_FFT*linspace(0,1/2,num_polar_freq)+1;       
            polars(:,index_deg,channel) = polars(:,index_deg,channel) + magnitude(index_freq);
        end
    end
end

polars = polars/5;
polars = mag2db(polars);
maxMag = max(polars(:));
polars = polars-maxMag*ones(size(polars));
flag = 1;
save("MatData/polars.mat", "polars");
disp("Job done");
%%
addpath("FromThomasDietzen\")
load("MatData/polars.mat");

%compensation coefficient for target direction range:
%elevation_deg-compen_width_deg/2:elevation_deg+compen_width_deg/2
compen_width_deg = 40;
compen_index = mod((elevation_deg-compen_width_deg/2:step:elevation_deg+compen_width_deg/2)+360,360)/step +1;
A_compen = -mean(polars(:,compen_index,:),2);
%discard compensation for f<200Hz(sinesweep not covered) 
A_compen(1:ceil(200*N_STFT/fs),:) =  ones(ceil(200*N_STFT/fs),numel(A_compen(1,:)));
%discard compensation for f>8000Hz(not accurate) 
A_compen(ceil(10000*N_STFT/fs):end,:) =  ones(size(A_compen(ceil(10000*N_STFT/fs):end,:)));
A_compen = db2mag(A_compen);

% input data time domain
setNr = 1;
if setNr == 1
    totalFileNr = 10;
else
    totalFileNr = 8;
end
for fileNr = 1:totalFileNr
    [input_t, fs_input] = audioread(sprintf("Temporary/toBeTested/set%d/Recording (%d).flac", setNr, fileNr));
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
    n_freq_bins = numel(input_stft(:,1,1));
    
    % calculate coefficients then apply
    for i = 1:n_freq_bins
        % delay and sum algorithm
        dasb_delay = s_pos*m_pos/norm(s_pos)/c; % to compensate the delay aka alignment: times "-" to a "-"
        d_dasb = exp(-1j*2*pi*(fs_input/N_STFT*n_freq_bins)*dasb_delay)/numel(m_pos(1,:));
        w_dasb =(d_dasb.*A_compen(i,:)).';
%         w_dasb = d_dasb.';
        %normalization
        w_dasb = w_dasb/sum(abs(w_dasb));
        output_stft(i,:,:) = squeeze(input_stft(i,:,:)) * w_dasb;
    end
    
    % plot output
    fig_out = figure;
    plotSpec(output_stft(:,:,1),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('Freq (kHz)');
    xlabel("Time");
    
    % output in time domain
    output_t = calc_ISTFT(output_stft, win, N_STFT, R_STFT, 'onesided');
    
    % write output signal
    audiowrite(sprintf("Temporary/toBeTested/set%d_out_DS/Recording (%d).flac", setNr, fileNr), output_t, fs_input);
    
    % save plots
    saveas(fig_in, sprintf("Temporary/Figures/set%d/DS_in_rec%d.svg", setNr, fileNr));
    saveas(fig_out,sprintf("Temporary/Figures/set%d/DS_out_rec%d.svg", setNr, fileNr));
    disp("Job done!");
end

%% calculate SNR
clc;
setNr = 1;
fileNr = 1;
input_t = audioread(sprintf("Temporary/toBeTested/set%d/Recording (%d).flac", fileNr));
N_in = audioread(sprintf("Temporary/toBeTested/set%d/Noise (%d).flac", fileNr));

output_t = audioread(sprintf("Temporary/toBeTested/set%d_out_DS/Recording (%d).flac", setNr, fileNr));
N_out = audioread(sprintf("Temporary/toBeTested/set%d_out_DS/Noise (%d).flac", setNr, fileNr));
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
