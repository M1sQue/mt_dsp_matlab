%% show figures
clear;
addpath("OldScripts\");
polars_cell = cell(1,6);
for i = 1:6
    load(sprintf("MatData/polars_average_channel_%d.mat", i));
    polars_cell{i} = polars;
end

f_sound = 1000; % sound frequency for beamforming
r = 0.057; % coordinate unit length
c = 343.3; % speed of sound

% microphone system parameters definition
m_pos = [r 0 0; r/2 -sqrt(3)/2*r 0; -r/2 -sqrt(3)/2*r 0; -r 0 0; -r/2 sqrt(3)/2*r 0; r/2 sqrt(3)/2*r 0]';
temp = -1;
for i=1:numel(polar_freq)
    if polar_freq(i) == f_sound
        temp = i;
    end
end
if temp == -1
    error("Frequency %dHz is not found in the directivity pattern.", f_sound);
end
sys = [polars_cell{1}(temp, :); polars_cell{2}(temp, :); polars_cell{3}(temp, :); polars_cell{4}(temp, :); polars_cell{5}(temp, :); polars_cell{6}(temp, :)]';
sys = db2mag(sys);

% sound source parameters definition
azimuth_deg = 0;
azimuth = deg2rad(azimuth_deg);

% mwf algorithm
[N, fs_noiseOnly] = audioread("Temporary/00N03_MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
% [X, fs_target]= audioread("Temporary/00X01_SNR_15.wav");
[Y, fs_targetPlusNoise]= audioread("Temporary/00Y02_MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
n_mics = numel(m_pos(1,:));

% convert data to GPU arrays
N_gpu = gpuArray(N);
% X_gpu = gpuArray(X);
Y_gpu = gpuArray(Y);

P_NN = calculateCPSD(N_gpu, N_gpu, n_mics, fs_noiseOnly, f_sound);
% P_XX = calculateCPSD(X_gpu, X_gpu, n_mics, fs_target, f_sound);
P_YY = calculateCPSD(Y_gpu, Y_gpu, n_mics, fs_targetPlusNoise, f_sound);
% w_mwf = P_XX/(P_XX+P_NN); % enhance the target signal while minimizing the effect of noise
w_mwf = (P_YY-P_NN)/P_YY; % attempt to remove the noise component from the total signal-plus-noise matrix, focusing on enhancing the clarity of the signal
w_mwf = gather(w_mwf);
for i = 1:numel(w_mwf(:,1))
    w_mwf(i, :) = w_mwf(i, :)/sum(abs(w_mwf(i, :)));
end

% simulation parameters
sound_delay_angles = deg2rad(0:5:360);
sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
sound_delays = -sound_delay_positions*m_pos/c; % to simulate the propagation: with "-"
simulations = w_mwf*(exp(-1j*2*pi*f_sound*sound_delays)).'; % to simulate signals from all directions: with "-"; do NOT use Hermitian

% simulation polar plot
figure;
for i = 1:numel(simulations(:,1))
    subplot(2, 3, i);
    H_threshold = max(mag2db(abs(simulations(i,:))));
    L_threshold = H_threshold - 20;
    for j = 1:numel(sound_delay_angles)
        if mag2db(abs(simulations(i, j))) < L_threshold
            simulations(i, j) = db2mag(L_threshold);
        end
    end
    polarplot(sound_delay_angles, mag2db(abs(simulations(i, :))));
    thetalim([0 360]);
    thetaticks(0:45:315);
    rlim([L_threshold H_threshold]);
    rticks(L_threshold:5:H_threshold);
    title(sprintf("channel %d, frequency %dHz", i, f_sound));
end

%% only calculate coefficients XN
clear;
addpath("OldScripts\");

f_sound_group = 100:100:8000; % sound frequency for beamforming
W_MWF = cell(1, length(f_sound_group));

% mwf algorithm
[N, fs_noiseOnly] = audioread("Temporary/00N04_MO202501-N798NLFT-20220716-235500-MULTICHANNEL.flac");
[X, fs_target]= audioread("Temporary/00X02_MO202501-N798NLFT-20220716-235500-MULTICHANNEL.flac");
n_mics = numel(N(1, :));

% convert data to GPU arrays
N_gpu = gpuArray(N);
X_gpu = gpuArray(X);

for k = 1:length(f_sound_group)
    f_sound = f_sound_group(k);
    P_NN = calculateCPSD(N_gpu, N_gpu, n_mics, fs_noiseOnly, f_sound);
    P_XX = calculateCPSD(X_gpu, X_gpu, n_mics, fs_target, f_sound);
    w_mwf = P_XX/(P_XX+P_NN); % enhance the target signal while minimizing the effect of noise
    for i = 1:numel(w_mwf(:,1))
        w_mwf(i, :) = w_mwf(i, :)/sum(abs(w_mwf(i, :)));
    end 
    W_MWF{k} = gather(w_mwf);
end

flag = 6;
save("MatData/w_mwf_00N01_00X01.mat", "W_MWF", "f_sound_group", "flag");
disp("Job done");

%% only calculate coefficients YN
clear;
addpath("OldScripts\");

f_sound_group = 100:100:8000; % sound frequency for beamforming
W_MWF = cell(1, length(f_sound_group));

% mwf algorithm
[N, fs_noiseOnly] = audioread("Temporary/00N05_MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
[Y, fs_targetPlusNoise]= audioread("Temporary/00Y04_MO201701-9R88J6JG-20221013-235500-MULTICHANNEL.flac");
n_mics = numel(N(1, :));

% convert data to GPU arrays
N_gpu = gpuArray(N);
Y_gpu = gpuArray(Y);

for k = 1:length(f_sound_group)
    f_sound = f_sound_group(k);
    P_NN = calculateCPSD(N_gpu, N_gpu, n_mics, fs_noiseOnly, f_sound);
    P_YY = calculateCPSD(Y_gpu, Y_gpu, n_mics, fs_targetPlusNoise, f_sound);
    w_mwf = (P_YY-P_NN)/P_YY; % attempt to remove the noise component from the total signal-plus-noise matrix, focusing on enhancing the clarity of the signal
    for i = 1:numel(w_mwf(:,1))
        w_mwf(i, :) = w_mwf(i, :)/sum(abs(w_mwf(i, :)));
    end 
    W_MWF{k} = gather(w_mwf);
end

flag = 6;
save("MatData/w_mwf_00N04_00Y03.mat", "W_MWF", "f_sound_group", "flag");
disp("Job done");
