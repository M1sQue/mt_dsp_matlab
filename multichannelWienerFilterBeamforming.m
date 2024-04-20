clear;
polars_cell = cell(1,6);
for i = 1:6
    load(sprintf("MatData/polars_average_channel_%d.mat", i));
    polars_cell{i} = polars;
end

f_sound = 4000; % sound frequency for beamforming
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
elevation_deg = -90;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = 50*r*[cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

% mwf algorithm
[noiseOnly, fs_noiseOnly] = audioread("Temporary/SNR_-15_pure_noise.wav", [1 200000]);
[targetPlusNoise, fs_targetPlusNoise]= audioread("Temporary/SNR_-15.wav", [1 200000]);
n_mics = numel(m_pos(1,:));
P_NN = ones(n_mics, n_mics);
P_XX = ones(n_mics, n_mics);
for i = 1:n_mics
    for j = i:n_mics  % Symmetric matrix, compute half and mirror
        [cpsd_ij, f] = cpsd(noiseOnly(:,i), noiseOnly(:,j), [], [], [], fs_noiseOnly);
        f_index = compareIndex(f_sound,fs_noiseOnly,f);
        P_NN(i,j) = cpsd_ij(f_index);
        [cpsd_ji, f] = cpsd(noiseOnly(:,j), noiseOnly(:,i), [], [], [], fs_noiseOnly);
        f_index = compareIndex(f_sound,fs_noiseOnly,f);
        P_NN(j,i) = cpsd_ji(f_index);
    end
end
for i = 1:n_mics
    for j = i:n_mics  % Symmetric matrix, compute half and mirror
        [cpsd_ij, f] = cpsd(targetPlusNoise(:,i), targetPlusNoise(:,j), [], [], [], fs_targetPlusNoise);
        f_index = compareIndex(f_sound,fs_targetPlusNoise,f);
        P_XX(i,j) = cpsd_ij(f_index);
        [cpsd_ji, f] = cpsd(targetPlusNoise(:,j), targetPlusNoise(:,i), [], [], [], fs_targetPlusNoise);
        f_index = compareIndex(f_sound,fs_targetPlusNoise,f);
        P_XX(j,i) = cpsd_ji(f_index);
    end
end
w_mwf = P_XX/(P_XX+P_NN);
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
    H_threshold = 0;
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
