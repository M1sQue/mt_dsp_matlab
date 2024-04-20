%% figure
clear;
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
elevation_deg = -90;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = 50*r*[cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

% delay and sum algorithm
dasb_delay = s_pos*m_pos/norm(s_pos)/c; % to compensate the delay aka alignment: times "-" to a "-"
d_dasb = exp(-1j*2*pi*f_sound*dasb_delay)/numel(m_pos(1,:));
w_dasb = d_dasb;

% mvdr algorithm
d_mvdr = d_dasb.';
n_mics = numel(m_pos(1,:));
Phi_NN = ones(n_mics, n_mics);
[noise, fs] = audioread("Temporary/SNR_-15_pure_noise.wav", [1 100000]);
for i = 1:n_mics
    for j = i:n_mics
        [cpsd_ij, f] = cpsd(noise(:,i), noise(:,j), [], [], [], fs);
        f_index = compareIndex(f_sound,fs,f);
        Phi_NN(i,j) = cpsd_ij(f_index);
        [cpsd_ji, f] = cpsd(noise(:,j), noise(:,i), [], [], [], fs);
        f_index = compareIndex(f_sound,fs,f);
        Phi_NN(j,i) = cpsd_ji(f_index);
    end
end
Phi_NN_inv = inv(Phi_NN);
w_mvdr = (Phi_NN_inv*d_mvdr/(d_mvdr'*Phi_NN_inv*d_mvdr)).';
w_mvdr = w_mvdr/sum(abs(w_mvdr));

% simulation parameters
sound_delay_angles = deg2rad(0:5:360);
sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
sound_delays = -sound_delay_positions*m_pos/c; % to simulate the propagation: with "-"
% simulations = w_mvdr*(exp(-1j*2*pi*f_sound*sound_delays).*sys).'; % to simulate signals from all directions: with "-"; do NOT use Hermitian
simulations = w_mvdr*(exp(-1j*2*pi*f_sound*sound_delays)).'; % to simulate signals from all directions: with "-"; do NOT use Hermitian

% simulation polar plot
H_threshold = max(mag2db(abs(simulations)));
L_threshold = H_threshold - 20;
for i = 1:numel(sound_delay_angles)
    if mag2db(abs(simulations(i))) < L_threshold
        simulations(i) = db2mag(L_threshold);
    end
end
polarplot(sound_delay_angles, mag2db(abs(simulations)));
thetalim([0 360]);
thetaticks(0:45:315);
rlim([L_threshold H_threshold]);
rticks(L_threshold:5:H_threshold);
title(sprintf("azimuth %d°, elevation %d°, frequency %dHz", azimuth_deg, elevation_deg, f_sound));

%% animation
clear;
% Load microphone array data
polars_cell = cell(1,6);
for i = 1:6
    load(sprintf("MatData/polars_average_channel_%d.mat", i));
    polars_cell{i} = polars;
end

% Sound and system parameters
f_sound = 1000; % sound frequency for beamforming
r = 0.057; % coordinate unit length
c = 343.3; % speed of sound

% Microphone positions
m_pos = [r 0 0; r/2 -sqrt(3)/2*r 0; -r/2 -sqrt(3)/2*r 0; -r 0 0; -r/2 sqrt(3)/2*r 0; r/2 sqrt(3)/2*r 0]';
% m_pos = [r 0 0; 2*r 0 0; 3*r 0 0; 4*r 0 0; 5*r 0 0; 6*r 0 0]';

% Find frequency index
temp = -1;
for i=1:numel(polars_cell{1}(:,1))
    if polar_freq(i) == f_sound
        temp = i;
        break;
    end
end
if temp == -1
    error("Frequency %dHz is not found in the directivity pattern.", f_sound);
end
sys = [polars_cell{1}(temp, :); polars_cell{2}(temp, :); polars_cell{3}(temp, :); polars_cell{4}(temp, :); polars_cell{5}(temp, :); polars_cell{6}(temp, :)]';
sys = db2mag(sys);

elevation_range = -180:10:180; % Range of elevation angles
% Initialize video writer
% v = VideoWriter('polar_animation.avi');
% open(v);

% Animation loop
for elevation_deg = elevation_range
    azimuth_deg = 0;
    azimuth = deg2rad(azimuth_deg);
    elevation = deg2rad(elevation_deg);
    s_pos = 50*r*[cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

    % Delay and sum beamforming
    dasb_delay = s_pos*m_pos/norm(s_pos)/c;
    d_dasb = exp(-1j*2*pi*f_sound*dasb_delay)/numel(m_pos(1,:));
    w_dasb = d_dasb;

    % MVDR beamforming
    d_mvdr = d_dasb.';
    n_mics = numel(m_pos(1,:));
    Phi_NN = ones(n_mics, n_mics);
    [noise, fs] = audioread("Temporary/SNR_-15_pure_noise.wav", [1 100000]);
    for i = 1:n_mics
        for j = i:n_mics
            [cpsd_ij, f] = cpsd(noise(:,i), noise(:,j), [], [], [], fs);
            f_index = compareIndex(f_sound,fs,f);
            Phi_NN(i,j) = cpsd_ij(f_index);
            [cpsd_ji, f] = cpsd(noise(:,j), noise(:,i), [], [], [], fs);
            f_index = compareIndex(f_sound,fs,f);
            Phi_NN(j,i) = cpsd_ji(f_index);
        end
    end
    Phi_NN_inv = inv(Phi_NN);
    w_mvdr = (Phi_NN_inv*d_mvdr/(d_mvdr'*Phi_NN_inv*d_mvdr)).';
    w_mvdr = w_mvdr/sum(abs(w_mvdr));

    % Simulate sound from all directions
    sound_delay_angles = deg2rad(0:5:360);
    sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
    sound_delays = -sound_delay_positions*m_pos/c;
    % simulations = w_mvdr*(exp(-1j*2*pi*f_sound*sound_delays).*sys).';
    simulations = w_mvdr*(exp(-1j*2*pi*f_sound*sound_delays)).';

    % Plotting
    H_threshold = max(mag2db(abs(simulations)));
    L_threshold = H_threshold - 10;
    for i = 1:numel(sound_delay_angles)
        if mag2db(abs(simulations(i))) < L_threshold
            simulations(i) = db2mag(L_threshold);
        end
    end
    steering_direction = L_threshold*ones(1,numel(sound_delay_angles));
    if elevation_deg <= 0   
        steering_index = mod(round(elevation_deg/5)+72,73)+1;
    else 
        steering_index = mod(round(elevation_deg/5)+73,73)+1;
    end
    steering_direction(steering_index) = 0; %direction of elevation_deg

    polarplot(sound_delay_angles, mag2db(abs(simulations)), 'LineWidth', 2);
    hold on;
    polarplot(sound_delay_angles, steering_direction, 'LineWidth', 2);
    hold off;
    thetalim([0 360]);
    thetaticks(0:45:315);
    rlim([L_threshold H_threshold]);
    rticks(L_threshold:10:H_threshold);
    title(sprintf("Azimuth %d°, Elevation %d°, Frequency %dHz", azimuth_deg, elevation_deg, f_sound));
    drawnow;

    % Capture the frame
    % frame = getframe(gcf);
    % writeVideo(v, frame);
end

% Close the video file
% close(v);
