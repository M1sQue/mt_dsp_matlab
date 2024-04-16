%% single frequency display
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
elevation_deg = -45;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = 50*r*[cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

% delay and sum algorithm
dasb_delay = s_pos*m_pos/norm(s_pos)/c; % to compensate the delay aka alignment: times "-" to a "-"
d_dasb = exp(-1j*2*pi*f_sound*dasb_delay)/numel(m_pos(1,:));
w_dasb = d_dasb;

% simulation parameters
sound_delay_angles = deg2rad(0:5:360);
sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
sound_delays = -sound_delay_positions*m_pos/c; % to simulate the propagation: with "-"
simulations = w_dasb*(exp(-1j*2*pi*f_sound*sound_delays).*sys).'; % to simulate signals from all directions: with "-"; do NOT use Hermitian

% simulation polar plot
threshold = -60;
for i = 1:numel(sound_delay_angles)
    if mag2db(abs(simulations(i))) < threshold
        simulations(i) = 10^(threshold/20);
    end
end
polarplot(sound_delay_angles, mag2db(abs(simulations)));
thetalim([0 360]);
thetaticks(0:45:315);
rlim([threshold 0]);
rticks(threshold:5:0);
title(sprintf("azimuth %d°, elevation %d°, frequency %dHz", azimuth_deg, elevation_deg, f_sound));

%% multiple frequencies display
clear;
polars_cell = cell(1,6);
for i = 1:6
    load(sprintf("MatData/polars_average_channel_%d.mat", i));
    polars_cell{i} = polars;
end

f_sound = polar_freq; % sound frequency for beamforming
r = 0.057; % coordinate unit length
c = 343.3; % speed of sound

% microphone system parameters definition
m_pos = [r 0 0; r/2 -sqrt(3)/2*r 0; -r/2 -sqrt(3)/2*r 0; -r 0 0; -r/2 sqrt(3)/2*r 0; r/2 sqrt(3)/2*r 0]';

% directivity patterns of target frequencies
numFreq = numel(polar_freq);
sys = cell(1,numFreq);

for i = 1:numFreq
    sys{i} = db2mag([polars_cell{1}(i, :); polars_cell{2}(i, :); polars_cell{3}(i, :); polars_cell{4}(i, :); polars_cell{5}(i, :); polars_cell{6}(i, :)]');
end

% sound source parameters definition
azimuth_deg = 0;
elevation_deg = -90;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = 50*r*[cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

% delay and sum algorithm
dasb_delay = s_pos*m_pos/norm(s_pos)/c;
d_dasb = exp(-1j*2*pi*f_sound'*dasb_delay)/numel(m_pos(1,:));
w_dasb = d_dasb;

% simulation parameters
sound_delay_angles = deg2rad(0:5:360);
sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
sound_delays = -sound_delay_positions*m_pos/c;
simulations = zeros(numFreq,numel(sound_delay_angles));
for i = 1:numFreq
    simulations(i,:) = w_dasb(i,:)*(exp(-1j*2*pi*f_sound(i)*sound_delays).*sys{i}).';
end

% simulation polar plot
threshold = -60;
simulations_dB = 20 * log10(abs(simulations));
under_threshold_indices = simulations_dB < threshold;
simulations_dB(under_threshold_indices) = threshold;

tbl = array2table([sound_delay_angles' simulations_dB']);
tbl = renamevars(tbl, "Var1", "Angles (rad)");
polar_freq_labels = string(polar_freq);
for i = 1:numFreq
    polar_freq_labels(i) = sprintf("%dHz",polar_freq(i));
    tbl = renamevars(tbl, sprintf("Var%d",i+1), sprintf("%dHz",polar_freq(i)));
end

polarplot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
legend;

thetalim([0 360]);
thetaticks(0:45:315);
rlim([threshold 0]);
rticks(threshold:5:0);
title(sprintf("Delay and Sum with 6 channel avg directivity patterns \n azimuth %d°, elevation %d°", azimuth_deg, elevation_deg));
