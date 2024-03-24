clear;
load("polars.mat");

f_sound = 1000; % sound frequency for beamforming
r = 0.1; % coordinate unit length
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
sys = db2mag(repmat(polars(temp,:), numel(m_pos(1,:)), 1)');

% sound source parameters definition
azimuth_deg = 0;
elevation_deg = -90;
azimuth = deg2rad(azimuth_deg);
elevation = deg2rad(elevation_deg);
s_pos = [cos(elevation)*cos(azimuth) cos(elevation)*sin(azimuth) sin(elevation)];

% delay and sum algorithm
dasb_delay = s_pos*m_pos/c;
dasb = exp(-1j*2*pi*f_sound*dasb_delay);

% simulation parameters
sound_delay_angles = deg2rad(0:5:360);
sound_delay_positions = [cos(sound_delay_angles')*cos(azimuth) cos(sound_delay_angles')*sin(azimuth) sin(sound_delay_angles')];
sound_delays = sound_delay_positions*m_pos/c;
simulations = dasb*(exp(1j*2*pi*f_sound*sound_delays).*sys)'/numel(m_pos(1,:));

% simulation polar plot
threshold = -30;
for i = 1:numel(sound_delay_angles)
    if 20*log10(abs(simulations(i))) < threshold
        simulations(i) = 10^(threshold/20);
    end
end
polarplot(sound_delay_angles, 20*log10(abs(simulations)));
thetalim([0 360]);
thetaticks(0:45:315);
rlim([threshold 10]);
rticks(threshold:5:10);
title(sprintf("azimuth %d°, elevation %d°, frequency %dHz", azimuth_deg, elevation_deg, f_sound));
