%% single directivity pattern
clear;

step = 5;
angles = 0:step:355;
num_angles = length(angles);
ir_data = cell(1, num_angles);
fft_data = cell(1, num_angles);

polar_freq = [500 1000 2000 4000 8000];
num_polar_freq = length(polar_freq);
polars = zeros(num_polar_freq, num_angles);
Zs = zeros(num_polar_freq, num_angles);

figure;
channel = 1;
hold on;
for i = 1:num_angles
    %filename = sprintf("Microphone_Impulse_Responses/ShureSM58_125cm_Normalised_IRs/IRs/ShureSM58_125cm_%dDeg.wav", angles(i));
    filename = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_monKD8N255G_%d_Channel_%d.wav", channel, angles(i), channel);
    [ir_data{i}, fs] = audioread(filename);
    N = length(ir_data{i});
    freq = (1:N)*(fs/N);
    freq = freq(1:N/2);
    fft_data{i} = fft(ir_data{i});
    Z = fft_data{1, i}(1:N/2);
    magnitude = mag2db(abs(Z));
    for j = 1:num_polar_freq
        index = 1+mod(i-1+270/step,360/step);
        polars(j, index) = magnitude(polar_freq(j)*N/fs);
        Zs(j, index) = Z(polar_freq(j)*N/fs);
    end
    plot(freq, magnitude);
end
hold off;
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");

title(sprintf("FFTs for Channel %d", channel));

figure;
angles = [angles 360];
polars = [polars polars(:,1)];
Zs = [Zs Zs(:,1)];
maxMag = max(polars(:));
polars = polars-maxMag*ones(numel(polars(:,1)), numel(polars(1,:)));
tbl = array2table([deg2rad(angles') polars']);
tbl = renamevars(tbl, "Var1", "Angles (rad)");
polar_freq_labels = string(polar_freq);
for i = 1:num_polar_freq
    polar_freq_labels(i) = sprintf("%dHz",polar_freq(i));
    tbl = renamevars(tbl, sprintf("Var%d",i+1), sprintf("%dHz",polar_freq(i)));
end

polarplot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
legend;
thetalim([0 360]);
thetaticks(0:45:315);
rlim([-60 0]);
rticks(-60:5:0);
title(sprintf("Directivity Pattern for Channel %d", channel));

save("MatData/polars.mat", "Zs", "polar_freq");

%% mean directivity pattern
clear;

step = 5;
angles = 0:step:355;
num_angles = length(angles);
ir_data = cell(1, num_angles);
fft_data = cell(1, num_angles);

polar_freq = [500 1000 2000 4000 8000];
num_polar_freq = length(polar_freq);
polars = zeros(num_polar_freq, num_angles);
Zs = zeros(num_polar_freq, num_angles);

monName = ["monCQGLL74L" "monG7SMCCBW" "monKD8N255G" "monGHHX" "mon99PJ"];
channel = 6;
for k = 1:5
    for i = 1:num_angles
        %filename = sprintf("Microphone_Impulse_Responses/ShureSM58_125cm_Normalised_IRs/IRs/ShureSM58_125cm_%dDeg.wav", angles(i));
        filename = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_%s_%d_Channel_%d.wav", channel, monName(k), angles(i), channel);
        [ir_data{i}, fs] = audioread(filename);
        N = length(ir_data{i});
        fft_data{i} = fft(ir_data{i});
        Z = fft_data{1, i}(1:N/2);
        magnitude = abs(Z);
        for j = 1:num_polar_freq
            index = 1+mod(i-1+270/step,360/step);
            polars(j, index) = polars(j, index)+magnitude(polar_freq(j)*N/fs);
            Zs(j, index) = Zs(j, index)+Z(polar_freq(j)*N/fs);
        end
    end
end
polars = polars/5;
Zs = Zs/5;
polars = mag2db(polars);
maxMag = max(polars(:));
polars = polars-maxMag*ones(numel(polars(:,1)), numel(polars(1,:)));

figure;
angles = [angles 360];
polars = [polars polars(:,1)];
Zs = [Zs Zs(:,1)];
tbl = array2table([deg2rad(angles') polars']);
tbl = renamevars(tbl, "Var1", "Angles (rad)");
polar_freq_labels = string(polar_freq);
for i = 1:num_polar_freq
    polar_freq_labels(i) = sprintf("%dHz",polar_freq(i));
    tbl = renamevars(tbl, sprintf("Var%d",i+1), sprintf("%dHz",polar_freq(i)));
end

polarplot(tbl, "Angles (rad)", polar_freq_labels, 'Linewidth', 1);
legend;
thetalim([0 360]);
thetaticks(0:45:315);
rlim([-60 0]);
rticks(-60:5:0);
title(sprintf("Average Directivity Pattern for Channel %d", channel));

save(sprintf("MatData/polars_average_channel_%d", channel), "Zs", "polar_freq");
