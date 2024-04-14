%% split measurements
for i = 1:6
    for j = 0:5:115
        fileName = sprintf("AnechoicRoomMeasurements/RecycleBin2/monKD8N255G_%d_%d.wav", j, i);
        [ss_data, fs] = audioread(fileName);
        ss_data_1 = ss_data((0*fs+1):1:17*fs);
        audiowrite(sprintf("AnechoicRoomMeasurements/Temp/monKD8N255G_%d_%d_seg.wav", j, i),ss_data_1, fs);
        ss_data_2 = ss_data((17*fs+1):1:34*fs);
        audiowrite(sprintf("AnechoicRoomMeasurements/Temp/monKD8N255G_%d_%d_seg.wav", (j+120), i),ss_data_2, fs);
        ss_data_3 = ss_data((34*fs+1):1:51*fs);
        audiowrite(sprintf("AnechoicRoomMeasurements/Temp/monKD8N255G_%d_%d_seg.wav", (j+240), i),ss_data_3, fs);
    end
end

%% calculate IRs with inverse sine sweeps
inverseSineSweep = audioread("SineSweep/inverse_zeineb.wav");
for i = 1:6
    for j = 0:5:355
        fileName = sprintf("AnechoicRoomMeasurements/Temp/monKD8N255G_%d_%d_seg.wav", j, i);
        [ss_seg_data, fs] = audioread(fileName);
        IR_temp = conv(ss_seg_data, inverseSineSweep);
        maxNum = max(IR_temp);
        IR = IR_temp/maxNum;
        audiowrite(sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Original/IR_monKD8N255G_%d_Channel_%d.wav", i, j, i), IR(15*fs:1:length(IR)), fs);
    end
end

%% apply Hanning windows to IRs
c = 343.3;
delay_diff = (sqrt(5)-1)/c;

for i = 1:6
    for j = 0:5:355
        fileName = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Original/IR_monKD8N255G_%d_Channel_%d.wav", i, j, i);
        [ir_data, fs] = audioread(fileName);
        [M, I] = max(ir_data);
        N_hann_half = round(fs*delay_diff);
        N_hann = 2*N_hann_half;
        w = hann(N_hann);
        N = length(ir_data);
        w_padd = [zeros(I-N_hann/2, 1);w;zeros(N-I-N_hann/2, 1)];
        ir_data_hann = ir_data.*w_padd;
        audiowrite(sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_monKD8N255G_%d_Channel_%d.wav", i, j, i), ir_data_hann, fs);
    end
end

%% fix rotation of monitor 99PJ
c = 343.3;
delay_diff = (sqrt(5)-1)/c;
fix = [3 3 3 -3 -3 -3];

for i = 1:6
    for j = 0:5:355
        fileName = sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Original/IR_mon99PJ_%d_Channel_%d.wav", i, j, i);
        [ir_data, fs] = audioread(fileName);
        [M, I] = max(ir_data);
        N_hann_half = round(fs*delay_diff);
        N_hann = 2*N_hann_half;
        w = hann(N_hann);
        N = length(ir_data);
        w_padd = [zeros(I-N_hann/2, 1);w;zeros(N-I-N_hann/2, 1)];
        ir_data_hann = ir_data.*w_padd;
        index = i + fix(i);
        audiowrite(sprintf("AnechoicRoomMeasurements/IRs_Channel_%d_Hann/IR_mon99PJ_%d_Channel_%d.wav", index, j, index), ir_data_hann, fs);
    end
end