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