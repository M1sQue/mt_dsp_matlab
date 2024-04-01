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