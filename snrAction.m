beamImprovements = zeros(4, 1);
avg = zeros(8, 1);
for setNr = 1:2
    if setNr == 1
        totalFileNr = 18;
    else
        totalFileNr = 8;
    end
    for fileNr = 1:totalFileNr
        input_t = audioread(sprintf("Temporary/toBeTested/set%d_Recording (%d).flac", setNr, fileNr));
        N_in = audioread(sprintf("Temporary/toBeTested/set%d_N (%d).flac", setNr, fileNr));
        
        beamformerTypes = ["DS" "DSC" "MVDR" "MWF"];
        for typeNr = 1:length(beamformerTypes)
            output_t = audioread(sprintf("Temporary/toBeTested/out_%s/set%d_Recording (%d).flac", beamformerTypes(typeNr), setNr, fileNr));
            N_out = audioread(sprintf("Temporary/toBeTested/out_%s/set%d_N (%d).flac", beamformerTypes(typeNr), setNr, fileNr));
            noise_in_psd = pwelch(N_in);
            noise_out_psd = pwelch(N_out);
            y_psd = pwelch(input_t);
            x_hat_psd = pwelch(output_t);
            
            signal_in_psd = mean(abs(y_psd)) - mean(abs(noise_in_psd));
            input_SNR = 10*log10(sum(signal_in_psd)/sum(mean(abs(noise_in_psd))));
            
            signal_psd_out = mean(abs(x_hat_psd)) - mean(abs(noise_out_psd));
            output_SNR = 10*log10(sum(signal_psd_out)/sum(mean(abs(noise_out_psd))));

            snr_improvement = output_SNR - input_SNR;
            beamImprovements(typeNr) = snr_improvement;
        end
        text = sprintf("& Rec.%d & %.2f & %.2f & %.2f & %.2f & %.2f \\", fileNr, input_SNR, beamImprovements(1), beamImprovements(2), beamImprovements(3), beamImprovements(4));
        disp(text);
        for j = 1+4*(setNr-1):4+4*(setNr-1)
            avg(j) = avg(j)+beamImprovements(j-4*(setNr-1));
        end
    end
    for j = 1+4*(setNr-1):4+4*(setNr-1)
        avg(j) = avg(j)/totalFileNr;
    end
end
%% plot averages
data_set = reshape(avg, 2, 4);
figure;
bar(1:size(data_set, 1), data_set(:, :), 'grouped');
title('Average SNR Improvement');
xlabel('Set number');
ylabel('SNR Improvement (dB)');
legend('DS', 'DS-C', 'MVDR', 'MWF');
