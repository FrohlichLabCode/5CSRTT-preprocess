function fig = AH_figure_eeg(freq, spec, regionNames)
fig = figure;
for iRegion = 1 : length(regionNames)
    subplot(2,2,iRegion);
    plot(freq,spec(iRegion*16-15:iRegion*16,:)');
    title(regionNames{iRegion});
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    xlim([0,150]);
end
end