% To create location map for 4 2x8 channels (from left to right, from word side)
% eg. if electrode face anterior, x+ is from left ear to midline, y+ is from anterior to posterior
% For EEG: x is towards the nose, y is towards the left ear, and z towards the vertex.

% cortical 2x8 electrode
coords.C16.x = [0 200 400 600 800 1000 1200 1400 0 200 400 600 800 1000 1200 1400]/1000; % um/1000 is mm
coords.C16.y = [0 0 0 0 0 0 0 0 200 200 200 200 200 200 200 200]/1000;
% Thalamic 2x8 electrode
coords.T16.x = [0 250 500 750 1000 1250 1500 1750 0 250 500 750 1000 1250 1500 1750]/1000; % um
coords.T16.y = [0 0 0 0 0 0 0 0 250 250 250 250 250 250 250 250]/1000;
% Thalamic circular optrode
coords.O16.x = [-250 -176.8 0 176.8 500 353.6 0 -353.6 -176.8 0 176.8 250 353.6 0 -353.6 -500]/1000; % um
coords.O16.y = [0 176.8 250 176.8 0 -353.6 -500 -353.6 -176.8 -250 -176.8 0 353.6 500 353.6 0]/1000;
coords.ID16  = [5 6 7 8 9 10 11 12 4 3 2 1 16 15 14 13];

animalCode  = '0201';
locs = table();
switch animalCode
    case '0201'
        regionNames = {'PMC','CLA','PPC'}; % facing: PPAA, + for facing posterior
        regionDistanceX  = [-1, -1, -0.3]; % +: right ear, 0:midline
        regionDistanceY  = [1, 0, -0.5]; % +: nose
        locs.channum = [1:48]';
        locs.x = [coords.C16.x' + regionDistanceX(1);...
             coords.O16.x' + regionDistanceX(2);...
             coords.C16.x' + regionDistanceX(3)];
        locs.y = [coords.C16.y' + regionDistanceY(1);...            
             coords.O16.y' + regionDistanceY(2);...
             coords.C16.y' + regionDistanceY(3)];
        locs.z = ones(48,1)*0.95; % to be able to show up on the surface
    case '0171'
        regionNames = {'PFC','LPl','PPC','VC'}; % facing: PPAA, + for facing posterior
        regionDistanceX  = [-1, 0.6, -0.3, -0.3]; % +: right ear
        regionDistanceY  = [1, 0, -0.2, -0.8]; % +: nose
        locs.channum = [1:64]';
        locs.x = [coords.C16.x' + regionDistanceX(1);...
             coords.O16.x' + regionDistanceX(2);...
            -coords.C16.x' + regionDistanceX(3);...
            -coords.C16.x' + regionDistanceX(4)];
        locs.y = [coords.C16.y' + regionDistanceY(1);...            
             coords.O16.y' + regionDistanceY(2);...
            -coords.C16.y' + regionDistanceY(3);...
            -coords.C16.y' + regionDistanceY(4)];
        locs.z = ones(64,1)*0.95; % to be able to show up on the surface
end
% add region labels
for i = 1:numel(regionNames)
    regionName = regionNames{i};
    for j = 1:16
        locs.label(16*(i-1)+j) = {[regionName num2str(coords.ID16(j))]};
    end
end

locsCell = table2cell(locs);
save(['ChnMap_' animalCode],'locs','locsCell','coords');
%writecell(locsCell,'ChnMap_0171.txt','Delimiter',' ');
%writetable(locs,'ChnMap_0171.txt','Delimiter',' ');

%% To plot channel location:
% copy and paste variable "locs" content into ChnMap_0171.xyz (open as txt), don't
% include title
loc = readlocs('ChnMap_0201.xyz','filetype','xyz');
figure;topoplot([],loc); % without label
figure;topoplot([],loc,'electrodes','ptslabels'); % with label
