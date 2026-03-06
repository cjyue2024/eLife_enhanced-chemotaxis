% Export tracks to XML file in imageJ
% Import and analyse Trackmate data in Matlab
clear
files = dir('Data04*.xml');

for II = 1:length(files)
    
    file_path_tracks = files(II).name;
    tracks = importTrackMateTracks(file_path_tracks);  % T in frame, X, Y, Z
    n_tracks = numel( tracks );
    fprintf('Found %d tracks in the file.\n', n_tracks)
    
    clipZ = true;    % Remove Z coordinates, if you know you can.
    scaleT = true;   % Use physical time for T, in seconds.
    tracks = importTrackMateTracks(file_path_tracks, clipZ, scaleT);
    
    [ tracks, md ] = importTrackMateTracks(file_path_tracks,clipZ, scaleT);
    md 
    save([files(II).name(1:end-4),'.mat'],'tracks','md')
    clearvars -except files II
end