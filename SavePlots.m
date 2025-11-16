% Get all open figure handles
figHandles = findall(0, 'Type', 'figure');
% Loop through each figure and save it as a JPEG
for i = 1:length(figHandles)
    figure(figHandles(i)); % Make the current figure active (optional, but good practice)
    filename = sprintf('figure_%d.jpg', i); % Create a unique filename
    saveas(figHandles(i), filename, 'jpg');
end