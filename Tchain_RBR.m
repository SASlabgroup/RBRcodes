% Bulk Read RBR Solo and Duo loggers on a T-chain
% 
% Michael James
% 06/27/2024
% 

clc, clear, close all
%% Manual inputs
% Please input your instrument information and estimated z locations in order by z, 
% and desired time start/stop. Z must be manually inputted since RBR solos
% do not have data on depth/position

t1 = 7.394247886921297e5; % time start
t2 = 7.394299614814814e5; % time end

% Input SN of RBR here IN ORDER of associated z below.
serial_list = [212832;212833;212834;212857;212858;212859;212860;212861;212862; ... 
    212863;212864;200488;212865;212866;212867;212868;212869;212870;...
    212871;212872;200492;212873;212874;212875;212876;212877;212878;...
    212879;212880; 212881 ;200495];

% Input Tchain lengths here. 
z = [0;0.25;0.5;0.75;1;1.25;1.5;1.75;2;...
    2.25;2.5;2.75;3.25;3.5;3.75;4;4.25;4.5;...
    4.75;5;5.25;5.75;6.25;6.75;7.25;7.75;8.25;8.75;9.25;9.75;10;];

% Distance from surface to first node. 
buoyoffset = 0.5; % This is an additional length to offset your first RBR in case you have slack between it and surface.
z = buoyoffset + z;



%% Loading loop

% Need RSK tools library from RBR

rawdatapath = [uigetdir()]; % Calls explorer to get directory.  
% rawdatapath = []; % Can also input raw dir here. 

files = dir(fullfile(rawdatapath, '*.rsk*')); % Use this to filter your files

for fi = 1:length(files)
    fid = RSKopen(fullfile(files(fi).folder, files(fi).name));
    rsk = RSKreaddata(fid,'t1',t1,'t2',t2);    
    idx = find(rsk.instruments.serialID == serial_list);
    Tchain.RSK{idx,1} = rsk;        
    Tchain.time(idx,:) = rsk.data.tstamp; %epoch
    Tchain.temp(idx,:) = rsk.data.values(:,1)'; %C
    try
        Tchain.P(idx,:) = rsk.data.values(:,2)'; %dbar
        disp('duo processed')
    catch 
        disp('solo processed')
        Tchain.P(idx,:) = nan(size(rsk.data.values(:,1)'));
    end; clear rsk;
end;clear fi fid idx;


%% Interpolate z to depth (Linear Fit)
if any(~isnan(Tchain.P))
    disp('Duos present, interpolating depth to P')
    % P1 Subtract absolute pressure as needed (Assuming we're in the top 10 m)
    [Pmin, r] = min(Tchain.P,[], "all"); % Calculate Min
    [r, ~] = find(Tchain.P == Tchain.P(r)); % Grab Row
    if abs(z(r)-Pmin) > 0.1*z(r) % examine offset % 
        Tchain.P  = Tchain.P -10; % Remove 1 bar atmospheric
    end
    
    %% P2 linear interp Assuming dbar = 1m depth
    [Prows ~] = find(~isnan(Tchain.P(:,1))) %Assume consistent
    depth = nan(size(Tchain.P));
    for i =1:length(Tchain.time)
        Stretch_ratio = Tchain.P(Prows,i) ./ z(Prows); % Linear Stretching between each P point
        % Stretch_ratio = Stretch_ratio(~isnan(Stretch_ratio));
        depth(1:Prows(1),i) = Stretch_ratio(1).*z(1:Prows(1));
        depth(Prows(1)+1:Prows(2),i) = Stretch_ratio(2).*z(Prows(1)+1:Prows(2));
        depth(Prows(2)+1:Prows(3),i) = Stretch_ratio(3).*z(Prows(2)+1:Prows(3));
    end;clear i;
else 
    disp('All Solos')
end;
    

%% Plot pcolor of Tchain
figure;

% z vs P to examine stretch
plot(mean(Tchain.P, 2), z, 'k.')
xlabel('Pressure [dbar]'); ylabel('Length Along tether unstretched [m]')
axis manual
l = [min(min(Tchain.P,[], "all"), min(z)), max(max(Tchain.P,[], "all"), max(z))
plot(l, l, 'k--')

figure;

limc = round([nanmean(Tchain.temp,'all') - 2*nanstd(Tchain.temp,[],'all'), ...
               nanmean(Tchain.temp,'all') + 2*nanstd(Tchain.temp,[],'all')]);


tiledlayout("vertical")
nexttile;
pcolor(Tchain.time, depth, Tchain.temp)
set(allchild(gca), 'EdgeAlpha', 0)
caxis(limc); colormap(cmocean('thermal'));ylabel(colorbar,'[deg C]'
ylabel('Pressure [dbar]'); xlabel('[UTC]'); datetick;axis tight;axis ij;

nexttile;
pcolor(Tchain.time, z, Tchain.temp)
set(allchild(gca), 'EdgeAlpha', 0)
caxis(limc); colormap(cmocean('thermal'));ylabel(colorbar,'[deg C]')
ylabel('Length Along tether unstretched [m]'); xlabel('[UTC]'); datetick;axis tight;axis ij;

linkaxes(findobj(allchild(gcf), 'Type', 'Axes'),'xy');

name = [rawdatapath , ' Tchain'];
sgtitle(name)

savefig(name)




