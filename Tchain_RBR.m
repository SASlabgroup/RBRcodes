% Bulk Read RBR Solo and Duo loggers on a T-chain
% 
% Michael James, 06/27/2024
% J. Thomson, modify for single pressure, 11/2025
% 

clc, clear, close all
%% Manual inputs
% Please input your instrument information and estimated z locations in order by z, 
% and desired time start/stop. Z must be manually inputted since RBR solos
% do not have data on depth/position

atmospres = 10;  % dB

t1 = datenum(2025,10,28,7,0,0); % time start
t2 = datenum(2025,11,5,18,0,0); % time end

% Input SN of RBR here IN ORDER of associated z below.
serial_list = [212832;212833;212834;212857;212858;212859;212860;212861;212862;240148; ... 
    212863;212864;212865;212866;212867;212868;212869;212870;...
    212871;212872;212873;];

% Input Tchain lengths here. 
z = [0:2:40];

% Distance from surface to first node. 
buoyoffset = 0; % This is an additional length to offset your first RBR in case you have slack between it and surface.
z = buoyoffset + z;

% force the last sensor to be on the bottom (depth same as anchor)
anchordepth = 20.5;


%% Loading loop

% Need RSK tools library from RBR

%rawdatapath = [uigetdir()]; % Calls explorer to get directory.  
rawdatapath = pwd % Can also input raw dir here. 

files = dir(fullfile(rawdatapath, '/*.rsk*')); % Use this to filter your files

for fi = 1:length(files)
    fid = RSKopen(fullfile(files(fi).folder, files(fi).name));
    rsk = RSKreaddata(fid,'t1',t1,'t2',t2);    
    idx = find(rsk.instruments.serialID == serial_list);
    Tchain.RSK{idx,1} = rsk;        
    Tchain.time(idx,:) = rsk.data.tstamp; %epoch
    Tchain.temp(idx,:) = rsk.data.values(:,1)'; %C
    try
        Tchain.P(idx,:) = rsk.data.values(:,2)'; %dbar
        disp('pressure data found and processed')
    catch 
        disp('solo temperature only processed')
        Tchain.P(idx,:) = nan(size(rsk.data.values(:,1)'));
    end; clear rsk;
end;clear fi fid idx;

% force last sensor to be on the bottom
Tchain.P(size(Tchain.P,1),:) = anchordepth + atmospres; 


%% Interpolate z to depth (Linear Fit)
if any(~isnan(Tchain.P))
    disp('Duos / Duets present, interpolating depth to P')
    % P1 Subtract absolute pressure as needed (Assuming we're in the top 10 m)
    [Pmin, r] = min(Tchain.P,[], "all"); % Calculate Min
    [r, ~] = find(Tchain.P == Tchain.P(r)); % Grab Row
    if abs(z(r)-Pmin) > 0.1*z(r) % examine offset % 
        Tchain.P  = Tchain.P - atmospres; % Remove  atmospheric
    end
    
    %% P2 linear interp Assuming dbar = 1m depth
    % edit this block for variable number of pressure sensors 

    [Prows ~] = find(~isnan(Tchain.P(:,1))) %Assume consistent
    depth = nan(size(Tchain.P));
    for i =1:length(Tchain.time)
        Stretch_ratio = Tchain.P(Prows,i) ./ z(Prows); % Linear Stretching between each P point
        % Stretch_ratio = Stretch_ratio(~isnan(Stretch_ratio));
        depth(1:Prows(1),i) = Stretch_ratio(1).*z(1:Prows(1));
        %depth((Prows(1)+1):length(z),i) = Stretch_ratio(2).*z((Prows(1)+1):end);
        depth(Prows(1)+1:Prows(2),i) = Stretch_ratio(2).*z(Prows(1)+1:Prows(2));
    end;clear i;
else 
    disp('All Solos')
end;
    

%% Plot pcolor of Tchain


% z vs P to examine stretch
figure(1)
plot(mean(Tchain.P, 2), z, 'k.')
xlabel('Pressure [dbar]'); ylabel('Length Along tether unstretched [m]')
axis manual
l = [min(min(Tchain.P,[], "all"), min(z)), max(max(Tchain.P,[], "all"), max(z))];
plot(l, l, 'k--')



limc = round([nanmean(Tchain.temp,'all') - 2*nanstd(Tchain.temp,[],'all'), ...
               nanmean(Tchain.temp,'all') + 2*nanstd(Tchain.temp,[],'all')]);

figure(2)
pcolor(Tchain.time, depth, Tchain.temp)
set(allchild(gca), 'EdgeAlpha', 0)
caxis(limc); colormap(cmocean('thermal')); ylabel(colorbar,'[deg C]')
ylabel('Pressure [dbar]'); xlabel('[UTC]'); datetick;axis tight;axis ij;

figure(3)
pcolor(Tchain.time, z, Tchain.temp)
set(allchild(gca), 'EdgeAlpha', 0)
caxis(limc); colormap(cmocean('thermal'));ylabel(colorbar,'[deg C]')
ylabel('Length Along tether unstretched [m]'); xlabel('[UTC]'); datetick;axis tight;axis ij;

%linkaxes(findobj(allchild(gcf), 'Type', 'Axes'),'xy');

%name = [rawdatapath , ' Tchain'];
%sgtitle(name)
%savefig(name)




