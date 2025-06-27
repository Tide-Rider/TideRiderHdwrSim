% Hardware-in-the-loop simulation of TideRider.
%
% Status as of 2025/04/14:
% * crude depth sim works; simulated TR does not track surface, however, so confusing.  Just have it "stick" to a band on the surface using logic rather than dynamics.
% * lat/lon sketched out but (kinematic) current not applied.
% * does not produce lat/lon fixes yet and TR code needs to be reloaded onto arduino.
%
%
% Revision History
% 2025-04-14    mvj    Created.
% 2025-04-18    mvj    added current; untested.  TR sleeps a lot; will need to update state on timeout, not just async.
% 2025-05-16    mvj    added ui window for sending passthru.




if exist('s','var')
    s.delete();
end

s = serialport('COM4',9600);
%s = serialport('/dev/ttyACM0',9600);
%s = serialport('/dev/ttyUSB0',9600,'Timeout',1);
configureTerminator(s,'CR/LF');

s.flush();

% Use to set phase of tide.
t0 = now*3600*24;

% environmental parameters
global MLW TIDE TIDAL_PERIOD LAT0 LON0 CURRENT TIDAL_PHASE
MLW = 0.5; % constant with tide.
TIDE = 2; % m
TIDAL_PERIOD = 12*3600;
LAT0 = 41;
LON0 = -70;
CURRENT = 0.3; % m/s Assumed to change on same tidal period.
TIDAL_PHASE = -2*pi*t0/TIDAL_PERIOD;

% state
[zw] = get_water_depth(LAT0,LON0,t0);
z = min(MLW-zw,MLW);
zt = 0;
easting = 0;
northing = 0;

% set up plots.
figure(1); clf reset; % depth versus time.

hbot = line(NaN,NaN,'linestyle','-','linewidth',2,'color',brown);
hsurf = line(NaN,NaN,'linestyle','-','linewidth',2,'color','b');
htr1 = line(NaN,NaN,'linestyle','-','linewidth',2,'color','r');
ydepth;

figure(2); clf reset; % n/e position.
htr2 = line(NaN,NaN,'linestyle','none','marker','o','markerfacecolor','r','color','r');
hwp = line(NaN,NaN,'linestyle','none','marker','o','markerfacecolor',[0 0.5 0],'color',[0 0.5 0]);
xlim(LON0 + 0.001*[-1 1]);
ylim(LAT0 + 0.001*[-1 1]);

% to insert commands
% $SMPAS,<cmd>*00  Note - generally TR cmds end with a ','
if exist('uifig','var'); delete uifig; end    
uifig = uifigure('Name', 'Passthru');

% Add a label
uilabel(uifig, ...
    'Text', 'Send:', ...
    'Position', [50, 120, 100, 30]);

% Add a text input field
textInput = uieditfield(uifig, 'text', ...
    'Position', [160, 120, 180, 30]);

% Add a button to display the input
uibutton(uifig, ...
    'Text', 'Send', ...
    'Position', [160, 70, 100, 30], ...
    'ButtonPushedFcn', @(btn, event) writeline(s,textInput.Value));

t_ = t0;
thrcmd = 1500;  % neutral.
while(1)

    % code is query/response so no need for threading.
    ln = s.readline();  % returns empty on timeout.
    
    if ~isempty(ln)
        if ln.startsWith("$TRSMQ")
            
            fprintf(1,'Got sim query.\n');
            
            if ln.startsWith("$TRSMQ,BRD");
                
                fprintf(1,'Got BRD request.\n');
                sprintf("$SMBRD,%.3f,m,0,C*00",z-zs)
                s.writeline(sprintf("$SMBRD,%.3f,m,0,C*00",z-zs)); % pressure
                
            elseif ln.startsWith("$TRSMQ,FIX");
                
                fprintf(1,'Got FIX request.\n');
                sprintf("$SMFIX,%.6f,%.6f,V*00",lat,lon)
                s.writeline(sprintf("$SMFIX,%.6f,%.6f,V*00",lat,lon));
                
            end
            
        elseif ln.startsWith("$TRTHR")
            
            % $TRTHR,1500*3c
            fprintf(1,'Got thruster command: %s\n',ln);
            thrcmd = sscanf(ln,"$TRTHR,%f*");
            
        else
            
            fprintf(1,'Unhandled: %s\n',ln);
            
        end
    end
        
    % simulate.  Happens on message or timeout.
    t = now*3600*24;
    dt = t-t_;
    t_ = t;
    
    % compute water depth
    [zw] = get_water_depth(LAT0,LON0,t);
    zb = MLW;
    zs = min(0,MLW-zw);

    % compute acceleration and integrate.
    ztt = dyn(z,zt,thrcmd);
    zt = zt + ztt*dt;
    z = z + zt*dt;

    % apply bounds.
    if z <= zs & ztt < 0
        z = zs;
        zt = 0;
    elseif z >= zb & ztt > 0
        z = zb;
        zt = 0;
    end

    % action with current is purely kinematic.
    if z == zb
        nt = 0;
    else
        nt = CURRENT*cos(2*pi/TIDAL_PERIOD*t + TIDAL_PHASE); % 90 out of phase with tide.
    end
    northing = northing + nt*dt;
    [lat,lon] = ned2geodetic(northing,easting,0,LAT0,LON0,0,wgs84Ellipsoid);
    
    % update plots
    figure(1);
    TWIN = [t-600,t+60];
    xlim(TWIN);
    gt = get(htr1,'xdata');
    gz = get(htr1,'ydata');
    ii = gt >= TWIN(1) & gt < TWIN(2);
    set(htr1,'xdata',[gt(ii) t],'ydata',[gz(ii) z]);
    gz = get(hbot,'ydata');
    set(hbot,'xdata',[gt(ii) t],'ydata',[gz(ii) zb]);
    gz = get(hsurf,'ydata');
    set(hsurf,'xdata',[gt(ii) t],'ydata',[gz(ii) zs]);

    figure(2);
    set(htr2,'xdata',lon,'ydata',lat);
    set(hwp,'xdata',LON0,'ydata',LAT0);  % no mechanism to synchronize WPTs between device and this code. (could be added).
    
end

s.delete();


%
% Dynamics etc.
%

function [zw] = get_water_depth(lat,lon,t)

    global MLW TIDE TIDAL_PERIOD TIDAL_PHASE
    zw = MLW + TIDE/2 + TIDE/2*sin(2*pi/TIDAL_PERIOD*t + TIDAL_PHASE);

end

function ztt = dyn(z,zt,cmd)

    
    % get force from thruster command
    Z = (cmd-1500)*10/400;  % no idea which way is positive.

    % compute drag
    Zww = 2;
    
    % compute buoyancy
    % this is somewhat tricky if want to sim action of scuttle.  is there some crude way to capture the essentials?
    % yes - do nothing.
    Zbuoy = 0;

    % compute nominal acceration
    m = 10;
    ztt = 1/m*(Z - Zww*zt*abs(zt) - Zbuoy);
    
end    
    
