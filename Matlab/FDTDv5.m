function [echo] = FDTDv5(scenario,f,saveRCS,saveSignal,viewMovie,saveMovie,speedFactor,fps)

%Constants
eps0 = 8.854e-12;
u0 = 4*pi*1e-7;
c0 = 1/sqrt(eps0*u0);
eta0 = sqrt(u0 / eps0);

%Simulation Inputs
d = 0.5; %[m] height of simulation
Tk = 273.15-10; %Temperature, Kelvin
snowDensity = 0.3;%[g/cm^3] Typical value 0.1 < x < 0.4.  Taken from figure 4-17, pg 141.

%Unused Sim Values
% scenario = 2;
% saveRCS = 0; %Save NCRS data in table
% viewMovie = 1;
% saveMovie = 0;
% saveSignal = 0;
% speedFactor = 5; %Only display every nth frame
% fps = 48; %default 12
% f = 1e9; %[Hz]

depths = [0.2, 0.3];
diels = ["snow", "ice"];
epsrSpecified = 2;

%Preset Scenarios
if scenario == 1 %Air
    depths = [d/2];
    diels = ["air"];
    Tks = ones(1,length(depths)) * Tk;

elseif scenario == 2 %Basic dielectric relative permittivity of 4
    d = 1;
    depths = [d/2];
    diels = ["0"];
    epsrSpecified = [4];
    Tks = ones(1,length(depths)) * Tk;

elseif scenario == 3 %Snow on ice
    depths = [0.2, 0.3];
    diels = ["snow", "ice"];
    Tks = ones(1,length(depths)) * Tk;

elseif scenario == 4 %Same as below, but without oil
    depths = [0.05];
    diels = ["ice"];
    Tks = ones(1,length(depths)) * Tk;

elseif scenario == 5 %Examining the Impact of a Crude Oil Spill on the Permittivity Profile and Normalized Radar Cross Section of Young Sea Ice
    d = 0.4;
    depths = (0:0.025:0.2) + 0.05;
    diels = string(zeros(1,length(depths)));
    erprime = [3.75, 4, 4.5, 5, 5.5, 5.5, 6.25, 6.25, 6.25];
    erDoublePrime = [0, 0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.7, 0.7];
    epsrSpecified = complex(erprime, erDoublePrime);
    Tks = ones(1,length(depths)) * Tk;
    
elseif scenario == 6 %Snow-ice-sea
    d = 3;
    depths = [0.5,1,1.5,1.6,2.5];
    diels = ["snow", "ice", "brine", "ice", "sea water"];
    Tks = ones(1,length(depths)) * Tk;
    
end


%depths = [0.25];
%diels = ["air"];

%Wave properties
omega = 2*pi*f;
Emag = 1; %[V/m]
EmagNorm = Emag / eta0;
k = omega / c0;
lambda0 = 2 * pi / k;

nlayers = length(diels);
[epsLayers,~,~,~] = generateGrids(nlayers + 1, nlayers + 1, (1:nlayers) - 0.5, diels, epsrSpecified, Tks, f, snowDensity);
epsrLayers = epsLayers / eps0;
epsrMax = max(epsrLayers);
lambdaMin = lambda0 / sqrt(epsrMax);

%Simulation Parameters
dz = lambdaMin / 20;
dt = dz / c0 / 2;
zs = 0:dz:d;
ncells = length(zs);

[eps, conds, epsImag, diels] = generateGrids(d, ncells, depths, diels, epsrSpecified, Tks, f, snowDensity);

lossTangent = epsImag ./ eps;
alpha = omega * sqrt(u0*eps/2.*(sqrt(1+lossTangent.^2)-1));
beta = omega * sqrt(u0*eps/2.*(sqrt(1+lossTangent.^2)+1));
vp = omega ./ beta;
simTime = sum(dz./vp) * 2;
nsteps = ceil(simTime / dt);
t = 0:dt:dt*(nsteps - 1);

%Source
q = 0.05; %damping factor, 0.2 for 250MHz, 0.05 for 1GHz
Esource = Emag * sinPulse(t,f);
Esource = Emag * rickerPulse(t,q,f);
%Esource = Emag * sin(omega*t);

%Initialize Fields
Ex = zeros(1, ncells);
Hy = zeros(1, ncells);

%Set up ABC
%Far Side
oldE2far = Ex(ncells-1);
oldE1far = oldE2far;
%Close Side
oldE2LHS = Ex(1);
oldE1LHS = oldE2LHS;
oldSource1 = 0;
oldSource2 = oldSource1;

%Constant Terms used in calculations
epsr = eps/eps0;
q1 = conds * dt ./ (2 * eps);
q2 = c0 * dt / dz;

k1 = (1 - q1) ./ (1 + q1);
k2 = q2 ./ epsr ./ (1 + q1);
k3 = 1;
k4 = q2;

%Reciever
echo = zeros(1,nsteps);

%Simulation Movie
j = 0;
if (viewMovie || saveMovie)
    M = moviein(floor(nsteps/speedFactor));
    p = simPlots(depths, diels, zs, -2, 2);
end

for i = 1:(nsteps)

    %Efield calculations
    Ex(1) = EmagNorm * Esource(i) + oldE2LHS;
    Ex(2:ncells-1) = k1(2:ncells-1) .* Ex(2:ncells-1) + k2(2:ncells-1) .* (Hy(1:ncells-2) - Hy(2:ncells-1));

    %Absorbing Boundary Conditions
    %Far Side
    Ex(ncells) = oldE2far;
    %oldE2far + (c0*dt - dz)/(c0*dt + dz) * (oldE2far) + 2*dz / (c0*dt + dz) * (oldE1far);
    %*exp(-alpha(ncells)*dz);
    
%     oldE2far = oldE1far;
%     oldE1far = Ex(ncells);
    
    oldE2far = oldE1far;
    oldE1far = Ex(ncells - 1);
    %Close Side
    oldE2LHS = oldE1LHS;
    oldE1LHS = Ex(2) - oldSource2;
    oldSource2 = oldSource1;
    oldSource1 = EmagNorm * Esource(i);

    %H Field calculations
    Hy(1:ncells-1) = k3 .* Hy(1:ncells-1) + k4 .* (Ex(1:ncells-1) - Ex(2:ncells));

    %Recievd Field
    echo(i) = Ex(1) * eta0;

    if(viewMovie || saveMovie)
        j = j + 1;

        %Plots
        if (j >= speedFactor)
            j = 0;
            ex = Ex * eta0;
            p.YData = ex;
            M(:,round(i/speedFactor)) = getframe;
        end
    end

end

if viewMovie
    movie(M,1,fps);
end

%Plot recieved signal
figure()
plot(t,echo)

if saveMovie
    name = "simulationScenario" + string(scenario) + " f = " + f/1e9 + "GHz";
    name = char(name);
    v = VideoWriter(name,'Archival');
    open(v);
    writeVideo(v,M);
    close(v);

end

if saveSignal
    name = "simulationScenario" + string(scenario) + " f = " + f/1e9 + "GHz";
    name = char(name);
    saveas(1,name,'png')
    csvwrite('echoSignal.csv',echo);
    csvwrite('timeStamps.csv',t);
    csvwrite('sentPulse.csv', Esource);
end

if (scenario > 0 && saveRCS)
    %Save
    nrcsIndex = csvread('nrcsIndex.csv');
    returnSignal = echo - Esource;
    echoPower = sum(returnSignal.^2);
    sentPower = sum(Esource.^2);
    nrcsIndex(scenario) = echoPower / sentPower;
    csvwrite('nrcsIndex.csv', nrcsIndex);
end
    
