% Assignment on link budget for satellite communications

clc
clear
close all

%% Setting up the scenario %%
startTime=datetime(2023,01,01,00,00,0);
stopTime=startTime+days(7);
sampleTime=60;
sc=satelliteScenario(startTime,stopTime,sampleTime);

%% Satellite object definition %%
hSat=500e3;
inclAngleSat=97.4;
EOSat.semiMajorAxis=earthRadius+hSat; % [m]
EOSat.inclination=inclAngleSat; % [deg]
EOSat.eccentricity=0;
EOSat.rightAscensionOfAscendingNode=0;
EOSat.argumentOfPeriapsis=0;
EOSat.trueAnomaly=0;
EOSat.name="EO Sat";
EOSat.satellite=satellite(sc,EOSat.semiMajorAxis,EOSat.eccentricity,EOSat.inclination, ...
    EOSat.rightAscensionOfAscendingNode,EOSat.argumentOfPeriapsis,EOSat.trueAnomaly,"Name",EOSat.name);

% Satellite gimbal definition
EOSat.gimbal=gimbal(EOSat.satellite);

% Satellite transmitter definition
fc=8200e6;      %%% Carrier frequency in 8 GHz - 8.4 GHz (X-band)
pTxSat=1.76;    %%% Tx power: 1.5 W (Enduro SAT) -> 1.76 dBW
sysLSat=1;      %%% System loss
m=2;            %%% Modulation: QPSK
ro=0.2;         %%% Roll-off factor of the shaping filter
Rbnet=64e6;     %%% Target bit rate: constrained by Rbgross < 150 Mbps
Rcode=0.5;      %%% Code rate: 1/2
Rs=Rbnet/(Rcode*m);
Bandwidth=((1+ro)*Rs)/1e6; %[MHz]
bitRate=Bandwidth; 
txSat=transmitter(EOSat.gimbal,Name="Sat Tx",Frequency=fc,Power=pTxSat,BitRate=bitRate,SystemLoss=sysLSat);

% Satellite Gaussian antenna definition
dishSat=0.183;
effSat=0.65;
gaussianAntenna(txSat,DishDiameter=dishSat,ApertureEfficiency=effSat);

%% Inuvik ground station %%
elAngle=10;
gs1=groundStation(sc,Name="Inuvik",Latitude=68.35,Longitude=-133.72,Altitude=15,MinElevationAngle=elAngle);

% Station gimbal definition
gimbgs1_clear=gimbal(gs1,MountingAngles=[0;180;0]);
gimbgs1_rain=gimbal(gs1,MountingAngles=[0;180;0]);

% Station receiver definition for clear sky condition
GT_clear=25;
sysLGS=1;              %%% Ground station loss
EbNoThresholdMu=2+2;   %%% 
Rbgross=Rbnet/Rcode;   %%% Max Rbgross: 150 Mbps
B=Bandwidth*1e6;
SpectEff=m*Rcode;
thSNR=EbNoThresholdMu+SpectEff-10*log10(1+ro);
rxgs1_clear=receiver(gimbgs1_clear,Name="Inuvik RX",GainToNoiseTemperatureRatio=GT_clear,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station receiver definition for rain condition
GrxGS=47.8; % Obtained from 10*log10(141.50+50)+G/T
dishGS=3.7;
effGS=0.65;
tempRx=141.5;
PLcfgP618=p618Config(Frequency=fc,ElevationAngle=elAngle,Latitude=68.35,Longitude=-133.72,TotalAnnualExceedance=0.1,AntennaDiameter=dishGS,AntennaEfficiency=effGS);
[PL,~,Tsky_rain]=p618PropagationLosses(PLcfgP618);
GT_rain=GrxGS-10*log10(tempRx+Tsky_rain);
rxgs1_rain=receiver(gimbgs1_rain,Name="Inuvik RX",GainToNoiseTemperatureRatio=GT_rain,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station antenna definition
gaussianAntenna(rxgs1_clear,DishDiameter=dishGS,ApertureEfficiency=effGS);
gaussianAntenna(rxgs1_rain,DishDiameter=dishGS,ApertureEfficiency=effGS);

%% Link budget evaluation for Inuvik ground station %%

% Gimbals orientation to each other
pointAt(EOSat.gimbal,gs1);
pointAt(gimbgs1_clear,EOSat.satellite);
pointAt(gimbgs1_rain,EOSat.satellite);

% Link creation for clear sky
SatLink_clear=link(txSat,rxgs1_clear);
intClear1=linkIntervals(SatLink_clear); % when satellite can communicate, SNR being up than the threshold
DailyLinkDurationClear1=retime(table2timetable([intClear1(:,4),intClear1(:,6)]),'daily','sum');
DailyCapacityClear1=DailyLinkDurationClear1.Duration*(Rs*SpectEff);
meanDailyCapacityClear1=sum(DailyCapacityClear1)/7;
totAvailabilityPercentage1=(sum(DailyLinkDurationClear1.Duration)*100)/(168*3600);
Tab1=table2timetable([intClear1(:,4),intClear1(:,5)]);
Tab1_start=Tab1.StartTime(2:end);
Tab1_end=Tab1.EndTime(1:end-1);
Latency1=between(Tab1_end,Tab1_start);
meanLatency1=mean(time(Latency1));
maxLatency1=max(time(Latency1));
[SNR_clear,timeSamples]=ebno(SatLink_clear);

% Link creation for rainy sky
SatLink_rain=link(txSat,rxgs1_rain);
intRain1=linkIntervals(SatLink_rain);
DailyLinkDurationRain1=retime(table2timetable([intRain1(:,4),intRain1(:,6)]),'daily','sum');
DailyCapacityRain1=DailyLinkDurationRain1.Duration*Rs;
SNR_rain=ebno(SatLink_rain)-PL.At;

% Check the LOS of the satellite and the station
acc_gs1=access(EOSat.satellite,gs1); % when the satellite can see the station
accInt1=accessIntervals(acc_gs1);
DailyAccesses1=retime(table2timetable(accInt1(:,3:4)),'daily','count');
meanDailyAccesses1=sum(DailyAccesses1.IntervalNumber)/7;

% Plotting data
figure;
hold on;
plot(timeSamples,SNR_clear,'b',timeSamples,SNR_rain,'r'),grid on;
yline(thSNR,'--',"Label","Threshold+Margin","Color",'k',"Interpreter","Latex","LabelHorizontalAlignment","left");
axx=xlabel("time");
set(axx,"Interpreter","Latex");
axy=ylabel("SNR [dB]");
set(axy,"Interpreter","Latex");
title("Link performance at Inuvik","Interpreter","Latex");
leg=legend("Clear Sky","Rainy Sky","Location","southeast");
set(leg,"Interpreter","Latex")
hold off;

%% Svalbard ground station %%
gs2=groundStation(sc,Name="Svalbard",Latitude=78.22,Longitude=15.38,Altitude=440,MinElevationAngle=elAngle);

% Station gimbal definition
gimbgs2_clear=gimbal(gs2,MountingAngles=[0;180;0]);
gimbgs2_rain=gimbal(gs2,MountingAngles=[0;180;0]);

% Station receiver definition for clear sky condition
rxgs2_clear=receiver(gimbgs2_clear,GainToNoiseTemperatureRatio=GT_clear,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station receiver definition for rain condition
PLcfgP618=p618Config(Frequency=fc,ElevationAngle=elAngle,Latitude=78.22,Longitude=15.38,TotalAnnualExceedance=0.1,AntennaDiameter=dishGS,AntennaEfficiency=effGS);
[PL,~,Tsky_rain]=p618PropagationLosses(PLcfgP618);
GT_rain=GrxGS-10*log10(tempRx+Tsky_rain);
rxgs2_rain=receiver(gimbgs2_rain,GainToNoiseTemperatureRatio=GT_rain,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station antenna definition
gaussianAntenna(rxgs2_clear,DishDiameter=dishGS,ApertureEfficiency=effGS);
gaussianAntenna(rxgs2_rain,DishDiameter=dishGS,ApertureEfficiency=effGS);

%% Link budget evaluation for Svalbard ground station %%

% Gimbals orientation to each other
pointAt(EOSat.gimbal,gs2);
pointAt(gimbgs2_clear,EOSat.satellite);
pointAt(gimbgs2_rain,EOSat.satellite);

% Link creation for clear sky
SatLink_clear=link(txSat,rxgs2_clear);
intClear2=linkIntervals(SatLink_clear); % when satellite can communicate, being greater than the threshold
DailyLinkDurationClear2=retime(table2timetable([intClear2(:,4),intClear2(:,6)]),'daily','sum');
DailyCapacityClear2=DailyLinkDurationClear2.Duration*(Rs*SpectEff);
meanDailyCapacityClear2=sum(DailyCapacityClear2)/7;
totAvailabilityPercentage2=(sum(DailyLinkDurationClear2.Duration)*100)/(168*3600);
Tab2=table2timetable([intClear2(:,4),intClear2(:,5)]);
Tab2_start=Tab2.StartTime(2:end);
Tab2_end=Tab2.EndTime(1:end-1);
Latency2=between(Tab2_end,Tab2_start);
meanLatency2=mean(time(Latency2));
maxLatency2=max(time(Latency2));
[SNR_clear,timeSamples]=ebno(SatLink_clear);

% Link creation for rainy sky
SatLink_rain=link(txSat,rxgs2_rain);
intRain2=linkIntervals(SatLink_rain);
DailyLinkDurationRain2=retime(table2timetable([intRain2(:,4),intRain2(:,6)]),'daily','sum');
DailyCapacityRain2=DailyLinkDurationRain2.Duration*Rs;
SNR_rain=ebno(SatLink_rain)-PL.At;

% Check the LOS of the satellite and the station
acc_gs2=access(EOSat.satellite,gs2); % when the satellite can see the station
accInt2=accessIntervals(acc_gs2);
DailyAccesses2=retime(table2timetable(accInt2(:,3:4)),'daily','count');
meanDailyAccesses2=sum(DailyAccesses2.IntervalNumber)/7;

% Plotting data
figure;
hold on;
plot(timeSamples,SNR_clear,'b',timeSamples,SNR_rain,'r'),grid on;
yline(thSNR,'--',"Label","Threshold+Margin","Color",'k',"Interpreter","Latex","LabelHorizontalAlignment","left");
axx=xlabel("time");
set(axx,"Interpreter","Latex");
axy=ylabel("SNR [dB]");
set(axy,"Interpreter","Latex");
title("Link performance at Svalbard","Interpreter","Latex");
leg=legend("Clear Sky","Rainy Sky","Location","southeast");
set(leg,"Interpreter","Latex")
hold off;

%% Awarua ground station %%
gs3=groundStation(sc,Name="Awarua",Latitude=-46.52,Longitude=168.48,Altitude=0,MinElevationAngle=elAngle);

% Station gimbal definition
gimbgs3_clear=gimbal(gs3,MountingAngles=[0;180;0]);
gimbgs3_rain=gimbal(gs3,MountingAngles=[0;180;0]);

% Station receiver definition for clear sky condition
rxgs3_clear=receiver(gimbgs3_clear,GainToNoiseTemperatureRatio=GT_clear,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station receiver definition for rain condition
PLcfgP618=p618Config(Frequency=fc,ElevationAngle=elAngle,Latitude=-46.52,Longitude=168.48,TotalAnnualExceedance=0.1,AntennaDiameter=dishGS,AntennaEfficiency=effGS);
[PL,~,Tsky_rain]=p618PropagationLosses(PLcfgP618);
GT_rain=GrxGS-10*log10(tempRx+Tsky_rain);
rxgs3_rain=receiver(gimbgs3_rain,GainToNoiseTemperatureRatio=GT_rain,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station antenna definition
gaussianAntenna(rxgs3_clear,DishDiameter=dishGS,ApertureEfficiency=effGS);
gaussianAntenna(rxgs3_rain,DishDiameter=dishGS,ApertureEfficiency=effGS);

%% Link budget evaluation for Awarua ground station %%

% Gimbals orientation to each other
pointAt(EOSat.gimbal,gs3);
pointAt(gimbgs3_clear,EOSat.satellite);
pointAt(gimbgs3_rain,EOSat.satellite);

% Link creation for clear sky
SatLink_clear=link(txSat,rxgs3_clear);
intClear3=linkIntervals(SatLink_clear); % when satellite can communicate, being greater than the threshold
DailyLinkDurationClear3=retime(table2timetable([intClear3(:,4),intClear3(:,6)]),'daily','sum');
DailyCapacityClear3=DailyLinkDurationClear3.Duration*(Rs*SpectEff);
meanDailyCapacityClear3=sum(DailyCapacityClear3)/7;
totAvailabilityPercentage3=(sum(DailyLinkDurationClear3.Duration)*100)/(168*3600);
Tab3=table2timetable([intClear3(:,4),intClear3(:,5)]);
Tab3_start=Tab3.StartTime(2:end);
Tab3_end=Tab3.EndTime(1:end-1);
Latency3=between(Tab3_end,Tab3_start);
meanLatency3=mean(time(Latency3));
maxLatency3=max(time(Latency3));
[SNR_clear,timeSamples]=ebno(SatLink_clear);

% Link creation for rainy sky
SatLink_rain=link(txSat,rxgs3_rain);
intRain3=linkIntervals(SatLink_rain);
DailyLinkDurationRain3=retime(table2timetable([intRain3(:,4),intRain3(:,6)]),'daily','sum');
DailyCapacityRain3=DailyLinkDurationRain3.Duration*Rs;
SNR_rain=ebno(SatLink_rain)-PL.At;

% Check the LOS of the satellite and the station
acc_gs3=access(EOSat.satellite,gs3); % when the satellite can see the station
accInt3=accessIntervals(acc_gs3);
DailyAccesses3=retime(table2timetable(accInt3(:,3:4)),'daily','count');
meanDailyAccesses3=sum(DailyAccesses3.IntervalNumber)/7;

% Plotting data
figure;
hold on;
plot(timeSamples,SNR_clear,'b',timeSamples,SNR_rain,'r'),grid on;
yline(thSNR,'--',"Label","Threshold+Margin","Color",'k',"Interpreter","Latex","LabelHorizontalAlignment","left");
axx=xlabel("time");
set(axx,"Interpreter","Latex");
axy=ylabel("SNR [dB]");
set(axy,"Interpreter","Latex");
title("Link performance at Awarua","Interpreter","Latex");
leg=legend("Clear Sky","Rainy Sky","Location","southeast");
set(leg,"Interpreter","Latex")
hold off;

%% Troll ground station %%
gs4=groundStation(sc,Name="Troll",Latitude=-72.01,Longitude=2.53,Altitude=1200,MinElevationAngle=elAngle);

% Station gimbal definition
gimbgs4_clear=gimbal(gs4,MountingAngles=[0;180;0]);
gimbgs4_rain=gimbal(gs4,MountingAngles=[0;180;0]);

% Station receiver definition for clear sky condition
rxgs4_clear=receiver(gimbgs4_clear,GainToNoiseTemperatureRatio=GT_clear,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station receiver definition for rain condition
PLcfgP618=p618Config(Frequency=fc,ElevationAngle=elAngle,Latitude=-72.01,Longitude=2.53,TotalAnnualExceedance=0.1,AntennaDiameter=dishGS,AntennaEfficiency=effGS);
[PL,~,Tsky_rain]=p618PropagationLosses(PLcfgP618);
GT_rain=GrxGS-10*log10(tempRx+Tsky_rain);
rxgs4_rain=receiver(gimbgs4_rain,GainToNoiseTemperatureRatio=GT_rain,SystemLoss=sysLGS,RequiredEbNo=thSNR,PreReceiverLoss=0);

% Station antenna definition
gaussianAntenna(rxgs4_clear,DishDiameter=dishGS,ApertureEfficiency=effGS);
gaussianAntenna(rxgs4_rain,DishDiameter=dishGS,ApertureEfficiency=effGS);

%% Link budget evaluation for Troll ground station %%

% Gimbals orientation to each other
pointAt(EOSat.gimbal,gs4);
pointAt(gimbgs4_clear,EOSat.satellite);
pointAt(gimbgs4_rain,EOSat.satellite);

% Link creation for clear sky
SatLink_clear=link(txSat,rxgs4_clear);
intClear4=linkIntervals(SatLink_clear); % when satellite can communicate, being greater than the threshold
DailyLinkDurationClear4=retime(table2timetable([intClear4(:,4),intClear4(:,6)]),'daily','sum');
DailyCapacityClear4=DailyLinkDurationClear4.Duration*(Rs*SpectEff);
meanDailyCapacityClear4=sum(DailyCapacityClear4)/7;
totAvailabilityPercentage4=(sum(DailyLinkDurationClear4.Duration)*100)/(168*3600);
Tab4=table2timetable([intClear4(:,4),intClear4(:,5)]);
Tab4_start=Tab4.StartTime(2:end);
Tab4_end=Tab4.EndTime(1:end-1);
Latency4=between(Tab4_end,Tab4_start);
meanLatency4=mean(time(Latency4));
maxLatency4=max(time(Latency4));
[SNR_clear,timeSamples]=ebno(SatLink_clear);

% Link creation for rainy sky
SatLink_rain=link(txSat,rxgs4_rain);
intRain4=linkIntervals(SatLink_rain);
DailyLinkDurationRain4=retime(table2timetable([intRain4(:,4),intRain4(:,6)]),'daily','sum');
DailyCapacityRain4=DailyLinkDurationRain4.Duration*Rs;
SNR_rain=ebno(SatLink_rain)-PL.At;

% Check the LOS of the satellite and the station
acc_gs4=access(EOSat.satellite,gs4); % when the satellite can see the station
accInt4=accessIntervals(acc_gs4);
DailyAccesses4=retime(table2timetable(accInt4(:,3:4)),'daily','count');
meanDailyAccesses4=sum(DailyAccesses4.IntervalNumber)/7;

% Plotting data
figure;
hold on;
plot(timeSamples,SNR_clear,'b',timeSamples,SNR_rain,'r'),grid on;
yline(thSNR,'--',"Label","Threshold+Margin","Color",'k',"Interpreter","Latex","LabelHorizontalAlignment","left");
axx=xlabel("time");
set(axx,"Interpreter","Latex");
axy=ylabel("SNR [dB]");
set(axy,"Interpreter","Latex");
title("Link performance at Troll","Interpreter","Latex");
leg=legend("Clear Sky","Rainy Sky","Location","southeast");
set(leg,"Interpreter","Latex")
hold off;

%% Plots per station and per network %%

stations=categorical({'Inuvik','Svalbard','Awarua','Troll'});

% Mean Daily Capacity [Tbit] 
DailyCapacityNetClear=horzcat(DailyCapacityClear1,DailyCapacityClear2,DailyCapacityClear3,DailyCapacityClear4);
meanDailyCapacityNet=sum(sum(DailyCapacityNetClear,2))/7
meanCapacityArray=[meanDailyCapacityClear1;meanDailyCapacityClear2;meanDailyCapacityClear3;meanDailyCapacityClear4];
figure;
bar(stations,meanCapacityArray./1e12),grid on;
tit=title('Mean daily capacity chart');
axy=ylabel('Tbit');
set(axy,"Interpreter","Latex");
set(tit,"Interpreter","Latex");

% Mean/Max Latencies [min] -> duration(x,'Format','hh:mm:ss') 
% for standard display of latency times
TabNet=sortrows(vertcat(Tab1,Tab2,Tab3,Tab4));
TabNet_start=TabNet.StartTime(2:end);
TabNet_end=TabNet.EndTime(1:end-1);
LatencyNet=between(TabNet_end,TabNet_start);
meanLatencyNet=duration(mean(time(LatencyNet)),'Format','m')
maxLatencyNet=duration(max(time(LatencyNet)),'Format','m')
meanLatencyArray=[duration(meanLatency1,'Format','m'),duration(meanLatency2,'Format','m'),duration(meanLatency3,'Format','m'),duration(meanLatency4,'Format','m')];
figure;
bar(stations,meanLatencyArray),grid on;
tit=title('Mean latency chart');
set(tit,"Interpreter","Latex");
maxLatencyArray=[duration(maxLatency1,'Format','m'),duration(maxLatency2,'Format','m'),duration(maxLatency3,'Format','m'),duration(maxLatency4,'Format','m')];
figure;
bar(stations,maxLatencyArray),grid on;
tit=title('Max latency chart');
set(tit,"Interpreter","Latex");

% Mean Daily Accesses
DailyAccessesNet=horzcat(DailyAccesses1.IntervalNumber,DailyAccesses2.IntervalNumber,DailyAccesses3.IntervalNumber,DailyAccesses4.IntervalNumber);
meanDailyAccessesNet=sum(sum(DailyAccessesNet,2))/7
meanAccessesArray=[meanDailyAccesses1,meanDailyAccesses2,meanDailyAccesses3,meanDailyAccesses4];
figure;
bar(stations,meanAccessesArray),grid on;
tit=title('Mean daily accesses chart');
set(tit,"Interpreter","Latex");

% Availability [%]
DailyDurationNet=sum(horzcat(DailyLinkDurationClear1.Duration,DailyLinkDurationClear2.Duration,DailyLinkDurationClear3.Duration,DailyLinkDurationClear4.Duration),2);
totAvailabilityPercentageNet=(sum(DailyDurationNet)*100)/(168*3600)
totAvailabilityPercentageArray=[totAvailabilityPercentage1,totAvailabilityPercentage2,totAvailabilityPercentage3,totAvailabilityPercentage4];
figure;
bar(stations,totAvailabilityPercentageArray),grid on;
tit=title('Availability chart');
axy=ylabel('Percentage');
set(axy,"Interpreter","Latex");
set(tit,"Interpreter","Latex");

%% LEO Data Relay System Performance %%

% Walker Constellation definition
LDRS.planes=3;
LDRS.satsPerPlane=6;
LDRS.numSats=LDRS.planes*LDRS.satsPerPlane;
LDRS.eccentricity=0;
LDRS.inclination=98.2;
LDRS.altitude=1000e3;
LDRS.semiMajorAxis=earthRadius+LDRS.altitude;
LDRS.argOfPeriapsis=0;
% Walker
LDRS.planePhase=360/LDRS.planes;
LDRS.satsPhase=360/LDRS.satsPerPlane;
LDRS.RAAN=zeros(1,LDRS.numSats);
LDRS.trueAnomaly=zeros(1,LDRS.numSats);
idx=0;
for plane=1:LDRS.planes
    for sat=1:LDRS.satsPerPlane
        idx=idx+1;
        LDRS.RAAN(idx)=plane*LDRS.planePhase;
        LDRS.trueAnomaly(idx)=sat*LDRS.satsPhase;
    end
end
LDRS.Satellites=satellite(sc,LDRS.semiMajorAxis*ones(1,LDRS.numSats),LDRS.eccentricity*ones(1,LDRS.numSats),LDRS.inclination*ones(1,LDRS.numSats), ...
                            LDRS.RAAN,LDRS.argOfPeriapsis*ones(1,LDRS.numSats),LDRS.trueAnomaly);

% Adding gimbal with conical sensor to EOSat
EOSat.LCTGimbal=gimbal(EOSat.satellite);
EOSat.LCT(1)=conicalSensor(EOSat.LCTGimbal,MaxViewAngle=100);

% Compute thes LOS between the satellite and each constellation's satellite
acc=access(EOSat.satellite,LDRS.Satellites);
accTable=accessStatus(acc);

% Check whether they're in the laser range
[azimuth,elevation,range] = aer(EOSat.satellite,LDRS.Satellites);
accTableOnRangeTable=accTable.*range;
LDRSVisibleSatsOverSimTime=0;
connectivityLEO=zeros(18,10081);
for i=1:18
    for j=1:10081
        if (accTableOnRangeTable(i,j)<0.4e+07) && (accTableOnRangeTable(i,j)>0)
            LDRSVisibleSatsOverSimTime=LDRSVisibleSatsOverSimTime+1;
            connectivityLEO(i,j)=1;
        end
    end
end
LDRSVisibleSatsOverSimTime;


%% Plots for LEO data relay system %%

% Connectivity, Latency, Availability
connectivityLEOArray=sum(connectivityLEO,1);
positions=zeros(10081,1);
for i=1:length(connectivityLEOArray)
    if connectivityLEOArray(i)~=0
        connectivityLEOArray(i)=1;
    end
end
connectivityCheck=all(connectivityLEOArray);
if connectivityCheck==1
    figure;
    plot(timeSamples,connectivityLEOArray,'.'),grid on;
    tit=title('LEO data relay system - Connectivity');
    set(tit,"Interpreter","Latex");
    LatencyLEO=0
    AvailabilityPercentageLEO=100
end