%function Network_N_C_direct_calc2(runID)

%% Network N C direct calc
% This program uses the Le Sueur River network/wetland complex to directly
% calculate N and C concentration without numerical simulation
% Set-up for Nitrogen and Carbon transport
% Set-up for the Le Sueur River Network/Wetland Complex

% Jon Czuba
% July 19, 2016

%% Import Network, Initialize Variables
clear all
close all
clc

% QQ=1:1:400;
% 
% NN(1:LinkNum,1:length(QQ))=NaN;
% CC(1:LinkNum,1:length(QQ))=NaN;
% CN(1:LinkNum,1:length(QQ))=NaN;
% DA(1:LinkNum,1:length(QQ))=NaN;
%OUT(1:length(QQ),1:7)=NaN;

% for qq=1:length(QQ)
%     clearvars -except NN CC CN DA qq QQ
%    clearvars -except OUT qq QQ
%     qq
     
% XX(1:LinkNum*2,1:103)=NaN;
% YYB1(1:LinkNum*2,1:103)=NaN;
% YYB2(1:LinkNum*2,1:103)=NaN;
% YYB3(1:LinkNum*2,1:103)=NaN;
% YYB4(1:LinkNum*2,1:103)=NaN;
% YYR1(1:LinkNum*2,1:103)=NaN;
% YYR2(1:LinkNum*2,1:103)=NaN;
% YYR3(1:LinkNum*2,1:103)=NaN;
% YYR4(1:LinkNum*2,1:103)=NaN;

% XX(1:LinkNum*2,1)=NaN;
% YYBn(1:LinkNum*2,1:4)=NaN;
% YYBc(1:LinkNum*2,1:4)=NaN;
% YYBr(1:LinkNum*2,1:4)=NaN;
% YYBd(1:LinkNum*2,1:4)=NaN;
% YYRn(1:LinkNum*2,1:4)=NaN;
% YYRc(1:LinkNum*2,1:4)=NaN;
% YYRr(1:LinkNum*2,1:4)=NaN;
% YYRd(1:LinkNum*2,1:4)=NaN;


%QQ=[1.7 6.4 54 160];

% for qq=1:length(QQ)
% %for xwet=541:LinkNum
% %clearvars -except XX YYB1 YYB2 YYB3 YYB4 YYR1 YYR2 YYR3 YYR4 xwet LinkNum
% clearvars -except XX YYBn YYBc YYBr YYBd YYRn YYRc YYRr YYRd xwet LinkNum QQ qq
% xwet=585;

% %sensitivity on Cin and Jleach
% %cval=(1:0.5:10);%Cin, mg/L
% cval=cat(2,1,10:10:200);%Cin for non ag, mg/L
% jval=cat(2,1,5:5:100);%Jleach, mg/m2/hr
% dcval=cat(2,0.1,0.5:0.5:10);%denitrif. coeff. for DOC
% %dcval=3:0.5:6;%denitrif. coeff. for DOC
% 
% cc=10;
% jj=18;
% dd=8;

% mxmcc=find(MCC==max(max(max(MCC))));
% %mxmcc=find(CRMSE==min(min(min(CRMSE))));
% CVAL=repmat(cval',1,length(jval),length(dcval));
% JVAL=repmat(jval,length(cval),1,length(dcval));
% DCVAL=ones(length(cval),length(jval),length(dcval));
% for i=1:length(dcval)
%     DCVAL(:,:,i)=DCVAL(:,:,i).*dcval(i);
% end

% NRMSE(1:length(cval),1:length(jval),1:length(dcval))=NaN;
% CRMSE(1:length(cval),1:length(jval),1:length(dcval))=NaN;
% ACC(1:length(cval),1:length(jval),1:length(dcval))=NaN;
% MCC(1:length(cval),1:length(jval),1:length(dcval))=NaN;
% for cc=1:length(cval)    
% for jj=1:length(jval)   
% for dd=1:length(dcval) 
%     disp([cc,jj,dd,length(cval),length(jval),length(dcval)]);
%     
%sensitivity on n
% nval=(0.01:0.001:0.1);
% RMSE(1:length(nval),1:4)=NaN;
% for nn=1:length(nval)

%VAL(1:108,1:3)=NaN;
%QVAL(1:76,1:4)=NaN;
%for dmyq=1:7
%for dmyq=7;
%clearvars -except VAL QVAL dmyq cval jval dcval NRMSE CRMSE ACC MCC cc jj dd
%clearvars -except VAL QVAL dmyq cval jval dcval NRMSE CRMSE cc jj dd
%clearvars -except VAL QVAL dmyq nval RMSE nn
%clearvars -except VAL QVAL dmyq
%clearvars -except runID
%clearvars -except Nsave
%clearvars -except network wetland
%close all
%clc

% Import a network
% Determine connectivity of links
load('LS_attributes3.mat');%load river network with lakes and attributes
%load('LS_isolated_wetlands.mat');
%load('LS_iwea.mat');
load('LS_fainN.mat');
load('LS_fainC.mat');
load('LS_site_data_v2.mat');
load('LS_Q_site_data.mat');
load('LS_QBrating.mat');
load('CALidx2.mat');

%load shapefiles for plotting
network = shaperead('Shapefiles\LS_channels2_utm.shp');
wetland = shaperead('Shapefiles\LS_wetlands2.shp');
boundary = shaperead('Shapefiles\LS_boundary1_utm.shp');
%areas = shaperead('Shapefiles\LS_areas1.shp');
%sites = shaperead('Shapefiles\LS_sites_edit.shp');

%% Remove an individual wetland and convert to channel
%xwet=560;
% Feature(xwet)=1;
% WA(xwet)=0;
% pEM(xwet)=0;

% Feature(582)=1;
% WA(582)=0;
% pEM(582)=0;

% %remove all in-channel lakes/wetlands
% Feature(:,1)=1;
% WA(:,1)=0;
% pEM(:,1)=0;
% 
% %remove all isolated wetlands
% fainN(:,1)=1;
%% Assign values into a shapefile structure and write files
% for i = 1:540
%       [network(i).paeff] = paeff(i);
%       [network(i).paw] = paw(i);
%       [network(i).pEM] = pEM(i);
%       [network(i).Feature] = Feature(i);
% end
% for i = 541:LinkNum
%       [wetland(i-540).paeff] = paeff(i);
%       [wetland(i-540).paw] = paw(i);
%       %[wetland(i-540).pEM] = pEM(i);
%       [wetland(i-540).Feature] = Feature(i);
% end
%shapewrite(network,'Shapefiles\LS_channels2_utm.shp');
%shapewrite(wetland,'Shapefiles\LS_wetlands2.shp');

% Nout(1:118,7)=0;
% for dmyq=1:7
%     idx=((NDay==gageDMYQ(dmyq,1))&(NMonth==gageDMYQ(dmyq,2))&(NYear==gageDMYQ(dmyq,3)));
%     sassign=sSiteID(idx);
%     vassign=NO3(idx);
%     
%     for ss=1:length(sassign)
%         Nout(find(sassign(ss)==site),dmyq)=vassign(ss);
%     end
% end
% Nout(isnan(Nout))=0;
% %%
% Qout(1:118,7)=0;
% for dmyq=1:7
%     idx2=((QDay==gageDMYQ(dmyq,1))&(QMonth==gageDMYQ(dmyq,2))&(QYear==gageDMYQ(dmyq,3)));
%     sassign=QsSiteID(idx2);
%     vassign=gQ(idx2);
%     
%     for ss=1:length(sassign)
%         Qout(find(sassign(ss)==site),dmyq)=vassign(ss);
%     end
% end
% Qout(isnan(Qout))=0;
% %%
% for i = 1:118
%       [sites(i).Site_ID2] = site(i,1);
%       [sites(i).N1] = Nout(i,1);
%       [sites(i).N2] = Nout(i,2);
%       [sites(i).N3] = Nout(i,3);
%       [sites(i).N4] = Nout(i,4);
%       [sites(i).N5] = Nout(i,5);
%       [sites(i).N6] = Nout(i,6);
%       [sites(i).N7] = Nout(i,7);
%       
%       [sites(i).Q1] = Qout(i,1);
%       [sites(i).Q2] = Qout(i,2);
%       [sites(i).Q3] = Qout(i,3);
%       [sites(i).Q4] = Qout(i,4);
%       [sites(i).Q5] = Qout(i,5);
%       [sites(i).Q6] = Qout(i,6);
%       [sites(i).Q7] = Qout(i,7);
% 
% end
% shapewrite(sites,'Shapefiles\LS_sites_edit3.shp');
%% Export data to shapefile for plotting in ArcGIS

% for i = 1:540
%       [network(i).N] = N_conc_ds(i);
%       [network(i).fto] = ftransout(i);
%       [network(i).CNlim] = CNlim(i);
%       [network(i).BQlim] = BQlim(i);
% end
% for i = 541:LinkNum
%       [wetland(i-540).N] = N_conc_ds(i);
%       [wetland(i-540).fto] = ftransout(i);
%       [wetland(i-540).CNlim] = CNlim(i);
%       [wetland(i-540).BQlim] = BQlim(i);
% end
% shapewrite(network,'Shapefiles\LS_channels2_utm_2OUT041917.shp');
% shapewrite(wetland,'Shapefiles\LS_wetlands2_2OUT041917.shp');

%% Write Da number

% for i = 1:540
%       [network(i).Da1] = DaOUT(i,1);%1.7m3/s
%       [network(i).Da2] = DaOUT(i,2);
%       [network(i).Da3] = DaOUT(i,3);
%       [network(i).Da4] = DaOUT(i,4);%160m3/s
% end
% for i = 541:LinkNum
%       [wetland(i-540).Da1] = DaOUT(i,1);%1.7m3/s
%       [wetland(i-540).Da2] = DaOUT(i,2);
%       [wetland(i-540).Da3] = DaOUT(i,3);
%       [wetland(i-540).Da4] = DaOUT(i,4);%160m3/s
% end
% shapewrite(network,'Shapefiles\LS_channels2_utm_Da.shp');
% shapewrite(wetland,'Shapefiles\LS_wetlands2_Da.shp');

%% write json; initially for setting file structure
% paeff=fainN.*100;
% paeff(isnan(paeff))=0;
% paw=WA./Area.*100;%percentage of area that is a wetland/lake
% paw(isnan(paw))=0;
% mods=cat(2,GridID,paeff,Feature,paw,pEM);
% data2json=struct('flow','L','modifications',mods);
% savejson('',data2json,strcat('modifications_',runID,'.json'));

%% read json
% if ~ischar(runID)
%     runID=num2str(runID);
% end
% dat=loadjson(strcat('modifications_',runID,'.json'));
% flow=dat.flow;
% fainN=dat.modifications(:,2)./100;
% Feature=dat.modifications(:,3);
% WA=dat.modifications(:,4).*Area./100;
% pEM=dat.modifications(:,5);

%% Assign Q
% switch flow
%     case 'H'
%         Qgage=400;%m3/s
%     case 'M'
%         Qgage=40;%m3/s
%     case 'L'
%         Qgage=4;%m3/s
%     otherwise
%         disp('flow undefined')
% end

%Qgage=QQ(qq);
Qgage=160;%160;%54;%6.4;%1.7;%m3/s
%Qgage=gageDMYQ(dmyq,4);
%Q=gageDMYQ(dmyq,4)./usarea(472,1).*usarea;
%Q=400./usarea(472,1).*usarea;
Q=Qgage./usarea(472,1).*usarea;

%% Assign Q to ensure continuity from tributaries
hw=(sum(Sub,1)==1);%identifies headwater, first-order links
Q(~hw)=NaN;
q=Q;

tocomp=sum(isnan(Q));
while tocomp>0
    for i=1:LinkNum

        if isnan(Q(i)) &&  ~isnan(sum(Q(Connect(:,2)==i)))
            %add incremenal Q via Area and upstream Qs
            q(i)=Qgage./usarea(472,1).*Area(i);
            Q(i)=q(i)+sum(Q(Connect(:,2)==i));
        end
        
    end
    tocomp=sum(isnan(Q));
end

%% Assign B
%assign B at gage based on rating curve
B(1:LinkNum,1)=NaN;
if Q(472)<Qbf
    B(472)=a1.*Q(472).^b1;
else
    B(472)=a2.*Q(472).^b2;
end

%scale B throughout the basin
B=(B(472)./(usarea(472).^(0.5))).*usarea.^(0.5);

%% Determine U, H, and wetland hydraulics
n=0.035; %Manning's roughness
%n=nval(nn);
g=9.81; %m/s2 - acceleration due to gravity
% not used for anything at the moment except for a check on conditions
U=(1/n.*(Q./B).^(2/3).*Slope.^(1/2)).^(3/5);
H=Q./U./B;
Fr=U./sqrt(g.*H);

%wetland volume, also not used for anything at the moment
WV(1:LinkNum,1)=0;%channel
WV(Feature==2,1)=2.1*WA(Feature==2,1);%m3 from Bevis MS thesis, lake
WV(Feature==3,1)=0.0032*WA(Feature==3,1).^1.47;%m3 Gleason 2007, wetland

%pEM(Feature==2)=20;
%pEM(Feature==3)=80;
%pEM(:,:)=100;

% Feature defined as channel=1, lake=2, wetland=3
% Lake wetland break based on palustrine fraction > 50%
% pP(1:LinkNum,1)=0;
% pP(541:LinkNum,1)=[wetland.pP].';
% clear Feature
% Feature(1:LinkNum,1)=1;%channel
% Feature(541:LinkNum,1)=2;%lake
% Feature(pP>50)=3;%wetland

Length(Feature>1,1)=sqrt(WA(Feature>1,1));%wetland length
B(Feature>1,1)=sqrt(WA(Feature>1,1));%wetland width
H(Feature>1,1)=WV(Feature>1,1)./WA(Feature>1,1);%wetland depth
U=Q./B./H;%compute velocity

%% N and C conc - initially for headwater links
%fainN(fainN<0)=0;
N_conc_ri(1:LinkNum,1)=30.*fainN;%mg/L
N_conc_ri(isnan(N_conc_ri))=0;
N_conc_us=N_conc_ri;

%C_conc_ri(1:LinkNum,1)=cval(cc);%mg/L
%Jleach=jval(jj)/3600; %mg/m2/s
Jleach=85/3600; %mg/m2/s
%C_conc_ri(1:LinkNum,1)=7;%mg/L
%C_conc_ri(1:LinkNum,1)=200.*fainC+3.2.*fainN;%mg/L
%C_conc_ri(1:LinkNum,1)=cval(cc).*fainC+4.5.*fainN;%mg/L
C_conc_ri(1:LinkNum,1)=90.*fainC+4.5.*fainN;%mg/L

%Jleach=60/3600; %mg/m2/s
%C_conc_us=C_conc_ri;
C_conc_us=C_conc_ri+(Jleach.*WA.*pEM./100./Q./1000);%C generation added here

CNrat=C_conc_us./N_conc_us;
CNrat(isnan(CNrat))=1;
%compute denitrification rate for C:N>=1.5
%nitrate limited
Jden(CNrat>=1,1)=(11.5.*N_conc_us(CNrat>=1).^(0.5))./3600;%mg/m2/s of N or C
%compute denitrification rate for C:N<1.5
%carbon limited
Jden(CNrat<1,1)=(3.5.*C_conc_us(CNrat<1))./3600;%mg/m2/s of N or C
%Jden(CNrat<1,1)=(dcval(dd).*C_conc_us(CNrat<1))./3600;%mg/m2/s of N or C
%Jden(CNrat<1,1)=(11.5.*C_conc_us(CNrat<1).^(0.5))./3600;%mg/m2/s of N or C

N_conc_ds=N_conc_us-(Jden.*B.*Length./Q./1000);%mg/L, nitrate conc at ds end of link
N_conc_ds(N_conc_ds<0)=0;

%C_conc_ds=C_conc_us-(Jden.*B.*Length./Q./1000)+(Jleach.*WA.*pEM./100./Q./1000);
C_conc_ds=C_conc_us-(Jden.*B.*Length./Q./1000); %mg/L, C conc at ds end of link
C_conc_ds(C_conc_ds<0)=0;

%% N and C conc for all links ensuring continuity at tributaries
hw=(sum(Sub,1)==1);%identifies headwater, first-order links
N_conc_us(~hw)=NaN;%mg/L
N_conc_ds(~hw)=NaN;%mg/L
Jden(~hw)=NaN;%mg/m2/s
C_conc_us(~hw)=NaN;%mg/L
C_conc_ds(~hw)=NaN;%mg/L

tocomp=sum(isnan(N_conc_ds));
while tocomp>0
    for i=1:LinkNum

        if isnan(N_conc_ds(i)) &&  ~isnan(sum(N_conc_ds(Connect(:,2)==i)))

            N_conc_us(i)=(N_conc_ri(i).*q(i)+...
                sum(N_conc_ds(Connect(:,2)==i).*Q(Connect(:,2)==i)))./Q(i); 
            C_conc_us(i)=(C_conc_ri(i).*q(i)+...
                sum(C_conc_ds(Connect(:,2)==i).*Q(Connect(:,2)==i)))./Q(i);
            %add C generation here
            C_conc_us(i)=C_conc_us(i)+(Jleach.*WA(i).*pEM(i)./100./Q(i)./1000);
            
            CNrat(i)=C_conc_us(i)./N_conc_us(i);
            if isnan(CNrat(i))
                CNrat(i)=1;
            end
            if CNrat(i)>=1
                Jden(i)=(11.5.*N_conc_us(i).^(0.5))./3600;%mg/m2/s of N or C
            else
                Jden(i)=(3.5.*C_conc_us(i))./3600;%mg/m2/s of N or C
                %Jden(i)=(dcval(dd).*C_conc_us(i))./3600;%mg/m2/s of N or C
                %Jden(i)=(11.5.*C_conc_us(i).^(0.5))./3600;%mg/m2/s of N or C
            end
            
            N_conc_ds(i)=N_conc_us(i)-(Jden(i).*B(i).*Length(i)./Q(i)./1000);
            N_conc_ds(N_conc_ds<0)=0;
            
            %C_conc_ds(i)=C_conc_us(i)-(Jden(i).*B(i).*Length(i)./Q(i)./1000)+...
            %    (Jleach.*WA(i).*pEM(i)./100./Q(i)./1000);
            C_conc_ds(i)=C_conc_us(i)-(Jden(i).*B(i).*Length(i)./Q(i)./1000);
            C_conc_ds(C_conc_ds<0)=0;

        end
        
    end
    tocomp=sum(isnan(N_conc_ds));
end

%% Compute additional values for output
N_load_ds=N_conc_ds.*Q.*1000./1000./1000;%kg/s
C_load_ds=C_conc_ds.*Q.*1000./1000./1000;%kg/s
%outlet_N_load=N_load_ds(392);%kg/s, outlet for Education extent
outlet_N_load=N_load_ds(OutletLinkID);%kg/s
netavg_N_conc=sum(N_conc_ds.*Length)./sum(Length);%mg/L
%%

%determine links DS of in-channel wetlands
DSwet(1:LinkNum,1)=0;
for i=1:LinkNum
    dsflag=sum(Connect(i,~isnan(Connect(i,:)))>540);
    if dsflag==0
        DSwet(i,1)=1;
    end
end
%determine links US of in-channel wetlands
USwet(1:LinkNum,1)=0;
USwet(DSwet==0,1)=1;
USwet(Feature>1,1)=0;

netavg_N_conc_wet=sum(N_conc_ds(Feature>1).*Length(Feature>1))./sum(Length(Feature>1));%mg/L
netavg_N_conc_DSwet=sum(N_conc_ds(DSwet>0).*Length(DSwet>0))./sum(Length(DSwet>0));%mg/L
netavg_N_conc_USwet=sum(N_conc_ds(USwet>0).*Length(USwet>0))./sum(Length(USwet>0));%mg/L

%% Additional values used in analyzing results
%TT=Length.*B.*H./Q;%seconds./60./60./24;%days , travel time
prem=(N_conc_us-N_conc_ds)./N_conc_us.*100;% percent N removal in a link
prem(isnan(prem))=0;
frem=prem./100;% fraction N removal in a link
ftrans=1-frem;% fraction N transported through a link
%%
ftransout(1:LinkNum,1)=NaN;
for i=1:LinkNum
    ftransout(i,1)=prod(ftrans(Connect(i,1:find(~isnan(Connect(i,:)),1,'last'))));
end
%%
%netavg_frem=sum(frem.*Length)./sum(Length);
%fremout=1-ftransout;
%netavg_fremout=sum(fremout.*Length)./sum(Length);
Nmassin=q.*N_conc_ri.*1000./1000./1000;%kg/s
total_Nmassin=sum(Nmassin);
totalfracNret=1-(outlet_N_load./total_Nmassin);

%%
%Nlim(1:LinkNum,1)=0;%0 C limited
%Nlim(CNrat>=1,1)=1;%1 N limited, 
%Nlim(prem<50,1)=2;%2 H limited

%% Damkohler Number
Da=-log(N_conc_ds./N_conc_us);
Da(isinf(Da))=100;%infinite
Da(isnan(Da))=-100;%undefined

%%
CNlim(1:LinkNum,1)=0;%0 Biogeochemically limited by C
CNlim(CNrat>=1,1)=1;%1 Biogeochemically limited by N

BQlim(1:LinkNum,1)=1;%1 Limited by flow
BQlim(Da>=1,1)=0;%0 Limited by biogeochemistry

CNQlim(1:LinkNum,1)=0;%0 Limited by C
CNQlim(CNrat>=1,1)=1;%1 Limited by N
CNQlim(Da<1,1)=2;%2 Limited by flow

%%
netavg_frac_CN_Clim=sum(Length(CNlim==0))./sum(Length);%fraction biogeochem. C lim
netavg_frac_BQ_Qlim=sum(Length(BQlim==1))./sum(Length);%fraction flow lim

netavg_frac_CNQ_Qlim=sum(Length(CNQlim==2))./sum(Length);%fraction flow lim
netavg_frac_CNQ_Clim=sum(Length(CNQlim==0))./sum(Length);%fraction C lim

%%

% OUT(qq,1)=netavg_N_conc;
% OUT(qq,2)=outlet_N_load;
% OUT(qq,3)=totalfracNret;
% OUT(qq,4)=netavg_frac_CN_Clim;
% OUT(qq,5)=netavg_frac_BQ_Qlim;
% OUT(qq,6)=netavg_frac_CNQ_Qlim;
% OUT(qq,7)=netavg_frac_CNQ_Clim;
% 
% end
%% Plot Basin Summary Figs
% i=5;
% 
% figure; hold on; box on
% set(gca,'XScale','Log')
% %set(gca,'YScale','Log')
% plot(QQ,OUT_Base(:,i),'k')
% plot(QQ,OUT_pNin(:,i),'-r')
% plot(QQ,OUT_mNin(:,i),'-b')
% plot(QQ,OUT_pCinag(:,i),'--r')
% plot(QQ,OUT_mCinag(:,i),'--b')
% plot(QQ,OUT_pCinnag(:,i),':r')
% plot(QQ,OUT_mCinnag(:,i),':b')
% plot(QQ,OUT_pJleach(:,i),'-.r')
% plot(QQ,OUT_mJleach(:,i),'-.b')
% plot(QQ,OUT_pCdc(:,i),'m')
% plot(QQ,OUT_mCdc(:,i),'g')
% ylim([0 1])
% %ylim([0 20])
% %ylim([10^-6 10])
% xlim([1 400])
% 
% plot([1.7 1.7],[0 20],'k')
% plot([6.4 6.4],[0 20],'k')
% plot([160 160],[0 20],'k')
% plot([1.7 1.7],[10^-6 10],'k')
% plot([6.4 6.4],[10^-6 10],'k')
% plot([160 160],[10^-6 10],'k')
% plot([1.7 1.7],[0 1],'k')
% plot([6.4 6.4],[0 1],'k')
% plot([160 160],[0 1],'k')
% 
% %%
% 
% figure; hold on; box on
% set(gca,'XScale','Log')
% plot(QQ,OUT_Base(:,6),'k')
% plot(QQ,OUT_Base(:,6)+OUT_Base(:,7),'b')
% 
% ylim([0 1])
% xlim([1 400])
% 
% plot([1.7 1.7],[0 1],'k')
% plot([160 160],[0 1],'k')

%%
% X=Jden.*B.*Length./Q./1000./N_conc_us*100;
% Y=prem;
% figure; hold on; box on
% plot(X,Y,'.k')
% plot(X(Nlim==1),Y(Nlim==1),'.g')
% %plot(X(Feature==1),Y(Feature==1),'.k')
% %plot(X(Feature==2),Y(Feature==2),'.b')
% %plot(X(Feature==3),Y(Feature==3),'.g')
% set(gca,'XScale','log','YScale','log')


%% write json
% clear data2json
% sres=cat(2,outlet_N_load,netavg_N_conc);
% mres=cat(2,GridID,N_conc_ds,C_conc_ds,N_load_ds,C_load_ds);
% data2json=struct('summaryresults',sres,'modelresults',mres);
% savejson('',data2json,strcat('results_',runID,'.json'));

%% Extract model results for comparison to data
% idx=((NDay==gageDMYQ(dmyq,1))&(NMonth==gageDMYQ(dmyq,2))&(NYear==gageDMYQ(dmyq,3)));
% VAL(idx,2)=N_conc_ds(sGridID(idx));
% VAL(idx,3)=C_conc_ds(sGridID(idx));
% 
% VAL(~CALidx,:)=NaN;%toss out non cal values for calibration
% %VAL(logical(CALidx),:)=NaN;%toss out non val values for validation
% 
% idx2=((QDay==gageDMYQ(dmyq,1))&(QMonth==gageDMYQ(dmyq,2))&(QYear==gageDMYQ(dmyq,3)));
% QVAL(idx2,1)=Q(QsGridID(idx2));
% QVAL(idx2,2)=U(QsGridID(idx2));
% QVAL(idx2,3)=B(QsGridID(idx2));
% QVAL(idx2,4)=H(QsGridID(idx2));
% 
% end

%%
%CALidx=round(rand(108,1));
% sum(~isnan(NO3(logical(CALidx))))
% sum(~isnan(NO3(logical(~CALidx))))
% sum(~isnan(DOC(logical(CALidx))))
% sum(~isnan(DOC(logical(~CALidx))))

%% C/N verification
% CN=DOC./NO3;
% VAL(:,1)=VAL(:,3)./VAL(:,2);
% 
% TPidx=and(CN>=1,VAL(:,1)>=1);
% TNidx=and(CN<1,VAL(:,1)<1);
% FNidx=and(CN>=1,VAL(:,1)<1);
% FPidx=and(CN<1,VAL(:,1)>=1);
% 
% TP=sum(TPidx);
% TN=sum(TNidx);
% FN=sum(FNidx);
% FP=sum(FPidx);
% 
% acc=(TP+TN)./(TP+TN+FP+FN).*100;% in percent
% 
% mcc=((TP.*TN)-(FP*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));

%% Sensitivity on n
% RMSE(nn,1)=sqrt(nanmean((gU-QVAL(:,2)).^2));
% RMSE(nn,2)=sqrt(nanmean((gH-QVAL(:,4)).^2));
% RMSE(nn,3)=sqrt(nanmean((NO3-VAL(:,2)).^2));
% RMSE(nn,4)=sqrt(nanmean((DOC-VAL(:,3)).^2));

%% Sensitivity on Cin and Jleach
% NRMSE(cc,jj,dd)=sqrt(nanmean((BVAL(:,2)-VAL(:,2)).^2));
% CRMSE(cc,jj,dd)=sqrt(nanmean((BVAL(:,3)-VAL(:,3)).^2));


% NRMSE(cc,jj,dd)=sqrt(nanmean((NO3-VAL(:,2)).^2));
% CRMSE(cc,jj,dd)=sqrt(nanmean((DOC-VAL(:,3)).^2));
% ACC(cc,jj,dd)=acc;
% MCC(cc,jj,dd)=mcc;
% end
% end
% end

%% Calibration variability
% %% MCC
% x=cval;
% y=jval;
% z=MCC(:,:,dd);
% figure; box on; hold on
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(cval(cc),jval(jj),'ok','MarkerFaceColor','r')
% xlabel('nominal C input from vegetated lands, mg/L')
% ylabel('DOC leaching rate, mg/m2/hr')
% colorbar
% caxis([nanmin(nanmin(nanmin(MCC))), nanmax(nanmax(nanmax(MCC)))])
% 
% x=cval;
% y=dcval;
% z=[];
% for i=1:size(MCC,3)
% z=cat(2,z,MCC(:,jj,i));
% end
% figure; box on; hold on
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(cval(cc),dcval(dd),'ok','MarkerFaceColor','r')
% xlabel('nominal C input from vegetated lands, mg/L')
% ylabel('C-limited denitrif. coeff.')
% colorbar
% caxis([nanmin(nanmin(nanmin(MCC))), nanmax(nanmax(nanmax(MCC)))])
% 
% x=dcval;
% y=jval;
% z=[];
% for i=1:size(MCC,3)
% z=cat(1,z,MCC(cc,:,i));
% end
% figure; box on; hold on
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(dcval(dd),jval(jj),'ok','MarkerFaceColor','r')
% xlabel('C-limited denitrif. coeff.')
% ylabel('DOC leaching rate, mg/m2/hr')
% colorbar
% caxis([nanmin(nanmin(nanmin(MCC))), nanmax(nanmax(nanmax(MCC)))])
% %% NRMSE
% x=cval;
% y=jval;
% z=NRMSE(:,:,dd);
% figure; box on; hold on
% colormap(flipud(parula))
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(cval(cc),jval(jj),'ok','MarkerFaceColor','r')
% xlabel('nominal C input from vegetated lands, mg/L')
% ylabel('DOC leaching rate, mg/m2/hr')
% colorbar
% caxis([nanmin(nanmin(nanmin(NRMSE))), nanmax(nanmax(nanmax(NRMSE)))])
% 
% x=cval;
% y=dcval;
% z=[];
% for i=1:size(NRMSE,3)
% z=cat(2,z,NRMSE(:,jj,i));
% end
% figure; box on; hold on
% colormap(flipud(parula))
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(cval(cc),dcval(dd),'ok','MarkerFaceColor','r')
% xlabel('nominal C input from vegetated lands, mg/L')
% ylabel('C-limited denitrif. coeff.')
% colorbar
% caxis([nanmin(nanmin(nanmin(NRMSE))), nanmax(nanmax(nanmax(NRMSE)))])
% 
% x=dcval;
% y=jval;
% z=[];
% for i=1:size(NRMSE,3)
% z=cat(1,z,NRMSE(cc,:,i));
% end
% figure; box on; hold on
% colormap(flipud(parula))
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(dcval(dd),jval(jj),'ok','MarkerFaceColor','r')
% xlabel('C-limited denitrif. coeff.')
% ylabel('DOC leaching rate, mg/m2/hr')
% colorbar
% caxis([nanmin(nanmin(nanmin(NRMSE))), nanmax(nanmax(nanmax(NRMSE)))])
% %% CRMSE
% x=cval;
% y=jval;
% z=CRMSE(:,:,dd);
% figure; box on; hold on
% colormap(flipud(parula))
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(cval(cc),jval(jj),'ok','MarkerFaceColor','r')
% xlabel('nominal C input from vegetated lands, mg/L')
% ylabel('DOC leaching rate, mg/m2/hr')
% colorbar
% caxis([nanmin(nanmin(nanmin(CRMSE))), nanmax(nanmax(nanmax(CRMSE)))])
% 
% x=cval;
% y=dcval;
% z=[];
% for i=1:size(CRMSE,3)
% z=cat(2,z,CRMSE(:,jj,i));
% end
% figure; box on; hold on
% colormap(flipud(parula))
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(cval(cc),dcval(dd),'ok','MarkerFaceColor','r')
% xlabel('nominal C input from vegetated lands, mg/L')
% ylabel('C-limited denitrif. coeff.')
% colorbar
% caxis([nanmin(nanmin(nanmin(CRMSE))), nanmax(nanmax(nanmax(CRMSE)))])
% 
% x=dcval;
% y=jval;
% z=[];
% for i=1:size(CRMSE,3)
% z=cat(1,z,CRMSE(cc,:,i));
% end
% figure; box on; hold on
% colormap(flipud(parula))
% contourf(repmat(x',1,length(y)),repmat(y,length(x),1),z,'LineStyle','none')
% plot(dcval(dd),jval(jj),'ok','MarkerFaceColor','r')
% xlabel('C-limited denitrif. coeff.')
% ylabel('DOC leaching rate, mg/m2/hr')
% colorbar
% caxis([nanmin(nanmin(nanmin(CRMSE))), nanmax(nanmax(nanmax(CRMSE)))])
%%

% figure; hold on; box on
% plot(nval,RMSE(:,1),'b')
% plot(nval,RMSE(:,2),'r')
% legend('U','H')

% figure; hold on; box on
% contourf(repmat(cval',1,length(jval)), repmat(jval,length(cval),1),CRMSE)

%% J value
% C=0:0.1:30;%mg/L
% JdenNlim=11.5.*C.^(0.5);%mg/m2/hr of N or C
% JdenClim=3.5.*C;%mg/m2/hr of N or C
% figure
% hold on
% plot(C,JdenNlim,'g')
% plot(C,JdenClim,'k')
% xlabel('Concentration, mg/L')
% ylabel('Denitrification rate, mg/m2/hr')
%%
% CN=DOC./NO3;
% JdenNlim(1:size(NO3,1),1:size(NO3,2))=NaN;
% JdenNlim(CN>=1.5)=11.5.*NO3(CN>=1.5).^(0.5);%mg/m2/hr of N or C
% JdenClim(1:size(DOC,1),1:size(DOC,2))=NaN;
% JdenClim(CN<1.5)=3.*DOC(CN<1.5);%mg/m2/hr of N or C
%%
% figure
% hold on
% plot(NO3,JdenNlim,'.g')
% plot(DOC,JdenClim,'.k')
%%
% sum(sum(~isnan(JdenClim)))
% nanmean(nanmean(JdenClim))
%%
% nanmean(nanmean(NO3(CN<1.5)))

%%
% for i=1:113
%     if ~isempty(sGridID(find(sS(i)==sSiteID,1,'first')))
%         sG(i,1)=sGridID(find(sS(i)==sSiteID,1,'first'));
%     end
% end


%%

% NN(:,qq)=N_conc_ds;
% CC(:,qq)=C_conc_ds;
% CN(:,qq)=CNrat;
% DA(:,qq)=Da;
% 
% end
%%
% % for wet=541:LinkNum
% %     close all
%     %wet=143;%channel
%     %wet=541;%wetland
%     wet=560;%lake
%     figure; hold on; box on
%     plot(QQ,NN(wet,:),'g')
%     plot(QQ,CC(wet,:),'k')
%     ylim([0 30])
%     xlim([1 400])
%     set(gca,'XScale','log')
%     %set(gcf,'Position',[141 183 560 420])
%     
%     %plot([1.7 1.7],[0 30],'k')
%     plot([6.4 6.4],[0 30],'k')
%     %plot([54 54],[0 30],'k')
%     plot([160 160],[0 30],'k')
%     
%     
%     figure; hold on; box on
%     plot(QQ,CN(wet,:),'g')
%     %plot([1 400],[1.5 1.5],'-.g')
%     plot(QQ,DA(wet,:),'r')
%     plot([1 400],[1 1],'-k')
%     set(gca,'YScale','log')
%     ylim([10^-4 10^2])
%     xlim([1 400])
%     set(gca,'XScale','log')
%     %set(gcf,'Position',[705 184 560 420])
%     %title(num2str(wet))
%     %pEM(wet)
% %     
% %     pause
% % end

%% plot long profile of given value
% %loc=563;%left
% loc=585;%middle
% %loc=124;%right
% %loc=xwet;
% 
% index=Connect(loc,~isnan(Connect(loc,:)));
% 
% xx=[];
% yy=[];
% x=cumsum(Length(index))./1000;%km
% y=C_conc_ds(index);%mg/L
% % y=usarea_km2(index);
% % yy(1)=(y(1)-0)./y(1).*100;
% % for i=2:length(y)
% %    yy(i)=(y(i)-y(i-1))./y(i-1).*100;
% % end
% %xx=cat(1,0,x);
% %yy=cat(1,C_conc_us(loc),y);
% xx=x;
% yy=y;
% %y=CNrat(index);
% %y=Da(index);
% %y=GridID(index);
% %y=Jden(index).*3600;%mg/m2/hr
% 
% % xx=cat(1,xx,0,x(1));
% % yy=cat(1,yy,y(1),y(1));
% % for i=2:length(x)
% %     xx=cat(1,xx,x(i-1:i,1));
% %     yy=cat(1,yy,y(i),y(i));
% % end
% 
% %figure; hold on; box on
% plot(xx,yy,':b');%'r')%'y')%'g')%'b')
% %ylim([0 30])
% %ylim([10^-2 10^3])
% %plot([0 140],[1 1],'k')
% 
% %XX(1:length(xx),xwet-540)=xx;
% %YYB1(1:length(yy),xwet-540)=yy;
% %YYB2(1:length(yy),xwet-540)=yy;
% %YYB3(1:length(yy),xwet-540)=yy;
% %YYB4(1:length(yy),xwet-540)=yy;
% %YYR1(1:length(yy),xwet-540)=yy;
% %YYR2(1:length(yy),xwet-540)=yy;
% %YYR3(1:length(yy),xwet-540)=yy;
% %YYR4(1:length(yy),xwet-540)=yy;
% 
% %XX(1:length(xx),qq)=xx;
% %YYBn(1:length(yy),qq)=yy;
% %YYBc(1:length(yy),qq)=yy;
% %YYBr(1:length(yy),qq)=yy;
% %YYBd(1:length(yy),qq)=yy;
% %YYRn(1:length(yy),qq)=yy;
% %YYRc(1:length(yy),qq)=yy;
% %YYRr(1:length(yy),qq)=yy;
% %YYRd(1:length(yy),qq)=yy;
% 
% %end

%%
% %close all
% xwet=580;
% figure; hold on; box on;
% plot(XX(:,xwet-540),YYB1(:,xwet-540),'-.b')
% plot(XX(:,xwet-540),YYB2(:,xwet-540),'-.c')
% plot(XX(:,xwet-540),YYB3(:,xwet-540),'-.g')
% plot(XX(:,xwet-540),YYB4(:,xwet-540),'-.r')
% plot(XX(:,xwet-540),YYR1(:,xwet-540),'b')
% plot(XX(:,xwet-540),YYR2(:,xwet-540),'c')
% plot(XX(:,xwet-540),YYR3(:,xwet-540),'g')
% plot(XX(:,xwet-540),YYR4(:,xwet-540),'r')
% %%
% D1=YYR1-YYB1;
% D2=YYR2-YYB2;
% D3=YYR3-YYB3;
% D4=YYR4-YYB4;
% %%
% % for xwet=541:LinkNum
% % close all
% % %xwet=613;
% figure; hold on; box on;
% plot(XX(:,xwet-540),D1(:,xwet-540),'b')
% plot(XX(:,xwet-540),D2(:,xwet-540),'c')
% plot(XX(:,xwet-540),D3(:,xwet-540),'g')
% plot(XX(:,xwet-540),D4(:,xwet-540),'r')
% title(num2str(xwet))
% % pause
% % end
% 
% %%
% %d2max=nanmax(D2,[],1);
% d2dist(1,1:103)=NaN;
% for xwet=541:LinkNum
%     idx=find(D2(:,xwet-540)>=0.1,1,'last');
%     if isempty(idx)
%         d2dist(1,xwet-540)=NaN;
%     else
%         d2dist(1,xwet-540)=XX(idx,xwet-540);
%     end
% end
% 
% figure; hold on; box on
% plot(WA(541:end)./usarea(541:end),d2dist,'.c')

%%

% figure; hold on; box on;
% plot(XX(:,1),YYBd(:,1),'-.b')
% plot(XX(:,2),YYBd(:,2),'-.c')
% plot(XX(:,3),YYBd(:,3),'-.g')
% plot(XX(:,4),YYBd(:,4),'-.r')
% plot(XX(:,1),YYRd(:,1),'b')
% plot(XX(:,2),YYRd(:,2),'c')
% plot(XX(:,3),YYRd(:,3),'g')
% plot(XX(:,4),YYRd(:,4),'r')
% 
% plot([0 140],[1 1],'k')
% %%
% Dn=YYRn-YYBn;
% Dc=YYRc-YYBc;
% Dr=YYRr-YYBr;
% Dd=YYRd-YYBd;
% %%
% % for xwet=541:LinkNum
% % close all
% % %xwet=613;
% figure; hold on; box on;
% plot(XX(:,1),Dc(:,1),'b')
% plot(XX(:,2),Dc(:,2),'c')
% plot(XX(:,3),Dc(:,3),'g')
% plot(XX(:,4),Dc(:,4),'r')
% 
% plot([0 140],[0 0],'k')

%% compare DOC vs wetlands

% wet = shaperead('intersec1.shp');
% %%
% em(1:11557,1)=0;
% wa(1:11557,1)=NaN;
% site(1:11557,1)=NaN;
% sua(1:11557,1)=NaN;
% 
% for i=1:11557
%    em(i,1)=strcmp(wet(i).Att3,'EM');
%    wa(i,1)=wet(i).WA_m2;
%    site(i,1)=str2num(wet(i).Name);
%    sua(i,1)=wet(i).Shape_Ar_1;    
% end
% %%
% pwa=wa./sua;%percent wetland area
% 
% pWetland(1:108,1)=NaN;
% pEM(1:108,1)=NaN;
% 
% for i=1:108
%     sid=(sSiteID(i)==site);
%     pWetland(i,1)=nansum(pwa(sid));
%     pEM(i,1)=nansum(pwa(and(sid,em)));
% end
% %%
% figure; hold on; box on
% plot(pWetland.*100,DOC,'.k')
% xlabel('%wetland')
% ylabel('DOC, mg/L')
% 
% figure; hold on; box on
% plot(pEM.*100,DOC,'.k')
% xlabel('%EM')
% ylabel('DOC, mg/L')
% %%
% figure; hold on; box on
% plot(pWetland.*100,C_conc_ds(sGridID),'.k')
% xlabel('%wetland')
% ylabel('DOC, mg/L')
% 
% figure; hold on; box on
% plot(pEM.*100,C_conc_ds(sGridID),'.k')
% xlabel('%em')
% ylabel('DOC, mg/L')
%% PLOT

%
% % 
%% Plot Spatial Data
%set timestep to plot
t = 1;

%concN(isnan(concN))=0;
%nttl = 'N concentration';
%val=1-(aggrem./usarea);
%val=1-(aggrem./usefflandarea);
%val=1-rem;
val=N_conc_ds;
%val=ftransout;
%val=Nlim;
%val=pEM;
%val=Da;

for i = 1:540
      [network(i).val] = val(i);
end
for i = 541:LinkNum
      [wetland(i-540).val] = val(i);
end

f1 = figure;
%set(f1,'Position',[213 50 938 632]);
%a1 = axes('Position',[0.08 0.1 0.4 0.8]);
a1 = axes;
%cbpos = [0.78 0.22 0.01 0.6];%mrb,hw,c
cbpos = [0.90 0.22 0.01 0.6];%be
%cbpos = [0.93 0.22 0.01 0.6];%w
box(a1);

%edge = [0 1 2 4 8 16 32 64 128 200];
%edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
%edge = [0 1 2 4 6 8 10 12 16 20];
%edge = [0 1 2 3 4 5 10 20 50 100];
%edge = [-100 0 0.01 0.1 0.5 1 1.5 2 5 100];
%edge = [0 10 30 40 50 60 70 80 90 100];
%edge = [0 .2 .3 .4 .5 .6 .7 .8 .9 1];
%edge = [0 .01 .02 .03 .04 .05 .1 .2 .5 1];
%edge = [0 1 10 50 100 500 1000 5000 10000 100000];
%edge = [0 1 10 50 70 100 200 300 500 1000];
%edge = [0 .5 .75 1 1.25 1.5 1.75 2 2.25 2.5];
%edge = cat(2,0,logspace(-2,0,9));
%edge = cat(2,0,round(exp(linspace(log(1),log(max([network.sed])),9))));
%edge = cat(2,0,(exp(linspace(log(min([network.val])),log(max([network.val])),9))));
edge = cat(2,0,(exp(linspace(log(1),log(max([network.val])),9))));

%edge = cat(2,0,((linspace((min([network.val])),(max([network.val])),9))));
%edge = cat(2,0,(exp(linspace(log(1),log(max(max(LCS))),9))));
%edge = cat(2,0,round(exp(linspace(log(1),log(max(max(numpar))),9))));
%edge = cat(2,0,round(linspace(mn,ceil(max(max(LCS))),9)));
%edge = cat(2,0,round(exp(linspace(log(mn),log(ceil(max(max(LCS)))),9))));

%edge = cat(2,0,(exp(linspace(log(0.1),log(max([network.Nconc])),9))));
%edge = cat(2,0,round((linspace((min([network.Nconc])),(max([network.Nconc])),9))));
%edge = cat(2,0,round((linspace((1),(max(max(numpar))),9))));
%edge = cat(2,0,round(quantile([network.Nconc],0.1:0.1:0.8)),37);%round(max([network.Nconc])));
%edge = cat(2,0,round(linspace(1,100,9)));
%edge = [0, 1, 5, 10, 15, 20, 23, 26, 30, 37];
%edge = [0, 1, 5, 10, 15, 10, 13, 15, 18, 21];
%edge = [0, .1, .25, .4, .5, .6, .7, .8, .9, 1];
%edge = [0, 3, 5, 8, 10, 12, 15, 18, 20, 25];

% ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];nttl};
%ttl={['Time = ',num2str(fix(time(t))),' days ',...
%    num2str(round(mod(time(t)*24,24))),' hours'];['Time step = ',num2str(t)];nttl};
%ttl=''; 
%[a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);

%color=cat(1,[0.85 0.85 0.85],winter(8));
%color=cat(1,[0 0 0],[0 1 0],[0 0 1],repmat([0.85 0.85 0.85],6,1));
color=cat(1,[0.85 0.85 0.85],jet(7),[0 0 0]);
%color=jet(9);
%color=lines(9);
colormap(color)
caxis([0 90])
% cb=colorbar('YTick', linspace(0,90,10),...
%     'YTickLabel',{'0','1','2','4','8','16','32','64','128','200'},...
%     'Position',[0.5 0.2 0.01 0.6]);
%edge=edge./10000;

cb=colorbar('YTick', linspace(0,90,10),...
    'YTickLabel',{num2str(edge(1)),num2str(edge(2)),num2str(edge(3)),...
    num2str(edge(4)),num2str(edge(5)),num2str(edge(6)),num2str(edge(7)),...
    num2str(edge(8)),num2str(edge(9)),num2str(edge(10))},...
    'Position',cbpos);

%convert m to km
% cb=colorbar('YTick', linspace(0,90,10),...
%     'YTickLabel',{num2str(round(edge(1)./1000)),num2str((edge(2)./1000)),num2str(round(edge(3)./1000)),...
%     num2str(round(edge(4)./1000)),num2str(round(edge(5)./1000)),num2str(round(edge(6)./1000)),num2str(round(edge(7)./1000)),...
%     num2str(round(edge(8)./1000)),num2str(round(edge(9)./1000)),num2str(round(edge(10)./1000))},...
%     'Position',cbpos);

% cb=colorbar('Location','North','XTickMode','manual','XTick', linspace(0,90,10),...
%     'XTickLabel',{'0','1','2','4','8','16','32','64','128','200'},'Position',[0.1 0.9 0.3 0.01]);
%edge=edge.*10000;

netspec = makesymbolspec('Line',...
    {'Default', 'Color','black'}, ...
    {'val',[edge(1) edge(2)], 'Color',color(1,:), 'LineWidth', 1.5},...
    {'val',[edge(2) edge(3)], 'Color',color(2,:), 'LineWidth', 1.5},...
    {'val',[edge(3) edge(4)], 'Color',color(3,:), 'LineWidth', 1.5},...
    {'val',[edge(4) edge(5)], 'Color',color(4,:), 'LineWidth', 1.5},...
    {'val',[edge(5) edge(6)], 'Color',color(5,:), 'LineWidth', 1.5},...
    {'val',[edge(6) edge(7)], 'Color',color(6,:), 'LineWidth', 1.5},...
    {'val',[edge(7) edge(8)], 'Color',color(7,:), 'LineWidth', 1.5},...
    {'val',[edge(8) edge(9)], 'Color',color(8,:), 'LineWidth', 1.5},...
    {'val',[edge(9) edge(10)], 'Color',color(9,:), 'LineWidth', 1.5});

wetspec = makesymbolspec('Polygon',...
    {'Default', 'FaceColor','black'}, ...
    {'val',[edge(1) edge(2)], 'FaceColor',color(1,:)},...
    {'val',[edge(2) edge(3)], 'FaceColor',color(2,:)},...
    {'val',[edge(3) edge(4)], 'FaceColor',color(3,:)},...
    {'val',[edge(4) edge(5)], 'FaceColor',color(4,:)},...
    {'val',[edge(5) edge(6)], 'FaceColor',color(5,:)},...
    {'val',[edge(6) edge(7)], 'FaceColor',color(6,:)},...
    {'val',[edge(7) edge(8)], 'FaceColor',color(7,:)},...
    {'val',[edge(8) edge(9)], 'FaceColor',color(8,:)},...
    {'val',[edge(9) edge(10)], 'FaceColor',color(9,:)},...
    {'Default', 'EdgeColor','black'}, ...
    {'val',[edge(1) edge(2)], 'EdgeColor',color(1,:)},...
    {'val',[edge(2) edge(3)], 'EdgeColor',color(2,:)},...
    {'val',[edge(3) edge(4)], 'EdgeColor',color(3,:)},...
    {'val',[edge(4) edge(5)], 'EdgeColor',color(4,:)},...
    {'val',[edge(5) edge(6)], 'EdgeColor',color(5,:)},...
    {'val',[edge(6) edge(7)], 'EdgeColor',color(6,:)},...
    {'val',[edge(7) edge(8)], 'EdgeColor',color(7,:)},...
    {'val',[edge(8) edge(9)], 'EdgeColor',color(8,:)},...
    {'val',[edge(9) edge(10)], 'EdgeColor',color(9,:)});

a1b = mapshow(a1,boundary, 'FaceColor', 'none', 'EdgeColor', [0.4 0.4 0.4]);

xlabel(a1,'Easting in meters','FontSize',12);
ylabel(a1,'Northing in meters','FontSize',12);
%xlim(a1,[1.5e5 5e5]);
%ylim(a1,[4.75e6 5.15e6]);
hold(a1);

a1n = mapshow(a1,network,'SymbolSpec',netspec);
a1w = mapshow(a1,wetland,'SymbolSpec',wetspec);

%title(a1,ttl,'FontSize',14);

%% Verification Plotting
% 
% 
% %% plot NO3 vs NO3 with 1:1 line
% X=VAL(:,2);
% Y=NO3;
% dlab='NO3';
% 
% figure; hold on; box on
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% %mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% %mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% mec=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% mfc=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% %mfc=[1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((NDay==gageDMYQ(dmyq,1))&(NMonth==gageDMYQ(dmyq,2))&(NYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthWest')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% %set(gca,'XScale','log','YScale','log')
% 
% whichstats = {'beta','rsquare'};
% stats = regstats(Y,X,'linear',whichstats);
% a=(stats.beta(1));
% b=stats.beta(2);
% x=linspace(min(X),max(X),100);
% plot(x,a+b.*x,'-k');
% 
% plot([0 30],[0 30],'--k')
% %xlim([0 30])
% 
% RMSE=sqrt(nanmean((Y-X).^2));
% title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)
% 
% %% DOC 1:1
% X=VAL(:,3);
% Y=DOC;
% dlab='DOC';
% 
% figure; hold on; box on
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% %mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% %mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% mec=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% mfc=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% %mfc=[1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((NDay==gageDMYQ(dmyq,1))&(NMonth==gageDMYQ(dmyq,2))&(NYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthWest')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% %set(gca,'XScale','log','YScale','log')
% 
% whichstats = {'beta','rsquare'};
% stats = regstats(Y,X,'linear',whichstats);
% a=(stats.beta(1));
% b=stats.beta(2);
% x=linspace(min(X),max(X),100);
% plot(x,a+b.*x,'-k');
% %title(sprintf('Y = %.2g + {%.2f} * X\nR^2 = %.2f',[a b stats.rsquare]),'FontSize',14)
% 
% plot([0 30],[0 30],'--k')
% plot([0 60],[0 30],'--k')
% plot([0 60],[0 20],'--k')
% %plot([0 90],[0 30],'--k')
% %xlim([0 70])
% 
% RMSE=sqrt(nanmean((Y-X).^2));
% title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)
% 
% %% C/N 1:1
% X=VAL(:,1);
% Y=CN;
% dlab='DOC/NO3';
% 
% figure; hold on; box on
% set(gca,'YScale','Log','XScale','Log')
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% %mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% %mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% mec=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% mfc=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% %mfc=[1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((NDay==gageDMYQ(dmyq,1))&(NMonth==gageDMYQ(dmyq,2))&(NYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthWest')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% %set(gca,'XScale','log','YScale','log')
% 
% % whichstats = {'beta','rsquare'};
% % stats = regstats(Y,X,'linear',whichstats);
% % a=(stats.beta(1));
% % b=stats.beta(2);
% % x=linspace(min(X),max(X),100);
% % plot(x,a+b.*x,'-k');
% %title(sprintf('Y = %.2g + {%.2f} * X\nR^2 = %.2f',[a b stats.rsquare]),'FontSize',14)
% 
% plot([10^-2 10^2],[1 1],'k')
% plot([1 1],[10^-2 10^2],'k')
% plot([10^-2 10^2],[10^-2 10^2],':k')
% xlim([10^-2 10^2])
% ylim([10^-2 10^2])
% 
% %RMSE=sqrt(nanmean((Y-X).^2));
% %title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)
% 
% %%
% % figure; hold on; box on
% % set(gca,'XScale','Log','YScale','Log')
% % plot(CN,VAL(:,1),'.k')
% % plot([10^-2 10^2],[1 1],'k')
% % plot([1 1],[10^-2 10^2],'k')
% % plot([10^-2 10^2],[10^-2 10^2],':k')
% % xlim([10^-2 10^2])
% % ylim([10^-2 10^2])
% % xlabel('Observed C/N ratio')
% % ylabel('Predicted C/N ratio')
% %% plot Q with 1:1 line
% X=QVAL(:,1);
% Y=gQ;
% dlab='Q';
% 
% figure; hold on; box on
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% %mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% %mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% mec=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% %mfc=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% mfc=[1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((QDay==gageDMYQ(dmyq,1))&(QMonth==gageDMYQ(dmyq,2))&(QYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthWest')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% set(gca,'XScale','log','YScale','log')
% 
% whichstats = {'beta','rsquare'};
% stats = regstats(Y,X,'linear',whichstats);
% a=(stats.beta(1));
% b=stats.beta(2);
% x=linspace(min(X),max(X),100);
% plot(x,a+b.*x,'-k');
% %title(sprintf('Y = %.2g + {%.2f} * X\nR^2 = %.2f',[a b stats.rsquare]),'FontSize',14)
% 
% plot([10^-2 10^3],[10^-2 10^3],'--k')
% %plot([0 400],[0 400],'--k')
% %ylim([0 400])
% 
% RMSE=sqrt(nanmean((Y-X).^2));
% title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)
% 
% %% plot B with 1:1 line
% X=QVAL(:,3);
% Y=gB;
% dlab='B';
% 
% figure; hold on; box on
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% %mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% %mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% mec=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% %mfc=[0,0,0; 0,0,0; 0.5,0.5,0.5; 0.5,0.5,0.5; 0.75,0.75,0.75; 0.75,0.75,0.75; 0.75,0.75,0.75;];
% mfc=[1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1; 1,1,1;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((QDay==gageDMYQ(dmyq,1))&(QMonth==gageDMYQ(dmyq,2))&(QYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthWest')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% %set(gca,'XScale','log','YScale','log')
% 
% whichstats = {'beta','rsquare'};
% stats = regstats(Y,X,'linear',whichstats);
% a=(stats.beta(1));
% b=stats.beta(2);
% x=linspace(min(X),max(X(~isnan(gB))),100);
% plot(x,a+b.*x,'-k');
% % title(sprintf('Y = %.2g + {%.2f} * X\nR^2 = %.2f',[a b stats.rsquare]),'FontSize',14)
% 
% plot([0 10],[0 10],'--k')
% xlim([0 10])
% 
% RMSE=sqrt(nanmean((Y-X).^2));
% title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)
%%
% %% plot U with 1:1 line
% X=QVAL(:,2);
% Y=gU;
% dlab='U';
% 
% figure; hold on; box on
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((QDay==gageDMYQ(dmyq,1))&(QMonth==gageDMYQ(dmyq,2))&(QYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthWest')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% %set(gca,'XScale','log','YScale','log')
% 
% whichstats = {'beta','rsquare'};
% stats = regstats(Y,X,'linear',whichstats);
% a=(stats.beta(1));
% b=stats.beta(2);
% x=linspace(min(X),max(X(~isnan(gU))),100);
% plot(x,a+b.*x,'-k');
% % title(sprintf('Y = %.2g + {%.2f} * X\nR^2 = %.2f',[a b stats.rsquare]),'FontSize',14)
% 
% plot([0 1],[0 1],'--k')
% xlim([0 1])
% 
% RMSE=sqrt(nanmean((Y-X).^2));
% title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)
% 
% %% plot H with 1:1 line
% X=QVAL(:,4);
% Y=gH;
% dlab='H';
% 
% figure; hold on; box on
% %plot(X,Y,'.k')
% 
% m=['^','o','^','o','^','o','s'];
% mec=[0,0,0; 0,0,0; 0,0,0; 0,0,0; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% mfc=[0,0,0; 0,0,0; 1,1,1; 1,1,1; 0.7,0.7,0.7; 0.7,0.7,0.7; 0.7,0.7,0.7;];
% 
% for dmyq=1:7
%     clear idx2
%     idx2=((QDay==gageDMYQ(dmyq,1))&(QMonth==gageDMYQ(dmyq,2))&(QYear==gageDMYQ(dmyq,3)));
%     plot(X(idx2),Y(idx2),'Marker',m(dmyq),'MarkerEdgeColor',mec(dmyq,:),'MarkerFaceColor',mfc(dmyq,:),'LineStyle','none')
% end
% legend('6/11/2013','6/12/2013','6/23/2014','6/24/2014','6/15/2015','6/16/2015','6/18/2015','Location','NorthEast')
% xlabel(cat(2,'predicted ',dlab))
% ylabel(cat(2,'observed ',dlab))
% %set(gca,'XScale','log','YScale','log')
% 
% whichstats = {'beta','rsquare'};
% stats = regstats(Y,X,'linear',whichstats);
% a=(stats.beta(1));
% b=stats.beta(2);
% x=linspace(min(X),max(X(~isnan(gB))),100);
% plot(x,a+b.*x,'-k');
% % title(sprintf('Y = %.2g + {%.2f} * X\nR^2 = %.2f',[a b stats.rsquare]),'FontSize',14)
% 
% plot([0 1],[0 1],'--k')
% xlim([0 1])
% 
% RMSE=sqrt(nanmean((Y-X).^2));
% title(sprintf('Y = %.2g + {%.2f} * X; R^2 = %.2f \n RMSE = %.2f',[a b stats.rsquare RMSE]),'FontSize',14)

%%
%end