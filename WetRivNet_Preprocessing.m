%% WetRivNet_Preprocessing
% This file contains useful snippets of code for preprocessing the 
% river-wetland network shapefiles for use in the model code.

% Jon Czuba
% January 16, 2018

%% load shapefiles
%network = shaperead('Map\Headwaters_Mo\Headwaters_prj.shp');
%network = shaperead('Map\MNR_network10_updateattributetable4.shp');
% %NHD MRB
% network = shaperead('Map\MRB_NHD\MRB_NHD_network2_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_diss_edit_prj.shp');
%NHD BE
% network = shaperead('Map\MRB_NHD\BE_NHD_network_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_diss_prj.shp');
%NHD HW
% network = shaperead('Map\MRB_NHD\HW_NHD_network_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_HW2_edit_diss2_prj.shp');
%NHD C
% network = shaperead('Map\MRB_NHD\C_NHD_network2_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_C3_diss_prj.shp');
%NHD LBE
% network = shaperead('Map\MRB_NHD\LBE_NHD_network2_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_LBE3_diss_prj.shp');
%NHD LS
% network = shaperead('Map\MRB_NHD\LS_NHD_network2_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_LS3_diss_prj.shp');
% points = shaperead('Map\MRB_NHD\MN_River_project_all_study_sites_2014_LS_snap.shp');
%NHD W
% network = shaperead('Map\MRB_NHD\W_NHD_network_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_W2_diss_prj.shp');
%NHD Cott
% network = shaperead('Map\MRB_NHD\Cott_NHD_network_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_Cott2_diss_prj.shp');

%node = shaperead('Map\NodeXYElev.shp');
%% Martin Lakes
network = shaperead('Map\MRB_NHD\MartinLakes\BE_NHD_network_MartinLake_ds_prj_att.shp');
%%
boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_diss_prj.shp');
%ravines = shaperead('Map\MRB_NHD\MartinLakes\GBER_ravines_merge_points_snapedit_BEnetMartin_noli.shp');
%bluffs = shaperead('Map\MRB_NHD\MartinLakes\GBER_bluffs_906_trim_points2_editsnap2_basinE_BEnetMartin_noli.shp');
lakes = shaperead('Map\MRB_NHD\MartinLakes\nhd_waterbody_netintsct_gt04km_poly2_prj.shp');
%% Martin Lakes LS
network = shaperead('Map\MRB_NHD\MartinLakes\LS_NHD_network_MartinLake_ds_prj_att.shp');
boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_LS3_diss_prj.shp');
lakes = shaperead('Map\MRB_NHD\MartinLakes\nhd_waterbody_netintsct_gt04km_poly2_LS.shp');
%% Kloiber wetland layer
wetlands = shaperead('Map\MRB_NHD\KloiberWetlands\SMN_NWI_UTM_LS.shp');
%%
wets = shaperead('Shapefiles\SMN_NWI_UTM_attbrk_isoandfp_intsct.shp');
iwea = shaperead('Shapefiles\wat_iswet_intsct_diss.shp');

%% LS Wetland/Network READ
network = shaperead('RivNet\N\Shapefiles\LS_NHD_network2_add_diss_feat2line2_slap_elevj.shp');
%boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_LS3_diss_prj.shp');
wetland = shaperead('RivNet\N\Shapefiles\SMN_NWI_UTM_attbrk_noFPorISf_edit_diss3.shp');
wetland_brk = shaperead('RivNet\N\Shapefiles\SMN_NWI_UTM_attbrk_noFPorISf_edit.shp');
%dca = shaperead('RivNet\N\Shapefiles\ep2_es_wat2_edit1_diss.shp');
dca = shaperead('RivNet\N\Shapefiles\ep2_es_wat2_edit1_diss_mergeduplands.shp');
%%
network = shaperead('RivNet\N\Shapefiles\LS_channels1.shp');
wetland = shaperead('RivNet\N\Shapefiles\LS_wetlands1.shp');
dca = shaperead('RivNet\N\Shapefiles\LS_areas1.shp');
%%
nlcd_area = shaperead('Shapefiles\LS_nlcd_areas_intsct_diss.shp');
nlcd_iwea = shaperead('Shapefiles\LS_nlcd_iwea_intsct_diss.shp');

%% LS Wetland/Network WRITE
shapewrite(network,'RivNet\N\Shapefiles\LS_channels1.shp');
shapewrite(wetland,'RivNet\N\Shapefiles\LS_wetlands1.shp');
%shapewrite(dca,'RivNet\N\Shapefiles\ep2_es_wat2_edit1_diss_mergeuplands.shp');
%shapewrite(dca,'RivNet\N\Shapefiles\LS_areas1.shp');

%% Tushar Patrick
%network = shaperead('Tushar_Patrick\NHDFlowline_Tushar2_nozm.shp');
for i = 1:LinkNum
    [network(i).CPI] = CPI(1,i);
end
%% Export shapefile
shapewrite(network,'Tushar.shp');


%%
for i=1:LinkNum
   [network(i).seddepth] = seddepth(6001,i);
   %[network(i).numfracocap] = numfracocap(i); 
   %[network(i).exceedHs] = exceedHs(i);
   %[network(i).wfpeak1] = lnkpl(i);
end
% Write Shapefile
%shapewrite(network,'Map\MRB_NHD\MartinLakes\BE_NHD_network_MartinLake_ds_prj_att_overcap.shp');
shapewrite(network,'Map\MRB_NHD\MartinLakes\BE_NHD_network_stor6001seddepth.shp');

%%
% filename='Map\MRB_NHD\mnr_30m_w_prj_r2a.txt';
% %read grid
% [Z,R] = arcgridread(filename);

%% plot network
figure
box on
mapshow(network);
mapshow(boundary, 'FaceColor', 'none', 'EdgeColor', 'black');
mapshow(wetland, 'FaceColor', 'blue', 'EdgeColor', 'none');
%mapshow(bluffs);
%mapshow(ravines);
xlabel('Easting in meters')
ylabel('Northing in meters')
%%

DSlink(1:1360,1)=0;
for i=1339:1360
DSlink(i)=FeatID(find(dslink(i)==GridID,1),1);
end

%%
LinkNum=length(network);
GridID(1:LinkNum,1)=NaN;
%FromNode(1:LinkNum,1)=NaN;
ToNode(1:LinkNum,1)=NaN;
Length(1:LinkNum,1)=NaN;
Area_km(1:LinkNum,1)=NaN;
Slope(1:LinkNum,1)=NaN;
%Node(1:LinkNum+1,1)=NaN;
%NodeElev(1:LinkNum+1,1)=NaN;
usarea_km(1:LinkNum,1)=NaN;
%maxelev(1:LinkNum,1)=NaN;
%minelev(1:LinkNum,1)=NaN;
Lake(1:LinkNum,1)=NaN;

% GridID = [network.GridID]';
% FromNode = [network.from_node]';
% ToNode = [network.to_node]';
% Length = [network.Shape_Leng]';
% Area = [network.Shape_Area]';
% %Slope = [network.Slope]';
% Node = [node.Node]';
% NodeElev = [node.Elev]';
% 
% for i=1:LinkNum
%    Slope(i,1)= (NodeElev(find(FromNode(i,1)==Node))-...
%        NodeElev(find(ToNode(i,1)==Node)))./Length(i,1);
% end

% GridID = [network.Hydroseq]';
% FromNode = [network.UpHydroseq]';
% TN = [network.DnHydroseq]';
% Length = [network.LENGTHKM]';
% usarea_km = [network.TotDASqKM]';
% Area_km = [network.AreaSqKM]';
% Slope = [network.SLOPE]';
% %Node = [node.Node]';
% %NodeElev = [node.Elev]';

% GridID = [network.GridID]';
% %FromNode = [network.UpHydroseq]';
% ToNode = [network.dslink]';
% Length = [network.length_m]';
% Length = [network.LengthKM]'.*1000;
% usarea_km = [network.usarea_km2]';
% Slope = [network.SLOPE]';
% maxelev = [network.maxelevsm_]';%m
% minelev = [network.minelevsm_]';%m

GridID = [network.FeatID]';
ToNode = [network.FeatDSlnk]';
Length = [network.Length_m]';
%usarea_km = [network.A_km2]';
Area_km = [network.a_km2]';
%Slope = [network.slope]';
mxelev = [network.mxelev_m]';%m
mnelev = [network.mnelev_m]';%m
Lake = [network.LakeID]';
%usarea=usarea_km.*10.^6;
Area=Area_km.*10.^6;

Slope=(mxelev-mnelev)./Length;
Slope(Slope<1e-4)=1e-4;

%%
GridID=(1:643)';

gid=[dca.GridID]';
a=[dca.a_m2]';

Area(1:643,1)=0;
for i=1:640
    Area(gid(i))=a(i);
end
Area_km2=Area./10.^6;

LinkNum=643;

ToNode=cat(1,[network.ToGridID]',[wetland.ToGridID]');
OutletLinkID=483;
ToNode(OutletLinkID)=NaN;

Length=cat(1,[network.Length_m]',sqrt([wetland.WA_m2])');
Slope(1:LinkNum,1)=1e-4;
Slope(1:540,1)=[network.Slope]';

WA(1:LinkNum,1)=0;
WA(541:LinkNum,1)=[wetland.WA_m2]';

Feature(1:LinkNum,1)='X';
Feature(1:540,1)='C';
Feature(541:LinkNum,1)='W';

%% Determine outlet link ID
% for i=1:LinkNum
%     %find element of ToNode with one unique number
%     if numel(find(ToNode(i)==ToNode))==1
%         OutletLinkID=find(ToNode(i)==ToNode);
%         break
%     end
% end
%%
%OutletLinkID=find(ToNode==0);
%OutletLinkID_mrb=1802;%Cottonwood%1065;%5745;
OutletLinkID_mrbbe=1296;%MartinLakeLS
%OutletLinkID=1338;%MartinLakeBE
OutletLinkID_old=800008270;
%% Reassign GridID, ToNode, OutletLinkID for subbasins
GridID_old=FN;
ToNode_old=TN;
%GridID=(1:LinkNum)';
OutletLinkID=find(GridID_old==OutletLinkID_old);
for i=1:LinkNum
    if i==OutletLinkID
        ToNode(i)=0;
    else
        ToNode(i)=find(ToNode_old(i)==GridID_old);
    end
end

%% Establish connectivity structure of each link to the outlet
% each row in Connect is the connectivity structure from that
% row to the outlet. the columns correspond to the subsequent
% links to connect to the outlet.
% Connect(1:LinkNum,1:LinkNum)=NaN;
% Connect(1:LinkNum,1)=1:LinkNum;
% for i=1:LinkNum
%     j=1;
%     while ToNode(Connect(i,j),1)~=ToNode(OutletLinkID)
%         Connect(i,j+1)=find(ToNode(Connect(i,j),1)==FromNode);
%         j=j+1;
%     end
% end
%%
Connect(1:LinkNum,1:LinkNum)=NaN;
Connect(1:LinkNum,1)=1:LinkNum;
for i=1:LinkNum
    j=1;
    %while ToNode(Connect(i,j),1)~=ToNode(OutletLinkID)
    while ~isnan(ToNode(Connect(i,j),1))
        Connect(i,j+1)=find(ToNode(Connect(i,j),1)==GridID);
        j=j+1;
    end
end
%% Compute width function for every link
% for every row of Connect i (every link in the basin),
% find which links connect through that point and mark these
% links j in the column of Sub
% diagonal = 1

% Sub(i,j)
% *each row i indicates which indicies that node connects
% through (including itself) to reach the outlet
% *each column j corresponds to all indicies upstream of
% and including j
Sub(1:LinkNum,1:LinkNum)=0;
for j=1:LinkNum
    for i=1:LinkNum
        if ~isempty(find(Connect(i,:)==j,1))
            Sub(i,j)=1;
        end
    end
end

%% Sum upstream area
usarea(1:LinkNum,1)=NaN;
for i=1:LinkNum
    usarea(i,1)=sum(Sub(:,i).*Area);
end
usarea_km=usarea./10.^6;
%give Sub NaN for now for plotting below
%Sub(Sub==0)=NaN;
%Sub(isnan(Sub))=0;

%% compute isolated wetland area
%iWA(1:LinkNum,1)=0;
iwea(1:LinkNum,1)=0;
for i=1:LinkNum
   %iWA(i,1)=sum(wa_m2_2(GridID_1==i)); 
   iwea(i,1)=sum(iwea_m2(GridID1==i)); 
end
%%
landArea=Area-iWA-WA;
plen=(1-landArea./Area).*100;
plen(isnan(plen))=0;
plandArea=100-plen;
%%
efflandArea=Area-iwea-WA;
efflandArea(efflandArea<0)=0;
(efflandArea./Area).*100;

%%
%uslandarea(1:LinkNum,1)=NaN;
usefflandarea(1:LinkNum,1)=NaN;
for i=1:LinkNum
    %uslandarea(i,1)=sum(Sub(:,i).*landArea);
    usefflandarea(i,1)=sum(Sub(:,i).*efflandArea);
end
usplen=(1-uslandarea./usarea).*100;

%% Compute directly contributing area contributing N
% determine area of any in-channel wetlands, isolated wetlands,
% the area contributed to those isolated wetlands, and, soon to
% be accounted for, the areas that are not row-crop agriculture

fainN(1:LinkNum,1)=NaN;

for i = 1:LinkNum
    
    if sum(A_GridID==i)==0
        continue
    end
    
    clear tota aga iaga
    tota=sum(A_a_m2(A_GridID==i));
    
    if sum(((A_GRIDCODE==82)+(A_GridID==i))==2)~=0  
        aga=A_a_m2(find(((A_GRIDCODE==82)+(A_GridID==i))==2));
    else
        aga=0;
    end
    
    if sum(((I_GRIDCODE==82)+(I_GridID==i))==2)~=0
        iaga=I_a_m2(find(((I_GRIDCODE==82)+(I_GridID==i))==2));
    else
        iaga=0;
    end
    
    fainN(i,1)=(aga-iaga)./tota;
    
end

fainN(isnan(fainN))=0;
fainN(fainN<0)=0;

%% Compute directly contributing area contributing C

fainC(1:LinkNum,1)=NaN;
A_a_m2_nag=A_a_m2;
I_a_m2_nag=I_a_m2;

A_a_m2_nag(A_GRIDCODE==82)=0;
A_a_m2_nag(A_GRIDCODE==31)=0;
A_a_m2_nag(A_GRIDCODE==21)=0;
A_a_m2_nag(A_GRIDCODE==22)=0;
A_a_m2_nag(A_GRIDCODE==23)=0;
A_a_m2_nag(A_GRIDCODE==24)=0;

I_a_m2_nag(I_GRIDCODE==82)=0;
I_a_m2_nag(I_GRIDCODE==31)=0;
I_a_m2_nag(I_GRIDCODE==21)=0;
I_a_m2_nag(I_GRIDCODE==22)=0;
I_a_m2_nag(I_GRIDCODE==23)=0;
I_a_m2_nag(I_GRIDCODE==24)=0;

for i = 1:LinkNum
    
    if sum(A_GridID==i)==0
        continue
    end
    
    clear tota naga niaga
    tota=sum(A_a_m2(A_GridID==i));
    
    naga=sum(A_a_m2_nag(A_GridID==i));
    inaga=sum(I_a_m2_nag(I_GridID==i));
      
    fainC(i,1)=(naga-inaga-WA(i))./tota;
    
end

fainC(isnan(fainC))=0;
fainC(fainC<0)=0;

%% Construct Adjaceny matrix
Adj(1:LinkNum,1:LinkNum)=0;
%establish DS connection
for i=1:LinkNum
    if isnan(Connect(i,2))
        continue
    end
    Adj(Connect(i,1),Connect(i,2))=1;
end
%establish US connection
%Adj=Adj+Adj';
%Adj(Adj==0)=NaN;
%Adj(isnan(Adj))=0;

%% Determine incremental area
Area_km(1:LinkNum,1)=NaN;
for i=1:LinkNum
Area_km(i,1)=usarea_km(i)-sum(usarea_km(Adj(:,i)>0));
end

%% Stream Order for 1st order links
%Ord1=find(nansum(Sub,1)==1)';

Ord1(1:LinkNum,:)=0;
for i=1:LinkNum
   if isempty(find(i==ToNode,1)) 
       Ord1(i)=1;
   end
end

%% Stream Order for all links
Ord(1:LinkNum,1)=NaN;
hw=(sum(Sub,1)==1);%identifies headwater, first-order links
Ord(hw)=1;

tocomp=sum(isnan(Ord));
while tocomp>0
    for i=1:LinkNum

        if isnan(Ord(i)) &&  ~isnan(sum(Ord(Connect(:,2)==i)))
            if length(Ord(Connect(:,2)==i))==1
                Ord(i)=Ord(Connect(:,2)==i);
            else 
                uo=sort(Ord(Connect(:,2)==i),'descend');
                if uo(1)>uo(2)
                    Ord(i)=uo(1);
                else if uo(1)==uo(2)
                        Ord(i)=uo(1)+1;
                    end
                end
            end
        end
       
    end
    tocomp=sum(isnan(Ord));
end

%%
Lke_ID=[lakes.LakeID]';
Lke_Larea_km=[lakes.Shape_Area]'./10^6;
%%
Lke_ToNode(1:length(Lke_ID),1)=NaN;
Lke_Warea_km(1:length(Lke_ID),1)=NaN;
for i=1:length(Lke_ID)
    %find GridID of links with same LakeID
    ntlk=find(Lake==Lke_ID(i));
    %find links those links connect to
    nlcon=ToNode(ntlk);
    for j=1:length(ntlk)
        %find link which connects to only 1 ds link
        if numel(find(nlcon(j)==nlcon))==1
            Lke_ToNode(i,1)=nlcon(j);
            Lke_Warea_km(i,1)=usarea_km(ntlk(j));
            break
        end
    end
end
%%
%% Break out wetland codes from Kloiber
WetNum=length(wetlands);
ridx(1:WetNum,1)=1;

LLWW1=[];
LLWW2=[];
LLWW3=[];

Att1=[];
Att2=[];
Att3=[];
Att4=[];
Att5=[];
Att6=[];

for i=1:WetNum
    LLWW1=cat(1,LLWW1,[wetlands(i,1).HGM_CODE(1:2)]);
    LLWW2=cat(1,LLWW2,[wetlands(i,1).HGM_CODE(3:4)]);
    LLWW3=cat(1,LLWW3,[wetlands(i,1).HGM_CODE(5:6)]);
    
    Att1=cat(1,Att1,[wetlands(i,1).ATTRIBUTE(ridx(i))]);
    ridx(i)=ridx(i)+1;
    
    if strcmp(Att1(i),'P')
        Att2=cat(1,Att2,'0');
    else
        Att2=cat(1,Att2,[wetlands(i,1).ATTRIBUTE(ridx(i))]);
        ridx(i)=ridx(i)+1;
    end
    
    %for those with multiple classifications
    if length([wetlands(i,1).ATTRIBUTE])>=5
        if strcmp([wetlands(i,1).ATTRIBUTE(5)],'/')
            Att3=cat(1,Att3,'EM');
            Att4=cat(1,Att4,'1');
            Att5=cat(1,Att5,[wetlands(i,1).ATTRIBUTE(9)]);
            
            if length([wetlands(i,1).ATTRIBUTE])<10
                Att6=cat(1,Att6,'0');
            else
                Att6=cat(1,Att6,[wetlands(i,1).ATTRIBUTE(10)]);
            end
            
            continue
        end
    end
    
    
    Att3=cat(1,Att3,[wetlands(i,1).ATTRIBUTE(ridx(i):ridx(i)+1)]);
    ridx(i)=ridx(i)+2;
    
    if or(strcmp(Att3(i,:),'SS'),...
            or(strcmp(Att3(i,:),'EM'),...
            strcmp(Att3(i,:),'FO')))
        Att4=cat(1,Att4,[wetlands(i,1).ATTRIBUTE(ridx(i))]);
        ridx(i)=ridx(i)+1;
    else
        Att4=cat(1,Att4,'0');
    end
    
    Att5=cat(1,Att5,[wetlands(i,1).ATTRIBUTE(ridx(i))]);
    ridx(i)=ridx(i)+1;
    
    if length([wetlands(i,1).ATTRIBUTE])<(ridx(i))
        Att6=cat(1,Att6,'0');
    else
        Att6=cat(1,Att6,[wetlands(i,1).ATTRIBUTE(ridx(i))]);
        ridx(i)=ridx(i)+1;
    end
        
end

%% Assign values to shapefile
for i = 1:WetNum
    [wetlands(i).LLWW1] = LLWW1(i,:);
    [wetlands(i).LLWW2] = LLWW2(i,:);
    [wetlands(i).LLWW3] = LLWW3(i,:);
    
    [wetlands(i).Att1] = Att1(i,:);
    [wetlands(i).Att2] = Att2(i,:);
    [wetlands(i).Att3] = Att3(i,:);
    [wetlands(i).Att4] = Att4(i,:);
    [wetlands(i).Att5] = Att5(i,:);
    [wetlands(i).Att6] = Att6(i,:);
end
%% Assign values to network
for i = 1:size(network,1)
    [network(i).Slope] = Slope(i,:);
end
%%
for i = 1:540
    [network(i).usarea_m2] = usarea(i,:);
    [network(i).usarea_km2] = usarea_km2(i,:);
end
%% Export shapefile
shapewrite(wetlands,'Map\MRB_NHD\KloiberWetlands\SMN_NWI_UTM_LS_attbrk.shp');

%% Determine merged wetland properties
% linked via ID field
EMarea(1:size(wetland,1),1)=0;
Parea(1:size(wetland,1),1)=0;
L1area(1:size(wetland,1),1)=0;
L2area(1:size(wetland,1),1)=0;
Warea(1:size(wetland,1),1)=0;
for i=1:size(wetland,1)
    clear idx
    %find entries composing wetland feature
    idx=find([wetland_brk.ID]==wetland(i).ID);
    for j=1:length(idx)
       
        Warea(i,1)=Warea(i,1)+wetland_brk(idx(j)).WA_m2;
        
        %if emergent
        if strcmp(wetland_brk(idx(j)).Att3,'EM')
            %add the area to emergent area
            EMarea(i,1)=EMarea(i,1)+wetland_brk(idx(j)).WA_m2;
        end
        
        switch wetland_brk(idx(j)).Att1
            case 'P'
                Parea(i,1)=Parea(i,1)+wetland_brk(idx(j)).WA_m2;
            case 'L'
                switch wetland_brk(idx(j)).Att2
                    case '1'
                        L1area(i,1)=L1area(i,1)+wetland_brk(idx(j)).WA_m2;
                    case '2'
                        L2area(i,1)=L2area(i,1)+wetland_brk(idx(j)).WA_m2;
                end
            case 'R'
                L2area(i,1)=L2area(i,1)+wetland_brk(idx(j)).WA_m2;
        end
 
    end
end

%WA=[wetland.WA_m2]';%entire wetland area
pEM=EMarea./Warea.*100;%percent emergent
pP=Parea./Warea.*100;
pL1=L1area./Warea.*100;
pL2=L2area./Warea.*100;
%% Assign values to wetland
for i = 1:size(wetland,1)
    [wetland(i).pEM] = pEM(i,:);
    [wetland(i).pP] = pP(i,:);
    [wetland(i).pL1] = pL1(i,:);
    [wetland(i).pL2] = pL2(i,:);
end
%%
for i = 541:643
    [wetland(i-540).usarea_m2] = usarea(i,:);
    [wetland(i-540).usarea_km2] = usarea_km2(i,:);
end
%%
GridID(GridID>=2000)=GridID(GridID>=2000)-2000;
 %% Assign values to dca
for i = 1:size(dca,1)
    [dca(i).GridID] = GridID(i,:);
end   


%% Divide links greater than 5km into multiple links
network = shaperead('RivNet\N\Shapefiles\LS_NHD_network2_add_diss_feat2line2.shp');
LinkNum=length(network);
len(1:LinkNum,1)=NaN;
%id those links greater than 5 km
for i=1:LinkNum
    %len(i,1)=network(i).L1_m;
    len(i,1)=network(i).LengthKM.*1000;
end

% compute incremental length through each link, see code below

%%
X=[];
Y=[];
slen=100;%m
%split=find(len>5000);%break for <5km
split=find(len>slen);%break for <100m
for i=1:length(split)
    i
    clear dist numlnks incdist
    dist=(network(split(i)).L.*len(split(i)))';%distance along link
    numlnks=ceil(len(split(i))./slen);%number of links to break into
    incdist=len(split(i))./numlnks;%break links into this length
    
    j=1;
    while j <= numlnks-1
        clear C I
        [C,I]=min(abs(dist-(j.*incdist)));
        X=cat(1,X,network(split(i)).X(I));
        Y=cat(1,Y,network(split(i)).Y(I));
        j=j+1;
    end
end

%% site values

Site_ID=str2num(char(Site_ID));
data(data==9999)=NaN;

[sSite_ID ix] = sort(Site_ID);
sGridID=GridID(ix);

%%
DOC(1:length(Site_ID),1:3)=NaN;
TDN(1:length(Site_ID),1:3)=NaN;
NO3(1:length(Site_ID),1:3)=NaN;
Year=[2013, 2014, 2015];

for i=1:length(data)
    row=find(data(i,2)==sSite_ID);
    col=find(data(i,1)==Year);
    
    DOC(row,col)=data(i,3);
    TDN(row,col)=data(i,4);
    NO3(row,col)=data(i,5);
end

%%
for i=34:78 
sGridID(i,1)=sGrid_ID(find(sSiteID(i)==sSite_ID));
end

%% QBHU rating at Le Sueur Gage

figure; hold on; box on;
set(gca,'XScale','log','YScale','log')
plot(Q,B,'.b')

%Fr=U./sqrt(9.81.*H);

Qbf=5;

whichstats = {'beta','rsquare'};
X1=Q(Q<Qbf);
Y1=B(Q<Qbf);
stats1 = regstats(log10(Y1),log10(X1),'linear',whichstats);
a1=10.^(stats1.beta(1));
b1=stats1.beta(2);
x1=linspace(min(X1),max(X1),100);
plot(x1,a1.*x1.^b1,'-k');

X2=Q(Q>=Qbf);
Y2=B(Q>=Qbf);
stats2 = regstats(log10(Y2),log10(X2),'linear',whichstats);
a2=10.^(stats2.beta(1));
b2=stats2.beta(2);
x2=linspace(min(X2),max(X2),100);
plot(x2,a2.*x2.^b2,'-k');

%% compute n
Slope=0.001113;%slope at gage = Slope(472)
n=(H.^(2/3))./U.*(Slope.^(1/2));

figure; hold on; box on
plot(Q,n,'.b')
set(gca,'XScale','log')
mean(n(Q>40))
%%

%%
%% Plot Spatial Data
%set timestep to plot
t = 30;%64;%83%%115,83,40,20 

%set n as # of parcels in each state
%n = hist(State(t,:),0:1:LinkNum);
%nttl = 'Number of parcels';
%nttl = 'Cluster size, km';
%nttl = 'Parcels removed';
nttl = 'N concentration';

%set n as # of inactive states
% n = hist(State(t,logical(Inactive(t,:))),0:1:LinkNum);
% nttl = 'Number of inactive parcels';

%LCS=LCS_spdist;
%LCS=LCS_parconc;

% out2zero=1850;
% LCS=nansum(LCS_parconc.*LCS_spdist,1)./1000./(out2zero./365./0.175.*4);%#*km/yr
%LCS=nansum(LCS_parconc.*(LCS_spdist./1000).*(dt./60./60./24./365./0.175.*4),1);%#*km*yr

%LCS=nansum(LCS_parconc.*(dt./60./60./24./365./0.175.*4),1);%./USlinks';%#*yr

%LCS=nansum(numpar);
%  LCS=nansum(LCS_parconc.*LCS_spdist,1)./1000./...
%      (sum(~isnan(LCS_parconc),1)./365./0.175.*4);%./B';%#*km/yr
% LCSB=nansum(LCS_parconc.*LCS_spdist,1)./1000./...
%      (sum(~isnan(LCS_parconc),1)./365./0.175.*4)./B';%#*km/yr

% LCS(LCS==0)=NaN;                    
% mn=nanmin(nanmin(LCS));
% LCS(isnan(LCS))=0;
%LCSB(isnan(LCSB))=0;
%LCS(LCS~=repmat(nanmax(LCS_spdist,[],2),1,LinkNum))=0;

% csource=LCS;
% csource(:,:)=0;
% loc=find(LCS(t,:)>0);
% for loci=1:length(loc)
% csource(t,Source{t,loc(loci)})=1;
% end
%LCSa(1,1:LinkNum)=0;
%LCSa=max(LCS,[],1);
% for i=1:LinkNum
% LCSa(1,LCS(i,:)>0)=LCS(i,LCS(i,:)>0);
% end
%LCS=loc594';
% mn=nanmin(nanmin(numpar(numpar>0)));

% n_zero = n(1);
% n = n(2:end);

for i = 1:LinkNum
    %[network(i).sed] = n(i);
    %[network(i).sed] = input(i,1);
     %[network(i).sed] = numpar(t,i);
     %[network(i).sed] = CapRate(i,1);
     %[network(i).sed] = exceedHs(i,1);
     %[network(i).sed] = lnkconc(t,i);
     [network(i).sed] = abs(herr(i,1));
     %[network(i).npar7] = numpar(7,i);
     %[network(i).sed64] = numpar(64,i);
     %[network(i).sed631] = numpar(631,i);
%     [network(i).sed115] = numpar(115,i);
    %[network(i).sed] = CsizeComp(t,i);
    %[network(i).sed] = csou(1,i);
  %  [network(i).sed] = LCS(1,i);
%     [network(i).lcs1] = LCS(1,i);
%     [network(i).lcs7] = LCS(7,i);
%     [network(i).lcs64] = LCS(64,i);
%     [network(i).lcs631] = LCS(631,i);
    %[network(i).lcsb] = LCSB(t,i);
    %[network(i).cs64] = csource(t,i);
    %[network(i).cs631] = csource(t,i);
    %[network(i).sed] = LCSa(1,i);
    %[network(i).sed] = TimeS(i)/60/60/24/365;
    %[network(i).pathway] = pathway(i);
    %[network(i).sed] = Dgm(i);
    %[network(i).sed] = max(numpar(:,i),[],1);%max in link
    %[network(i).sed] = tleng(i);
    %[network(i).sed] = VelocityS(i).*1000000;
    %[network(i).sed] = Sub(i,j);
    %[network(i).tau] = tau(i);
    %[network(i).sed] = Length(i);
    %[network(i).sed] = input(i,1);
    %[network(i).sub] = subbasin(i,1);
    %[network(i).Dist] = Dist(i,1)./1000;
    %[network(i).TimeS] = TimeS(i,1)./60./60./24./365./0.175;
    %[network(i).slopeML] = Slope(i,1);
end
%
%
f1 = figure;
set(f1,'Position',[213 50 938 632]);
%a1 = axes('Position',[0.08 0.1 0.4 0.8]);
a1 = axes;
%cbpos = [0.78 0.22 0.01 0.6];%mrb,hw,c
cbpos = [0.90 0.22 0.01 0.6];%be
%cbpos = [0.93 0.22 0.01 0.6];%w
box(a1);

%edge = [0 1 2 4 8 16 32 64 128 200];
%edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
%edge = [0 1 2 4 6 8 10 12 16 20];
edge = [0 1 2 3 4 5 10 20 50 100];
%edge = [0 1 10 50 100 500 1000 5000 10000 100000];
%edge = [0 1 10 50 70 100 200 300 500 1000];
%edge = [0 .5 .75 1 1.25 1.5 1.75 2 2.25 2.5];
%edge = cat(2,0,logspace(-2,0,9));
%edge = cat(2,0,round(exp(linspace(log(1),log(max([network.sed])),9))));
%edge = cat(2,0,(exp(linspace(log(min([network.sed])),log(max([network.sed])),9))));
%edge = cat(2,0,(exp(linspace(log(1),log(max(max(LCS))),9))));
%edge = cat(2,0,round(exp(linspace(log(1),log(max(max(numpar))),9))));
%edge = cat(2,0,round(linspace(mn,ceil(max(max(LCS))),9)));
%edge = cat(2,0,round(exp(linspace(log(mn),log(ceil(max(max(LCS)))),9))));
%edge = cat(2,0,round((linspace((min([network.sed])),(max([network.sed])),9))));
%edge = cat(2,0,round((linspace((1),(max(max(numpar))),9))));

% ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];nttl};
ttl={['Time = ',num2str(fix(time(t))),' days ',...
    num2str(round(mod(time(t)*24,24))),' hours'];['Time step = ',num2str(t)];nttl};
%ttl=''; 
[a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);

%% Create Movie

%set timestep to plot
t = 1;
nttl = 'Number of parcels';
%nttl = 'Comp. Cluster Size';
%nttl = 'Cluster Size, m';

% LCS=LCS_spdist;
% mn=nanmin(nanmin(LCS));
% LCS(isnan(LCS_spdist))=0;
% LCS(LCS~=repmat(nanmax(LCS_spdist,[],2),1,LinkNum))=0;

for i = 1:LinkNum
    [network(i).sed] = numpar(t,i);
    %[network(i).sed] = LCS(t,i);
end
%
%
f1 = figure;
set(f1,'Position',[213 50 938 632]);
%a1 = axes('Position',[0.08 0.1 0.4 0.8]);
a1 = axes;
%cbpos = [0.78 0.22 0.01 0.6];%mrb,hw,c
cbpos = [0.90 0.22 0.01 0.6];%be
%cbpos = [0.93 0.22 0.01 0.6];%w
box(a1);

%edge = [0 1 2 4 8 16 32 64 128 200];
%edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
%edge = cat(2,0,round(exp(linspace(log(1),log(max([network.sed])),9))));
%edge = cat(2,0,round(exp(linspace(log(1),log(max(max(CsizeComp))),9))));
%edge = cat(2,0,round(exp(linspace(log(1),log(max(max(LCS))),9))));
%edge = cat(2,0,round(linspace(mn,max(max(LCS)),9)));
edge = cat(2,0,round(exp(linspace(log(1),log(max(max(numpar))),9))));
edge(4)=3;
% ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];nttl};
% ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];['Time step = ',num2str(t)];nttl};
%sand 0.4 mm
ttl={['Time = ',num2str(fix(t/365/0.175*4)),' years ',...
     num2str(round(mod(t/0.175*4,365))),' days'];['Time step = ',num2str(t)];nttl};
%uniform 1 m/s
% ttl={['Time = ',num2str(fix(t.*144./60./60./24)),' days ',...
%      num2str(round(mod(t.*144./60./60,24))),' hours ',...
%      num2str(round(mod(t.*144./60,60))),' minutes'];['Time step = ',num2str(t)];nttl}; 
 
[a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);
%
writerObj = VideoWriter('ex9','MPEG-4');
writerObj.FrameRate = 15;
writerObj.Quality =100;
open(writerObj);

set(gcf,'Renderer','Painters');

frame = getframe(gcf);
writeVideo(writerObj,frame);
cla(a1);
colorbar('delete');

for t=2:1:1782%LS-932%W-1034
    for i = 1:LinkNum
        [network(i).sed] = numpar(t,i);
        %[network(i).sed] = CsizeComp(t,i);
        %[network(i).sed] = LCS(t,i);
    end
    
%     ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];nttl};
%     ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%         num2str(round(mod(t/0.175,365))),' days'];['Time step = ',num2str(t)];nttl};
%sand 0.4 mm
    ttl={['Time = ',num2str(fix(t/365/0.175*4)),' years ',...
        num2str(round(mod(t/0.175*4,365))),' days'];['Time step = ',num2str(t)];nttl};
%uniform 1 m/s
%     ttl={['Time = ',num2str(fix(t.*144./60./60./24)),' days ',...
%         num2str(round(mod(t.*144./60./60,24))),' hours ',...
%         num2str(round(mod(t.*144./60,60))),' minutes'];['Time step = ',num2str(t)];nttl};
    
    [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);

    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    cla(a1);
    colorbar('delete');
    
end

close(writerObj);

%% Determine incremental length spatially though links 2-D
%length through link
for i=1:LinkNum
    network(i).L=network(i).X;
end
for i=1:LinkNum
    network(i).L(1)=0;
    for j=2:length(network(i).X)-1
        network(i).L(j)=sqrt((network(i).X(j)-network(i).X(j-1)).^2+...
            (network(i).Y(j)-network(i).Y(j-1)).^2);
    end
end
for i=1:LinkNum
    for j=3:length(network(i).X)-1
        network(i).L(j)=network(i).L(j)+network(i).L(j-1);
    end
    network(i).L=network(i).L./nanmax(network(i).L);
end

%% Determine distance between 2 points
su=50;%42  %54;%46%site id of us point
lu=1119;%686;%link of us point
sd=51;%43  %53;%45%site id of ds point
ld=1113;%1054;%link of ds point

iu=find(lu==GridID_mrb);
id=find(ld==GridID_mrb);

d=sum(Length(Connect(iu,2:find(Connect(iu,:)==id)-1)))+...
(1-network(iu).L(find(points(42).X==network(iu).X))).*Length(iu)+...%us point
network(id).L(93).*Length(id);%ds point
%network(id).L(find(points(43).X==network(id).X)).*Length(id);%ds point

%% Associate subbasins with individual links
%% load shapefiles
network = shaperead('Map\MRB_NHD\MartinLakes\BE_NHD_network_MartinLake_ds_prj_att.shp');
%basins = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_prj.shp');
basins = shaperead('Map\MRB_NHD\MartinLakes\gsmsoilmu_a_mn_ia_prj_BEclip_mat_prj2.shp');


%% 

numvert(1:1360,1:131)=0;

for j=1:131
    %for each basin
    Xpolygon=basins(j).X;
    Ypolygon=basins(j).Y;
    
    for i=1:1360
        %check every link to see if any verticies are within that basin
        Xpoint=network(i).X;
        Ypoint=network(i).Y;
        
        IN=inpolygon(Xpoint(~isnan(Xpoint)),Ypoint(~isnan(Ypoint)),...
            Xpolygon(~isnan(Xpolygon)),Ypolygon(~isnan(Ypolygon)));
        %determine the number of verticies from that link in each basin
        numvert(i,j)=sum(IN);
        
        clear IN Xpoint Ypoint
    end
    clear Xpolygon Ypolygon
end
%%
%linkforbasin(1:1674,1)=NaN;
numlinks(1:1674,1)=NaN;
for j=1:1674
    numlinks(j,1)=sum(numvert(:,j)>0,1);
    %[c linkforbasin(j,1)]=max(numvert(:,j)>0);
end
%%
[C linkforbasin]=max(numvert,[],1);
linkforbasin=linkforbasin';
linkforbasin(numlinks==0)=0;
%%
[D basinforlink]=max(numvert,[],2);
%for each link what is the dominant underlying basin
%%
%linkmat=char(1360,1);

for i=1:1360
   linkmat(i,1)=basins(basinforlink(i)).type(1);
end

%%
Up_sand_frac(1:1360,1)=NaN;

for i=1:1360
    
    switch linkmat(i)
        case 't'
            Up_sand_frac(i,1)=0.35;
        case 'l'
            Up_sand_frac(i,1)=0.10;
        case 'a'
            Up_sand_frac(i,1)=0.50;
        case 'o'
            Up_sand_frac(i,1)=0.50;
    end
    
end


