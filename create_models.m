%% Create models for single file
clc
clear all
close all

%filedirectory = fullfile("..","Tapes Data","SuperOx","SuperOx GdBCO 2G HTS","All data.xlsx");
%filedirectory = fullfile("..","Tapes Data","SuperPower","SuperPower Advanced Pinning 2G HTS","All data.xlsx");
filedirectory = fullfile("..","Tapes Data","AMSC","AMSC Amperium® 2G HTS","All data.xlsx");

%[Table,TableRelevant] = read_alldata(filedirectory,-1,-1,-1,-1,-1,-1,0);
[Table,TableRelevant] = read_alldata_opt(filedirectory,-1,-1,-1,-1,-1,-1,0);
%[ParameterTable] = compute_kim_param(TableRelevant);
%[ParameterTable] = compute_kim_param_weighted(TableRelevant);
%[ParameterTable] = compute_kim_param_interp(TableRelevant);
%% Construct Samples Table
clc
clear all
close all

%Samples = create_sample_table();
%Samples = add_data_sample_table(Samples,-1,-1,0,3,-1,-1);
%save AllSamplesBase3T.mat
%load AllSamplesBase3T.mat
%Samples = add_kim_sample_table(Samples);
%save AllSamplesKim3T.mat
%load AllSamplesKim3T.mat
%Samples = add_custom_sample_table(Samples);
%save AllSamplesCustom3T.mat
load AllSamplesCustom3T.mat
Samples = add_advanced_custom_sample_table(Samples);
save AllSamplesAllModels3T.mat
%% Functions to create model
function [Samples] = create_sample_table()
filedirectory = [fullfile("..","Tapes Data","SuperPower","SuperPower Advanced Pinning 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","SuperOx","SuperOx GdBCO 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","SuperOx","SuperOx YBCO 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","AMSC","AMSC Amperium® 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","Sumitomo","Sumitomo new type H DI-BSCCO®","All data.xlsx");...
    fullfile("..","Tapes Data","Fujikura","Fujikura FESC 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","Fujikura","Fujikura FYSC 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","SST","Shanghai Superconductor High Field Low Temperature 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","SST","Shanghai Superconductor Low Field High Temperature 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","AMSC","AMSC Amperium® Type 8502-350 coil formulation 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","AMSC","AMSC Amperium® Type 8700 cable formulation 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","AMSC","AMSC Amperium® Type 8702 non-magnetic cable formulation 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","Samri","Samri 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","STI","STI Conductus® 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","THEVA","THEVA Pro-Line 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","suNAM","SuNAM HAN04200 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","suNAM","SuNAM SAN04200 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","SCSC","Shanghai Creative Superconductor Technologies 2G HTS","All data.xlsx");...
    fullfile("..","Tapes Data","Nexans","Nexans Bi-2212","All data.xlsx");...
    fullfile("..","Tapes Data","InnoST","InnoST 1G HTS","All data.xlsx")];

company = ["SuperPower";...
    "SuperOx";...
    "SuperOx";...
    "AMSC";...
    "Sumitomo";...
    "Fujikura";...
    "Fujikura";...
    "SST";...
    "SST";...
    "AMSC";...
    "AMSC";...
    "AMSC";...
    "Samri";...
    "STI";...
    "THEVA";...
    "suNAM";...
    "suNAM";...
    "SCSC";...
    "Nexans";...
    "InnoST"];

label = ["SuperPowerAP";...
    "SuperOxGdBCO";...
    "SuperOxYBCO";...
    "AMSC1";...
    "Sumitomo";...
    "FujikuraFESC";...
    "FujikuraFYSC";...
    "SSTHFLT";...
    "SSTLFHT";...
    "AMSC2";...
    "AMSC3";...
    "AMSC4";...
    "Samri";...
    "STI";...
    "THEVA";...
    "suNAMHAN";...
    "suNAMSAN";...
    "SCSC";...
    "Nexans";...
    "InnoST"];

composition = ["YBCO";...
    "GdBCO";...
    "YBCO";...
    "???";...
    "Bi-2223";...
    "EuBCO";...
    "GdBCO";...
    "???";...
    "GdBCO";...
    "???";...
    "???";...
    "???";...
    "YBCO";...
    "SmBCO";...
    "GdBCO";...
    "???";...
    "GdBCO";...
    "???";...
    "Bi-2212";...
    "Bi-2223"];

width = [0.1;...
    0.4;...
    0.1;...
    0.44;...
    0.43;...
    0.05;...
    0.05;...
    0.05;...
    0.05;...
    1.2;...
    0.04;...
    0.04;...
    0.2;...
    0.05;...
    0.4;...
    0.4;...
    0.4;...
    0.5;...
    0.4;...
    0.4];

datestrings = {'29-Jan-2021';...
    '20-Apr-2018';...
    '20-Jan-2021';...
    '23-Jul-2014';...
    '18-Dec-2013';...
    '14-Jul-2021';...
    '11-Sep-2021';...
    '20-Jan-2022';...
    '15-Dec-2021';...
    '04-Oct-2017';...
    '05-Jul-2017';...
    '08-Aug-2017';...
    '18-Apr-2019';...
    '06-Aug-2014';...
    '05-Feb-2019';...
    '31-Jul-2015';...
    '23-May-2017';...
    '23-Sep-2016';...
    '05-Nov-2015';...
    '25-Sep-2019'};
datevar = datetime(datestrings,'InputFormat','dd-MMM-yyyy')

Samples = table(company,label,composition,width,datevar,filedirectory);
Samples.Properties.VariableNames = ["Company","Label","Composition","WidthEff","DateReceived","FileDirectory"];
end

function [Table,TableRelevant] = read_alldata(filedirectory,T1,T2,B1,B2,A1,A2,width)
%% Reads table
disp("Reading from "+filedirectory)
opts = detectImportOptions(filedirectory);
opts = setvartype(opts,'double');
opts = setvaropts(opts,'TreatAsMissing',{'.','NA','','--','NaN','-99'});
Table = readtable(filedirectory,opts);
%% Cleans table
disp('Cleaning')
Table = rmmissing(Table);
%% Creates Critical Current Sheet Density column if non existent
if string(Table.Properties.VariableNames(4)) == 'Ic_A_'
    Table.Ic_w_A_cm_ = Table.Ic_A_ / width;%creates the new column
    Table = movevars(Table, 'Ic_w_A_cm_', 'After', 'Angle___');
    disp('Fixed no Critical Current Sheet Density column')
end
%% Formats into type 1 or type 2 table (type 1: 13 columns, type 2: 10 columns), and more converting into numeric
Names = {'Temperature'; 'Field'; 'Angle'; 'CriticalCurrentSheetDensity'; 'CriticalCurrent'; 'IdealityFactor';... 
    'V0'; 'V1'; 'HallField'; 'HallAngle'; 'SetTemperature'; 'SetField'; 'SetAngle'};

n_columns = numel(Table(1,:));
disp('Converting into correct type')
if n_columns == 13
    disp('type 1 table')
    table_type=1;
    Table.Properties.VariableNames = Names;
elseif n_columns == 10
    disp('type 2 table')
    table_type=2;
    Table.Properties.VariableNames = Names(1:10);
end
%% Converts type 2 table into type 1
if table_type == 2
    Table = addvars(Table,zeros(numel(Table.(1)),1),zeros(numel(Table.(1)),1),zeros(numel(Table.(1)),1));
    Table.Properties.VariableNames(11:13) = Names(11:13);

    Table.SetTemperature = round(Table.Temperature,1);
    Table.SetAngle = round(Table.Angle);
    roundTargets = [0;0.01;0.015;0.02;0.03;0.05;0.07;0.1;0.15;0.2;0.3;0.5;0.7;1;1.5;2;3;4;5;6;7;8];
    Table.SetField = interp1(roundTargets,roundTargets,Table.Field,'nearest','extrap');
end
%% Rounds errors in SetTemperature
idx = abs(Table.SetTemperature - round(Table.SetTemperature,2,"significant")) < 0.5;
Table.SetTemperature(idx) = round(Table.SetTemperature(idx),2,"significant");
%% Removes rows with multiple 0 cells (errors)
n_rows = numel(Table.(1));
j=1;
cenas = 0;
for i = 1:n_rows
    if ~any(Table{j,4:10})
        Table(j,:) = [];
        cenas = cenas + 1;
    else
        j=j+1;
    end
end
if cenas
    disp("Removed "+num2str(cenas)+" lines with full 0")
end
%% Sorts and removes duplicate rows
Table = sortrows(Table, ["SetTemperature","SetField","SetAngle"]);
n_rows = numel(Table.CriticalCurrent);
i=1;
j=0;
while i <= n_rows - j - 1
    if(Table.SetTemperature(i) == Table.SetTemperature(i+1) && Table.SetField(i) == Table.SetField(i+1) && Table.SetAngle(i) == Table.SetAngle(i+1))
        Table(i,:) = [];
        j = j + 1;
        i = i - 1;
    end
    i = i + 1;
end
%% Filters for min and max temperatures, field and angles
[Table] = filter_table(Table,T1,T2,B1,B2,A1,A2);
%% Filters for relevant temperatures and relevant angles
TableRelevant = Table;
[all_temps_counts,all_temps] = groupcounts(TableRelevant.SetTemperature);
[all_angles_counts,all_angles] = groupcounts(TableRelevant.SetAngle);
for i=1:numel(all_temps)
    if (all_temps_counts(i) < 5)
        TableRelevant(ismember(TableRelevant.SetTemperature,all_temps(i)),:)=[];
    end
end

for i=1:numel(all_angles)
    if (all_angles_counts(i) < 5)
        TableRelevant(ismember(TableRelevant.SetAngle,all_angles(i)),:)=[];
    end
end
%% Removes curves without B=0
TableRelevant = sortrows(TableRelevant, ["SetAngle","SetTemperature","SetField"]);
n_rows_relevant = numel(TableRelevant.CriticalCurrent);
toberemoved = zeros(n_rows_relevant,1);
i=1;
while(i<=n_rows_relevant)%keeps only curves that dont have B=0 data point
    temp_buffer = TableRelevant.SetTemperature(i);
    if (TableRelevant.SetField(i) == 0)
        while (TableRelevant.SetTemperature(i) == temp_buffer)
            i=i+1;
            if (i>n_rows_relevant)
                break
            end
        end
        continue
    else
        while (TableRelevant.SetTemperature(i) == temp_buffer)
            toberemoved(i)=1;
            i=i+1;
            if (i>n_rows_relevant)
                break
            end
        end
        continue
    end
end

n_rows_relevant = numel(TableRelevant.CriticalCurrent);
i=1;
j=1;
while (i<=n_rows_relevant)
    if (toberemoved(i) == 1)
        TableRelevant(j,:) = [];
    else
        j=j+1;
    end
    i=i+1;
end
TableRelevant = sortrows(TableRelevant, ["SetTemperature","SetField","SetAngle"]);
end

function [Table] = filter_table(Table,T1,T2,B1,B2,A1,A2)
if T1 ~= -1
    Table = Table(and(Table.SetTemperature >= T1,Table.SetTemperature <= T2),:);
end

if B1 ~= -1
    Table = Table(and(Table.SetField >= B1,Table.SetField <= B2),:);
end

if A1 ~= -1
    Table = Table(and(Table.SetAngle >= A1,Table.SetAngle <= A2),:);
end
end

function Samples = add_data_sample_table(Samples,T1,T2,B1,B2,A1,A2)
Samples = addvars(Samples,zeros(numel(Samples.(1)),1),zeros(numel(Samples.(1)),1));
Samples.Properties.VariableNames(end-1:end) = ["Data","DataRelevant"];
Samples.Data = num2cell(Samples.Data);
Samples.DataRelevant = num2cell(Samples.DataRelevant);

for i = 1:numel(Samples.(1))
    [Table,TableRelevant] = read_alldata(Samples.FileDirectory(i),T1,T2,B1,B2,A1,A2,Samples.WidthEff(i));
    Samples.Data{i} = Table;
    Samples.DataRelevant{i} = TableRelevant;
end
end

function Samples = add_kim_sample_table(Samples)
Samples = addvars(Samples,zeros(numel(Samples.(1)),1));
Samples.Properties.VariableNames(end) = ["KimParam"];
Samples.KimParam = num2cell(Samples.KimParam);

for i = 1:numel(Samples.(1))
    disp("Adding Kim from "+Samples.FileDirectory(i))
    ParameterTable = compute_kim_param_interp(Samples.DataRelevant{i});
    Samples.KimParam{i} = ParameterTable;
end
end

function Samples = add_custom_sample_table(Samples)
Samples = addvars(Samples,zeros(numel(Samples.(1)),1));
Samples.Properties.VariableNames(end) = "CustomParam";
Samples.CustomParam = num2cell(Samples.CustomParam);

for i = 1:numel(Samples.(1))
    disp("Adding custom from "+Samples.FileDirectory(i))
    ParameterTable = compute_custom_param_interp_weighted(Samples.DataRelevant{i});
    Samples.CustomParam{i} = ParameterTable;
end
end

function Samples = add_advanced_custom_sample_table(Samples)
Samples = addvars(Samples,zeros(numel(Samples.(1)),1));
Samples.Properties.VariableNames(end) = "AdvancedCustomParam";
Samples.AdvancedCustomParam = num2cell(Samples.AdvancedCustomParam);

for i = 1:numel(Samples.(1))
    disp("Adding advanced custom from "+Samples.FileDirectory(i))
    ParameterTable = compute_advanced_custom_param_interp_weighted(Samples.DataRelevant{i});
    Samples.AdvancedCustomParam{i} = ParameterTable;
end
end

function [ParameterTable] = compute_kim_param_interp(TableRelevant)
%% Creates Parameter table
ParameterTable = unique(TableRelevant(:,[11,13]));
for i = 1:8
    ParameterTable = addvars(ParameterTable,zeros(numel(ParameterTable.(1)),1),'After',2);
end
ParameterTable.Properties.VariableNames(3:10) = ["IcB0","KimB0","KimPlusB0","KimPlusbeta","curveKim","curveKimPlus","gofKim","gofKimPlus"];
ParameterTable.curveKim = num2cell(ParameterTable.curveKim);
ParameterTable.gofKim = num2cell(ParameterTable.gofKim);
ParameterTable.curveKimPlus = num2cell(ParameterTable.curveKimPlus);
ParameterTable.gofKimPlus = num2cell(ParameterTable.gofKimPlus);
%% Computes and adds parameters to parameter table
max_field = max(TableRelevant.SetField);
relevant_temps = unique(TableRelevant.SetTemperature);
relevant_angles = unique(TableRelevant.SetAngle);

fit_kim_o = fitoptions('Method','NonlinearLeastSquares',...
    ...%         B0
    'Lower',     0,...
    'Upper',     30,...
    'StartPoint',max_field/2);
fit_kim_t = fittype('IcB0/(1+B/B0)','problem','IcB0','independent','B','options',fit_kim_o);

fit_kim_plus_o = fitoptions('Method','NonlinearLeastSquares',...
    ...%          B0           beta
    'Lower',[     0,           0.5],...
    'Upper',[     30,          1.5],...
    'StartPoint',[max_field/2, 1]);
fit_kim_plus_t = fittype('IcB0/(1+B/B0)^beta','problem','IcB0','independent','B','options',fit_kim_plus_o);

for i = 1:numel(relevant_temps)
    for j = 1:numel(relevant_angles)
        i_pt = find(and(ParameterTable.SetTemperature == relevant_temps(i),ParameterTable.SetAngle == relevant_angles(j)));
        if isempty(i_pt)
            continue
        end
        Table_current = TableRelevant(and(TableRelevant.SetTemperature == relevant_temps(i),TableRelevant.SetAngle == relevant_angles(j)),:);
        SetField_interp = min(Table_current.SetField):Table_current.SetField(2)-Table_current.SetField(1):max(Table_current.SetField);
        Current_interp = interp1(Table_current.SetField,Table_current.CriticalCurrent,SetField_interp);

        ParameterTable.IcB0(i_pt) = Table_current.CriticalCurrent(Table_current.SetField == 0);

        [curve_fitted_buffer,gof_buffer] = fit(SetField_interp',Current_interp',fit_kim_t,'problem',ParameterTable.IcB0(i_pt));
        parameters_buffer = coeffvalues(curve_fitted_buffer);
        ParameterTable.KimB0(i_pt) = parameters_buffer;
        ParameterTable.curveKim{i_pt} = curve_fitted_buffer;
        ParameterTable.gofKim{i_pt} = gof_buffer;

        [curve_fitted_buffer,gof_buffer] = fit(SetField_interp',Current_interp',fit_kim_plus_t,'problem',ParameterTable.IcB0(i_pt));
        parameters_buffer = coeffvalues(curve_fitted_buffer);
        ParameterTable.KimPlusB0(i_pt) = parameters_buffer(1);
        ParameterTable.KimPlusbeta(i_pt) = parameters_buffer(2);
        ParameterTable.curveKimPlus{i_pt} = curve_fitted_buffer;
        ParameterTable.gofKimPlus{i_pt} = gof_buffer;
    end
end
end

function [ParameterTable] = compute_kim_aniso_param(TableRelevant)
%% Creates Parameter table
ParameterTable = table(unique(TableRelevant.SetTemperature));
for i = 1:6
    ParameterTable = addvars(ParameterTable,zeros(numel(ParameterTable.(1)),1),'After',1);
end
ParameterTable.Properties.VariableNames = ["SetTemperature","IcB0","KimAnisoB0","KimAnisobeta","KimAnisoka","curveKimAniso","gofKimAniso"];
ParameterTable.curveKimAniso = num2cell(ParameterTable.curveKimAniso);
ParameterTable.gofKimAniso = num2cell(ParameterTable.gofKimAniso);
%% Computes and adds parameters to parameter table
max_field = max(TableRelevant.SetField);
relevant_temps = unique(TableRelevant.SetTemperature);
fit_kim_aniso_o = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[inf,10,10],...
               'StartPoint',[max_field / 2,1,1]);
fit_kim_aniso_t = fittype('IcB0/(1+B*sqrt(cos(angle)^2+ka^2*sin(angle)^2)/B0)^beta','problem','IcB0','independent',{'B','angle'},...
'options',fit_kim_aniso_o);

for i = 1:numel(relevant_temps)%computes parameters for every angle
    Table_current = TableRelevant(TableRelevant.SetTemperature == relevant_temps(i),:);
    
    ParameterTable.IcB0(i) = TableRelevant.CriticalCurrent(and(and(TableRelevant.SetField == 0,TableRelevant.SetAngle == 0),TableRelevant.SetTemperature == relevant_temps(i)));
    
    [curve_fitted_buffer,gof_buffer] = fit([Table_current.SetField, Table_current.SetAngle * pi / 180],...
        Table_current.CriticalCurrent, fit_kim_aniso_t, 'problem', ParameterTable.IcB0(i));
    parameters_buffer = coeffvalues(curve_fitted_buffer);
    ParameterTable.KimAnisoB0(i) = parameters_buffer(1);
    ParameterTable.KimAnisobeta(i) = parameters_buffer(2);
    ParameterTable.KimAnisoka(i) = parameters_buffer(3);
    ParameterTable.curveKimAniso{i} = curve_fitted_buffer;
    ParameterTable.gofKimAniso{i} = gof_buffer;
end
end

function [ParameterTable] = compute_custom_param_interp_weighted(TableRelevant)
%% Creates Parameter table
ParameterTable = table(unique(TableRelevant.SetTemperature));
for i = 1:11
    ParameterTable = addvars(ParameterTable,zeros(numel(ParameterTable.(1)),1));
end
ParameterTable.Properties.VariableNames = ...
    ["SetTemperature","IcB0","B0","beta","SFa","SFw","SFp","IRa","IRw","IRp","curve","gof"];
ParameterTable.curve = num2cell(ParameterTable.curve);
ParameterTable.gof = num2cell(ParameterTable.gof);
%% Computes and adds parameters to parameter table
relevant_temps = unique(TableRelevant.SetTemperature);
fit_kim_custom_o = fitoptions('Method','NonlinearLeastSquares',...
               ...%          B0,  IRa, IRp,   IRw,  SFa, SFp,        SFw,  beta
               'Lower',     [0,   0,   -pi/6, 0.4,  0,   -15*pi/180, 1,    0],...
               'Upper',     [30,  0.9, pi/6,  1000, 0.9, 15*pi/180,  1000, 10],...
               'StartPoint',[1.5, 0,   0,     20,   0,   0,          20,   1]);
fit_kim_custom_t = ...
    fittype('IcB0/(1+B*(1-IRa*exp(IRw*(cos(2*(angle-IRp))-1))-SFa*exp(SFw*(cos(2*(angle-pi/2-SFp))-1)))/B0)^beta',...
    'problem','IcB0','independent',{'B','angle'},'options',fit_kim_custom_o);

for i = 1:numel(relevant_temps)%computes parameters for every angle
    Table_current = TableRelevant(TableRelevant.SetTemperature == relevant_temps(i),:);
    relevant_fields = unique(Table_current.SetField);
    relevant_angles = unique(Table_current.SetAngle);
    i_pt = find(ParameterTable.SetTemperature == relevant_temps(i));
    if any(numel(relevant_angles) == [1 2 3])
        disp("Skipped"+num2str(relevant_temps(i)))
        ParameterTable(i_pt,:) = [];
        continue
    end
    
    ParameterTable.IcB0(i_pt) = ...
        TableRelevant.CriticalCurrent(and(and(TableRelevant.SetField == 0,TableRelevant.SetAngle == 0),...
        TableRelevant.SetTemperature == relevant_temps(i)));

    [X,Y] = meshgrid(relevant_fields,relevant_angles);
    clearvars Z;
    for j = 1:numel(relevant_fields)
        for k = 1:numel(relevant_angles)
            i_tr = find(and(Table_current.SetTemperature == relevant_temps(i),...
                and(Table_current.SetField == relevant_fields(j),Table_current.SetAngle == relevant_angles(k))));
            if isempty(i_tr)
                Z(k,j) = NaN;
                continue
            end
            Z(k,j) = Table_current.CriticalCurrent(i_tr);
        end
    end

    el_max = numel(Z);
    el_nan = sum(transpose(sum(isnan(Z))));
    if el_nan > 0
        disp(num2str(relevant_temps(i)))
        disp(num2str(el_max))
        disp(num2str(el_nan))
        Z = fillmissing(Z,'linear',2,'EndValues','extrap');
    end

    relevant_fields_interp = min(relevant_fields):relevant_fields(2)-relevant_fields(1):max(relevant_fields);
    relevant_angles_interp = min(relevant_angles):1:max(relevant_angles);
    [X_interp,Y_interp] = meshgrid(relevant_fields_interp,relevant_angles_interp);
    Z_interp = interp2(X,Y,Z,X_interp,Y_interp);
    W_interp = ones(numel(X_interp(:,1)),numel(X_interp(1,:)));

    gap0 = 30;
    gap90 = 30;
    magnet = 75;
    high = 0.2;
    for j = 1:numel(relevant_fields_interp)
        idx_angles_peaks(1,:) = Y_interp(:,j) <= gap0;
        idx_angles_peaks(2,:) = (Y_interp(:,j) >= 90 - gap90) & (Y_interp(:,j) <= 90 + gap90);
        idx_angles_peaks(3,:) = (Y_interp(:,j) >= 180 - gap0) & (Y_interp(:,j) <= 180 + gap0);
        idx_angles_middle(1) = find(Y_interp(:,j) == 45);
        idx_angles_middle(2) = find(Y_interp(:,j) == 135);

        max_current(1) = max(Z_interp(idx_angles_peaks(1,:),j));
        max_current(2) = max(Z_interp(idx_angles_peaks(2,:),j));
        max_current(3) = max(Z_interp(idx_angles_peaks(3,:),j));
        max_current(4) = Z_interp(idx_angles_middle(1),j);
        max_current(5) = Z_interp(idx_angles_middle(2),j);

        if Y_interp(find(Z_interp(:,j) == max_current(1)),j) == gap0
            max_current(1) = Z_interp(find(relevant_angles_interp == 0),j);
            max_current(3) = Z_interp(find(relevant_angles_interp == 180),j);
        end

        W_interp(find(sum(Z_interp(:,j) == max_current,2)),j) = magnet;
    end
    W_interp(:,find(X_interp(1,:) <= high)) = 1;

    SetField_interp = reshape(X_interp,[],1);
    SetAngle_interp = reshape(Y_interp,[],1);
    Current_interp = reshape(Z_interp,[],1);
    Weights_interp = reshape(W_interp,[],1);

    [curve_fitted_buffer,gof_buffer] = fit([SetField_interp, SetAngle_interp * pi / 180],...
        Current_interp, fit_kim_custom_t, 'problem', ParameterTable.IcB0(i_pt), 'Weights', Weights_interp);
    parameters_buffer = coeffvalues(curve_fitted_buffer);
    ParameterTable.B0(i_pt) = parameters_buffer(1);
    ParameterTable.IRa(i_pt) = parameters_buffer(2);
    ParameterTable.IRp(i_pt) = parameters_buffer(3);
    ParameterTable.IRw(i_pt) = parameters_buffer(4);
    ParameterTable.SFa(i_pt) = parameters_buffer(5);
    ParameterTable.SFp(i_pt) = parameters_buffer(6);
    ParameterTable.SFw(i_pt) = parameters_buffer(7);
    ParameterTable.beta(i_pt) = parameters_buffer(8);
    ParameterTable.curve{i_pt} = curve_fitted_buffer;
    ParameterTable.gof{i_pt} = gof_buffer;
end
end

function [ParameterTable] = compute_advanced_custom_param_interp_weighted(TableRelevant)
%% Creates Parameter table
ParameterTable = table(unique(TableRelevant.SetTemperature));
for i = 1:16
    ParameterTable = addvars(ParameterTable,zeros(numel(ParameterTable.(1)),1));
end
ParameterTable.Properties.VariableNames = ...
    ["SetTemperature","IcB0","B0","beta","SFa","SFwi","SFwf","SFws","SFpi","SFpf","SFps","IRaf","IRas","IRw","IRp","curve","gof"];
ParameterTable.curve = num2cell(ParameterTable.curve);
ParameterTable.gof = num2cell(ParameterTable.gof);
%% Computes and adds parameters to parameter table
relevant_temps = unique(TableRelevant.SetTemperature);
fit_kim_custom_o = fitoptions('Method','NonlinearLeastSquares',...
               ...%          B0,  IRaf, IRas, IRp,   IRw,  SFa, SFpf,      SFpi,      SFps, SFwf, SFwi, SFws, beta
               'Lower',     [0,   0,    0.5,  -pi/6, 1,    0,   -5*pi/180, 0,         1,    0,    0.6,  0.1,  0],...
               'Upper',     [30,  0.9,  1000, pi/6,  1000, 0.9, 5*pi/180,  15*pi/180, inf,  inf,  1000, 100,  10],...
               'StartPoint',[1.5, 0,    1,    0,     20,   0,   0,         0,         100,  5,    20,   2,    1]);
fit_kim_custom_t = ...
    fittype('IcB0/(1+B*(1-tanh(IRas*B)*IRaf*exp(IRw*(cos(2*(angle-IRp))-1))-SFa*exp((SFwi+SFwf*tanh(SFws*B))*(cos(2*(angle-pi/2-(-SFpi*exp(-SFps*B)+SFpf)))-1)))/B0)^beta',...
    'problem','IcB0','independent',{'B','angle'},'options',fit_kim_custom_o);

for i = 1:numel(relevant_temps)%computes parameters for every angle
    Table_current = TableRelevant(TableRelevant.SetTemperature == relevant_temps(i),:);
    relevant_fields = unique(Table_current.SetField);
    relevant_angles = unique(Table_current.SetAngle);
    i_pt = find(ParameterTable.SetTemperature == relevant_temps(i));
    if any(numel(relevant_angles) == [1 2 3])
        disp("Skipped"+num2str(relevant_temps(i)))
        ParameterTable(i_pt,:) = [];
        continue
    end
    
    ParameterTable.IcB0(i_pt) = ...
        TableRelevant.CriticalCurrent(and(and(TableRelevant.SetField == 0,TableRelevant.SetAngle == 0),...
        TableRelevant.SetTemperature == relevant_temps(i)));

    [X,Y] = meshgrid(relevant_fields,relevant_angles);
    clearvars Z;
    for j = 1:numel(relevant_fields)
        for k = 1:numel(relevant_angles)
            i_tr = find(and(Table_current.SetTemperature == relevant_temps(i),...
                and(Table_current.SetField == relevant_fields(j),Table_current.SetAngle == relevant_angles(k))));
            if isempty(i_tr)
                Z(k,j) = NaN;
                continue
            end
            Z(k,j) = Table_current.CriticalCurrent(i_tr);
        end
    end

    el_max = numel(Z);
    el_nan = sum(transpose(sum(isnan(Z))));
    if el_nan > 0
        disp(num2str(relevant_temps(i)))
        disp(num2str(el_max))
        disp(num2str(el_nan))
        Z = fillmissing(Z,'linear',2,'EndValues','extrap');
    end

    relevant_fields_interp = min(relevant_fields):relevant_fields(2)-relevant_fields(1):max(relevant_fields);
    relevant_angles_interp = min(relevant_angles):1:max(relevant_angles);
    [X_interp,Y_interp] = meshgrid(relevant_fields_interp,relevant_angles_interp);
    Z_interp = interp2(X,Y,Z,X_interp,Y_interp);
    W_interp = ones(numel(X_interp(:,1)),numel(X_interp(1,:)));

    gap0 = 30;
    gap90 = 30;
    magnet = 75;
    high = 0.2;
    for j = 1:numel(relevant_fields_interp)
        idx_angles_peaks(1,:) = Y_interp(:,j) <= gap0;
        idx_angles_peaks(2,:) = (Y_interp(:,j) >= 90 - gap90) & (Y_interp(:,j) <= 90 + gap90);
        idx_angles_peaks(3,:) = (Y_interp(:,j) >= 180 - gap0) & (Y_interp(:,j) <= 180 + gap0);
        idx_angles_middle(1) = find(Y_interp(:,j) == 45);
        idx_angles_middle(2) = find(Y_interp(:,j) == 135);

        max_current(1) = max(Z_interp(idx_angles_peaks(1,:),j));
        max_current(2) = max(Z_interp(idx_angles_peaks(2,:),j));
        max_current(3) = max(Z_interp(idx_angles_peaks(3,:),j));
        max_current(4) = Z_interp(idx_angles_middle(1),j);
        max_current(5) = Z_interp(idx_angles_middle(2),j);

        if Y_interp(find(Z_interp(:,j) == max_current(1)),j) == gap0
            max_current(1) = Z_interp(find(relevant_angles_interp == 0),j);
            max_current(3) = Z_interp(find(relevant_angles_interp == 180),j);
        end

        W_interp(find(sum(Z_interp(:,j) == max_current,2)),j) = magnet;
    end
    W_interp(:,find(X_interp(1,:) <= high)) = 1;

    SetField_interp = reshape(X_interp,[],1);
    SetAngle_interp = reshape(Y_interp,[],1);
    Current_interp = reshape(Z_interp,[],1);
    Weights_interp = reshape(W_interp,[],1);

    [curve_fitted_buffer,gof_buffer] = fit([SetField_interp, SetAngle_interp * pi / 180],...
        Current_interp, fit_kim_custom_t, 'problem', ParameterTable.IcB0(i_pt), 'Weights', Weights_interp);
    parameters_buffer = coeffvalues(curve_fitted_buffer);
    ParameterTable.B0(i_pt) = parameters_buffer(1);
    ParameterTable.IRaf(i_pt) = parameters_buffer(2);
    ParameterTable.IRas(i_pt) = parameters_buffer(3);
    ParameterTable.IRp(i_pt) = parameters_buffer(4);
    ParameterTable.IRw(i_pt) = parameters_buffer(5);
    ParameterTable.SFa(i_pt) = parameters_buffer(6);
    ParameterTable.SFpf(i_pt) = parameters_buffer(7);
    ParameterTable.SFpi(i_pt) = parameters_buffer(8);
    ParameterTable.SFps(i_pt) = parameters_buffer(9);
    ParameterTable.SFwf(i_pt) = parameters_buffer(10);
    ParameterTable.SFwi(i_pt) = parameters_buffer(11);
    ParameterTable.SFws(i_pt) = parameters_buffer(12);
    ParameterTable.beta(i_pt) = parameters_buffer(13);
    ParameterTable.curve{i_pt} = curve_fitted_buffer;
    ParameterTable.gof{i_pt} = gof_buffer;
end
end