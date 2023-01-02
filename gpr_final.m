%% Create GPR model using training data
clc
clear all
close all

disp('Loading Data')
load AllSamplesAllModels3T.mat
disp('Data Loaded')

Table = Samples.Data{1};

X_total = [Table.SetTemperature Table.SetField Table.SetAngle];
Y_total = Table.CriticalCurrent;

idx_train = mod(abs(Table.SetTemperature), 10) == 0;
X_train = [Table.SetTemperature(idx_train) Table.SetField(idx_train) Table.SetAngle(idx_train)];
Y_train = Table.CriticalCurrent(idx_train);

disp('Started training model')

% gpr_model = fitrgp(X_train,Y_train,...
%     'PredictMethod','exact',...
%     'FitMethod','sd',...
%     'Standardize',1,...
%     'OptimizeHyperparameters',{'BasisFunction','KernelFunction','KernelScale','Sigma'});

% gpr_model = fitrgp(X_train,Y_train,...
%     'PredictMethod','exact',...
%     'FitMethod','sd',...
%     'Standardize',1,...
%     'KernelFunction','exponential',...
%     'OptimizeHyperparameters',{'BasisFunction','KernelScale','Sigma'});

% gpr_model = fitrgp(X_train,Y_train,...
%     'PredictMethod','exact',...
%     'FitMethod','sd',...
%     'Standardize',1,...
%     'KernelFunction','ardrationalquadratic',...
%     'OptimizeHyperparameters',{'BasisFunction','Sigma'});

gpr_model = fitrgp(X_train,Y_train,...
    'PredictMethod','exact',...
    'FitMethod','sd',...
    'Standardize',0,...
    'KernelFunction','exponential');

disp('Model is trained')
save('NoStandardizationEstimatedExponentialSuperPower3T','gpr_model')
%% Plots
clc
clear all
close all

disp('Loading Data and Model')
load AllSamplesAllModels3T.mat

%load AutoAllSuperPower3T.mat
%load AutoAllSuperOxGdBCO3T.mat
%load AutoAllSuperOxYBCO3T.mat

%load ExponentialAutoRestSuperPower3T.mat
%load Matern32AutoRestSuperPower3T.mat
%load Matern52AutoRestSuperPower3T.mat
%load SquaredExponentialAutoRestSuperPower3T.mat
%load RationalQuadraticAutoRestSuperPower3T.mat

%load ARDExponentialAutoRestSuperPower3T.mat
%load ARDMatern32AutoRestSuperPower3T.mat
%load ARDMatern52AutoRestSuperPower3T.mat
%load ARDSquaredExponentialAutoRestSuperPower3T.mat
load ARDRationalQuadraticAutoRestSuperPower3T.mat

%load NonEstimatedExponentialSuperPower3T.mat
%load EstimatedExponentialSuperPower3T.mat

%load('./PIC2/ExponentialGeneric.mat')
%gpr_model = ExponentialGeneric;

%load('./PIC2/ExponentialSD.mat')
%gpr_model = ExponentialSD;

%load NoStandardizationEstimatedExponentialSuperPower3T.mat

disp('Data and Model Loaded')
idx = 1;
label_ml = gpr_model.KernelFunction;

simple_string = label_ml;
if simple_string(1,1) == 'A'
    simple_string = label_ml(4:end);
end
switch simple_string
    case 'Exponential'
        colour = "#4DBEEE";
    case 'Matern32'
        colour = "#7E2F8E";
    case 'Matern52'
        colour = "#EDB120";
    case 'SquaredExponential'
        colour = "#0072BD";
    case 'RationalQuadratic'
        colour = "#77AC30";
end

plot_inputoutput_gpr(Samples.Data{idx},Samples.Label{idx},gpr_model,label_ml,colour)
plot_curve_gpr(Samples.Data{idx},Samples.Label(idx),Samples.WidthEff(idx),gpr_model,label_ml,[20 45 60 70 77.5],[0 45 90],"log");
plot_anglecurve_gpr(Samples.Data{idx},Samples.Label{idx},Samples.WidthEff(idx),gpr_model,label_ml,[20 45 60 77.5],[0.2 0.5 1 2 3])
plot_relative_error_gpr(Samples.DataRelevant{idx},Samples.Label{idx},gpr_model,label_ml,[0 45 90])
disp('Finished saving figures')
%% Problem plots
clc
clear all
close all

disp('Loading Data')
load AllSamplesAllModels3T.mat
load('./PIC2/ExponentialGeneric.mat')
load('./PIC2/ExponentialSD.mat')
disp('Data Loaded')

Table = Samples.Data{1};

X_total = [Table.SetTemperature Table.SetField Table.SetAngle];
Y_total = Table.CriticalCurrent;

idx_train = mod(abs(Table.SetTemperature), 10) == 0;
X_train = [Table.SetTemperature(idx_train) Table.SetField(idx_train) Table.SetAngle(idx_train)];
Y_train = Table.CriticalCurrent(idx_train);

figure
hold on
set(gca,'FontSize',12)
plot(Table.CriticalCurrent,'ro','markersize',10)
plot(predict(ExponentialGeneric,X_total),'g.','markersize',12)
plot(predict(ExponentialSD,X_total),'b.','markersize',12)
legend('Database','Non-Estimated Params','Estimated Params','Location','best')
xlabel('Input #')
ylabel('Critical Current [A]')
grid on
hold off

temps_toplot = 45;
fields_toplot = 2;
points_toplot = [linspace(0,0,500)'+temps_toplot, linspace(0,0,500)'+fields_toplot, linspace(0,240, 500)'];
idx_tr = Table.SetTemperature == temps_toplot & Table.SetField == fields_toplot;

[CurrentPrediction,~,CurrentInterval] = predict(ExponentialSD,points_toplot,'Alpha',0.05);

figure
hold on
set(groot,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',12)
plot(Table.SetAngle(idx_tr),Table.CriticalCurrent(idx_tr),'r*')
plot(points_toplot(:,3),CurrentPrediction,'blue')
patch([points_toplot(:,3);flipud(points_toplot(:,3))],[CurrentInterval(:,1);flipud(CurrentInterval(:,2))],'blue','FaceAlpha',0.1);
legend('Database','Standardized/Estimated Params','Location','best')
xlabel('Angle [º]')
ylabel('Critical Current [A]')
grid on
hold off
%% Functions
function plot_inputoutput_gpr(Table,label_fita,gpr_model,label_ml,colour)
X_total = [Table.SetTemperature Table.SetField Table.SetAngle];

figure
hold on
set(groot,'DefaultLineLineWidth',1.5)
%title(label_fita+" input for "+label_ml)
plot(Table.CriticalCurrent,'red')
plot(predict(gpr_model,X_total),'Color',colour)
legend('Database',label_ml,'Location','best')
xlabel('Input #')
ylabel('Critical Current [A]')
grid on
set(gca,'FontSize',12)
savestring = label_fita+" inputoutput "+label_ml+".png";
exportgraphics(gcf,"./Plots ML/"+savestring)
end

function plot_curve_gpr(Table,label_fita,width,gpr_model,label_ml,temps_toplot,angles_toplot,scale)
max_field = max(Table.SetField);
B = linspace(0,max_field,500);

for i=1:numel(angles_toplot)
    figure
    hold on
    %title(label_fita+"'s "+width+"cm width:"+label_ml+"Ic/w vs B at "+num2str(angles_toplot(i))+"º angle")
    set(groot,'DefaultLineLineWidth',1.5)
    if scale == "log"
        set(gca,'YScale','log')
    end
    for j=1:numel(temps_toplot)
        points_toplot = [linspace(0,0,500)'+temps_toplot(j), B', linspace(0,0, 500)'+angles_toplot(i)];

        table_buffer = Table(Table.SetTemperature == temps_toplot(j) & Table.SetAngle == angles_toplot(i),:);
        plot(table_buffer.SetField, table_buffer.CriticalCurrent/width,"red")
        prediction = predict(gpr_model,points_toplot)/width;
        plot(B,prediction,"blue")
        height = predict(gpr_model,[temps_toplot(j) 4 angles_toplot(i)])/width+0.4*predict(gpr_model,[temps_toplot(j) 4 angles_toplot(i)])/width;
        text(max_field/2,height,num2str(temps_toplot(j))+"K")
    end
    legend('Database',label_ml,'Location','best')
    xlabel('Magnetic Field [T]')
    ylabel('Critical Current Sheet Density [A/cm]')
    grid on
    set(gca,'FontSize',12)
    savestring = label_fita+" curve "+label_ml+" "+num2str(angles_toplot(i))+".png";
    exportgraphics(gcf,"./Plots ML/"+savestring)
end
end

function plot_anglecurve_gpr(Table,label_fita,width,gpr_model,label_ml,temps_toplot,fields_toplot)
for i=1:numel(temps_toplot)
    for j=1:numel(fields_toplot)
        points_toplot = [linspace(0,0,500)'+temps_toplot(i), linspace(0,0,500)'+fields_toplot(j), linspace(0,240, 500)'];
        idx_tr = Table.SetTemperature == temps_toplot(i) & Table.SetField == fields_toplot(j);

        figure
        hold on
        %title(label_fita+" angle dependency for "+label_ml+" "+num2str(temps_toplot(i))+"K and "+fields_toplot(j)+"T")
        plot(Table.SetAngle(idx_tr),Table.CriticalCurrent(idx_tr)/width,'red')
        plot(points_toplot(:,3),predict(gpr_model,points_toplot)/width,'blue')
        legend('Database',label_ml,'Location','best')
        xlabel('Angle [º]')
        ylabel('Critical Current [A]')
        grid on
        set(gca,'FontSize',12)
        savestring = label_fita+" angle curve "+label_ml+" "+num2str(temps_toplot(i))+"K "+num2str(fields_toplot(j))+"T.png";
        exportgraphics(gcf,"./Plots ML/"+savestring)
    end
end
end

function plot_relative_error_gpr(TableRelevant,label_fita,gpr_model,label_ml,angles_toplot)
%% Compute and filter relative error
create_figure = [0,1];
filt_limit = 10;

relevant_temps = unique(TableRelevant.SetTemperature);
if numel(relevant_temps) == 1
    create_figure = [1,0];
end
relevant_fields = unique(TableRelevant.SetField);
[X,Y] = meshgrid(relevant_temps,relevant_fields);

for k = 1:numel(angles_toplot)
    for i = 1:numel(relevant_temps)
        for j = 1:numel(relevant_fields)
            idx_tr = find(and(TableRelevant.SetTemperature == relevant_temps(i),...
                and(TableRelevant.SetField == relevant_fields(j),TableRelevant.SetAngle == angles_toplot(k))));
            if isempty(idx_tr)
                Z(j,i,k) = NaN;
                continue
            end
            Z(j,i,k) = 100 * abs((TableRelevant.CriticalCurrent(idx_tr) - predict(gpr_model,[relevant_temps(i),relevant_fields(j),angles_toplot(k)]))...
                / TableRelevant.CriticalCurrent(idx_tr));
        end
    end
end
Z_filt = Z;
Z_filt(Z_filt > filt_limit) = NaN;
if create_figure(1) == 1
    for k = 1:numel(angles_toplot)
        if numel(relevant_temps) == 1
            figure
            plot(Y,Z_filt(:,1,k))
            %title("GPR <10% relative error "+num2str(angles_toplot(k))+"º "+relevant_temps(1)+"K")
            ylabel('Relative Error [%]')
            xlabel('Magnetic Field [T]')
            grid on
            set(gca,'FontSize',12)
        else
            figure
            surf(X,Y,Z_filt(:,:,k))
            %title("GPR <10% relative error "+num2str(angles_toplot(k))+"º")
            zlabel('Relative Error [%]')
            ylabel('Magnetic Field [T]')
            xlabel('Temperature [K]')
            grid on
            c = colorbar;   c.Label.String = 'Relative Error [%]';  c.Label.FontSize = 12;
            colormap cool
            set(gca,'FontSize',12)
            view(gca,[0.144782608695627 90]);
        end
    end
end
%% Interpolate and filter relative error
if create_figure(2) == 0
    return
end
interp_res = 6;
relevant_temps_interp = min(relevant_temps):(relevant_temps(2)-relevant_temps(1))/interp_res:max(relevant_temps);
n_fields = numel(relevant_fields);
relevant_fields_interp = min(relevant_fields):(relevant_fields(n_fields)-relevant_fields(n_fields-1))/2/interp_res:max(relevant_fields);
[X_interp,Y_interp] = meshgrid(relevant_temps_interp,relevant_fields_interp);

for k=1:numel(angles_toplot)
    Z_interp = interp2(X,Y,Z(:,:,k),X_interp,Y_interp);
    Z_interp_filt = Z_interp;
    Z_interp_filt(Z_interp_filt >= filt_limit) = NaN;

    if create_figure(2) == 1
        figure
        surf(X_interp,Y_interp,Z_interp_filt)
        %title(label_fita+" <10% relative error "+label_ml+" "+num2str(angles_toplot(k))+"º interpolated")
        zlabel('Relative Error [%]')
        ylabel('Magnetic Field [T]')
        xlabel('Temperature [K]')
        xlim([min(TableRelevant.SetTemperature) max(TableRelevant.SetTemperature)])
        ylim([0 3])
        grid on
        c = colorbar;   c.Label.String = 'Relative Error [%]';  c.Label.FontSize = 12;
        colormap cool
        set(gca,'FontSize',12)
        view(gca,[0.144782608695627 90]);
        savestring = label_fita+" Rel Error "+label_ml+" "+num2str(angles_toplot(k))+".png";
        exportgraphics(gcf,"./Plots ML/"+savestring)
    end
end
end