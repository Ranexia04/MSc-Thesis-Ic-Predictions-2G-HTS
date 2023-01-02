%% Analyse a single model
clc
clear all
close all

load AllSamplesAllModels3T.mat
i = 1;

plot_relative_error_kim(Samples.DataRelevant{i},Samples.KimParam{i},Samples.Label{i},[0 90],"Kim Plus");
plot_param_kim(Samples.DataRelevant{i},Samples.KimParam{i},Samples.Label{i});
plot_kim_curve(Samples.Data{i},Samples.KimParam{i},Samples.WidthEff(i),Samples.Label{i},[77.5], [0 90],"Kim Plus","normal");
%% Analyse all Kim sequentially
clc
clear all
close all

load AllSamplesAllModels3T.mat

for i = 1:3%numel(Samples.(1))
    close all
    disp("Displaying Kim plots for "+Samples.Label(i))
    plot_relative_error_kim(Samples.DataRelevant{i},Samples.KimParam{i},Samples.Label(i),[0 45 90],["Kim","Kim Plus"]);
    plot_param_kim(Samples.DataRelevant{i},Samples.KimParam{i},Samples.Label(i),["Kim","Kim Plus"]);
    plot_kim_curve(Samples.Data{i},Samples.KimParam{i},Samples.Label(i),Samples.WidthEff(i),[20 40 60 70 77.5],[0 45 90],["Kim","Kim Plus"],"log");
    disp('Waiting for user input.')
    %pause
end
%% Analyse all Custom sequentially
clc
clear all
close all

load AllSamplesAllModels3T.mat
for i = 1:3%numel(Samples.(1))
    close all
    disp("Displaying Custom plots for "+Samples.Label(i))
    plot_relative_error_custom(Samples.DataRelevant{i},Samples.CustomParam{i},Samples.Label(i),"Custom");
    plot_param_custom(Samples.CustomParam{i},Samples.Label(i));
    %plot_custom_curve(Samples.Data{i},Samples.CustomParam{i},Samples.Label(i),Samples.WidthEff(i));
    plot_custom_angle(Samples.DataRelevant{i},Samples.CustomParam{i},Samples.Label(i),Samples.WidthEff(i),[20 40 50 77.5],[0.2 0.5 1 2 3],"Custom");
    disp('Waiting for user input.')
    %pause
end
%% Analyse all Advanced Custom sequentially
clc
clear all
close all

load AllSamplesAllModels3T.mat
%load AllSamplesAllModelsIRw13T.mat
for i = 1:1%numel(Samples.(1))
    close all
    disp("Displaying Advanced Custom plots for "+Samples.Label(i))
    plot_relative_error_custom(Samples.DataRelevant{i},Samples.AdvancedCustomParam{i},Samples.Label(i),"Advanced Custom");
    %plot_param_advanced_custom(Samples.AdvancedCustomParam{i},Samples.Label(i));
    %plot_custom_curve(Samples.Data{i},Samples.AdvancedCustomParam{i},Samples.Label(i),Samples.WidthEff(i));
    %plot_custom_angle(Samples.DataRelevant{i},Samples.AdvancedCustomParam{i},Samples.Label(i),Samples.WidthEff(i),[20 40 50 77.5],[0.2 0.5 1 2 3],"Advanced Custom");
    %plot_advanced_custom_bdepend(Samples.AdvancedCustomParam{i},Samples.Label(i),[20 40 50 77.5]);
    disp('Waiting for user input.')
    %pause
end
%% Analyse angle dependency 1G
clc
clear all
close all

load AllSamplesAllModels3T.mat
Table = Samples.Data{5};

temp_toplot = 77.5;
field_toplot = 1;
Table_current = Table(and(Table.SetTemperature == temp_toplot,Table.SetField == field_toplot),:);

relevant_angles = Table_current.SetAngle;
for i = 1:numel(relevant_angles)
    i_pt = find(and(Samples.KimParam{1}.SetTemperature == temp_toplot,Samples.KimParam{1}.SetAngle == relevant_angles(i)));
    z_plus(i,1) = feval(Samples.KimParam{1}.curveKimPlus{i_pt},field_toplot);
end

icb0 = Table.CriticalCurrent(and(and(Table.SetTemperature == temp_toplot,Table.SetField == 0),Table.SetAngle == 0));
B0 = 0.58;
beta = 1;
ka = 18;
B_perp = cos(Table_current.SetAngle * pi / 180);
B_para = sin(Table_current.SetAngle * pi / 180);
z_aniso = icb0 ./ (1 + field_toplot * sqrt(B_para.^2 + ka^2 * B_perp.^2) / B0).^beta;

figure
hold on
plot3(B_perp,B_para,Table_current.CriticalCurrent)
plot3(B_perp,B_para,z_aniso)
xlabel('Bperp [pu]')
ylabel('Bpara [pu]')
zlabel('Critical Current [A]')
grid on
legend(["Database","Aniso"])
view(gca,[37.1657458563536 33.7415887850467]);
set(gca,'FontSize',15)

figure
hold on
plot(Table_current.SetAngle,Table_current.CriticalCurrent)
plot(Table_current.SetAngle,z_aniso)
xlabel('Angle [º]')
ylabel('Critical Current [A]')
grid on
legend(["Database","Aniso"])
set(gca,'FontSize',15)
%% Analyse angle dependency 2G
clc
clear all
close all

load AllSamplesAllModels3T.mat
Table = Samples.Data{1};

temp_toplot = 77.5;
field_toplot = 2;
Table_current = Table(and(Table.SetTemperature == temp_toplot,Table.SetField == field_toplot),:);

relevant_angles = Table_current.SetAngle;
for i = 1:numel(relevant_angles)
    i_pt = find(and(Samples.KimParam{1}.SetTemperature == temp_toplot,Samples.KimParam{1}.SetAngle == relevant_angles(i)));
    z_plus(i,1) = feval(Samples.KimParam{1}.curveKimPlus{i_pt},field_toplot);
end

icb0 = Table.CriticalCurrent(and(and(Table.SetTemperature == temp_toplot,Table.SetField == 0),Table.SetAngle == 0));
B0 = 0.0517207406709136;
beta = 0.687079719075761;
ka = 18;
B_perp = cos(Table_current.SetAngle * pi / 180);
B_para = sin(Table_current.SetAngle * pi / 180);
z_aniso = icb0 ./ (1 + field_toplot * sqrt(B_para.^2 + ka^2 * B_perp.^2) / B0).^beta;

B_custom_mult = 1.3 * sqrt(1.8 + cos(2 * (Table_current.SetAngle - 100) * pi / 180)) - 1.2 * exp(-1.5 * ((Table_current.SetAngle - 90) * pi / 180).^2 .* ((Table_current.SetAngle - 270) * pi / 180).^2);
z_custom = icb0 ./ (1 + field_toplot * B_custom_mult / B0).^beta;
i_pt2 = Samples.CustomParam{1}.SetTemperature == temp_toplot;
z_custom_fit = feval(Samples.CustomParam{1}.curve{i_pt2},field_toplot, Table_current.SetAngle * pi / 180);

wB = ((icb0./Table_current.CriticalCurrent).^(1/beta)-1).*B0/field_toplot;

picos = -wB + 1;

set(groot,'DefaultLineLineWidth',1.5)
%set(gca,'FontSize',12)%,'YScale','log'
figure
hold on
plot3(B_perp,B_para,Table_current.CriticalCurrent);
plot3(B_perp,B_para,z_plus)
%plot3(B_perp,B_para,z_aniso)
%plot3(B_perp,B_para,z_custom)
%plot3(B_perp,B_para,z_custom_fit)
xlabel('Bperp [pu]')
ylabel('Bpara [pu]')
zlabel('Critical Current [A]')
zlim([0 12])
grid on
legend(["Database","KimPlus"])
set(gca,'FontSize',15)
view(gca,[41.2419889502762 31.3373831775701]);

figure
hold on
plot(Table_current.SetAngle,Table_current.CriticalCurrent)
plot(Table_current.SetAngle,z_plus)
%plot(Table_current.SetAngle,z_aniso)
%plot(Table_current.SetAngle,z_custom)
%plot(Table_current.SetAngle,z_custom_fit)
xlabel('Angle [º]')
ylabel('Critical Current [A]')
ylim([0 12])
grid on
legend(["Database","KimPlus"])
set(gca,'FontSize',15)

figure
hold on
plot(Table_current.SetAngle,wB)
%plot(Table_current.SetAngle, B_mult)
xlabel('Angle \phi_B [º]')
ylabel('Magnetic Field weight w_B')
grid on
ylim([0 1])
set(gca,'FontSize',15)

figure
hold on
plot(Table_current.SetAngle,picos)
%plot(Table_current.SetAngle, B_mult)
xlabel('Angle \phi_B [º]')
ylabel('Picos')
grid on
ylim([0 1])
set(gca,'FontSize',15)
%% SuperPower Website vs DB
clc
clear all
close all

load AllSamplesInterp3T.mat

filedirectory = fullfile('..','superpower_data_figure_current_field_temp.txt');
Names = {'Temperature'; 'Bperp'; 'CriticalCurrent'};
Website = readtable(filedirectory);
Website.Properties.VariableNames = Names;

temps_toplot = [20 30 40 50];
angles_toplot = [0];
models = ["Kim Plus"];

max_field = max(Samples.DataRelevant{1}.SetField);

for i=1:numel(models)
    for j=1:numel(angles_toplot)
        figure
        hold on
        %title(num2str(angles_toplot(j))+"º angle : Experimental vs"+" "+models(i))
        set(groot,'DefaultLineLineWidth',1.5)
        set(gca,'FontSize',12)%,'YScale','log'
        for k=1:numel(temps_toplot)
            i_pt = find(and(ismember(Samples.KimParam{1}.SetTemperature,temps_toplot(k)),ismember(Samples.KimParam{1}.SetAngle,angles_toplot(j))));
            if isempty(i_pt)
                continue
            end
            Table_current = Samples.Data{1}(and(ismember(Samples.Data{1}.SetTemperature,temps_toplot(k)),ismember(Samples.Data{1}.SetAngle,angles_toplot(j))),:);
            plot(Table_current.SetField, Table_current.CriticalCurrent * 4,"red");
            Website_buffer = Website(Website.Temperature == temps_toplot(k),:);
            plot(Website_buffer.Bperp,Website_buffer.CriticalCurrent,'x')
        end
        legend('Experimental','Website','Location','northeast')
        xlabel('Magnetic Field [T]')
        ylabel('Critical Current [A]')
        grid on
        hold off
        xlim([0 3.2])
    end
end
%% Compare B0 and beta for all samples
clc
clear all
close all

load AllSamplesInterp3T.mat
temp_toanalyze = 77.5;
angle_toanalyze = 0;

for i = 1:numel(Samples.(1))
    i_pt = find(and(Samples.KimParam{i}.SetTemperature == temp_toanalyze,Samples.KimParam{i}.SetAngle == angle_toanalyze));
    if isempty(i_pt)
        B0(i) = NaN;
        beta(i) = NaN;
        disp('not available')
        continue
    end
    B0(i) = Samples.KimParam{i}.KimPlusB0(i_pt);
    beta(i) = Samples.KimParam{i}.KimPlusbeta(i_pt);
    labels(i) = Samples.Label(i);
end

figure
bar(B0)
grid on
xticks(linspace(1,numel(labels),numel(labels)));
xticklabels(labels);
ylabel('B0')

figure
bar(beta)
grid on
xticks(linspace(1,numel(labels),numel(labels)));
xticklabels(labels);
ylabel('beta')
%% Functions to analyse the models
function plot_relative_error_kim(Table_relevant,ParameterTable,label,angles_toplot,models_toplot)
%% Compute and filter relative error
create_figure = [0,1];
filt_limit = 10;

relevant_temps = unique(ParameterTable.SetTemperature);
relevant_fields = unique(Table_relevant.SetField);
[X,Y] = meshgrid(relevant_temps,relevant_fields);

for m = 1:numel(models_toplot)
    for k = 1:numel(angles_toplot)
        for i = 1:numel(relevant_temps)
            for j = 1:numel(relevant_fields)
                i_tr = find(and(Table_relevant.SetTemperature == relevant_temps(i),...
                    and(Table_relevant.SetField == relevant_fields(j),Table_relevant.SetAngle == angles_toplot(k))));
                i_pt = find(and(ParameterTable.SetTemperature == relevant_temps(i),ParameterTable.SetAngle == angles_toplot(k)));
                if isempty(i_tr)
                    Z(j,i,k,m) = NaN;
                    continue
                end
                switch models_toplot(m)
                    case "Kim"
                        Z(j,i,k,m) = 100*abs((Table_relevant.CriticalCurrent(i_tr) - feval(ParameterTable.curveKim{i_pt},Table_relevant.SetField(i_tr))) / Table_relevant.CriticalCurrent(i_tr));
                    case "Kim Plus"
                        Z(j,i,k,m) = 100*abs((Table_relevant.CriticalCurrent(i_tr) - feval(ParameterTable.curveKimPlus{i_pt},Table_relevant.SetField(i_tr))) / Table_relevant.CriticalCurrent(i_tr));
                end
            end
        end
    end
end
Z_filt = Z;
Z_filt(Z_filt > filt_limit) = NaN;
for m = 1:numel(models_toplot)
    for k = 1:numel(angles_toplot)
        if create_figure(1) == 1
            figure
            surf(X,Y,Z_filt(:,:,k,m))
            %title(label+" "+models_toplot(m)+" <10% relative error "+num2str(angles_toplot(k))+"º B")
            zlabel('Relative Error [%]')
            ylabel('Magnetic Field [T]')
            xlabel('Temperature [K]')
            grid on
            c = colorbar;   c.Label.String = 'Relative Error [%]';  c.Label.FontSize = 12;
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

for m=1:numel(models_toplot)
    for k=1:numel(angles_toplot)
        Z_interp = interp2(X,Y,Z(:,:,k,m),X_interp,Y_interp);
        Z_interp_filt = Z_interp;
        Z_interp_filt(Z_interp_filt >= filt_limit) = NaN;

        if create_figure(2) == 1
            figure
            surf(X_interp,Y_interp,Z_interp_filt)
            %title(label+" "+models_toplot(m)+" <10% relative error "+num2str(angles_toplot(k))+"º interpolated")
            zlabel('Relative Error [%]')
            ylabel('Magnetic Field [T]')
            xlabel('Temperature [K]')
            xlim([min(ParameterTable.SetTemperature) max(ParameterTable.SetTemperature)])
            ylim([0 3])
            grid on
            c = colorbar;   c.Label.String = 'Relative Error [%]';  c.Label.FontSize = 12;
            set(gca,'FontSize',12)
            view(gca,[0.144782608695627 90]);
            savestring = label+" Rel Error "+models_toplot(m)+" "+num2str(angles_toplot(k))+".png";
            exportgraphics(gcf,"./Plots EM/"+savestring)
        end
    end
end
end

function plot_param_kim(Table_relevant,ParameterTable,label,models_toplot)
relevant_temps = unique(Table_relevant.SetTemperature);
relevant_angles = unique(Table_relevant.SetAngle);
[X,Y] = meshgrid(relevant_temps,relevant_angles);

for m=1:numel(models_toplot)
    clearvars B0 beta
    for i=1:numel(relevant_temps)
        for j=1:numel(relevant_angles)
            i_pt = find(and(ParameterTable.SetTemperature == relevant_temps(i),ParameterTable.SetAngle == relevant_angles(j)));
            if isempty(i_pt)
                B0(j,i) = NaN;
                continue
            end
            switch models_toplot(m)
                case "Kim"
                    B0(j,i) = ParameterTable.KimB0(i_pt);
                case "Kim Plus"
                    B0(j,i) = ParameterTable.KimPlusB0(i_pt);
                    beta(j,i) = ParameterTable.KimPlusbeta(i_pt);
            end
        end
    end
    
    if numel(relevant_angles) == 1
        figure
        plot(relevant_temps,B0(1,:))
        %title(label+" "+models_toplot(m)+" B0")
        ylabel('B_0 [T]')
        xlabel('Temperature [K]')
        grid on
        set(gca,'FontSize',12)
        savestring = label+" B0 "+models_toplot(m)+".png";
        exportgraphics(gcf,"./Plots EM/"+savestring)
        if models_toplot(m) == "Kim Plus"
            figure
            plot(relevant_temps,beta(1,:))
            %title(label+" "+models_toplot(m)+" beta")
            ylabel('\beta')
            xlabel('Temperature [K]')
            grid on
            set(gca,'FontSize',12)
            savestring = label+" beta "+models_toplot(m)+".png";
            exportgraphics(gcf,"./Plots EM/"+savestring)
        end
    else
        figure
        surf(X,Y,B0)
        %title(label+" "+models_toplot(m)+" B0")
        zlabel('B_0 [T]')
        ylabel('Angle [º]')
        xlabel('Temperature [K]')
        grid on
        c = colorbar;   c.Label.String = 'B_0 [T]'; c.Label.FontSize = 12;
        set(gca,'FontSize',12)
        view(gca,[89.8105716251899 90]);
        savestring = label+" B0 "+models_toplot(m)+".png";
        exportgraphics(gcf,"./Plots EM/"+savestring)
        if models_toplot(m) == "Kim Plus"
            figure
            surf(X,Y,beta)
            %title(label+" "+models_toplot(m)+" beta")
            zlabel('\beta')
            ylabel('Angle [º]')
            xlabel('Temperature [K]')
            grid on
            c = colorbar;   c.Label.String = '\beta';   c.Label.FontSize = 12;
            set(gca,'FontSize',12)
            view(gca,[89.8105716251899 90]);
            savestring = label+" beta "+models_toplot(m)+".png";
            exportgraphics(gcf,"./Plots EM/"+savestring)
        end
    end
end
end

function plot_kim_curve(Table,ParameterTable,label,width,temps_toplot,angles_toplot,models_toplot,scale)
max_field = max(Table.SetField);
B = 0:0.01:max_field;

for i=1:numel(models_toplot)
    for j=1:numel(angles_toplot)
        figure
        hold on
        %title(label+"'s "+width+"cm width: Ic vs B at "+num2str(angles_toplot(j))+"º angle")
        set(groot,'DefaultLineLineWidth',1.5)
        if scale == "log"
            ylim([10E0,10E3])
            set(gca,'YScale','log')
        end
        for k=1:numel(temps_toplot)
            index = find(and(ismember(ParameterTable.SetTemperature,temps_toplot(k)),ismember(ParameterTable.SetAngle,angles_toplot(j))));
            if isempty(index)
                continue
            end
            table_buffer = Table(ismember(Table.SetTemperature,temps_toplot(k)),:);
            table_buffer = table_buffer(ismember(table_buffer.SetAngle,angles_toplot(j)),:);
            plot(table_buffer.SetField, table_buffer.CriticalCurrent/width,"red")
            switch models_toplot(i)
                case "Kim"
                    plot(B,feval(ParameterTable.curveKim{index},B)/width,"blue")
                    height = feval(ParameterTable.curveKim{index},max_field/2)/width+0.4*feval(ParameterTable.curveKim{index},max_field/2)/width;
                    text(max_field/2,height,"T="+num2str(temps_toplot(k))+" K B_0="+round(ParameterTable.KimB0(index),2),...
                        'HorizontalAlignment','center')
                case "Kim Plus"
                    plot(B,feval(ParameterTable.curveKimPlus{index},B)/width,"blue")
                    height = feval(ParameterTable.curveKimPlus{index},max_field/2)/width+0.4*feval(ParameterTable.curveKimPlus{index},max_field/2)/width;
                    text(max_field/2,height,"T="+num2str(temps_toplot(k))+" K, B_0="+round(ParameterTable.KimPlusB0(index),2)+", \beta="+...
                        round(ParameterTable.KimPlusbeta(index),2),'HorizontalAlignment','center')
            end
        end
        legend('Database',models_toplot(i),'Location','best')
        xlabel('Magnetic Field [T]')
        ylabel('Critical Current Sheet Density [A/cm]')
        grid on
        set(gca,'FontSize',12)
        savestring = label+" curve "+models_toplot(i)+" "+num2str(angles_toplot(j))+".png";
        exportgraphics(gcf,"./Plots EM/"+savestring)
    end
end
end

function plot_relative_error_custom(Table_relevant,ParameterTable,label,model_toplot)
if isempty(ParameterTable)
    disp("No custom relative error available for "+label)
    return
end
%% Compute and filter relative error
angles_toplot = [0 45 90];
create_figure = [0,1];
filt_limit = 10;

relevant_temps = unique(ParameterTable.SetTemperature);
if numel(relevant_temps) == 1
    create_figure = [1,0];
end
relevant_fields = unique(Table_relevant.SetField);
[X,Y] = meshgrid(relevant_temps,relevant_fields);

for k = 1:numel(angles_toplot)
    for i = 1:numel(relevant_temps)
        for j = 1:numel(relevant_fields)
            i_tr = find(and(Table_relevant.SetTemperature == relevant_temps(i),...
                and(Table_relevant.SetField == relevant_fields(j),Table_relevant.SetAngle == angles_toplot(k))));
            i_pt = find(ParameterTable.SetTemperature == relevant_temps(i));
            if isempty(i_tr)
                Z(j,i,k) = NaN;
                continue
            end
            Z(j,i,k) = 100 * abs((Table_relevant.CriticalCurrent(i_tr) - feval(ParameterTable.curve{i_pt},...
                Table_relevant.SetField(i_tr),Table_relevant.SetAngle(i_tr)*pi/180)) / Table_relevant.CriticalCurrent(i_tr));
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
            %title(label+" Custom <10% relative error "+num2str(angles_toplot(k))+"º "+relevant_temps(1)+"K")
            ylabel('Relative Error [%]')
            xlabel('Magnetic Field [T]')
            grid on
        else
            figure
            surf(X,Y,Z_filt(:,:,k))
            %title(label+" Custom <10% relative error "+num2str(angles_toplot(k))+"º")
            zlabel('Relative Error [%]')
            ylabel('Magnetic Field [T]')
            xlabel('Temperature [K]')
            grid on
            c = colorbar;   c.Label.String = 'Relative Error [%]';  c.Label.FontSize = 12;
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
        %title(label+" Custom <10% relative error "+num2str(angles_toplot(k))+"º interpolated")
        zlabel('Relative Error [%]')
        ylabel('Magnetic Field [T]')
        xlabel('Temperature [K]')
        xlim([min(ParameterTable.SetTemperature) max(ParameterTable.SetTemperature)])
        ylim([0 3])
        grid on
        c = colorbar;   c.Label.String = 'Relative Error [%]';  c.Label.FontSize = 12;
        set(gca,'FontSize',12)
        view(gca,[0.144782608695627 90]);
        savestring = label+" Rel Error "+model_toplot+" "+num2str(angles_toplot(k))+".png";
        exportgraphics(gcf,"./Plots EM/"+savestring)
    end
end
end

function plot_param_custom(ParameterTable,label)
if isempty(ParameterTable)
    disp("No custom params available for "+label)
    return
end

relevant_temps = unique(ParameterTable.SetTemperature);
figure_size = [100 300 800 300];
axes_size = [0.1 0.22 0.88 0.72];
names = ["B0","beta","SFa","SFw","SFp","IRa","IRw","IRp"];
names_plot = ["B_0 [T]","\beta","SF_a","SF_w","SF_p [º]","IR_a","IR_w","IR_p [º]"];

for i = 3:10
    figure('Position',figure_size)
    axes('Position',axes_size)
    set(groot,'DefaultLineLineWidth',1.5)
    if names(i-2) == "IRp" || names(i-2) == "SFp"
        plot(relevant_temps,ParameterTable.(i) * 180 / pi,'red')
    else
        plot(relevant_temps,ParameterTable.(i),'red')
    end
    %title(label+" Custom "+names(i-2))
    ylabel(names_plot(i-2))
    xlabel('Temperature [K]')
    %ylim([min(ParameterTable.(i)) max(ParameterTable.(i))])
    grid on
    set(gca,'FontSize',12)
    savestring = label+" "+names(i-2)+" Custom"+".png";
    exportgraphics(gcf,"./Plots EM/"+savestring)
end
end

function plot_param_advanced_custom(ParameterTable,label)
if isempty(ParameterTable)
    disp("No advanced custom params available for "+label)
    return
end

relevant_temps = unique(ParameterTable.SetTemperature);
figure_size = [100 300 800 300];
axes_size = [0.1 0.22 0.88 0.72];
names = ["B0","beta","SFa","SFwi","SFwf","SFws","SFpi","SFpf","SFps","IRaf","IRas","IRw","IRp"];
names_plot = ["B_0 [T]","\beta","SF_a","SF_{wi}","SF_{wf}","SF_{ws}","SF_{pi} [º]","SF_{pf} [º]","SF_{ps}","IR_{af}","IR_{as}","IR_w","IR_p [º]"];

for i = 3:15
    figure('Position',figure_size)
    axes('Position',axes_size)
    set(groot,'DefaultLineLineWidth',1.5)
    if names(i-2) == "IRp" || names(i-2) == "SFpi" || names(i-2) == "SFpf"
        plot(relevant_temps,ParameterTable.(i) * 180 / pi,'red')
    else
        plot(relevant_temps,ParameterTable.(i),'red')
    end
    %title(label+" Custom "+names(i-2))
    ylabel(names_plot(i-2))
    xlabel('Temperature [K]')
    %ylim([min(ParameterTable.(i)) max(ParameterTable.(i))])
    grid on
    set(gca,'FontSize',12)
    savestring = label+" "+names(i-2)+" Advanced Custom"+".png";
    exportgraphics(gcf,"./Plots EM/"+savestring)
end
end

function plot_advanced_custom_bdepend(ParameterTable,label,temps_toplot)
if isempty(ParameterTable)
    disp("No B dependent advanced custom params available for "+label)
    return
end

figure_size = [100 300 800 300];
axes_size = [0.075 0.2 0.9 0.72];

B = 0:0.01:3;
names = ["IRa","SFw","SFp"];
names_plot = ["IR_a","SF_w","SF_p [º]"];
for i = 1:numel(names)
    figure('Position',figure_size)
    axes('Position',axes_size)
    set(groot,'DefaultLineLineWidth',1.5)
    hold on
    for j = 1:numel(temps_toplot)
        i_pt = find(ParameterTable.SetTemperature == temps_toplot(j));
        if isempty(i_pt)
            continue
        end
        clear data
        if names(i) == "IRa"
            data = tanh(ParameterTable.IRas(i_pt)*B)*ParameterTable.IRaf(i_pt);
        elseif names(i) == "SFw"
            data = ParameterTable.SFwi(i_pt)+ParameterTable.SFwf(i_pt)*tanh(ParameterTable.SFws(i_pt)*B);
        elseif names(i) == "SFp"
            data = -ParameterTable.SFpi(i_pt)*exp(-ParameterTable.SFps(i_pt)*B)+ParameterTable.SFpf(i_pt);
            data = data * 180 / pi;
        end
        plot(B,data)
    end
    %title(label+" Advanced Custom "+names(i))
    ylabel(names_plot(i))
    xlabel('Magnetic Field [T]')
    legend(num2str(temps_toplot')+" K",'Location','northoutside','NumColumns',numel(temps_toplot))
    grid on
    set(gca,'FontSize',12)
    savestring = label+" "+names(i)+" Advanced Custom"+".png";
    exportgraphics(gcf,"./Plots Abstract/"+savestring)
end
end

function plot_custom_curve(Table,ParameterTable,label,width)
temps_toplot = [20 45 65 77.5 85];
angles_toplot = [0 45 90];

max_field = max(Table.SetField);

for i = 1:numel(angles_toplot)
    figure
    hold on
    %title(label+num2str(angles_toplot(i))+"º angle : Experimental vs Custom")
    set(groot,'DefaultLineLineWidth',1.5)
    set(gca,'FontSize',12,'YScale','log')
    for j = 1:numel(temps_toplot)
        i_pt = find(ParameterTable.SetTemperature == temps_toplot(j));
        if isempty(i_pt)
            continue
        end

        Table_current = Table(and(Table.SetTemperature == temps_toplot(j),Table.SetAngle == angles_toplot(i)),:);
        relevant_fields = unique(Table_current.SetField);
        plot(relevant_fields, Table_current.CriticalCurrent/width,"red")
        plot(relevant_fields,feval(ParameterTable.curve{i_pt},relevant_fields,ones(1,numel(relevant_fields))*angles_toplot(i)*pi/180)/width,"blue")
        height = feval(ParameterTable.curve{i_pt},max_field/2,angles_toplot(i))/width+0.4*feval(ParameterTable.curve{i_pt},max_field/2,angles_toplot(i))/width;

        text(max_field/2,height,num2str(temps_toplot(j))+"K")
    end
    legend('Experimental','Fitted','Location','best')
    xlabel('Magnetic Field [T]')
    ylabel('Critical Current Sheet Density [A/cm]')
    grid on
    hold off
end
end

function plot_custom_angle(Table_Relevant,ParameterTable,label,width,temps_toplot,fields_toplot,model_toplot)
for i=1:numel(temps_toplot)
    for j=1:numel(fields_toplot)
        Table_buffer = Table_Relevant(and(Table_Relevant.SetTemperature == temps_toplot(i),Table_Relevant.SetField == fields_toplot(j)),:);
        relevant_angles = unique(Table_buffer.SetAngle);
        i_pt = find(ParameterTable.SetTemperature == temps_toplot(i));
        if isempty(i_pt | any(numel(relevant_angles) == [1 2 3]))
            disp("Cannot see angle dependence.")
            continue
        end

        figure
        hold on
        %title(label+" angle dependency for "+num2str(temps_toplot(i))+"K and "+fields_toplot(j)+"T")
        set(groot,'DefaultLineLineWidth',1.5)
        set(gca,'FontSize',12)
        plot(relevant_angles,Table_buffer.CriticalCurrent/width,"red")
        plot(relevant_angles,feval(ParameterTable.curve{i_pt},ones(numel(relevant_angles),1)*fields_toplot(j),relevant_angles*pi/180)/width,"blue")
        legend('Database',model_toplot,'Location','best')
        xlabel('Angle [º]')
        ylabel('Critical Current Sheet Density [A/cm]')
        grid on
        set(gca,'FontSize',12)
        savestring = label+" angle curve "+model_toplot+" "+num2str(temps_toplot(i))+"K "+num2str(fields_toplot(j))+"T.png";
        exportgraphics(gcf,"./Plots EM/"+savestring)
    end
end
end

function plot_Ic_vs_Bphy(Table_relevant,temp_toplot)
relevant_fields = unique(Table_relevant.SetField);
relevant_angles = unique(Table_relevant.SetAngle);
[X,Y] = meshgrid(relevant_fields,relevant_angles);

for i = 1:numel(relevant_fields)
    for j = 1:numel(relevant_angles)
        i_tr = find(and(ismember(Table_relevant.SetTemperature,temp_toplot),and(ismember(Table_relevant.SetField,relevant_fields(i)),...
            ismember(Table_relevant.SetAngle,relevant_angles(j)))));
        if isempty(i_tr)
            Z(j,i) = NaN;
            continue
        end
        Z(j,i) = Table_relevant.CriticalCurrent(i_tr);
    end
end
surf(X,Y,Z)
%title("Ic vs (B,angle)")
zlabel('Critical Current [A]')
ylabel('Angle [º]')
xlabel('Field [T]')
grid on
end

function plot_Ic_vs_TB(Table_relevant,angle_toplot)
relevant_temps = unique(Table_relevant.SetTemperature);
relevant_fields = unique(Table_relevant.SetField);
[X,Y] = meshgrid(relevant_temps,relevant_fields);

for i=1:numel(relevant_temps)
    for j=1:numel(relevant_fields)
        index_tr = find(and(ismember(Table_relevant.SetTemperature,relevant_temps(i)),and(ismember(Table_relevant.SetField,relevant_fields(j)),ismember(Table_relevant.SetAngle,angle_toplot))));
        if isempty(index_tr)
            Z(j,i) = NaN;
            continue
        end
        Z(j,i) = Table_relevant.CriticalCurrent(index_tr);
    end
end
figure
surf(X,Y,Z)
%title("Critical Current")
zlabel('Critical Current [A]')
ylabel('Magnetic Field [T]')
xlabel('Temperature [K]')
grid on
view(gca,[142.515652173913 38.55]);
end
