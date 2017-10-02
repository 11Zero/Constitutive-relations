close all;
clear;
clc;
grade = [20 25 30 35 40 45 50 55 60 65 70 75 80];
%grade代表C20到C80的混凝土等级标号
%单位N/mm2,即MPa
f_ck = [13.4 16.7 20.1 23.4 26.8 29.6 32.4 35.5 38.5 41.5 44.5 47.4 50.2];%抗压强度标准值
f_tk = [1.54 1.78 2.01 2.20 2.39 2.51 2.64 2.74 2.85 2.93 2.99 3.05 3.11];%抗拉强度标准值
f_c = [9.6 11.9 14.3 16.7 19.1 21.1 23.1 25.3 27.5 29.7 31.8 33.8 35.9];%抗压强度设计值
f_t = [1.10 1.27 1.43 1.57 1.71 1.80 1.89 1.96 2.04 2.09 2.14 2.18 2.22];%抗拉强度设计值
Ec = [2.55 2.80 3.00 3.15 3.25 3.35 3.45 3.55 3.60 3.65 3.70 3.75 3.80].*10000;%弹性模量N/mm2
delta_c = [23.3 20.6 18.9 17.2 16.4 15.6 15.6 14.9 14.1 14 14 14 14]./100;%强度变异系数
%deta_c中的后四个14为自己的取值，按14取，如有实测数据再进行修改
f_cm = f_ck./(1-1.645.*delta_c);%抗压强度平均值
f_tm = f_tk./(1-1.645.*delta_c);%抗拉强度平均值
f_tr = f_t;%规范规定f_tr应根据实际需要取f_t f_tk f_tm之一
f_cr = f_c;%规范规定f_tr应根据实际需要取f_c f_ck f_cm之一
parameter_tr = [1.0 1.5 2.0 2.5 3.0 3.5 4.0;[65 81 95 107 118 128 137]/1000000;0.31 0.70 1.25 1.95 2.81 3.82 5];%参数插值表
rou_t = zeros(size(grade));
epsilon_tr = zeros(size(grade));
alpha_t = zeros(size(grade));
epsilon_length = 5;
epsilon_step = 20;
x = [0:epsilon_length*epsilon_step];
x = x./epsilon_step;
dt = zeros(numel(grade),epsilon_length*epsilon_step+1);
for i=1:numel(f_tr)
    temp_val = floor(f_tr(i)/0.5)*0.5;
    if(temp_val<1)
        temp_val = 1;
    elseif(temp_val >= 4)
        temp_val = 3.5;
    end
    [row fpos] = find(parameter_tr(1,:)==temp_val);
    rate = (f_tr(i)-parameter_tr(1,fpos))/(parameter_tr(1,fpos+1)-parameter_tr(1,fpos));
    epsilon_tr(i) = parameter_tr(2,fpos)+rate*(parameter_tr(2,fpos+1)-parameter_tr(2,fpos));
    alpha_t(i) = parameter_tr(3,fpos)+rate*(parameter_tr(3,fpos+1)-parameter_tr(3,fpos));
    rou_t(i) = f_tr(i)/Ec(i)/epsilon_tr(i);
    
    for j=1:numel(x)
        if(x(j)>1)
            dt(i,j) = 1-rou_t(i)/(alpha_t(i)*power(x(j)-1,1.7)+x(j));
        else
            dt(i,j) = 1-rou_t(i)*(1.2-0.2*power(x(j),5));
        end
    end
end
figure;
hold on;
title('受拉');
for i=1:13
    color = rand(1,3);
    selected = i;
    epsilon = x.*epsilon_tr(selected);
    rou = (1-dt(selected,:)).*Ec(selected).*x.*epsilon_tr(selected);
    xlswrite('concrete_t.xls',[-epsilon' -rou'],strcat('C',num2str(grade(i))));
    plot(x.*epsilon_tr(selected),(1-dt(selected,:)).*Ec(selected).*x.*epsilon_tr(selected),'Color',color);
end
hold off;

parameter_cr = [[20:5:80];...
    [1470 1560 1640 1720 1790 1850 1920 1980 2030 2080 2130 2190 2240]/1000000;...
    0.74 1.06 1.36 1.65 1.94 2.21 2.48 2.74 3.00 3.25 3.50 3.75 3.99;...
    3.0 2.6 2.3 2.1 2.0 1.9 1.9 1.8 1.8 1.7 1.7 1.7 1.6];%参数插值表
rou_c = zeros(size(grade));
epsilon_cr = zeros(size(grade));
alpha_c = zeros(size(grade));
epsilon_cu_epsilon_cr = zeros(size(grade));
epsilon_length = 2;
epsilon_step = 20;
x = [0:epsilon_length*epsilon_step];
x = x./epsilon_step;
dc = zeros(numel(grade),epsilon_length*epsilon_step+1);
for i=1:numel(f_cr)
    temp_val = floor(f_cr(i)/5)*5;
    if(temp_val<=20)
        temp_val = 20;
    elseif(temp_val >= 80)
        temp_val = 80;
    end
    [row fpos] = find(parameter_cr(1,:)==temp_val);
    rate = (f_cr(i)-parameter_cr(1,fpos))/(parameter_cr(1,fpos+1)-parameter_cr(1,fpos));
    epsilon_cr(i) = parameter_cr(2,fpos)+rate*(parameter_cr(2,fpos+1)-parameter_cr(2,fpos));
    alpha_c(i) = parameter_cr(3,fpos)+rate*(parameter_cr(3,fpos+1)-parameter_cr(3,fpos));
    epsilon_cu_epsilon_cr(i) = parameter_cr(4,fpos)+rate*(parameter_cr(4,fpos+1)-parameter_cr(4,fpos));
    rou_c(i) = f_cr(i)/Ec(i)/epsilon_cr(i);
    n = Ec(i)*epsilon_cr(i)/(Ec(i)*epsilon_cr(i)-f_cr(i));
    for j=1:numel(x)
        if(x(j)>1)
            dc(i,j) = 1-rou_c(i)/(alpha_c(i)*power(x(j)-1,2)+x(j));
        else
            dc(i,j) = 1-rou_c(i)*n/(n-1+power(x(j),n));
        end
    end
end
figure;
hold on;
title('受压');
for i=1:13
    color = 0.8*rand(1,3);
    selected = i;
    epsilon = x.*epsilon_cr(selected);
    rou = (1-dc(selected,:)).*Ec(selected).*x.*epsilon_cr(selected);
    xlswrite('concrete_c.xls',[epsilon' rou'],strcat('C',num2str(grade(i))));
    plot(epsilon,rou,'Color',color);
    plot(epsilon,rou,'x');
    plot(epsilon_cr(selected)*epsilon_cu_epsilon_cr(selected).*ones(1,2),[0 f_cr(selected)],'Color',color);
    plot([0 epsilon_cr(selected)*epsilon_cu_epsilon_cr(selected)],0.5*f_cr(selected).*ones(1,2),'Color',color);
end
hold off;






