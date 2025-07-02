%% Covergence Test for 2D example FEM
% 计算 -∇· σ(u1 u2) = (f1 f2) on \Omega
% u1 = 0; u2 = 0 on \Pratial \Omega;
% f1 = -(lamda + 3mu) * (-pi^2 * sin(pi * x) * sin(pi * y)) - (lamda + mu) * ((2x - 1) * (2y - 1));
% f2 = -(lamda + 2mu) * (2 * x * (x - 1)) - (lamda + mu) * (pi^2 * cos(pi* x) * cos(pi * y)) - mu * (2 * y * (y - 1));
% lamda = 1; mu = 2;
% real solution: [u1 u2] = [sin(pi * x) * sin(pi * y),x * (x - 1) * y * (y - 1)]
% \Omega = [0,1]^2
%% 初始化
clc;
clear;
N = 8;
result = zeros(N,6);
poolobj = gcp('nocreate');
if isempty(poolobj)
    threads = 32; % 设定机器线程数
    fprintf('Initialize parallel pool first, with threads %d.\n',threads)
    parpool(threads);
    fprintf('Now starting...\n');
else
    fprintf('Already initialized parallel pool, with threads %d.\n',poolobj.NumWorkers);
    fprintf('Now starting...\n');
end
%% 计算
t_s = tic;
for ii = 1:N
    [~,result(ii,:)] = FEM_2D_Linear_Dirichlet([1/2^(ii),1/2^(ii)],false,true,["max","L^inf","L^2","H^1(semi)"]);
    time = toc(t_s);
    fprintf('Completed task num: %d, time passed: %.8f s, Total: %d\n',ii,time,N);
end
result_order = log(result(1:end - 1,3:end)./result(2:end,3:end)) ./ log(result(1:end - 1,1)./result(2:end,1));
result_order = table(result_order(:,1),result_order(:,2),result_order(:,3),result_order(:,4),'VariableNames',{'Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
draw(result,'(Linear 2D FEM)');
result = table(sym(result(:,1)),sym(result(:,2)),result(:,3),result(:,4),result(:,5),result(:,6),...
               'VariableNames',{'hx','hy','Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
disp(result);
fprintf("Convergence Order: \n");
disp(result_order);
%% 绘图函数
function draw(data,txt)
f_1 = @(x)x;
f_2 = @(x)x.^2;
x = data(:,1);
% 绘制 O(*) 的参考线
figure;
loglog(x,f_1(x),'-o','LineWidth',2,'Color','r');
hold on;
loglog(x,f_2(x),'-o','LineWidth',2,'Color','g');
% 绘制实际的误差图
loglog(x,data(:,3), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,4), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,5), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,6), '--','LineWidth',2,'Marker',"pentagram");
grid on;
grid minor;
xlabel('log10(h)');
ylabel('log10(error)');
title(string(txt) + ' Error vs Mesh Width');
legend('$O(h)$ reference','$O(h^2)$ reference',...
        'max abs error for all finite nodes','$L^{\infty}$ error','$L^{2}$ error','$H^1(semi)$ error',...
        'Interpreter','latex','Location','best');
end