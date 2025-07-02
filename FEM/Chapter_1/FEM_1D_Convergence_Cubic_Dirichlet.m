%% Local Assembly Covergence Test for 1D example FEM
% 计算 -(exp(x)u')' = -exp(x)(cos(x)-2sin(x)-xcos(x)-xsin(x))
% 0 <= x <= 1
% analytic solution : u = xcos(x)
%% 初始化
clc;
clear;
N = 7;
result = zeros(N,5);
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
%% 串行计算
t_s = tic;
for ii = 1:N
    [~,result(ii,:)] = FEM_1D_Cubic_Dirichlet(1/2^(ii),false,true,["max","L^inf","L^2","H^1(semi)"]);
    time = toc(t_s);
    fprintf('Completed task num: %d, time passed: %.8f s, Total: %d\n',ii,time,N);
end
result_order = log(result(1:end - 1,2:end)./result(2:end,2:end)) ./ log(result(1:end - 1,1)./result(2:end,1));
result_order = table(result_order(:,1),result_order(:,2),result_order(:,3),result_order(:,4),'VariableNames',{'Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
draw(result,'(Cubic 1D FEM)');
result = table(sym(result(:,1)),result(:,2),result(:,3),result(:,4),result(:,5),'VariableNames',{'h','Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
disp(result);
fprintf("Convergence Order: \n");
disp(result_order);
%% 绘图函数
function draw(data,txt)
f_1 = @(x)x.^3;
f_2 = @(x)x.^4;
x = data(:,1);
% 绘制 O(*) 的参考线
figure;
loglog(x,f_1(x),'-o','LineWidth',2,'Color','r');
hold on;
loglog(x,f_2(x),'-o','LineWidth',2,'Color','g');
% 绘制实际的误差图
loglog(x,data(:,2), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,3), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,4), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,5), '--','LineWidth',2,'Marker',"pentagram");
grid on;
grid minor;
xlabel('log10(h)');
ylabel('log10(error)');
title(string(txt) + ' Error vs Mesh Width');
legend('$O(h^3)$ reference','$O(h^4)$ reference',...
        'max abs error for all finite nodes','$L^{\infty}$ error','$L^{2}$ error','$H^1(semi)$ error',...
        'Interpreter','latex','Location','best');
end