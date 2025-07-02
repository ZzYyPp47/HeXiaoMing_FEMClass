%% Covergence Test for 2D example FEM 
% 计算 ut-∇·(c∇u) = -3exp(x+y+t)
% c = 2;
% u = exp(x+y) t = 0 & in \Omega;u = exp(y+t) x = 0;
% u = exp(2+y+t) x = 2;u = exp(x+t) y = 0;u = exp(x+1+t) y = 1;
% real solution: exp(x+y+t)
% \omega x T = [0,2] x [0,1] x [0,1]
%% 初始化
clc;
clear;
N = 6;
result = zeros(N,7);
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
    [~,result(ii,:)] = FEM_2D_Linear_Dirichlet([1/2^(ii),1/2^(ii),1/2^(ii),1/2],"end",false,true,["max","L^inf","L^2","H^1(semi)"]);
    time = toc(t_s);
    fprintf('Completed task num: %d, time passed: %.8f s, Total: %d\n',ii,time,N);
end
result_order = log(result(1:end - 1,4:end)./result(2:end,4:end)) ./ log(result(1:end - 1,1)./result(2:end,1));
result_order = table(result_order(:,1),result_order(:,2),result_order(:,3),result_order(:,4),'VariableNames',{'Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
draw(result,'(Quadratic 2D FEM)');
result = table(sym(result(:,1)),sym(result(:,2)),sym(result(:,3)),result(:,4),result(:,5),result(:,6),result(:,7),...
               'VariableNames',{'hx','hy','ht','Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
disp(result);
fprintf("Convergence Order: \n");
disp(result_order);
%% 绘图函数
function draw(data,txt)
f_1 = @(x)x.^2;
f_2 = @(x)x.^3;
x = data(:,1);
% 绘制 O(*) 的参考线
figure;
loglog(x,f_1(x),'-o','LineWidth',2,'Color','r');
hold on;
loglog(x,f_2(x),'-o','LineWidth',2,'Color','g');
% 绘制实际的误差图
loglog(x,data(:,4), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,5), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,6), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,7), '--','LineWidth',2,'Marker',"pentagram");
grid on;
grid minor;
xlabel('log10(h)');
ylabel('log10(error)');
title(string(txt) + ' Error vs Mesh Width');
legend('$O(h^2)$ reference','$O(h^3)$ reference',...
        'max abs error for all finite nodes','$L^{\infty}$ error','$L^{2}$ error','$H^1(semi)$ error',...
        'Interpreter','latex','Location','best');
end