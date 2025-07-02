% 可视化函数
function show_eval(option,real_u,u,X,Y,Time,txt)
level = 100;
if option == "all"
    abs_data = zeros([size(X),size(Time,2)]);
    fem_data = zeros([size(X),size(Time,2)]);
    parfor ii = 1:1:size(Time,2)
        fem_data(:,:,ii) = reshape(u(:,ii),size(X));
        abs_data(:,:,ii) = abs(real_u(X,Y,Time(ii),0,0) - fem_data(:,:,ii));
    end
    for ii = 1:1:size(u,2)
        % 比较
        tiledlayout(2,1);
        cla(tiledlayout);% 清除之前的绘图

        nexttile;
        contourf(X,Y,fem_data(:,:,ii),level,'EdgeAlpha',0)
        title("FEM at time = " + num2str(Time(ii)));
        xlabel('x');
        ylabel('y');
        colorbar;
        if ii > 1
            min_old = min(fem_data(:,:,ii - 1),[],"all");
            min_now = min(fem_data(:,:,ii),[],"all");
            max_old = max(fem_data(:,:,ii - 1),[],"all");
            max_now = max(fem_data(:,:,ii),[],"all");
            clim([min(min_old,min_now) max(max_old,max_now)]);
        end
        axis equal;

        % 误差分析
        nexttile;
        contourf(X,Y,abs_data(:,:,ii),level,'EdgeAlpha',0)
        title("Absolute Error (" + txt + ") " + "at time = " + num2str(Time(ii)));
        xlabel('x');
        ylabel('y');
        colorbar;
        if ii > 1
            min_old = min(abs_data(:,:,ii - 1),[],"all");
            min_now = min(abs_data(:,:,ii),[],"all");
            max_old = max(abs_data(:,:,ii - 1),[],"all");
            max_now = max(abs_data(:,:,ii),[],"all");
            clim([min(min_old,min_now) max(max_old,max_now)]);
        end
        axis equal;

        drawnow;
    end
else
    % 比较
    tiledlayout(2,1);
    real = real_u(X,Y,Time(end));
    fem = u(:,:);

    nexttile;
    contourf(X,Y,fem,level,'EdgeAlpha',0)
    title("FEM at time = " + num2str(Time(end)));
    xlabel('x');
    ylabel('y');
    colorbar;
    axis equal;

    % 误差分析
    abs_error = abs(real - fem);
    nexttile;
    contourf(X,Y,abs_error,level,'EdgeAlpha',0)
    title("Absolute Error (" + txt + ") " + "at time = " + num2str(Time(end)));
    xlabel('x');
    ylabel('y');
    colorbar;
    axis equal;
end
end