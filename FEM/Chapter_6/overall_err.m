function result = overall_err(err_data,err_types)
result = err_data(1,1:2);
for ii = 1:1:length(err_types)
    switch err_types(ii)
        case "max"% 节点max误差
            temp = max(err_data(:,2 + ii));
        case "L^inf"% L^infty误差
            temp = max(err_data(:,2 + ii));
        case "L^2"% L^2误差
            temp = sqrt(err_data(1,2 + ii).^2 + err_data(2,2 + ii).^2);
        case "H^1(semi)"% H^1(semi)误差
            temp = sqrt(err_data(1,2 + ii).^2 + err_data(2,2 + ii).^2);
        case "H^1"% H^1误差
            temp = sqrt(err_data(1,2 + ii).^2 + err_data(2,2 + ii).^2);
        case "all"
            temp = overall_err(err_data,["max","L^inf","L^2","H^1","H^1(semi)"]);
            temp = temp(3:end);% 去除h
        otherwise
            error(['Unsupported error type: ' err_types(ii)]);
    end
    result = [result temp];
end
end