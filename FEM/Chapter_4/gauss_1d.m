% func need to be a function handle which can feed in x.
function result = gauss_1d(func,vertices)
[xmax,xmin] = deal(max(vertices),min(vertices));
gauss_points = [0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425];
gauss_weights = [0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834];
mapped_gauss_points = (xmax - xmin) ./ 2 .* gauss_points + (xmax + xmin) ./ 2;  
result = (xmax - xmin) ./ 2 * gauss_weights * func(mapped_gauss_points)';
end