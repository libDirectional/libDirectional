function x_mean = sphMeanShift(x, w)
% @author Kailai Li kailai.li@kit.edu
% @date 2018
x_mean = x(:, 1);
while 1
    x(:, x_mean'*x < 0) = -x(:, x_mean'*x < 0);
    x_t = sphLog(x_mean, x);
    x_mean_t = sum(x_t.*w, 2);
    if norm(x_mean_t) < 1E-6
        break
    end
    x_mean = sphExp(x_mean, x_mean_t);
end
end

function x_t = sphLog(x_center, x)
% spherical logarithm map
%
dot_prod = x_center' * x;
alpha = acos(dot_prod);
x_t = (x - dot_prod .* x_center) ./ sinc(alpha/pi);
end

function x = sphExp(x_center, x_t)
% spherical exponential map
%
norm_t = vecnorm(x_t); % sqrt(sum(x_t.^2,1));
x = cos(norm_t) .* x_center + x_t .* sinc(norm_t/pi);
end