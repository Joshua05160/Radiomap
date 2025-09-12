function [Distance,Distance_xy] = calculate_distance(ALPHA, Height, BETA)
M = size(ALPHA,2);
T = size(ALPHA,3);
N = size(BETA,2);
    for t = 1:T
        for m = 1:M
            for n = 1:N
                horizontal_distance_squared = sum((ALPHA(:, m, t) - BETA(:, n, t)).^2);
                vertical_distance_squared = Height(m, t)^2;
                Distance(m, n, t) = sqrt(horizontal_distance_squared + vertical_distance_squared);
                Distance_xy(m,m,t) = sqrt(horizontal_distance_squared);
            end
        end
    end
end
