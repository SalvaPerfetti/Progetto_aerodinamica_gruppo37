function M = pendenza_per_tratto(x,y)
M = [];
for i = 2 : size(x)
    m = (y(i)-y(i-1))/(x(i)-x(i-1));
    M = [M;m];
end
end

