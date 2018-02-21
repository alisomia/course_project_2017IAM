function u = bilinear_interpolate(x,y,x0,y0,h,U)
[m,n] = size(U);
j = floor((x - x0)/h)+1;
i = floor((y - y0)/h)+1;
x_res = x - x0 - (j-1)*h;
y_res = y - y0 - (i-1)*h;
if i < m && j < n
u = (U(i,j)*(h-x_res)*(h-y_res)+U(i,j+1)*x_res*(h-y_res)...
    +U(i+1,j)*(h-x_res)*y_res+U(i+1,j+1)*x_res*y_res)/(h^2);
elseif i < m
    u = (U(i,j)*(h-y_res)+U(i+1,j)*y_res)/h;
elseif j < m
    u = (U(i,j)*(h-x_res)+U(i,j+1)*x_res)/h;
else
    u = U(i,j);
end
end