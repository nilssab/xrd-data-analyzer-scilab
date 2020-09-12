function [y] = y_func(x, x_c, width, Area, y0_a, y0_b)
    Kalpha2 = 0.179285  //Co wavelengths
    Kalpha1 = 0.1788965
    mu = 1.5            //shape factor (1.5 is used in JADE 5)
    
    y = x.*y0_a+y0_b
    for i = 1:length(x_c)
        x_c_2(i) = 360.*asin(Kalpha2./(2.*(Kalpha1./(2.*sin(%pi*x_c(i)./(360))))))./(%pi)   //Kalpha 2 peak position
        Area_2(i) = Area(i)/2                                                               //alpha 2 area
        width_2(i) = width(i)                                                               //alpha 2 width
        
        y = y + Area(i)*((2*gamma(mu)*sqrt(2^(1/mu)-1))/(sqrt(%pi)*width(i)*gamma(mu-1/2)))*(1+4*((2^(1/mu)-1)/width(i)^2).*(x-x_c(i)).^2).^(-mu) + Area_2(i)*((2*gamma(mu)*sqrt(2^(1/mu)-1))/(sqrt(%pi)*width_2(i)*gamma(mu-1/2)))*(1+4*((2^(1/mu)-1)/width_2(i)^2).*(x-x_c_2(i)).^2).^(-mu)
    end
    
endfunction
