function y = singleGauss(std,center,xs,prop)


if prop==1
    y=1/(sqrt(2*pi*std^2))*exp(-(xs-center).^2/(2*std^2));
else
    y=exp(-(xs-center).^2/(2*std^2));
end
