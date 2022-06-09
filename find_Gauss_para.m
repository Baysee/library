function [s,c,vals]=find_Gauss_para(p,q,GivVals)

q_odd=mod(q,2); %returns 1 if odd
p_odd=mod(p,2);
pq_odd=mod(q*p,2);
%put mod 2q for s


if q_odd && p_odd==0
    s=mod(p*(f_mod_inv(p,q))^2,2*q);
    c=(q-1)/4+(1-JacobiSymbol(p,q))/2;
end

if p_odd && q_odd==0
    s=mod(p*(f_mod_inv(p,q))^2,2*q);
    c=-p/4-(1-JacobiSymbol(q,p))/2;
end

if p_odd && q_odd
    s=mod(8*p*f_mod_inv(2,q)*(f_mod_inv(2*p,q))^2,2*q);
    c=(q-1)/4+(1-JacobiSymbol(p,q))/2;
end



if GivVals==1
    vals=zeros(q,1);
for i=1:q
    
   vals(i)=mod(pi*s/q*(i-1)^2,2*pi);%mod(s/q*pi*(i-1)^2+c*pi,2*pi); 

end
else
    vals=0;
end

end
