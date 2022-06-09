% ps=1:10;
% qs=2:11;
% 
% s1=zeros(numel(ps),numel(qs));
% s2=s1;
% 
% for p_i=1:numel(ps)
%     for q_i=1:numel(qs)
%     p = ps(p_i) ; %assign current loop val
%     q = qs(q_i) ;
%     
%     if gcd(p,q)==1
%           s1(p_i,q_i)=generateSparameter2(p,q)   
%     else
%         s1(p_i,q_i)=0;
%     end
%         
%     end
% end
% 
function s=generateSparameter(p,q)

gcdpq=gcd(p,q);
if gcdpq~=1
   error('p and q must be coprime!!')
end

if mod(q,2)==1
    s=2*mulinv(2*p,q); %% sp=1+qe_q (mod 2q) -> s(2p)=1 mod(q)  (first comment from some time ago)
%     s=2*mulinv(2*p,q); %% sp=1+qe_q (mod 2q) -> s(2p)=2+2q (mod(2q)) =2 mod (2q)
else
    s= mulinv(p,2*q) ;
end

end


function y=mulinv(x,p)
% 
% if ~isprime(p)
%     disp('The field order is not a prime number');
%     return
% elseif sum(x>=p)
%     disp('All or some of the numbers do not belong to the field');
%     return
% elseif sum(~x)
%     disp('0 does not have a multiplicative inverse');
%     return
% end

k=zeros(size(x));   %set the counter array
m=mod(k*p+1,x);     %find reminders
while sum(m)        %are all reminders zero?
    k=k+sign(m);    %update the counter array
    m=mod(k*p+1,x); %caculate new reminders 
end
y=(k*p+1)./x;       %find the multiplicative inverses of X
end

