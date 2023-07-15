clear all
clf
p = [0, 0; 30, 70; 35, 72.5; 40, 75; 80, 75; 81, 75];
n = size(p, 1);
n1 = n - 1;
for i=0:1:n1
    sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values
end
l=[];
UB=[];
arySize = 4096;
xAxis = linspace(0, 1, arySize);
for u = 1 : arySize
    for d=1:n
        UB(d)=sigma(d)*((1-xAxis(u))^(n-d))*(xAxis(u)^(d-1));
    end
    l=cat(1,l,UB);
end
P=l*p;
plot(P(:,1),P(:,2))