adresseQ='C:\Users\Perle\Desktop\VIA\ddpls.txt';
Q = textread(adresseQ);
[n,p] = size(Q);

Q1 = Q(1:3:n-2);
Q2 = Q(2:3:n-1);
Q3 = Q(3:3:n);

xeff = zeros(1,n/3);
yeff = zeros(1,n/3);

for i = 1:n/3
    xeff(i) = cos(Q1(i)+Q2(i)+Q3(i)) + cos(Q1(i)+Q2(i)) + cos(Q1(i));
    yeff(i) = sin(Q1(i)+Q2(i)+Q3(i)) + sin(Q1(i)+Q2(i)) + sin(Q1(i));
end;

plot(xeff,'g');
hold on;
plot(yeff,'b');
%legend('xeff','yeff');
   

adresseQ='C:\Users\Perle\Desktop\VIA\q1.txt';
Q = textread(adresseQ);
[n,p] = size(Q);

Q1 = Q(1:3:n-2);
Q2 = Q(2:3:n-1);
Q3 = Q(3:3:n);

xeff = zeros(1,n/3);
yeff = zeros(1,n/3);

for i = 1:n/3
    xeff(i) = cos(Q1(i)+Q2(i)+Q3(i)) + cos(Q1(i)+Q2(i)) + cos(Q1(i));
    yeff(i) = sin(Q1(i)+Q2(i)+Q3(i)) + sin(Q1(i)+Q2(i)) + sin(Q1(i));
end;

plot(xeff,'g--');
hold on;
plot(yeff,'b--');
legend('xeff ddp LS','yeff ddp LS','xeff ddp classique','yeff ddp classique');