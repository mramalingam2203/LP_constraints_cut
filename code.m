a = [1,3; 2,1; 2,3; 1,1; 4,5];
b = [15;10;18;7;40];
c = [3;4];

brows = rows(b)
crows = rows(c)

lambda = zeros(brows,1);
zl = zeros(brows,1);
zb = zeros(brows,1);
bixj = zeros(brows,1);

sum =0;

for j = 1:brows
  for i = 1:crows
    sum = sum + a(j,i);
  endfor
  #disp(sum);
  lambda(j) = b(j)/sum;
  #disp(lambda(j));
  sum = 0;
endfor


[minval, idx] = min(lambda, [], 1);

for j = 1:brows
    sum = 0;
    for i = 1:crows
      sum = sum+c(i)*lambda(j);
    endfor
    zl(j) = sum;
endfor

#disp(z)

#for j=1:brows
  for j=1:crows
    bixj(j) = b(j)/a(j,1);
    disp(a(j,1))
  endfor
#endfor

bixj = b(:,1)./(a(1:brows))'

zb = (c(1,:))'.*bixj

r2 = (zb(2)-zl(2))/(bixj(2)-lambda(2))

r1 = (zb(1)-zl(1))/(bixj(1)-lambda(1))

#disp(bixj)



