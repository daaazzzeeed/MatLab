a = 0.1;
b = 1.7;

eps = 1*10^(-6);
n = [10,30,200,0];
fid = fopen('C:\Users\Student\Desktop\hist.txt', 'w+');
fclose(fid);
fPointer = @lab02_func; 
diary('C:\Users\Student\Desktop\hist.txt');

k = 1;
while(k < 4 )
    if k == 1
        t1 = zeros(n(k)+1, 1);
        u = 1;
        for l = a : (b-a)/n(k) : b
            t1(u) = l;  
            u = u + 1;
        end
         k = k+1;
    end
    if k == 2
        t2 = zeros(n(k)+1,1);
                u = 1;
         for l = a : (b-a)/n(k) : b
             t2(u) = l;
            u = u + 1;
         end
         k = k+1;
    end
    if k == 3
        t3 = zeros(n(k)+1,1);
                        u = 1;
         for l = a : (b-a)/n(k) : b
            t3(u) = l;
            u = u + 1;
         end
          k = k+1;
    end
    
end

step =4.8828e-06;
t4 = (a : step : b)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funcDiscrT1 = zeros(length(t1),1);  

 for y = 1 : 1 : length(t1)
     funcDiscrT1(y) = lab02_func(t1(y));
 end
 
 
 funcDiscrT2 = zeros(length(t2),1);
 for y = 1 : 1 : length(t2)
     funcDiscrT2(y) = lab02_func(t2(y));
 end
 
  funcDiscrT3 = zeros(length(t3),1);
  
  for y = 1 : 1 : length(t3)
     funcDiscrT3(y) = lab02_func(t3(y));
  end
  
  funcDiscrT4 = zeros(length(t4),1);
  
  for y = 1 : 1 : length(t4)
     funcDiscrT4(y) = lab02_func(t4(y));
  end
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

teta = [ 1/3; 1/15];

N0 = 10;

% global massive;
% massive = zeros(16,1);

d2nArr1 = zeros(1,1);
d2nArr2 = zeros(1,1);
d2nArr3 = zeros(1,16);
d2nArr4 = zeros(1,16);
d2nArr5 = zeros(1,16);
d2nArr6 = zeros(1,1);

[variable1, variable2] = divisions(N0,'leftRect',a,b,teta(1),eps,fPointer,d2nArr1);
% d2nArr1 = massive;
d2nArr1 = variable2;
[variable3, variable4]=divisions(N0,'rightRect',a,b,teta(1),eps,fPointer,d2nArr2);
% d2nArr2 = massive;
d2nArr2 = variable4;
[variable5, variable6]=divisions(N0,'midRect',a,b,teta(1),eps,fPointer,d2nArr3);
% d2nArr3 = massive;
d2nArr3 = variable6;
[variable7, variable8]=divisions(N0,'trapeze',a,b,teta(1),eps,fPointer,d2nArr4);
% d2nArr4 = massive;
d2nArr4 = variable8;
[variable9, variable10]=divisions(N0,'parabolic',a,b,teta(2),eps,fPointer,d2nArr5);
% d2nArr5 = massive;
d2nArr5 = variable10;

fprintf('\n');

[variable11, variable12] = divisions(N0,'leftRect',a,b,teta(1),eps,fPointer,d2nArr6);
n(4) = variable11;
% d2nArr6  = massive;
d2nArr6  = variable12;


 integralDiscr1 = lab02_int_digit( t1, funcDiscrT1 , 'rightRect');
 integralDiscr2 = lab02_int_digit( t2, funcDiscrT2 , 'rightRect');
 integralDiscr3 = lab02_int_digit( t3, funcDiscrT3 , 'rightRect');
 integralDiscr3_ = lab02_int_digit( t4, funcDiscrT4 , 'rightRect');
 
 integralDiscr4 = lab02_int_digit( t1, funcDiscrT1 , 'leftRect');
 integralDiscr5 = lab02_int_digit( t2, funcDiscrT2 , 'leftRect');
 integralDiscr6 = lab02_int_digit( t3, funcDiscrT3 , 'leftRect');
  integralDiscr6_ = lab02_int_digit( t4, funcDiscrT4 , 'leftRect');

 integralDiscr7 = lab02_int_digit( t1, funcDiscrT1 , 'trapeze');
 integralDiscr8 = lab02_int_digit( t2, funcDiscrT2 , 'trapeze');
 integralDiscr9 = lab02_int_digit( t3, funcDiscrT3 , 'trapeze');
  integralDiscr9_ = lab02_int_digit( t4, funcDiscrT4 , 'trapeze');

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 integralAnalog1 = lab02_int_analog(a,b,n(1),fPointer,'rightRect');
 integralAnalog2 = lab02_int_analog(a,b,n(2),fPointer,'rightRect');
 integralAnalog3 = lab02_int_analog(a,b,n(3),fPointer,'rightRect');
  integralAnalog3_ = lab02_int_analog(a,b,n(4),fPointer,'rightRect');

 
 integralAnalog4 = lab02_int_analog(a,b,n(1),fPointer,'leftRect');
 integralAnalog5 = lab02_int_analog(a,b,n(2),fPointer,'leftRect');
 integralAnalog6 = lab02_int_analog(a,b,n(3),fPointer,'leftRect');
  integralAnalog6_ = lab02_int_analog(a,b,n(4),fPointer,'leftRect');

 
 integralAnalog7 = lab02_int_analog(a,b,n(1),fPointer,'midRect');
 integralAnalog8 = lab02_int_analog(a,b,n(2),fPointer,'midRect');
 integralAnalog9 = lab02_int_analog(a,b,n(3),fPointer,'midRect');
  integralAnalog9_ = lab02_int_analog(a,b,n(4),fPointer,'midRect');

 
 integralAnalog10 = lab02_int_analog(a,b,n(1),fPointer,'trapeze');
 integralAnalog11 = lab02_int_analog(a,b,n(2),fPointer,'trapeze');
 integralAnalog12 = lab02_int_analog(a,b,n(3),fPointer,'trapeze');
  integralAnalog12_ = lab02_int_analog(a,b,n(4),fPointer,'trapeze');

 
 integralAnalog13 = lab02_int_analog(a,b,n(1),fPointer,'parabolic');
 integralAnalog14 = lab02_int_analog(a,b,n(2),fPointer,'parabolic');
 integralAnalog15 = lab02_int_analog(a,b,n(3),fPointer,'parabolic');
  integralAnalog15_ = lab02_int_analog(a,b,n(4),fPointer,'parabolic');


 integralToString1 = horzcat(integralDiscr1, horzcat(integralDiscr4, integralDiscr7));
 integralATS1 = horzcat(integralAnalog1, horzcat(integralAnalog4, horzcat(integralAnalog7, horzcat(integralAnalog10,integralAnalog13))));

 integralToString2 = horzcat(integralDiscr2, horzcat(integralDiscr5, integralDiscr8));
  integralATS2 = horzcat(integralAnalog2, horzcat(integralAnalog5, horzcat(integralAnalog8, horzcat(integralAnalog11,integralAnalog14))));

 integralToString3 = horzcat(integralDiscr3, horzcat(integralDiscr6, integralDiscr9));
 integralATS3 = horzcat(integralAnalog3, horzcat(integralAnalog6, horzcat(integralAnalog9, horzcat(integralAnalog12,integralAnalog15))));
 
 integralATS4 = horzcat(integralAnalog3_,integralAnalog6_,integralAnalog9_,integralAnalog12_,integralAnalog15_);
 integralToString4 = horzcat(integralDiscr3_,integralDiscr6_, integralDiscr9_);

 
integralValuesMatr = vertcat(integralToString1, vertcat(integralToString2,integralToString3,integralToString4));
S1 = integralValuesMatr;
integralValuesMatrAnalog = vertcat(integralATS1,integralATS2,integralATS3,integralATS4);
S2 = integralValuesMatrAnalog;

T1 = array2table(S1, 'VariableNames', {'rightRect','leftRect','trapeze'}, 'RowNames', {'N = 10','N = 30','N = 200','N = 327680'});
T2 = array2table(S2, 'VariableNames', {'rightRect','leftRect','midRect','trapeze', 'parabolic'}, 'RowNames', {'N = 10','N = 30','N = 200','N = 327680'});
 

d2nArr1_ = [0.0167,0.0093,0.0049 ,0.0025 ,0.0013 , 6.3753e-04 ,  3.1968e-04,  1.6007e-04,  8.0091e-05 ,  4.0060e-05 , 2.0033e-05 ,  1.0018e-05, 5.0090e-06 , 2.5046e-06, 1.2523e-06, 6.2615e-07];
fxx=subplot(1,1,1);
figure(1);
plot(t1,funcDiscrT1);
hold on;
plot(t2,funcDiscrT2);
hold on;
plot(t3,funcDiscrT3);
hold on;
plot(t4,funcDiscrT4);
hold off;
fxx.Title.String = 'график подынтегральной функции';
fxx.XLabel.String = 'x';
fxx.YLabel.String = 'f(x)';
legend('n = 10', 'n = 30','n = 200','n = 327680');

figure(2);
n_= (0:1:15);
fx =subplot(1,1,1);
plot(n_,d2nArr1_);
hold on;
plot(n_,d2nArr2);
hold on;
d2nArr3 = horzcat(d2nArr3);
for i=8:1:16
    d2nArr3(i,1) = 0;
    d2nArr4(i,1) = 0;
    d2nArr5(i,1) = 0;
end
plot(n_,d2nArr3);
hold on;
plot(n_,d2nArr4);
hold on;
plot(n_,d2nArr5);
hold off;
legend('rightRect','leftRect','midRect','trapeze','parabolic');


fx.Title.String = 'график зависимости d2n от N' ;
fx.XLabel.String = 'N, степени двойки';
fx.YLabel.String = 'delta2n';


barData = [S1(1,1),S1(1,2),S1(1,3);
   S1(2,1),S1(2,2),S1(2,3);
   S1(3,1),S1(3,2),S1(3,3);
   S1(4,1),S1(4,2),S1(4,3)];

figure(3);
b = bar(barData);
ax = gca;
xticks([1 2 3 4])
xticklabels({'N = 10','N = 30','N = 200', 'N = 327680'})
legend('rightRect','leftRect','midRect','trapeze','parabolic');
ax.YLabel.String = 'Value';
ax.Title.String = 'S1 Bar Graph';

legend('rightRect','leftRect','trapeze');

barData2 = [S2(1,1),S2(1,2),S2(1,3),S2(1,4),S2(1,5);
   S2(2,1),S2(2,2),S2(2,3),S2(2,4),S2(2,5);
   S2(3,1),S2(3,2),S2(3,3),S2(3,4),S2(3,5);
   S2(4,1),S2(4,2),S2(4,3),S2(4,4),S2(4,5)];

figure(4);
b2 = bar(barData2);
ax2 = gca;
ax2.YLabel.String = 'Value';
ax2.Title.String = 'S2 Bar Graph';
xticks([1 2 3 4])
xticklabels({'N = 10','N = 30','N = 200', 'N = 327680'})
legend('rightRect','leftRect','midRect','trapeze','parabolic');
fprintf('\n');
fprintf('S1');
fprintf('\n');
disp(S1);
fprintf('\n');
fprintf('S2');
fprintf('\n');
disp(S2);
fprintf('\n');
fprintf('T1');
fprintf('\n');
disp(T1);
fprintf('\n');
fprintf('T2');
fprintf('\n');
disp(T2);

diary off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [divs, array_] = divisions(N0_, methodToString, a_,b_,teta_,eps_,fPointer_, array)
% global massive;
% 
integral1 = lab02_int_analog(a_,b_,N0_,fPointer_,methodToString);
N2 = 2*N0_;
integral2 = lab02_int_analog(a_,b_,N2,fPointer_,methodToString);

acc = accuracy(teta_,integral1,integral2);
disp(acc); 
n = N0_;
i = 1;
array(i,1) = acc;
massive(i,1) = array(i,1);
array_(i,1) = acc;
i = i + 1;
while(acc > eps_)
    
   n = n*2;
   n2 = 2*n;
   integral1 = lab02_int_analog(a_,b_,n,fPointer_,methodToString);
   integral2 = lab02_int_analog(a_,b_,n2,fPointer_,methodToString);
   acc = accuracy(teta_,integral1,integral2);
    
   array(i,1) = acc;
   array_(i,1) = acc;
   massive(i,1) = array(i,1);
   i = i + 1;
   disp(acc);

end
divs = n;
 fprintf(methodToString);
 fprintf('\n');
 fprintf('N = %2f', n);
 fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function delta2n = accuracy(teta_, I1, I2)
delta2n = teta_*abs(I2-I1);
end

%  test
% if(((n(1)) + 1 ) == length(funcDiscrT1))
%     disp(length(funcDiscrT1));
% end

%  test
%  integralSum = 0;
%  for i = 2 : 1 : (200)
%    integralSum = integralSum + funcDiscrT3(i)*(t3(i) - t3(i-1));
%  end
%  disp(integralSum);


function f  = lab02_func(t)
w1 = 2*pi*1; 
w2 = 2*pi*10; 
f = sin(w1*t) + 0.2*cos(w2*t);
end


function analogInt = lab02_int_analog(a,b,n,fPtr,methodToStr)
integralSum = 0;
integralSum1 = 0;
integralSum2 = 0;

switch methodToStr
    case ('rightRect')
        for i = 1 : 1 : (n)
            integralSum = integralSum + fPtr( (a+((i-1+1)*((b-a)/n)) ))*( (a+((i-1+1)*((b-a)/n)) ) - ( (a+((i-1)*((b-a)/n)) ) ));
         end
    case ('midRect')
        for i = 1 : 1 : (n)
            integralSum = integralSum + fPtr( 0.5*(a+((i-1)*((b-a)/n)) + (a+((i-1+1)*((b-a)/n)))))*((a+((i-1+1)*((b-a)/n)) - (a+((i-1)*((b-a)/n)))));
        end
    case ('leftRect')
        for i = 0 : 1 : (n-1)
                    integralSum = integralSum + fPtr( (a+((i)*((b-a)/n)) ))*( (a+((i+1)*((b-a)/n)) ) - ( (a+((i)*((b-a)/n)) ) ));    
        end
    case ('trapeze')
        for i = 0 : 1 : (n-1)
            integralSum = integralSum +  0.5*( fPtr(a+((i)*((b-a)/n))) +fPtr(a+((i+1)*((b-a)/n))) )*((a+((i+1)*((b-a)/n))) - (a+((i)*((b-a)/n))) );
        end
        
    case ('parabolic')
        h = (b-a)/2/n;
%         N = 2*n;
%         ((b-a)/N)
         for j = 1 : 1 : (n-1)
         integralSum1 = integralSum1 + fPtr((a+(2*(j-1+1)*h)));
         end
        
       for j = 1 : 1 : n
        integralSum2 = integralSum2 + fPtr((a+((2*j-1)*h)));
       end
       
        integralTempSum = (h/3)*(fPtr(a) + fPtr(b) + 2*integralSum1 + 4*integralSum2);
        integralSum = integralTempSum;
       
    otherwise 
        disp('incorrect parameters!');
end
    analogInt = integralSum ;
end


function integralDiscr = lab02_int_digit( t, f , method)
integralSum = 0;
N = length(f) - 1;

 switch(method) 
     
           case ('leftRect')   
        for i = 0 : 1 : (N-1)
            integralSum = integralSum + f(i+1)*(t(i+1+1) - t(i+1));
        end
        
         case ('rightRect')   
         for i = 1 : 1 : (N)
            integralSum = integralSum + f(i+1)*(t(i+1) - t(i));
           
         end
         
      case ('trapeze') 
         for i = 0 : 1 : (N-1)
            integralSum = integralSum + 0.5*(f((i+1)) + f(i+1+1))*(t(i+1+1)-t(i+1));
         end
         
      otherwise 
        disp('incorrect parameters!');
 end
        integralDiscr = integralSum;

end
