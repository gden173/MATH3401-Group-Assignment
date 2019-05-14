% Math 3401 Complex Analysis Assignment 5
% Group Project
% Authors : Gabriel Dennis
%------------------------------------------
%------------------------------------------
%Question 1.b
% Partition \bar{\Om} into five concentric circles of radius |z| = N/5, 
% N =  1,2,..,5 and 8 radial lines at Arg(z) = -3pi/4,...,pi
%-----------------------------------------
clear all;  clc; clf;  % Clear all saved variables and figures
N = 24;                % set number of line segments 
% plot all five circles on the same graph, blue lines
plt_ax = [-1.1,1.1, -1.1,1.1];  % plot axis

coords = pltcords(N);

% plot on figure using blue for circles and red dotted for radial lines 
f = figure;

hold on 
for i = 1:5    
   plot(real(coords{i}), imag(coords{i}), 'b', 'Linewidth', 2)  % plot circles 
end
for i = 6:13    
   plot(real(coords{i}), imag(coords{i}), 'r.', 'Linewidth', 2)  % plot circles 
end

axis equal; grid on;
hold off
axis(plt_ax)
xlabel('Re', 'Fontsize', 20); ylabel('Im', 'Fontsize', 20)  % x,y labels
title('$$z$$','Interpreter','latex', 'Fontsize', 20)  % title 

%%
%-------------------------------------------------------------------
% Question 1.c
% Make two subplots with the pre image and image on the same plot
% We already have the pre image, now just creat the image
% function is f(z) = (1+z)^2/4

mapper = @(z)  (1 + z).^2/4;

for i = 1:13
    map{i} = mapper(coords{i});
end
subplot(1,2,1)
hold on 
for i = 1:5    
   plot(real(coords{i}), imag(coords{i}), 'b', 'Linewidth', 2)  % plot circles 
end
for i = 6:13    
   plot(real(coords{i}), imag(coords{i}), 'r.', 'Linewidth', 2)  % plot circles 
end
hold off; grid on 
axis equal; axis(plt_ax)
subplot(1,2,2)
hold on 
for i = 1:5    
   plot(real(map{i}), imag(map{i}), 'b', 'Linewidth', 2)  % plot circles 
end
for i = 6:13    
   plot(real(map{i}), imag(map{i}), 'r.', 'Linewidth', 2)  % plot circles 
end

hold off
grid on; axis equal;
axis(plt_ax);
title('$$f(z)$$', 'Interpreter','latex', 'Fontsize', 20)
xlabel('Re', 'Fontsize', 20)
ylabel('Im', 'Fontsize', 20)

%% Bonus question
%-------------------------------------------------------------------
% define a new function mapb, f(z) = (i + iz)/(1-z)
% To get different values for N it must be changed in 1.B.
% Set number of steps
  % use length of z,z1 as a guide 
  %set N
  N = 10;
 coords = pltcords(N);
mapper = @(z)  (1 + 1i*z)./(1-z);

for i = 1:13
    map{i} = mapper(coords{i});
end

hold on 
for i = 1:5    
   plot(real(map{i}), imag(map{i}), 'b', 'Linewidth', 2)  % plot circles 
end
for i = 7:13    
   plot(real(map{i}), imag(map{i}), 'r.', 'Linewidth', 2)  % plot circles 
end
 plot(real(map{6}), imag(map{6}), 'k.', 'Linewidth', 2)  % plot circles
hold off
grid on; axis equal;
axis([-10,10, -10,10]);
title(['$$f(z)$$, N = ',num2str(N)], 'Interpreter','latex', 'Fontsize', 20)
xlabel('Re', 'Fontsize', 20)
ylabel('Im', 'Fontsize', 20)

%%
%--------------------------------------------------------------------------
%Q2:  First plot the curve f(gamma), |z| = 2
% set x and y
num = linspace(0,2*pi, 100);
x = 2*cos(num);
y = 2*sin(num);
z = complex(x,y);  % create the complex coordinates
func = @(z) ((1 + z).^2/4);  % function
figure
plot(real(func(z)),imag(func(z)),'b', 'Linewidth',3)
grid on 
axis equal 
axis([-1.1, 3.1, -2.1, 2.1])
xlabel('Re','Fontsize', 20)
ylabel('Im','Fontsize', 20)
title('Plot of $$f(z) = \frac{(1 + z)^2}{4}$$ on the curve $$ \gamma = \{z:|z| = 2\}$$',...
    'Interpreter','latex', 'Fontsize',20)
print('gamma_curve','-dpng')


%% Contour integration function
% Q2.a and Q2.b
% Only works for contour integrals on a positvely orientated closed circle

%%
%-------------------------------------------------------------------------
%Q2.a
func = @(z) ((1 + z).^2/4);  % function
N = 100;  %steps
r = 2; % radius
[int1] = trap(func, r,N);  % int
fprintf('The Value of the Contour integral is %4.2f   \n', int1)
%% 
%--------------------------------------------------------------------------
%Q2.b

func = @(z) norm((1 + z).^2/4).^2 ;  % norm squared function 
N = 100;
r = 2;
int2 = trap(func, r,N);
fprintf('The Value of the Contour integral is %4.2f   \n', int2/(pi*1i)) 
%% 
%--------------------------------------------------------------------------
% Bonus question winding number
func = @(z) (0.5*(1 + z))/((1/4)*(1 + z).^2 + 1/2);
N = 100;
r = 2;
int3 = trap(func,r,N);  
fprintf('To get the winding number divide by 2*pi*i which results in a value of  %4.2f  \n', int3/(2*pi*1i))
%%  Contour Integration Function


function [coords] = pltcords(N)

theta = linspace(0,7*pi/4,8);
theta1 = linspace(0,2*pi,N);
pts = linspace(0,1,N);
rad = linspace(1,1/5,5);


for j = 1:5
    x1 = rad(j)*cos(theta1);
    y1 = rad(j)*sin(theta1);
    z = complex(x1,y1);
    coords{j} = z;
end
k = 6;
for j = 1:8
    x1 = pts*cos(theta(j));
    y1 = pts*sin(theta(j));
    z = complex(x1,y1);
    coords{k} = z;
    k = k+1;
end

end




function [int] = trap(func, r,N)
% [int] = trap(func, r,N) takes a function handle FUNC
% a radius r and number of steps N
% it then uses the trapexoidal rule to numerically approximate the contour
% integral on a closed positively orientated contour centred at the origin
steps = linspace(0,2*pi, N);  % number of intervals
y = r*sin(steps);  % polar x,y coordinates
x = r*cos(steps);
z = complex(x,y);  % complex vector
int = 0;
for k = 1:length(z)-1
    int = int + 0.5*(func(z(k+1)) + func(z(k)))*(z(k+1)-z(k));  % Trap rule
    
end
% Make it easier to see if the integral is zero
if(norm(int)) < 10^(-5)
    int = 0 ;
end

end

