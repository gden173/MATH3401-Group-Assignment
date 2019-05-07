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
N = 100;                % set number of line segments 
% plot all five circles on the same graph, blue lines
plt_ax = [-1.1,1.1, -1.1,1.1];  % plot axis
x = linspace(0,2*pi,N);  % theta from 0 to 2pi
z1 = linspace(-1/sqrt(2),1/sqrt(2),N*2); z = linspace(-1,1,N*2); % coordinates for  diagonal lines   


circles = zeros(10,length(x)); % Empty matrix to place coordinates in 
rad = [1/5,2/5, 3/5, 4/5, 1];   % vector of circl radii

for i = 1:2:10  % plot x,y coordinates using polar coordinates x = rcos(theta), y = rsin(theta)
    circles(i,:) = rad(ceil(i/2))*cos(x);
    circles(i+1,:) = rad(ceil(i/2))*sin(x);
end

% plot on figure using blue for circles and red dotted for radial lines 
f = figure;
subplot(1,2,1)
hold on 
for i = 1:2:10    
   plot(circles(i,:), circles(i+1,:), 'b', 'Linewidth', 2)  % plot circles 
end
axis equal; grid on;
plot(z1, z1, 'r.', z ,zeros(1,length(N)),'r.',zeros(1,length(N)), z, 'r.', z1, -z1, 'r.', 'MarkerSize',6); % plot lines 
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

map = @(x,y) [(1 + 2*x + x.^2 - y.^2)/4 ; (y + x.*y)/2];  % map = [ u(x,y), v(x,y)]
mat = zeros(10,N);   % empty matrix 
line1 = map(z1,z1);line2 = map(z1,-z1); line3 = map(z, zeros(1,length(N))); line4 =  map(zeros(1,length(N)),z);  % map lines
for i = 1:2:10
    mat(i:i+1,:) = map(circles(i,:), circles(i+1,:));  % map circles
end
subplot(1,2,2)
hold on 
for i = 1:2:10    
   plot(mat(i,:), mat(i+1,:), 'b', 'Linewidth', 2)  % plot circles 
end
plot(line1(1,:), line1(2,:),'r.',line2(1,:), line2(2,:),'r.',line3(1,:), line3(2,:),'r.',...
    line4(1,:), line4(2,:),'r.','MarkerSize',6)   % plot lines 
hold off
grid on; axis equal;
axis([-1.1,1.1, -1.1,1.1]);
title('$$f(z)$$', 'Interpreter','latex', 'Fontsize', 20)
xlabel('Re', 'Fontsize', 20)
ylabel('Im', 'Fontsize', 20)

%% Bonus question
%-------------------------------------------------------------------
% define a new function mapb, f(z) = (i + iz)/(1-z)
% To get different values for N it must be changed in 1.B.

mapb = @(x,y)[-2*y./((1-x).^2 + y.^2); (1-x.^2 - y.^2)./((1-x).^2 + y.^2)];  % f(z) = u(x,y) + iv(x,y)

mat = zeros(10,N); % redefine mat
line1 = mapb(z1,z1);line2 = mapb(z1,-z1); line3 = mapb(z, zeros(1,length(N))); line4 =  mapb(zeros(1,length(N)),z); % map lines
for i = 1:2:10
    mat(i:i+1,:) = mapb(circles(i,:), circles(i+1,:));  % map circles
end
figure   % plot co domain 
hold on 
for i = 1:2:10    
   plot(mat(i,:), mat(i+1,:), 'b', 'Linewidth', 2)
end
plot(line1(1,:), line1(2,:),'r.',line2(1,:), line2(2,:),'r.',line3(1,:), line3(2,:),'r.',...
    line4(1,:), line4(2,:),'r.','MarkerSize',6)
hold off
grid on; axis equal;
axis([-20,20,-20,20])
title('$$f(z)$$', 'Interpreter','latex', 'Fontsize', 20)
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
func = @(z) ((1 + z).^2/4);
N = 100;  %steps
r = 2; % radius
[int1] = trap(func, r,N);  % int
disp(int1)  % Integral is equal to zero as expected 
%% 
%--------------------------------------------------------------------------
%Q2.b

func = @(z) norm((1 + z).^2/4).^2 ;  % norm squared function 
N = 100;
r = 2;
int2 = trap(func, r,N);
disp(round(int2/(1i*pi),1))  % integral is equal to 5 which is what my analytic integration found
%|f(z)|^2 is not analytic at zero so integral is not equal to zero


%%  Contour Integration Function

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

