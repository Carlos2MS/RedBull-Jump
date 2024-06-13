%Red Bull Data plotted.
close all
clear 

% Read data from file
% First column is time in secs
% Second column is velocity
% Third column is theoretical terminal velocity at that time
jumpdata = csvread('RedBullJumpData.csv');
t_redbull = jumpdata(:,1);
v_redbull = jumpdata(:,2);
terminal_velocity = jumpdata(:,3); % You need to use this in the last question
N_timestamps = length(t_redbull);

%Calculate freefall velocity vector here
g = 9.81;
v_freefall = g * t_redbull;
%rand(size(t_redbull)) * 20 + 100;

% Part 1
% We are giving you this answer for free
figure(1);
h_part1 = plot(t_redbull, v_redbull, 'r-x', 'linewidth', 2.0);
shg;
hold on;

% Part 2
% This plot does not have the right linewidth. You fix it.
h_part2 = plot(t_redbull, v_freefall, 'k--','linewidth',2.0); 
shg;
% This is how to put on a grid 
grid on;
% This is how to fix an axis to a desired size
% axis([something goes in here]);
axis([0 180 0 400]);
% Set the fontsize and label the graph here!!
ylabel('Velocity m/s', 'fontsize', 24,'Color','b');
xlabel('Time (secs)', 'fontsize', 24,'Color','b');

% Calculate when he hits the atmosphere
% Part  3
% Need some stuff here ... or read it off from the graph
v_err = v_redbull - v_freefall;

z= abs(v_err)./v_freefall;
n= find((z) >= 0.02);

hit_instant = t_redbull(min(n));
fprintf('Mr. B hits the earth''s atmoshpere at %f secs after he jumps\n',...
       hit_instant);

% Part 4
% Now starting from the velocity and time = 56 secs 
% .. let's update and calculate v
g = 9.81;
drag_constant = 3/60;
v_numerical_1 = zeros(1, N_timestamps);
start = find(t_redbull == 56);
% Starting from this time instant, calculate the velocity required

v_numerical_1(1 : start) = v_redbull(1 : start);

%v_numerical_1(start) = (drag_constant^-1)*g*(1-exp((-drag_constant)*start));; % This is wrong .. for you to fix

for i = start : N_timestamps- 1
    h = t_redbull(i+1) - t_redbull(i);
    dvdt = g - (drag_constant*v_numerical_1(i));
    v_numerical_1(i+1) = v_numerical_1(i) + dvdt*h;
end


% Plot using the dashed green line with (+) markers
h_part4 = plot(t_redbull(start:N_timestamps), v_numerical_1(start:N_timestamps), 'g--o','linewidth',2.0,'MarkerSize',2.5);shg

% Part 5 
% Calculate the percentage error as required
t1= find(t_redbull == 69);
t2 = find(t_redbull == 180);

diff_1= abs(v_numerical_1(t1)-v_redbull(t1));
per_error(1) = (diff_1/v_redbull(t1))*100;% This is just some random number .. which is wrong

diff_2 = abs(v_numerical_1(t2)-v_redbull(t2));
per_error(2) = (diff_2/v_redbull(t2))*100;

fprintf('The percentage error at 69 and 180 secs is %1.1f and\n', per_error(1));
fprintf('%3.1f  respectively \n', per_error(2));

% Part 6 
% You'll need to modify your euler loop here to do heun instead
% Also update the drag constant at every timestamp and change the update
% calculation to allow for the new v^2(t) term
% A hint here that now you have to calculate the velocity using the new
% differental equation
% constant .. put it in v_numerical_2

v_numerical_2 = zeros(1, N_timestamps);
drag_mass = zeros(1,N_timestamps);

v_numerical_2(1 : start) = v_redbull(1 : start);

end_t = find(t_redbull == 100);

a=1/2;

for i = start : end_t-1
    h = t_redbull(i+1) - t_redbull(i);

    %drag measurements are shifted by one time interval
    drag_mass_1= g/(terminal_velocity(i+1)^2);
    drag_mass_2= g/(terminal_velocity(i+2)^2);

    
    k1= g -(drag_mass_1*(v_numerical_2(i)^2));
    k2= g - (drag_mass_2*(v_numerical_2(i) + k1*h)^2);
    
    v_numerical_2(i+1) = v_numerical_2(i) + a*(k1+k2)*h;
end



% This is the handle plot for part 6. You have to plot the right stuff not
% this stuff.
% Note that the plot linewidth and colour are wrong. Fix it.
h_part6 = plot(t_redbull(start:end_t), v_numerical_2(start:end_t), 'k--+','lineWidth',2.0);
shg

%Error at t=100, end_t has the index of t=100

diff_1= abs(v_numerical_2(end_t)-v_redbull(end_t));
est_error = (diff_1/v_redbull(end_t))*100;
fprintf('The error at t = 100 secs using my estimated drag information is %f\n', est_error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE THIS IS TO MAKE SURE YOU HAVE USED THE
% VARIABLES THAT WE ASKED FOR
% Check for existence of variables
if (~exist('v_freefall', 'var'))
  error('The variable v_freefall does not exist.')
end;
if (~exist('hit_instant', 'var'))
  error('The variable hit_instant does not exist.')
end;
if (~exist('per_error', 'var'))
  error('The variable per_error does not exist.')
end;
if (exist('per_error', 'var'))
  l = size(per_error);
  if ( sum(l - [1 2]) ~= 0)
    error('per_error is not a 2 element vector. Please make it so.')
  end;
end;
if (~exist('v_numerical_1', 'var'))
  error('The variable v_numerical_1 does not exist.')
end;  
if (~exist('est_error', 'var'))
  error('The variable est_error does not exist.')
end;  
if (~exist('h_part1', 'var'))
  error('The plot handle h_part1 is missing. Please create it as instructed.')
end;
if (~exist('h_part2', 'var'))
  error('The plot handle h_part2 is missing. Please create it as instructed.')
end;
if (~exist('h_part4', 'var'))
  error('The plot handle h_part4 is missing. Please create it as instructed.')
end;
if (~exist('h_part6', 'var'))
  error('The plot handle h_part6 is missing. Please create it as instructed.')
end;
