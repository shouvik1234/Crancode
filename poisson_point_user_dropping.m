function [ x, y ] = poisson_point_user_dropping( lambda, side_length )
% Inputs: user density, side length of square
% Outputs: x and y coordinates of all users

% We consider a square area of given side length, with corners (0, 0) and
% (side_length)*(1, 1)

numusers = poissrnd(lambda*(side_length^2), [1, 1]);
if ~numusers
    numusers = numusers + 1;
end
% generate the total number of users as a Poisson random variable

x = rand(1, numusers); 
y = rand(1, numusers);
% given the total number of users, the conditional distribution of the 
% users' location is simply uniform over the whole area

end

