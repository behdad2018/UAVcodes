clc;clear all;

tspan = [1 5];
ic = [0 0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) odefun(t,y), tspan, ic, opts);