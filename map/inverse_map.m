function [x_hat,t_hat] = inverse_map(x,t,up_front,lower_front,delta)
x_hat = x;
taui = up_front(x);
tuai_1 = lower_front(x);
t_hat = abs((t-tuai_1)./(delta));
end