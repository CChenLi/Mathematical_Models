# Mathematical Modeling

Mathematical examples from many areas, including population growth, epidemiology, cell biology, and traffic modeling. The models take the form of ordinary differential equations (ODEs), difference (recurrence) equations, and partial differential equations (PDEs).

## Before you start

You are suggested to have some basic math background before start this project. Some prerequisites are:
* Linear algebra
* Differential equations
* Recurrence relation
* Calculus
* Proobability theory
* Stocastic process


### An example on Modeling Cov-19 with SIR model

The following code is a SIR mathematical model that fitted the first 14 days of Cov-19 spreading data in U.S. 
from Center for Disease Control and Prevention. 
The code is based on R. 
Default assumption is enforce Social distancing policy 10 days after outbreak.
Explaination of parameters are in comment, feel free to adjust these parameters and check how these parameter
affect feature of spreading.

```
#Code for simulating COVID-19 outbreak

# Parameters
dt = 0.01; # time step (days)
tvecShort = seq(0, 14, length.out = dt); # vector for fitting CDC data
tsteps = length(tvecShort);

i = 1; #initial index
tvec = 0; #initial time

day_imposed = 20; # number of days after first case that regulations imposed (default is 20)

cutback = 0.2; # proportion of usual contacts you have with stay-at-home rules

infectivity_i = 0.12;  # prob. of infected transmitting infection per contact [1/pop] (Burke et al.)
infectivity_e = infectivity_i;  # prob. of exposed transmitting infection per contact [1/pop]

contact_i = 2.4;  # number of people infecteds encounter per day [pop/day] (Burke et al.)
contact_e = 5.8;  # number of people exposeds encounter per day [pop/day] (fit to CDC data)

# infection rates [1/day]
beta_i = infectivity_i*contact_i;  # infected
beta_e = infectivity_e*contact_e;  # exposed

eta = 0.2;  # latency [1/time] (Lauer et al.)
p = 0.25;  # probability of being asymptomatic [dimensionless]
mu = 1/14;  # rate of leaving infective state [1/day]
delta = 1/200;  # probability of death from infection [dimensionless]

# Data for parameter fitting and comparison
raw_data = c(1, 0, 8, 6, 23, 25, 20, 66, 47, 64, 147, 225, 290, 278, 414);  # first 14 days of new cases in US, starting 2/27 (CDC)
cdc_data = cumsum(raw_data);  # total number of cases for first 14 days
data = cdc_data/0.06;  # 6% of actual cases are detected (worldwide) Verity et al.

# Initial data
N = 3.28*10^8; # total population of U.S.
E = p*data[1]; # # of exposed 
I = (1-p)*data[1]; # # of infected
S = N-I-E; # # of susceptibles
R = 0; # # of recovered
D = 0; # # of dead

# Simulations

while(tail(I, 1)>10){ # run loop until disease has run its course
    # stay-at-home rules imposed
    if(i == round(day_imposed/dt)+1){
        beta_e = cutback*beta_e;
        beta_i = cutback*beta_i;
    }
        
    # recurrence equations
    S[i+1] = S[i] + dt*(-S[i]*(beta_i*I[i]+beta_e*E[i])/N);
    E[i+1] = E[i] + dt*(S[i]*(beta_i*I[i]+beta_e*E[i])/N - eta*E[i]);
    I[i+1] = I[i] + dt*(eta*(1-p)*E[i] - mu*I[i]);
    R[i+1] = R[i] + dt*(mu*(1-delta)*I[i] + eta*p*E[i]);
    D[i+1] = D[i] + dt*mu*delta*I[i];
    tvec[i+1] = tvec[i] + dt; #update time vector
    i = i+1; #increment index
}

print('Proportion Affected:')
print(((tail(R, 1)+tail(D, 1))/N))
# Solution plots
par(mfrow = c(2, 2))

#plot S, I, and R over time
plot(tvec, S, lwd = 1, type='l',xlim=c(0, tail(tvec, 1)), ylim = c(0, S[1]))
lines(tvec, R, lwd = 2)
legend('topright',c('S','R')) 

plot(tvec, E, type = 'l', lwd = 1, xlim=c(0, min(tail(tvec,1),730)), ylim = range(I))
lines(tvec, I)
legend('topright',c('E','I')) #add legend

plot(tvec, D, type='l',lwd = 1, xlim = range(tvec))

#set(gcf,'Position',[100 300 1200 400])

plot(seq_along(data), data, pch = '*', xlab = 'Days since February 27', ylab= 'Number of cases')
points(tvecShort, I[1:tsteps]+E[1:tsteps], lwd = 1)
```

## Author

* **Chen Li** 


## Acknowledgments

* This project is based on Math142 UCLA.
