# get bounds for PID circuit components from Routh-Hurwitz conditions

g= 9.81

# Power supply gain:
gamma=1.0 #20 #[A/V]

# L-coil time constant:
tau_L = 0.01 # seconds

# eddies delay
tau_eddies = 0.01 # seconds

# Convenience:
f_L= 1.0/tau_L
f_eddies = 1.0/tau_eddies

# alpha (from h-linearization)
alpha= 8.0 #11.4

# equilibrium L-coil current
IL0 = 40.0 #20

# laser sensor response
G5 = 20.0 #[V/m]

a0 = 1.0/(g*alpha*f_L*f_eddies)

a1 = (1.0/(f_L*g*alpha) + 1.0/(f_eddies*g*alpha))

a2 = (-1.0/(f_L*f_eddies) + 1.0/(g*alpha))

#a3 = (-1.0/tau_L - 1.0/tau_eddies + gamma*D / (IL0 * alpha / G5))

#a4 = gamma* P / (IL0 * alpha /G5)

cond1a = (IL0* alpha / (G5*gamma)) *(1.0/f_L +1.0/f_eddies)

cond1b = (IL0*alpha)/(G5*gamma)

cond2 = (IL0*alpha)/(G5*gamma) * (1.0/f_L + 1.0/f_eddies + a1*a2/a0)

cond3a = (IL0*alpha/(a1**2 * gamma* G5)) * (a0* (G5*gamma/(IL0*alpha))**2)

cond3b =  (IL0*alpha/(a1**2 * gamma* G5)) * ( (a1*a2*G5*gamma)/(IL0*alpha)+ 2*a0 *(1.0/f_L +1.0/f_eddies)*(gamma*G5/(IL0*alpha)))

cond3c = (IL0*alpha/(a1**2 * gamma* G5)) *( a1**2 - a1*a2*(1.0/f_L +1.0/f_eddies))


#######
print cond1a, " < D < ", cond2
print "P > ", cond1b
print "P < ", cond3a, " D^2 + ", cond3b, " D + ", cond3c

# Max variable resistance (knob)
Rmax=200e3

# Recall: P=R_43 / R_41 
R41= 1000 # Ohm

# Recall: D = R_63 * C_61
C61 = 9.4e-6 # Henry



#rho = P / Rmax
