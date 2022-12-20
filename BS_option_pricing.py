import numpy as np
from scipy.stats import norm

# Interest rate (risk-free)
r = 0.01
# Underlying asset price
S = 41.72
# Strike price
K = 32.0
# Option duration
T = 7.0 / 365.0
# Implied volatility (sigma)
sigma = 0.09

option_type = input("Enter 'c' or 'p' (then return) for a Call/Put:  ")


""" *********************  Function Definitions  ********************* """

# Calculate option price
def blackScholes(r, S, K, T, sigma, type):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma**np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    try:
        if type == "c":
            print("Option Type = Call")
            price = round(S*norm.cdf(d1, 0, 1) - K*np.exp(-r*T)*norm.cdf(d2, 0, 1),3)
        elif type == "p":
            print("Option Type = Put")
            price = round(K*np.exp(-r*T)*norm.cdf(-d2, 0, 1) - S*norm.cdf(-d1, 0, 1),3)
        return price
    except:
        print("Wrong input : Please confirm option type - 'c' for a Call, or 'p' for a Put!")

# Calculate Delta
def delta_calc(r, S, K, T, sigma, type):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma**np.sqrt(T))
    try:
        if type == "c":
            delta_calc = norm.cdf(d1, 0, 1)
        elif type == "p":
            delta_calc = norm.cdf(-d1, 0, 1)
        return delta_calc
    except:
        return

# Calculate Gamma
def gamma_calc(r, S, K, T, sigma, type):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma**np.sqrt(T))
    try:
        if type == "c":
            gamma_calc = norm.pdf(d1, 0, 1) / (S*sigma*np.sqrt(T))
        elif type == "p":
            gamma_calc = norm.pdf(d1, 0, 1) / (S*sigma*np.sqrt(T))
        return gamma_calc
    except:
        return

# Calculate Vega
def vega_calc(r, S, K, T, sigma, type):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma**np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    try:
        if type == "c":
            vega_calc = S*norm.pdf(d1)*np.sqrt(T)
        elif type == "p":
            vega_calc = S*norm.pdf(d1)*np.sqrt(T)
        return vega_calc
    except:
        return

# Calculate Rho
def rho_calc(r, S, K, T, sigma, type):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma**np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    try:
        if type == "c":
            rho_calc = K*T*np.exp(-r*T)*norm.cdf(d2)
        elif type == "p":
            rho_calc = -K*T*np.exp(-r*T)*norm.cdf(-d2)
        return rho_calc
    except:
        return

# Calculate Theta
def theta_calc(r, S, K, T, sigma, type):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T) / (sigma**np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    try:
        if type == "c":
            theta_calc = -S*sigma/2*np.sqrt(T)*norm.cdf(d1, 0, 1) + r*K*np.exp(-r*T)*norm.cdf(-d2)
        elif type == "p":
            theta_calc = -S*sigma/2*np.sqrt(T)*norm.cdf(d1, 0, 1) - r*K*np.exp(-r*T)*norm.cdf(d2)
        return theta_calc
    except:
        return

print("Option Price: ", (blackScholes(r, S, K, T, sigma, option_type)))
print("Delta: ", (delta_calc(r, S, K, T, sigma, option_type)))
print("Gamma: ", (gamma_calc(r, S, K, T, sigma, option_type)))
print("Vega: ", (vega_calc(r, S, K, T, sigma, option_type)))
print("Theta: ", (theta_calc(r, S, K, T, sigma, option_type)))
print("Rho = ", (rho_calc(r, S, K, T, sigma, option_type)))
