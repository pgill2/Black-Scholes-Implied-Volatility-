# implementing black scholes model. Assumptions: market consists of at least one risky asset (i.e. stock)
# and one (essentially) risk free asset (i.e. money market fund, cash, or governmnent bond)

# Assumes three properties of the two assets, and four of the market itself:

# Assumptions about the asset in the market include 
#       1. The rate of return on the risk free asset is constant meaning you can think of it as an interest rate
#       2. The instantaneous log return of the risky asset's price is assumed to behave as an infinitesimal random walk with constant drift and volatility, 
#       3. Risky asset does not pay a dividend

# you are welcome to plug in any data of your own (i.e. SPY CALLS)



# imports 

from math import *

# add whatever data


sigma = 3

T = 0.5
S = 187.85
r = 0.05
K = 80
# Calculations for the solution to BSM equation
dplus = (1/(sigma*sqrt(T)))*((log(S/K))+(r+(sigma**2)/2)*T)
dminus = (1/(sigma*sqrt(T)))*((log(S/K))+(r-(sigma**2)/2)*T)

# Calculating price of Call and Put
Call = S*norm.cdf(dplus) - K*exp(-r*T)*norm.cdf(dminus)
Put = K*exp(-r*T)*norm.cdf(-dminus)-S*norm.cdf(-dplus)

# Printing the values of call an put options
print("The Price of the Call option is %s" % round(Call, 2))
print("The Price of the Put option is %s" % round(Put,  2))




# Calculate Implied Volatlity with the Bisection Method (Bisection algorithm)


# Some stuff to think about when looking into the Bisection Algorithm


# c(lower, upper) the lower and upper bounds of the implied volatility range to search 
# underlying asset price
# strike price
# risk free rate
# continuous dividend yield rate for options on stocks and stock indicies paying a dividend, also the foreign risk free rate for options on currencies 
# time to maturity 
# market price 
# type of option you will be using: call or put
# tolerance used for stopping criteria 
# max number of iterations








from math import *
from scipy.stats import norm
import openpyxl

# Opening the data file


# my own data
# plug in your spreadsheet 
# 
# 
#  
wb = openpyxl.load_workbook('spy-options-exp-2023-01-27-weekly-near-the-money-stacked-12-15-2022.xlsx')
sheet1 = wb.get_sheet_by_name('SPY')

# Initializing arrays for variables
K = [0]*500  # Array of strike prices
BS = [0]*500  # Array of BSM prices
EU = [0]*500
AM = [0]*500


# Initializing constants:

T = (24/365)  # Time to maturity
S = 214.95  # Underlying asset price
r = 0.4  # interest rate
tol = 10**-6  # tolerance
A = 107  # starting row in excel
sigma = 0.3 # Implied volatility   Call
N = 100
option = 'Put'

# Calculating option value using BSM:
def CBS(S, K, T, r, sigma, option):
    # Calculations for the solution to BSM equation
    dplus = (1 / (sigma * sqrt(T))) * ((log(S / K)) + (r + (sigma ** 2) / 2) * T)
    dminus = (1 / (sigma * sqrt(T))) * ((log(S / K)) + (r - (sigma ** 2) / 2) * T)

    # Calculating price of Call and Put
    if option == 'Call':
        return S * norm.cdf(dplus) - K * exp(-r * T) * norm.cdf(dminus)
    elif option == 'Put':
        return K * exp(-r * T) * norm.cdf(-dminus) - S * norm.cdf(-dplus)

# Calculating option value using binomial for european:
def euro(S, r, sigma, K, N, option):
            # Compute more constants
            v = r - 0.5 * sigma ** 2
            deltaT = T / N
            disc = exp(-r * deltaT)
            dxu = sqrt((sigma ** 2) * deltaT + (v * deltaT) ** 2)
            dxd = -dxu
            #  p = (exp(r*deltaT)-d)/(u-d) # probability of increasing
        
            pu = 0.5 + 0.5 * ((v * deltaT) / dxu)
            pd = 1 - pu

            # Initialize arrays
            St = [0] * (N + 1)
            C = [0] * (N + 1)

            # Initialize asset prices at maturity N
            St[0] = S * exp(N * dxd)
            for j in range(1, N):
                St[j] = St[j - 1] * exp(dxu - dxd)

            # Initialize option values at maturity
            for j in range(0, N):
                if option == 'Call':
                    C[j] = max(0.0, St[j] - K)
                elif option == 'Put':
                    C[j] = max(0.0, K - St[j])

            # Step back through the tree
            for i in range(N - 1, 0):
                for j in range(0, i):
                    C[j] = disc * (pu * C[j + 1] + pd * C[j])
            if option == 'Call':
                return C[51]
            elif option == 'Put':
                return C[0]

def american(S, r, sigma, T, K, N, option):
    dt = T / N
    u = exp(sigma * sqrt(dt))
    d = 1 / u
    p = (exp(r * dt) - d) / (u - d)
    disc = exp(-r * dt)

    St = [0] * (N + 1)
    C = [0] * (N + 1)

    St[0] = S * d ** N
    for j in range(1, N):
        St[j] = St[j - 1] * (u / d)

    for j in range(0, N):
        if option == 'Put':
            C[j] = max(0.0, K - St[j])
        elif option == 'Call':
            C[j] = max(0.0, St[j] - K)
        else:
            break

    for i in range(N - 1, 0):
        for j in range(0, i):
            C[j] = disc * (p * C[j + 1] + (1 - p) * C[j])
            St[j] = (St[j]) / d
            if option == 'Put':
                C[j] = max(C[j], K - St[j])
            elif option == 'Call':
                C[j] = max(C[j], St[j] - K)
            next(j)
        next(i)
    if option == 'Call':
        return C[52]
    elif option == 'Put':
        return C[0]


## Calculating option value using Binomial Theorem:

# Loop for finding values for each row
while A < 211:
    A=A+1

    K[A] = float(sheet1.cell(row=(A), column=1).value)
    BS[A]= CBS(S,K[A],T,r,sigma,option)
    EU[A]= euro(S, r, sigma, K[A], N, option)
    AM[A]= american(S, r, sigma, T, K[A], N, option)


print(BS)
print(EU)
print(AM)





