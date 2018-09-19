# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 14:37:07 2018

@author: rstrepparava
"""

import numpy as np

def simulateOUrough(S0, mu, sigma, m_lambda, deltat, t):
    """ Approximate Ornstein-Uhlenbeck Generator. A more accurate version is preferred
         and available : simulateOU
         
         Args:
             S0         (double):   initial price level of the variable     
             mu         (double):   long term mean of OU process
             sigma      (double):   volatility of OU process
             m_lambda   (double):   mean reversion speed of OU process
             deltat     (double):   incremental time interval
             t          (double):   time span from time = 0
                 
         Returns:
             numpy.array:           (Tx1) array of dynamics of the simulated prices   
    """
    periods = int(t // deltat)
    S = np.zeros(periods)
    S[0] = S0
    #Precalculating all dWt's rather than one-per loop makes this function approx 50% faster.
    dWt = np.sqrt(deltat) * np.random.randn(periods) 
    for i in range(1, periods):
        dS = m_lambda*(mu-S[i-1])*deltat + sigma*dWt[i]
        S[i] = S[i-1] + dS
    return S

def simulateOU(S0, mu, sigma, m_lambda, deltat, t):
    """ Simulate an ornstein uhlenbeck process.
        Looks more complicated than expected, because if we don't include the
        exp() terms, we are not accurate as deltat becomes large.
         
         Args:
             S0         (double):   initial price level of the variable     
             mu         (double):   long term mean of OU process
             sigma      (double):   volatility of OU process
             m_lambda   (double):   mean reversion speed of OU process
             deltat     (double):   incremental time interval
             t          (double):   time span from time = 0
                 
         Returns:
             numpy.array:           (Tx1) array of dynamics of the simulated prices   
    """
    periods = int(t // deltat)
    S = np.zeros(periods)
    S[0] = S0
    exp_minus_lambda_deltat = np.exp(-m_lambda*deltat)
    # Calculate the random term. 
    if m_lambda == 0:  # Handle the case of lambda = 0 i.e. no mean reversion. 
        dWt = np.sqrt(deltat) * np.random.randn(periods) 
    else:
        dWt = np.sqrt((1-np.exp(-2*m_lambda*deltat))/(2*m_lambda))*np.random.randn(periods)
    #And iterate through time calculating each price. 
    for i in range(1, periods): 
        S[i] = S[i-1]*exp_minus_lambda_deltat + mu*(1-exp_minus_lambda_deltat) + sigma*dWt[i]
        
    #OPTIM Note : Precalculating all dWt's rather than one-per loop makes this function approx 50% faster. Useful for Monte-Carlo simulations. 
    #OPTIM Note : calculating exp(-lambda*deltat) makes it roughly 50% faster again. 
    #OPTIM Note : this is only about 25% slower than the rough calculation without the exp correction. end
    return S

def simulateABM(S0, mu, sigma, deltat, t):
    """ Simulate  an Arithmetic Brownian Motion (ABM).
         
         Args:
             S0         (double):                           initial price level of the variable     
             mu         (double)/numpy.array of doubles:    drift of the ABM process
             sigma      (double):                           volatility of the ABM process
             deltat     (double):                           incremental time interval
             t          (double):                           time span from time = 0
                 
         Returns:
             numpy.array:                                   (Tx1) array of dynamics of the simulated prices   
    """
    periods = int(t // deltat)
    S = np.zeros(periods)
    S[0] = S0
    #Precalculating all dWt's rather than one-per loop makes this function approx 50% faster.
    dWt = np.sqrt(deltat) * np.random.randn(periods)
    #Broadcast the drift
    if np.isscalar(mu):
        mu = mu*np.ones(periods)
    dS = mu*deltat + sigma*dWt
    for i in range(1, periods):
        S[i] = S[i-1] + dS[i]
    return S

def simulateGBM(S0, mu, sigma, deltat, t):
    """ Simulate  a Geometric Brownian Motion (GBM).
         
         Args:
             S0         (double):                           initial price level of the variable     
             mu         (double)/numpy.array of doubles:    drift of the GBM process
             sigma      (double):                           volatility of the GBM process
             deltat     (double):                           incremental time interval
             t          (double):                           time span from time = 0
                 
         Returns:
             numpy.array:                                   (Tx1) array of dynamics of the simulated prices   
    """
    periods = int(t // deltat)
    S = np.zeros(periods)
    S[0] = S0
    #Precalculating all dWt's rather than one-per loop makes this function approx 50% faster.
    dWt = np.sqrt(deltat) * np.random.randn(periods)
    #Broadcast the drift
    if np.isscalar(mu):
        mu = mu*np.ones(periods)
    dS = np.exp( (mu-sigma*sigma/2)*deltat + sigma*dWt ) 
    for i in range(1, periods):
        S[i] = S[i-1] * dS[i]
    return S


def calibrateOUregress(S, deltat):
    """ Calibrate an OU process by a simple discrete time regression.
        Does not properly take the reversion into account, meaning this will become inaccurate for large deltat.
        
        Args:
            S       (numpy.array):  array of data points (prices) to calibrate OU against 
            deltat  (double):       incremental time interval
        
        Returns:
            double:                 mu, long term mean of OU process
            double:                 sigma, volatility of OU process
            double:                 m_lambda, mean reversion speed of OU process         
    """
    #Regress S(t)-S(t-1) against S(t-1)
    t = len(S)*deltat
    res = np.linalg.lstsq(np.vstack([S[:-1], np.ones(len(S[:-1]))]).T, S[1:]-S[:-1])
    m_lambda = -res[0][0]/deltat
    mu = -res[0][1]/res[0][0]
    sigma = np.sqrt(res[1][0]/t)
    return mu, sigma, m_lambda
    
def calibrateOULS(S, deltat):
    """ Calibrate an OU process by least squares.
        Reference. % Based on the logic described at
        http://sitmo.com/doc/Calibrating_the_Ornstein-Uhlenbeck_model
        
        Args:
            S       (numpy.array):  array of data points (prices) to calibrate OU against 
            deltat  (double):       incremental time interval
        
        Returns:
            double:                 mu, long term mean of OU process
            double:                 sigma, volatility of OU process
            double:                 m_lambda, mean reversion speed of OU process         
    """
    #Regress S(t)-S(t-1) against S(t-1)
    t = len(S)*deltat
    res = np.linalg.lstsq(np.vstack([S[:-1], np.ones(len(S[:-1]))]).T, S[1:])
    m_lambda = -np.log(res[0][0])/deltat
    mu = res[0][1]/(1-res[0][0])
    sigma = np.sqrt(res[1][0]/t * 2*m_lambda*deltat/(1-res[0][0]**2))
    return mu, sigma, m_lambda
    
def calibrateOUMLE(S, deltat):
    """ Calibrate an OU process by maximum likelihood.
        Based on the algorithm and software described at :
        http://www.sitmo.com/doc/Calibrating_the_Ornstein-Uhlenbeck_model
        
        Args:
            S       (numpy.array):  array of data points (prices) to calibrate OU against 
            deltat  (double):       incremental time interval
        
        Returns:
            double:                 mu, long term mean of OU process
            double:                 sigma, volatility of OU process
            double:                 m_lambda, mean reversion speed of OU process         
    """
    n = len(S)-1
    Sx = np.sum(S[:-1])
    Sy = np.sum(S[1:])
    Sxx = np.sum(S[:-1]**2)
    Sxy = np.sum(S[:-1]*S[1:])
    Syy = np.sum(S[1:]**2)
    
    mu = (Sy*Sxx - Sx*Sxy) / (n*(Sxx - Sxy) - (Sx**2 - Sx*Sy))
    m_lambda = -(1/deltat)*np.log((Sxy - mu*Sx - mu*Sy + n*mu**2) / (Sxx -2*mu*Sx + n*mu**2))
    alpha = np.exp(- m_lambda*deltat)
    alpha2 = np.exp(-2*m_lambda*deltat)
    sigmahat2 = (1/n)*(Syy - 2*alpha*Sxy + alpha2*Sxx - 2*mu*(1-alpha)*(Sy - alpha*Sx) + n*mu**2*(1-alpha)**2)
    sigma = np.sqrt(sigmahat2*2*m_lambda/(1-alpha2))
    return mu, sigma, m_lambda

def calibrateOUMLEJackknife(S, deltat, m):
    """ Calibrate an OU process by maximum likelihood.
        Since the basic ML calibration has a bias (resulting in frequent estimates 
        of m_lambda which are much too high), we perform a 'jackknife' operation to reduce the bias.
        Reference.
        To get a less biased lambda, we just the jackknife procedure described in
        Phillips, Peter C. B., and Jun Yu. 2005. 'Jackknifing Bond Option Prices.'
        The Review of Financial Studies 18, no. 2 (Summer): 707-742.
        http://www.jstor.org/stable/3598050
        
        Args:
            S       (numpy.array):  array of data points (prices) to calibrate OU against 
            deltat  (double):       incremental time interval
            m       (int):          number of jackknife partitions
            
        Returns:
            double:                 mu, long term mean of OU process
            double:                 sigma, volatility of OU process
            double:                 m_lambda, mean reversion speed of OU process         
    """
    partlength = len(S) // m 
    Spart = np.zeros((m,partlength))
    for i in range(m): 
        Spart[i,:] = S[partlength*i:partlength*(i+1)]
    
    # Calculate for entire partition. 
    muT, sigmaT, m_lambdaT = calibrateOUMLE(S, deltat)
    # Calculate the individual partitions. 
    mupart          = np.zeros(m) 
    sigmapart       = np.zeros(m)
    m_lambdapart    = np.zeros(m)
    for i in range(m): 
        mupart[i], sigmapart[i], m_lambdapart[i] = calibrateOUMLE(Spart[i,:], deltat)
    # Now the jacknife calculation. 
    m_lambda = (m/(m-1))*m_lambdaT - (np.sum(m_lambdapart))/(m**2-m)
    # mu and sigma are not biased, so there's no real need for the jackknife
    # But we do it anyway for demonstration purposes. 
    mu = (m/(m-1))*muT - (np.sum(mupart))/(m**2-m)
    sigma = (m/(m-1))*sigmaT - (np.sum(sigmapart))/(m**2-m)    
    return mu, sigma, m_lambda

## EOF  