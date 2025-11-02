
import scipy.odr
import scipy.optimize

import numpy as np

import math 

import scipy
from scipy import optimize

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



def test1(x, z, y):
    return 3 * x + z

def plot_2d(func, x_min, x_max, x_amount, *args):
    plt.style.use('_mpl-gallery')
    i = args.index(None)
    
    x = np.linspace(x_min, x_max, x_amount)
    
    args = list(args)
    args[i] = x
    
    y = func(*args)
    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=2.0)
    ax.set_title("Automatic ticks")
    ax.axis('scaled')
    plt.show()

#plot_2d(test1, 0, 20, 20, 0, 20, 3, None, 4)


def plot_3d(func, x_min, x_max, x_amount, y_min, y_max, y_amount, *args):
    plt.style.use('_mpl-gallery')
    
    i_first = args.index(None)
    try:
        i_second = args.index(None, i_first + 1)
    except:
        i_second = None
    print(i_second)
    # Generate grid of x and y
    x = np.linspace(x_min, x_max, x_amount)
    y = np.linspace(y_min, y_max, y_amount)
    
    args = list(args)
    
    X, Y = np.meshgrid(x, y)
    
    args[i_first] = X
    #if i_second
    if [i_second] != None:
        args[i_second] = Y
    # Apply the function to calculate z-values
    Z = func(*args)
    
    # Plotting
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
    
    # Customize labels
    ax.set_title("Automatic ticks")
    
    plt.show()

def test2(x, y):
    return np.sin(x + y)

#plot_3d(test2, 0, 300, 50, 0, 20, 50, None, None)



