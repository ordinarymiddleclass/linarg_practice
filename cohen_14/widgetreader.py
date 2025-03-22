"""
    csvwidget.py
    ~~~~~~~~~~~~
    code for reading widget_data.csv and do lstsq, for use in Cohen's linear algebra code challenge, chapter 14.
"""
import csv
import numpy as np 
import matplotlib.pyplot as plt
#import 3d plot module
from mpl_toolkits.mplot3d import Axes3D

def read_widget_data(filename='widget_data.csv'):   
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        #read and convert to list of tuples of floats
        widget_data = [(float(row[0]), float(row[1]), float(row[2])) for row in reader]
    return widget_data

def widget_data_to_coeff_matrix(widget_data):
    #convert the first two columns of widget_data to a matrix
    A = np.array([[row[0], row[1], 1] for row in widget_data])
    #convert the third column of widget_data to a vector 
    b = np.array([[row[2]] for row in widget_data])
    return A, b

def read_to_coeff_matrix(filename='widget_data.csv'):
    widget_data = read_widget_data(filename)
    return widget_data_to_coeff_matrix(widget_data)

def solve_widget_data(filename='widget_data.csv'):
    widget_data = read_widget_data(filename)
    A, b = widget_data_to_coeff_matrix(widget_data)
    x, residuals, rank, s = np.linalg.lstsq(A, b)
    variance = np.var(b) * len(b)
    calculated_r_squared = 1 - residuals[0] / variance
    return x, calculated_r_squared

def plot_widget_data(filename='widget_data.csv'):
    widget_data = read_widget_data(filename)
    #convert the list of tuples to a list of lists
    widget_data = [[row[0], row[1], row[2]] for row in widget_data]
    #convert the list of lists to a matrix
    widget_data = np.array(widget_data)
    #plot the data
    #create a 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plot the data
    ax.scatter(widget_data[:,0], widget_data[:,1], widget_data[:,2])
    #label the axes
    ax.set_xlabel('time')
    ax.set_ylabel('age')
    ax.set_zlabel('widgets sold')
    #show the plot
    plt.show()

if __name__ == "__main__":    
    print("running the solve_widget_data function")
    x, calculated_r_squared = solve_widget_data()
    print(f"The model is: y = {x[0]} * time + {x[1]} * age + {x[2]}")    
    print(f"R square = {calculated_r_squared}")
    print("running the plot_widget_data function")
    plot_widget_data()
