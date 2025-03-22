"""
    main.py
    ~~~~~~~
    main program for widget data analysis for Cohen's linear algebra code challenge, chapter 14.
"""
from widgetreader import solve_widget_data, plot_widget_data

if __name__ == "__main__":    
    #print("running the plot_widget_data function")
    #plot_widget_data()
    print("running the solve_widget_data function")
    x, calculated_r_squared = solve_widget_data()
    print(f"The model is: y = {x[0]} * time + {x[1]} * age + {x[2]}")    
    print(f"R square = {calculated_r_squared}")