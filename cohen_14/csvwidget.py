"""
    csvwidget.py
    ~~~~~~~~~~~~
    code for reading widget_data.csv and returning a list of tuples, for use in Cohen's linear algebra code challenge, chapter 14.
"""
import csv
def read_widget_data(filename='widget_data.csv'):   
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        #read and convert to list of tuples of floats
        widget_data = [(float(row[0]), float(row[1]), float(row[2])) for row in reader]
    return widget_data

if __name__ == "__main__":
    print(read_widget_data())