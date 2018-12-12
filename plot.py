import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('data.csv', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

fig = plt.figure()
plt.scatter(x, y, marker='o', lw=0, s=(72./fig.dpi)**2)
plt.savefig('HofButterfly.eps', format='eps', dpi=1024)
