import shlex
from subprocess import Popen, PIPE
import numpy as np
import matplotlib.pyplot as plt

cmd_list = ['./test a', './test b', './test c', './test c 1', './test c 2', './test c 3', './test c 4', './test c 5', './test c 6', './test c 7', './test c 8', './test c 9', './test c 10']
test_times = 10

def get_time(cmd):
    process = Popen(shlex.split(cmd), stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
    return float(output)

def get_avg_time(cmd):
    time_list = [get_time(cmd) for i in range(test_times)]
    return sum(time_list)/float(len(time_list))

""" for cmd in cmd_list:
    print('{0}\t{1:.2f}'.format(cmd, get_avg_time(cmd))) """

time_list = [5.40, 2.20, 1.79, 3.01, 2.93, 3.61, 3.39, 3.77, 3.75, 4.53]
base_time = 5.44
group_list = [i+1 for i in range(10)]
plt.figure()
plt.plot(group_list, np.array(time_list)/base_time)
plt.xlabel("groups")
plt.ylabel("used time/base time")
plt.grid()
plt.show()
