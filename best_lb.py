import argparse
import matplotlib.pyplot as plt
import numpy as np
import os


def best_lb_vs_time(file_path):
    best_lb = []
    time_domain = []
    line_count = 0
    with open(file_path, mode='r') as in_csv:
        for line in in_csv:  # skip first line
            if line_count == 0:
                line_count += 1
                continue
            line = line.split(',')
            best_lb.append(float(line[0]))
            #time_domain.append(float(line[1][:-2]))
            time_domain.append(float(line[1]))
            # if len(time_domain) > 100000:
            #     break
            line_count += 1

    fig, ax = plt.subplots()

    ax.plot(np.array(time_domain), np.array(best_lb))

    ax.set(xlabel='Time [s]', ylabel='Best lower bound',
           title='Simulating-Annealing '+file_path)

    fig.savefig(os.path.join(file_path[:-3] + "png"))
    plt.show()


if __name__ == '__main__':
    # parse input
    parser = argparse.ArgumentParser(description='Process input param')
    parser.add_argument('--file', '-f', type=str, help='input csv file')

    args = parser.parse_args()

    best_lb_vs_time(args.file)

    i = 0.27
    k = 0.27
    j = 2
    #f = open("time.txt", "x")

    while j <= 100000:
      #f.write("%s \n" % k)
      k = i * j
      k = round(k, 2)
      j += 1
    #f.close()
