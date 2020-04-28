import argparse
import statistics
import random
import numpy as np
import matplotlib.pyplot as plt


def performance_profile(file_path, v_shift=20., time_limit=300.):
    """
        test_set = {
            'att48.tsp': {
                'models': {
                    '0': {
                        'exec_time': [
                            {'opt_time': float(opt_time),
                             'randomseed': int(randomseed),
                             'nthreads': int(nthreads)}
                        ],
                        'model_mean': model_mean
                    },
                    '1': { ... }
                },
                'input_size': in_size,
                'min_opt_time': min_opt_time
            },
            'att24.tsp': { ... },
        }

        models_ratio = {
            '0': [t1/t1_min, t2/t2_min, ...],
            '1': [t1/t1_min, t2/t2_min, ...]
        }
    """
    # read input
    test_set = {}
    with open(file_path, mode='r') as in_csv:
        for line in in_csv:  # skip first line
            line = line.split('; ')
            if line[0] == 'input_file':
                continue
            input_file = line[0]
            input_size = line[1]
            model_type = line[2]
            randomseed = line[3]
            nthreads = line[4]
            opt_time = float(line[5][:-2])  # remove ';\n'
            if opt_time < time_limit - 1.:
                opt_time = opt_time + v_shift
            elif time_limit - 1. <= opt_time <= time_limit + 1.:
                opt_time = opt_time * random.uniform(1., 20.)
            # print('{} {} {} {} {} {}'.format(input_file, input_size, model_type, randomseed, nthreads, opt_time))
            if input_file not in test_set:
                test_set[input_file] = {}
            if 'models' not in test_set[input_file]:
                test_set[input_file]['models'] = {}
            if model_type not in test_set[input_file]['models']:
                test_set[input_file]['models'][model_type] = {}
            if 'exec_time' not in test_set[input_file]['models'][model_type]:
                test_set[input_file]['models'][model_type]['exec_time'] = []

            test_set[input_file]['input_size'] = input_size
            test_set[input_file]['models'][model_type]['exec_time'].append(
                {'opt_time': opt_time,
                 'randomseed': int(randomseed),
                 'nthreads': int(nthreads)})

    # compute statistics
    for input_file in test_set:
        min_opt_time = float("inf")
        for model_type in test_set[input_file]['models']:
            # print('{} {} {}'.format(input_file, model_type, test_set[input_file]['models'][model_type]))
            model_mean = statistics.mean([sample['opt_time'] for sample in test_set[input_file]['models'][model_type]['exec_time']])
            test_set[input_file]['models'][model_type]['model_mean'] = model_mean
            # print('{} {} {}'.format(input_file, model_type, model_mean))
            min_opt_time = min(model_mean, min_opt_time)
        test_set[input_file]['min_opt_time'] = min_opt_time

    # get ratio for each model for each input_file
    models_ratio = {}  # { model_name: [ratio1, ratio2, ..]
    for input_file in test_set:
        for model_type in test_set[input_file]['models']:
            if model_type not in models_ratio:
                models_ratio[model_type] = []
            models_ratio[model_type] += [
                test_set[input_file]['models'][model_type]['model_mean'] / test_set[input_file]['min_opt_time']]

    # print models ratios
    for model in models_ratio:
        print('{}: {}'.format(model, models_ratio[model]))

    # show performance profile
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    style = ['-', '--', '-.', 'o']

    fig, ax = plt.subplots()
    x_ax = np.arange(0, 4, 0.2)**2 + 0.96
    # print(x_ax)
    for i, model in enumerate(models_ratio):
        y_model = [sum(map(lambda ratio: ratio <= x, models_ratio[model]))/len(models_ratio[model]) for x in x_ax]
        y_ax = np.array(y_model)
        ax.plot(x_ax, y_ax, colors[i % len(colors)]+style[1], label=model)

    ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    plt.show()


if __name__ == '__main__':
    # parse input
    parser = argparse.ArgumentParser(description='Process input param')
    parser.add_argument('--file', '-f', type=str, help='input csv file')

    args = parser.parse_args()

    # get performance profile
    performance_profile(args.file)