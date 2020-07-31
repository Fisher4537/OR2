import argparse
import statistics
import random
import numpy as np
import matplotlib.pyplot as plt


def get_testset(file_path, v_shift=20., time_limit=300., model_filter=None):
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
            opt_time = float(line[5])
            best_lb = float(line[6][:-2]) # remove ';\n'

            # filter model
            if model_filter is not None:
                if model_type not in model_filter:
                    continue

            # correct time
            if 0. < opt_time < time_limit - 1.:
                opt_time = opt_time + v_shift
            elif time_limit - 1. <= opt_time <= time_limit + 1.:
                opt_time = opt_time * random.uniform(1., 20.)
            else:  # invalid parameter: opt_time <= 0
                continue
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
                 'nthreads': int(nthreads),
                 'best_lb': best_lb})
    return test_set


def performance_profile(file_path, v_shift=20., time_limit=300., model_filter=None):
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
                        'model_time_mean': model_time_mean
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
    test_set = get_testset(file_path, v_shift=v_shift, time_limit=time_limit, model_filter=model_filter)


    # compute statistics
    for input_file in test_set:
        min_opt_time = float("inf")
        min_best_lb = float("inf")
        for model_type in test_set[input_file]['models']:
            # print('{} {} {}'.format(input_file, model_type, test_set[input_file]['models'][model_type]))
            model_time_mean = statistics.mean([sample['opt_time'] for sample in test_set[input_file]['models'][model_type]['exec_time']])
            model_best_lb_mean = statistics.mean([sample['best_lb'] for sample in test_set[input_file]['models'][model_type]['exec_time']])

            test_set[input_file]['models'][model_type]['model_time_mean'] = model_time_mean
            test_set[input_file]['models'][model_type]['model_best_lb_mean'] = model_best_lb_mean

            # print('{} {} {}'.format(input_file, model_type, model_time_mean))
            min_opt_time = min(model_time_mean, min_opt_time)
            min_best_lb = min(model_best_lb_mean, min_best_lb)

        test_set[input_file]['min_opt_time'] = min_opt_time
        test_set[input_file]['min_best_lb'] = min_best_lb

    # get ratio for each model for each input_file
    models_ratio = {}  # { model_name: [ratio1, ratio2, ..]
    models_best_lb_ratio = {}
    for input_file in test_set:
        for model_type in test_set[input_file]['models']:
            if model_type not in models_ratio:
                models_ratio[model_type] = []
            if model_type not in models_best_lb_ratio:
                models_best_lb_ratio[model_type] = []

            models_ratio[model_type] += [
                test_set[input_file]['models'][model_type]['model_time_mean'] / test_set[input_file]['min_opt_time']]
            models_best_lb_ratio[model_type] += [
                test_set[input_file]['models'][model_type]['model_best_lb_mean'] / test_set[input_file]['min_best_lb']]

    # print models ratios
    for model in models_ratio:
        print('{}: {}'.format(model, models_ratio[model]))
        print('{}: {}'.format(model, models_best_lb_ratio[model]))

    # show performance profile
    plot_pp(models_ratio, domain='time')
    plot_pp(models_best_lb_ratio, domain='lb')


def plot_pp(models_ratio, domain='time', name='res'):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    style = ['-', '--', '-.', 'o', '*', '+', 'x']

    fig, ax = plt.subplots(1)
    
    x_ax = np.arange(0, 4, 0.2)**2 + 0.96 if domain == 'time' else np.arange(0, 2, 0.1)**2 + 0.96
    # print(x_ax)
    for i, model in enumerate(models_ratio):
        y_model = [sum(map(lambda ratio: ratio <= x, models_ratio[model]))/len(models_ratio[model]) for x in x_ax]
        y_ax = np.array(y_model)
        ax.plot(x_ax, y_ax, colors[i % len(colors)]+style[i % len(style)], label=model)

    if domain == 'time':
        ax.legend(loc='lower right', shadow=True, fontsize='x-large')
        ax.set_title('Time measure performance profile')
    elif domain == 'lb':
        ax.legend(loc='lower right', shadow=True, fontsize='x-large')
        ax.set_title('Best_lb measure performance profile')

    plt.xlabel('$ k $')
    plt.ylabel('$ \\% $')
    plt.savefig(name + '_' + domain + '.png')
    plt.show()


if __name__ == '__main__':
    # parse input
    parser = argparse.ArgumentParser(description='Process input param')
    parser.add_argument('--file', '-f', type=str, help='input csv file')
    parser.add_argument('--model_filter', type=str, help='string with the name of the filter')

    args = parser.parse_args()

    # get performance profile
    model_filter = args.model_filter.split(",") if "model_filter" in args else None
    performance_profile(args.file, model_filter=model_filter)
