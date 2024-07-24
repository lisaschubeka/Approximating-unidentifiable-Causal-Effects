from autobounds.causalProblem import causalProblem
from autobounds.DAG import DAG

import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# input: graph and U's, ate, data, min_sample, max_sample, step_sample
DAG_graph = "Uz -> Z, Z -> X, X -> Y, Uxy -> X, Uxy -> Y"
unobserved_variables = "Uz, Uxy"
ate = "X,Y"
data = "raw_iv_synth.csv"
min_sample = 100
max_sample = 1000
step_sample = 50
data_points_per_num_samples = 50

def draw_graph(DAG_graph,
               unobserved_variables,
               ate,
               data,
               min_sample,
               max_sample,
               step_sample,
               data_points_per_num_samples=50):

    dag = DAG()
    dag.from_structure(edges=DAG_graph, unob=unobserved_variables)

    x_axis_ate = []
    y_axis_ate = []
    colour_of_data_point = []

    raw_data = pd.read_csv(data)

    pattern = r'\b[A-Z]+\b'
    variables = re.findall(pattern, DAG_graph)

    for i in range(0, int(max_sample/step_sample)):


        num_samples = min_sample + step_sample*i

        if num_samples > max_sample:
            break
        print(f"======== reading {num_samples} rows =======")

        for j in range(0, data_points_per_num_samples):

            sample_size = raw_data.sample(min_sample + step_sample * i)
            dat = pd.DataFrame(sample_size.groupby(['Z', 'X', 'Y']).value_counts().reset_index())
            dat['prob'] = dat['count'] / dat['count'].sum()
            dat = dat.drop(columns='count', axis=0)

            problem = causalProblem(dag)
            problem.load_data(dat)
            problem.add_prob_constraints()
            problem.set_ate(ind=ate.split(",")[0], dep=ate.split(",")[1])
            prog_ate = problem.write_program()
            try:
                prog_ate_optim = prog_ate.run_scip("results.csv")
                x_axis_ate.append(num_samples)
                y_axis_ate.append(np.round(prog_ate_optim[0]['dual'], 3))
                colour_of_data_point.append("lightcoral")

                x_axis_ate.append(num_samples)
                y_axis_ate.append(np.round(prog_ate_optim[1]['dual'], 3))
                colour_of_data_point.append("darkturquoise")
            except:
                shortest_length = min(len(x_axis_ate), len(y_axis_ate), len(colour_of_data_point))
                x_axis_ate = x_axis_ate[:shortest_length]
                y_axis_ate = y_axis_ate[:shortest_length]
                colour_of_data_point = colour_of_data_point[:shortest_length]

    fig, scatter = plt.subplots(nrows=1, ncols=1, figsize=(12,7))

    if len(x_axis_ate) < len(y_axis_ate):
        y_axis_ate = y_axis_ate[0:len(x_axis_ate)]
    elif len(x_axis_ate) > len(y_axis_ate):
        x_axis_ate = x_axis_ate[0:len(y_axis_ate)]

    scatter.scatter(
        x=x_axis_ate,
        y=y_axis_ate,
        s=0.6,
        c=colour_of_data_point,
    )

    scatter.set_xlabel("Number of Samples from Observational Data")
    scatter.set_ylabel("Upper and Lower Bound of the Causal Effect")

    plt.xticks(rotation="vertical")

    plt.show()

if "__main__" == __name__:
    draw_graph(DAG_graph,
               unobserved_variables,
               ate,
               data,
               min_sample,
               max_sample,
               step_sample,
               data_points_per_num_samples)