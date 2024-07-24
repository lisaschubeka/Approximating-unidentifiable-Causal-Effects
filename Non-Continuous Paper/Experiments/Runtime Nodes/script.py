import csv
import time
from autobounds.causalProblem import causalProblem
from matplotlib.ticker import FormatStrFormatter
from autobounds.DAG import DAG
import pandas as pd
import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math

def generate_random_dag(num_nodes):
    dag = nx.DiGraph()

    edge_prob = 0.8 * math.exp(-0.3 * num_nodes)

    dag.add_nodes_from(range(num_nodes))

    possible_edges = [(i, j) for i in range(num_nodes) for j in range(i + 1, num_nodes)]

    random.shuffle(possible_edges)

    for edge in possible_edges:
        if random.random() < edge_prob:
            dag.add_edge(edge[0], edge[1])
            # Check if adding this edge forms a cycle
            if not nx.is_directed_acyclic_graph(dag):
                dag.remove_edge(edge[0], edge[1])

    while not nx.is_weakly_connected(dag):
        weakly_connected_components = list(nx.weakly_connected_components(dag))
        for i in range(len(weakly_connected_components) - 1):

            node_i = random.choice(list(weakly_connected_components[i]))
            j = random.randint(i + 1, len(weakly_connected_components) - 1)
            node_j = random.choice(list(weakly_connected_components[j]))

            if not nx.has_path(dag, node_i, node_j):
                dag.add_edge(node_i, node_j)
                if not nx.is_directed_acyclic_graph(dag):
                    dag.remove_edge(node_i, node_j)
    return dag

def graph_to_string(graph):
    edges_list = []
    for edge in graph.edges():
        edges_list.append(f"{edge[0]} -> {edge[1]}")
    return ", ".join(edges_list)

def define_scm_functions(dag):
    functions = {}
    for node in dag.nodes():
        parents = list(dag.predecessors(node))
        if not parents:
            # No parents, exogenous variable
            functions[node] = lambda n_samples: np.random.binomial(1, 0.5, size=n_samples)
        else:
            # Define a function based on parents
            def func(n_samples, parents=parents):
                data = np.zeros(n_samples, dtype=int)
                for parent in parents:
                    data += np.random.binomial(1, random.random(), size=n_samples)  # Random binary data
                data = (data + np.random.binomial(1, random.random(), size=n_samples)) % 2  # Logical XOR with noise
                return data

            functions[node] = func
    return functions


def generate_data(dag, scm_functions, n_samples):
    data = {}
    for node in nx.topological_sort(dag):
        if not list(dag.predecessors(node)):
            data[node] = scm_functions[node](n_samples)
        else:
            parent_data = np.zeros(n_samples, dtype=int)
            for parent in dag.predecessors(node):
                parent_data += data[parent]
            parent_data = parent_data % 2  # Logical XOR for combining parent data
            data[node] = scm_functions[node](n_samples) ^ parent_data  # Logical XOR with parent data
    return data

def visualize_dag(dag, file_path):
    pos = nx.spring_layout(dag)
    plt.figure(figsize=(8, 6))
    nx.draw(dag, pos, with_labels=True, node_size=700, node_color='skyblue',
            arrowsize=20, font_size=15, font_weight='bold')
    plt.title('Random DAG')
    plt.savefig(file_path)
    plt.close()

def generate_dag_and_csv(num_nodes, file_path):
    n_samples = num_nodes * num_nodes * 1000
    dag = generate_random_dag(num_nodes)
    scm_functions = define_scm_functions(dag)
    data = generate_data(dag, scm_functions, n_samples)

    with open(file_path, 'w') as f:
        writer = csv.writer(f)
        headers = [i for i in data.keys()]
        writer.writerow(headers)
        num_lines = len(data[1])
        for i in range(num_lines):
            row = [data[key][i] for key in sorted(data.keys())]
            writer.writerow(row)
    return dag


def measure_time_for_i_variables(start_node, end_node, file_path, num_iterations):

    x_axis = []
    y_axis = []
    x_axis_exception = []
    y_axis_exception = []
    colour_of_data_point = []

    for i in range(start_node, end_node+1):
        for j in range(num_iterations):
            dag_nx = generate_dag_and_csv(i, file_path)
            input_dag = graph_to_string(dag_nx)
            csv_data = pd.read_csv(file_path)

            visualize_dag(dag_nx, f"/home/lisa/Documents/TUM/Semester6/thesis/non_continous/"
                                  f"code/autobounds-main/Experiments/graphs/end_node_{i}_iter_{j}.png")

            dep_node = None
            for node in dag_nx.nodes():
                if not list(dag_nx.successors(node)):
                    dep_node = node
                    break

            while True:
                ind_node = random.choice(list(dag_nx.nodes()))
                if list(dag_nx.predecessors(ind_node)):
                    break

            unobserved = ",".join([str(node) for node in dag_nx.nodes() if dag_nx.in_degree(node) == 0])
            observed = sorted([str(node) for node in dag_nx.nodes() if dag_nx.in_degree(node) > 0])

            dat = csv_data.loc[:, observed]
            dat = pd.DataFrame(dat.groupby(observed).value_counts().reset_index())
            dat['prob'] = dat['count'] / dat['count'].sum()
            dat = dat.drop(columns='count', axis=0)


            dag = DAG()
            dag.from_structure(edges=input_dag, unob=unobserved)
            problem = causalProblem(dag)
            problem.load_data(dat)
            problem.add_prob_constraints()
            problem.set_ate(ind=ind_node, dep=dep_node)
            prog_ate = problem.write_program()

            prog_ate_optim = prog_ate.run_scip(
            "/home/lisa/Documents/TUM/Semester6/thesis/non_continous/code/autobounds-main/Experiments/results.csv",
                maxtime=1000)

            x_axis.append(i)
            y_axis.append(prog_ate_optim[4])
            colour_of_data_point.append("indigo")
    fig, ax = plt.subplots()
    ax.scatter(x_axis, y_axis, color="indigo",s=20, alpha=0.3, marker='o' )
    ax.scatter(x_axis_exception, y_axis_exception, color="maroon",s=20, alpha=0.3, marker='x' )
    ax.set_xlabel('Number of variables from a randomly generated DAG')
    ax.set_ylabel('Time taken to calculate the causal bounds')
    ax.set_yscale('log')
    ax.set_yticks([0.001, 0.01, 0.1, 1, 10, 100, 1000])
    ax.set_xticks([3,4,5,6])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    plt.tight_layout()
    plt.show()
    plt.savefig("plot_bounds_against_#nodes.png")

measure_time_for_i_variables(3,
                             6,
                             "/Experiments/plot_bounds_against_#nodes.csv",
                             100)