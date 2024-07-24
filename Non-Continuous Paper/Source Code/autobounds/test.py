from DAG import DAG
from causalProblem import causalProblem

dag = DAG()
dag.from_structure(edges="U1 -> X1, U1 -> X2,"
                         "U2 -> X2, U2 -> X3,"
                         "U3 -> X2, U3 -> X4,"
                         "U4 -> X3, U4 -> X4,"
                         "X1 -> X5,"
                         "X3 -> X2,"
                         "X4 -> X3, X4 -> X5,"
                         "X5 -> X2, X5 -> X3",
                         unob="U1, U2, U3, U4",)
# Algo 1, step 2 (App. B.4.2): initialize program
# this automatically constructs generalized principal strata parameters
problem = causalProblem(dag) # Z, X, Y binary by default
# Algo 1, step 3 (App. B.4.3): state & automatically polynomialize estimand
problem.set_ate("X", "Y")
# Algo 1, step 4 (App. B.4.4): add & automatically polynomialize constraints
# this includes both empiricial evidence and probability axiom constraints
problem.load_data("/home/lisa/Documents/TUM/Semester6/thesis/non_continous/code/autobounds-main/autobounds/test.csv", optimize = False) # no Sec. 5 simplifications
problem.add_prob_constraints()
# Algo 2 (Sec. 6): compile program, implement primal-dual optimization
program = problem.write_program()
program.run_couenne() # epsilon = 0.01 by default