import operator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import copy
import string
import numpy as np
import igraph
import scipy as scipy
import scipy.stats
n = 100
mink = 1
maxk = n-1
gamma = -2.5
seed_count = 1
def heatmaps_deterministic(results):
  groups = ["final", "rzero_avg", "rzero_initial"]
  titles = {"final": "Total Infected (Deterministic Contagion)", "rzero_avg": "Rolling R0 (Deterministic Contagion)", "rzero_initial": "Initial R0 (Deterministic Contagion)"}
  clabels = {"final": "Total Infected %", "rzero_avg": "Rolling R0", "rzero_initial": "Initial R0"}
  for group in groups:
    print group
    sorted_set = []
    if "avg" in group or "initial" in group:
      for r in flatten_results(results):
        if r != 0:
          sorted_set.append(r[group])
      sorted_set = sorted(sorted_set)
    else:
      for r in flatten_results(results):
        if r != 0:
          sorted_set.append(r[group]/100.0)
      sorted_set = sorted(sorted_set)
    mapped_range = np.arange(0.01,1,0.01).tolist()
    heatmap = np.zeros((len(mapped_range), len(mapped_range))).tolist()
    betas = mapped_range
    mus = mapped_range
    for beta in betas:
      for mu in mus:
        if "avg" in group or "initial" in group:
          if results[betas.index(beta)][mus.index(mu)] != 0:
            heatmap[betas.index(beta)][mus.index(mu)] = results[betas.index(beta)][mus.index(mu)][group]
        else:
          if results[betas.index(beta)][mus.index(mu)] != 0:
            heatmap[betas.index(beta)][mus.index(mu)] = results[betas.index(beta)][mus.index(mu)][group]/100.0
    plt.clf()
    plt.imshow(np.matrix(heatmap), origin='lower', interpolation='none', vmax=1, extent=[0, 1, 0, 1])
    plt.ylabel(r'$\beta$')
    plt.xlabel(r'$\mu$')
    cb = plt.colorbar()
    cb.set_label(clabels[group])
    plt.clim(sorted_set[0],sorted_set[-1])
    plt.title(titles[group])
    plt.savefig(group)

def run_many_deterministic():
  mapped_range = np.arange(0.01,1,0.01).tolist()
  results = np.zeros((len(mapped_range), len(mapped_range))).tolist()
  betas = mapped_range
  mus = mapped_range
  for beta in betas:
    for mu in mus:
      print beta
      print mu
      final, rzero_avg, rzero_initial = sir_deterministic(n, beta, mu)
      results[betas.index(beta)][mus.index(mu)] = {"beta": beta, "mu": mu, "final": final, "rzero_avg": rzero_avg, "rzero_initial": rzero_initial}
    heatmaps_deterministic(results)
  return results
  
def s_to_i(current_status, beta):
  return (beta*current_status['susceptible']*float(current_status['infected']))/sum(current_status.values())

def i_to_r(current_status, mu):
  return (mu*current_status['infected'])

def sir_deterministic(n, beta, mu):
  statuses = [{'susceptible': n-1, 'infected': 1, 'recovered': 0}]
  while statuses[-1]['infected'] > 0.0001:
    s_t_i = s_to_i(statuses[-1], beta)
    i_t_r = i_to_r(statuses[-1], mu)
    statuses.append({'susceptible': statuses[-1]['susceptible']-s_t_i, 'infected': statuses[-1]['infected']+s_t_i-i_t_r, 'recovered': statuses[-1]['recovered']+i_t_r})
  return [statuses[-1]['recovered'], r_zeros([s['infected'] for s in statuses]), statuses[1]['infected']/1]

def edges_to_graph():
  graph = igraph.Graph.Read_Ncol('workfile')
  graph.vs['status'] = 'susceptible'
  return graph

def generate_graph(n, gamma, mink, maxk, net_type):
  if net_type == "barbell":
    first = generate_edges(n/2, gamma, mink, maxk/2, leading_char="A")
    second = generate_edges(n/2, gamma, mink, maxk/2, leading_char="B")
    for i in range(3):
      n1 = "A"+str(int((n/2)*np.random.random()))
      n2 = "B"+str(int((n/2)*np.random.random()))
      first.append([n1, n2])
    for edge in second:
      first.append(edge)
    f = open('workfile', 'w')
    f.write(str.join("\n", [str.join(" ", edge) for edge in first]))
    f.close()
  else:
    first = generate_edges(n, gamma, mink, maxk, leading_char="A")
    f = open('workfile', 'w')
    f.write(str.join("\n", [str.join(" ", edge) for edge in first]))
    f.close()
  return edges_to_graph()

def generate_edges(n, gamma, mink, maxk, leading_char="A"):
  edges = []
  degree_dist = get_degree_dist(n, gamma, mink, maxk)
  nodes = range(n)
  for i,degree in enumerate(degree_dist):
    first = i
    for d in range(degree):
      second = int(len(nodes)*np.random.random())
      while first == second:
        second = int(len(nodes)*np.random.random())
      edges.append([leading_char+str(first), leading_char+str(second)])
  return edges

def get_degree_dist(n, power, mink, maxk):
  return np.array(((maxk**(power+1)-mink**(power+1) )*scipy.random.random(n)+mink**(power+1))**(1/(power + 1))).astype(int).tolist()

def r_zeros(infected_arr):
  rzeros = []
  for i in range(len(infected_arr)-1):
    if infected_arr[i] != 0:
      rzeros.append(infected_arr[i+1]/infected_arr[i])
  if len(rzeros) == 0:
    return 0
  else:
    return np.mean(rzeros)

def run_both(n, mink, maxk, gamma, beta, mu, seed_count):
  graph = generate_graph(n, gamma, mink, maxk, "basic")
  graph_barbell = generate_graph(n, gamma, mink, maxk, "barbell")
  seeds = []
  for i in range(seed_count+1):
    seed = int(np.random.random()*n)
    if seed not in seeds:
      graph.vs(seed)[0]['status'] = 'infected'
      seeds.append(seed)
  initial_state = graph.vs['status']
  full_result = []
  for i in range(10):    
    graph.vs['status'] = initial_state
    final_simple, rzero_avg_simple, rzero_initial_simple = run_simple(graph, seeds, beta, mu)
    graph.vs['status'] = initial_state
    final_complex, rzero_avg_complex, rzero_initial_complex = run_complex(graph, seeds, beta, mu)
    graph_barbell.vs['status'] = initial_state
    final_simple_b, rzero_avg_simple_b, rzero_initial_simple_b = run_simple(graph_barbell, seeds, beta, mu)
    graph_barbell.vs['status'] = initial_state
    final_complex_b, rzero_avg_complex_b, rzero_initial_complex_b = run_complex(graph_barbell, seeds, beta, mu)
    full_result.append([final_simple, final_complex, r_zeros(rzero_avg_simple), r_zeros(rzero_avg_complex), final_simple_b, final_complex_b, r_zeros(rzero_avg_simple_b), r_zeros(rzero_avg_complex_b), rzero_initial_simple, rzero_initial_complex, rzero_initial_simple_b, rzero_initial_complex_b])
  return [np.mean(el) for el in np.matrix(full_result).transpose().tolist()]

def run_simple(graph, seeds, beta, mu):
  infected_arr = [graph.vs['status'].count('infected')]
  while 'infected' in graph.vs['status']:
    count_infected = 0
    next_statuses = []
    for v in graph.vs():
      if v['status'] == 'susceptible':
        neighbor_statuses = [graph.vs()[n]['status'] for n in graph.neighbors(v)]
        if np.random.random() < (1-(1-beta)**int(neighbor_statuses.count('infected'))):
          next_statuses.append('infected')
          count_infected += 1
        else:
          next_statuses.append('susceptible')
      elif v['status'] == 'infected':
        if np.random.random() < mu:
          next_statuses.append('recovered')
        else:
          next_statuses.append('infected')
      elif v['status'] == 'recovered':
        next_statuses.append('recovered')
    infected_arr.append(count_infected)
    graph.vs['status'] = next_statuses
  initial_r0 = 0
  if len(infected_arr) > 1:
    initial_r0 = infected_arr[1]/float(infected_arr[0])
  return graph.vs['status'].count("recovered"), infected_arr, initial_r0

def run_complex(graph, seeds, beta, mu):
  infected_arr = [graph.vs['status'].count('infected')]
  while 'infected' in graph.vs['status']:
    count_infected = 0
    next_statuses = []
    for v in graph.vs():
      if v['status'] == 'susceptible':
        neighbor_statuses = [graph.vs()[n]['status'] for n in graph.neighbors(v)]
        if neighbor_statuses.count('infected')/float(len(neighbor_statuses)) > beta or len(neighbor_statuses) == neighbor_statuses.count('infected'):
          next_statuses.append('infected')
          count_infected += 1
        else:
          next_statuses.append('susceptible')
      elif v['status'] == 'infected':
        if np.random.random() < mu:
          next_statuses.append('recovered')
        else:
          next_statuses.append('infected')
      elif v['status'] == 'recovered':
        next_statuses.append('recovered')
    infected_arr.append(count_infected)
    graph.vs['status'] = next_statuses
  initial_r0 = 0
  if len(infected_arr) > 1:
    initial_r0 = infected_arr[1]/float(infected_arr[0])
  return graph.vs['status'].count("recovered"), infected_arr, initial_r0

def flatten_results(results):
  return [item for sublist in results for item in sublist]

def heatmaps(results):
  groups = ["final_simple", "final_complex", "final_simple_b", "final_complex_b", "rzero_avg_simple", "rzero_avg_complex", "rzero_avg_simple_b", "rzero_avg_complex_b", "rzero_initial_simple", "rzero_initial_complex", "rzero_initial_simple_b", "rzero_initial_complex_b"]
  titles = {"final_simple": "Total Infected (Simple Contagion)", "final_complex": "Total Infected  (Complex Contagion)", "final_simple_b": "Total Infected (Barbell Simple Contagion)", "final_complex_b": "Total Infected (Barbell Complex Contagion)", "rzero_avg_simple": "Rolling R0 (Simple Contagion)", "rzero_avg_complex": "Rolling R0 (Complex Contagion)", "rzero_avg_simple_b": "Rolling R0 (Barbell Simple Contagion)", "rzero_avg_complex_b": "Rolling R0 (Barbell Complex Contagion)", "rzero_initial_simple": "Initial R0 (Simple Contagion)", "rzero_initial_complex": "Initial R0 (Complex Contagion)", "rzero_initial_simple_b": "Initial R0 (Barbell Simple Contagion)", "rzero_initial_complex_b": "Initial R0 (Barbell Complex Contagion)"}
  clabels = {"final_simple": "Total Infected %", "final_complex": "Total Infected %", "final_simple_b": "Total Infected %", "final_complex_b": "Total Infected %", "rzero_avg_simple": "Rolling R0", "rzero_avg_complex": "Rolling R0", "rzero_avg_simple_b": "Rolling R0", "rzero_avg_complex_b": "Rolling R0", "rzero_initial_simple": "Initial R0", "rzero_initial_complex": "Initial R0", "rzero_initial_simple_b": "Initial R0", "rzero_initial_complex_b": "Initial R0"}
  for group in groups:
    print group
    sorted_set = []
    if "avg" in group or "initial" in group:
      for r in flatten_results(results):
        if r != 0:
          sorted_set.append(r[group])
      sorted_set = sorted(sorted_set)
    else:
      for r in flatten_results(results):
        if r != 0:
          sorted_set.append(r[group]/100.0)
      sorted_set = sorted(sorted_set)
    mapped_range = np.arange(0.01,1,0.01).tolist()
    heatmap = np.zeros((len(mapped_range), len(mapped_range))).tolist()
    betas = mapped_range
    mus = mapped_range
    for beta in betas:
      for mu in mus:
        if "avg" in group or "initial" in group:
          if results[betas.index(beta)][mus.index(mu)] != 0:
            heatmap[betas.index(beta)][mus.index(mu)] = results[betas.index(beta)][mus.index(mu)][group]
        else:
          if results[betas.index(beta)][mus.index(mu)] != 0:
            heatmap[betas.index(beta)][mus.index(mu)] = results[betas.index(beta)][mus.index(mu)][group]/100.0
    plt.clf()
    plt.imshow(np.matrix(heatmap), origin='lower', interpolation='none', vmax=1, extent=[0, 1, 0, 1])
    plt.ylabel(r'$\beta$')
    plt.xlabel(r'$\mu$')
    cb = plt.colorbar()
    cb.set_label(clabels[group])
    plt.clim(sorted_set[0],sorted_set[-1])
    plt.title(titles[group])
    plt.savefig(group)

def run_many():
  mapped_range = np.arange(0.01,1,0.01).tolist()
  results = np.zeros((len(mapped_range), len(mapped_range))).tolist()
  betas = mapped_range
  mus = mapped_range
  for beta in betas:
    for mu in mus:
      print beta
      print mu
      final_simple, final_complex, rzero_avg_simple, rzero_avg_complex, final_simple_b, final_complex_b, rzero_avg_simple_b, rzero_avg_complex_b, rzero_initial_simple, rzero_initial_complex, rzero_initial_simple_b, rzero_initial_complex_b = run_both(n, mink, maxk, gamma, beta, mu, seed_count)
      results[betas.index(beta)][mus.index(mu)] = {"beta": beta, "mu": mu, "final_simple": final_simple, "final_complex": final_complex, "rzero_avg_simple": rzero_avg_simple, "rzero_avg_complex": rzero_avg_complex, "final_simple_b": final_simple_b, "final_complex_b": final_complex_b, "rzero_avg_simple_b": rzero_avg_simple_b, "rzero_avg_complex_b": rzero_avg_complex_b, "rzero_initial_simple": rzero_initial_simple, "rzero_initial_complex": rzero_initial_complex, "rzero_initial_simple_b": rzero_initial_simple_b, "rzero_initial_complex_b": rzero_initial_complex_b}
    heatmaps(results)
  return results

