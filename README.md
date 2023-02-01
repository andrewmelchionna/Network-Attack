# Network-Attack
This code builds random networks, and then tries to disconnect the network by removing nodes.
The nodes are chosen for removal based on different metrics of importance, and are calculated with noisy information about the network.

# NetworkAttack
This code builds large random networks (Erdos-Renyi and scale free), then creates a noisy copy by removing or adding some edges. 
Nodes are ranked according to different centrality measures, including degree, betweenness (number of geodesics on which the given node lies), and dynamical importance
(change in largest eigenvalue of the adjacency matrix upon removing the given node). These measures are calculated using the noisy copy of the network. We then remove these
nodes in order of a given centrality measure, and watch how the size of the giant connected component decreases.

We found that for Erdos-Renyi networks, attack effectiveness is more robust to the presence of false edges as opposed to missing edges, and that for small amounts of noise,
betweenness generally out-performs the other two measures, but for moderate amounts of noise, dynamical importance and betweenness are comparably effective.
