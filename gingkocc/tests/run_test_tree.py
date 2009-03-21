# /usr/bin/env python

import sys
from dendropy import datasets

def to_parent_array(tree, include_labels=False, include_edge_lens=False):
    pa = []
    edge_lens = []
    node_pa = {}
    for node in tree.preorder_node_iter():
        node_pa[node] = len(pa)
        if node.parent_node is not None:
            if node.taxon is not None:
                pa.append("%s:%s" % (node_pa[node.parent_node], node.taxon.label))
            else:
                pa.append(str(node_pa[node.parent_node]))
        else:
            pa.append(str(-1))
        if include_edge_lens:
            assert node.edge is not None
            edge_lens.append(node.edge.length)
    return pa, edge_lens

d = datasets.Dataset()
tree = d.trees_from_string("(a:2, (b:1, c:2):1):1", format="newick")[0]
pa, edge_lens = to_parent_array(tree, True)
print " ".join([str(i) for i in pa[1:]])