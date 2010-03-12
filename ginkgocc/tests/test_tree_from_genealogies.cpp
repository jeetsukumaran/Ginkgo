///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
//
// Copyright 2009 Jeet Sukumaran and Mark T. Holder.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include "../organism.hpp"
#include "../population.hpp"
#include "../species.hpp"
#include "../tree.hpp"
#include <iostream>
#include <string>
#include <sstream>

/*
                   0
                  / \
                 1   2
                / \   \
               3   4   5
              / \    \   \
             6    7   8   9
            / \  / \   \   \
           a   b c d    e   f


          ((((a:1,b:1):1,(c:1,d:1):1):1,e:3):1,f:4):1;

*/

std::string build_tree1() {
    ginkgo::GenealogyNode n0;
    ginkgo::GenealogyNode n1;
    ginkgo::GenealogyNode n2;
    n1.link(&n0);
    n2.link(&n0);
    ginkgo::GenealogyNode n3;
    ginkgo::GenealogyNode n4;
    ginkgo::GenealogyNode n5;
    n3.link(&n1);
    n4.link(&n1);
    n5.link(&n2);
    ginkgo::GenealogyNode n6;
    ginkgo::GenealogyNode n7;
    ginkgo::GenealogyNode n8;
    ginkgo::GenealogyNode n9;
    n6.link(&n3);
    n7.link(&n3);
    n8.link(&n4);
    n9.link(&n5);
    ginkgo::GenealogyNode na;
    ginkgo::GenealogyNode nb;
    ginkgo::GenealogyNode nc;
    ginkgo::GenealogyNode nd;
    ginkgo::GenealogyNode ne;
    ginkgo::GenealogyNode nf;
    na.link(&n6);
    nb.link(&n6);
    nc.link(&n7);
    nd.link(&n7);
    ne.link(&n8);
    nf.link(&n9);

    ginkgo::Tree tree;

    std::string label_a("a");
    tree.add_leaf(&na, &label_a);
    std::string label_b("b");
    tree.add_leaf(&nb, &label_b);
    std::string label_c("c");
    tree.add_leaf(&nc, &label_c);
    std::string label_d("d");
    tree.add_leaf(&nd, &label_d);
    std::string label_e("e");
    tree.add_leaf(&ne, &label_e);
    std::string label_f("f");
    tree.add_leaf(&nf, &label_f);
    std::ostringstream s;
    tree.write_newick_tree(s);
    return s.str();
}

std::string build_tree2() {
    ginkgo::GenealogyNode * n0 = NULL;
    ginkgo::GenealogyNode * n1 = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * n2 = new ginkgo::GenealogyNode();
    n1->link(n0);
    n2->link(n0);
    ginkgo::GenealogyNode * n3 = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * n4 = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * n5 = new ginkgo::GenealogyNode();
    n3->link(n1);
    n4->link(n1);
    n5->link(n2);
    ginkgo::GenealogyNode * n6 = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * n7 = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * n8 = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * n9 = new ginkgo::GenealogyNode();
    n6->link(n3);
    n7->link(n3);
    n8->link(n4);
    n9->link(n5);
    ginkgo::GenealogyNode * na = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * nb = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * nc = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * nd = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * ne = new ginkgo::GenealogyNode();
    ginkgo::GenealogyNode * nf = new ginkgo::GenealogyNode();
    na->link(n6);
    nb->link(n6);
    nc->link(n7);
    nd->link(n7);
    ne->link(n8);
    nf->link(n9);

    ginkgo::Tree tree;

    std::string label_a("a");
    tree.add_leaf(na, &label_a);
    std::string label_b("b");
    tree.add_leaf(nb, &label_b);
    std::string label_c("c");
    tree.add_leaf(nc, &label_c);
    std::string label_d("d");
    tree.add_leaf(nd, &label_d);
    std::string label_e("e");
    tree.add_leaf(ne, &label_e);
    std::string label_f("f");
    tree.add_leaf(nf, &label_f);
    std::ostringstream s;
    tree.write_newick_tree(s);

    delete nf;
    delete ne;
    delete nd;
    delete nc;
    delete nb;
    delete na;
    delete n9;
    delete n8;
    delete n7;
    delete n6;
    delete n5;
    delete n4;
    delete n3;
    delete n2;
    delete n1;

    return s.str();
}

std::string build_tree3() {
    ginkgo::HaploidMarker h0;
    ginkgo::HaploidMarker h1;
    ginkgo::HaploidMarker h2;
    h1.inherit(h0);
    h2.inherit(h0);
    ginkgo::HaploidMarker h3;
    ginkgo::HaploidMarker h4;
    ginkgo::HaploidMarker h5;
    h3.inherit(h1);
    h4.inherit(h1);
    h5.inherit(h2);
    ginkgo::HaploidMarker h6;
    ginkgo::HaploidMarker h7;
    ginkgo::HaploidMarker h8;
    ginkgo::HaploidMarker h9;
    h6.inherit(h3);
    h7.inherit(h3);
    h8.inherit(h4);
    h9.inherit(h5);
    ginkgo::HaploidMarker ha;
    ginkgo::HaploidMarker hb;
    ginkgo::HaploidMarker hc;
    ginkgo::HaploidMarker hd;
    ginkgo::HaploidMarker he;
    ginkgo::HaploidMarker hf;
    ha.inherit(h6);
    hb.inherit(h6);
    hc.inherit(h7);
    hd.inherit(h7);
    he.inherit(h8);
    hf.inherit(h9);

    ginkgo::Tree tree;

    std::string label_a("a");
    tree.add_leaf(ha.node(), &label_a);
    std::string label_b("b");
    tree.add_leaf(hb.node(), &label_b);
    std::string label_c("c");
    tree.add_leaf(hc.node(), &label_c);
    std::string label_d("d");
    tree.add_leaf(hd.node(), &label_d);
    std::string label_e("e");
    tree.add_leaf(he.node(), &label_e);
    std::string label_f("f");
    tree.add_leaf(hf.node(), &label_f);
    std::ostringstream s;
    tree.write_newick_tree(s);
    return s.str();
}

int main(int, char * []) {

    std::string expected1 = "((((a:1, b:1):1, (c:1, d:1):1):1, e:3):1, f:4):1";
    std::string expected2 = "((((a:1, b:1):1, (c:1, d:1):1):1, e:3):1, f:4):9999";

    std::string tree1 = build_tree1();
    if (true) {
        //std::cout << "FAIL" << std::endl;
        std::cerr << "\nExpecting:" << std::endl;
        std::cerr << expected1 << std::endl;
        std::cerr << "Observed:" << std::endl;
        std::cerr << tree1 << std::endl;
        //exit(1);
    }
}
