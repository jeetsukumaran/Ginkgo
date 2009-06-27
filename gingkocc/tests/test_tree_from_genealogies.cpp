///////////////////////////////////////////////////////////////////////////////
//
// GINGKO Biogeographical Evolution Simulator.
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

#include "../biosys.hpp"
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
    gingko::GenealogyNode n0;
    gingko::GenealogyNode n1;
    gingko::GenealogyNode n2;   
    n1.link(&n0);
    n2.link(&n0);
    gingko::GenealogyNode n3;
    gingko::GenealogyNode n4;
    gingko::GenealogyNode n5;
    n3.link(&n1);
    n4.link(&n1);
    n5.link(&n2);
    gingko::GenealogyNode n6;
    gingko::GenealogyNode n7;
    gingko::GenealogyNode n8;
    gingko::GenealogyNode n9;
    n6.link(&n3);
    n7.link(&n3);
    n8.link(&n4);
    n9.link(&n5);
    gingko::GenealogyNode na;
    gingko::GenealogyNode nb;
    gingko::GenealogyNode nc;
    gingko::GenealogyNode nd;
    gingko::GenealogyNode ne;
    gingko::GenealogyNode nf;    
    na.link(&n6);
    nb.link(&n6);
    nc.link(&n7);
    nd.link(&n7);
    ne.link(&n8);
    nf.link(&n9);
    
    gingko::Tree tree;
    
    std::string label_a("a");
    tree.add_node(&na, &label_a);    
    std::string label_b("b");
    tree.add_node(&nb, &label_b);
    std::string label_c("c");
    tree.add_node(&nc, &label_c);
    std::string label_d("d");
    tree.add_node(&nd, &label_d);
    std::string label_e("e");
    tree.add_node(&ne, &label_e);
    std::string label_f("f");
    tree.add_node(&nf, &label_f);
    std::ostringstream s;
    tree.write_newick_tree(s);
    return s.str();
}

std::string build_tree2() {
    gingko::GenealogyNode * n0 = NULL;
    gingko::GenealogyNode * n1 = new gingko::GenealogyNode();
    gingko::GenealogyNode * n2 = new gingko::GenealogyNode();   
    n1->link(n0);
    n2->link(n0);
    gingko::GenealogyNode * n3 = new gingko::GenealogyNode();
    gingko::GenealogyNode * n4 = new gingko::GenealogyNode();
    gingko::GenealogyNode * n5 = new gingko::GenealogyNode();
    n3->link(n1);
    n4->link(n1);
    n5->link(n2);
    gingko::GenealogyNode * n6 = new gingko::GenealogyNode();
    gingko::GenealogyNode * n7 = new gingko::GenealogyNode();
    gingko::GenealogyNode * n8 = new gingko::GenealogyNode();
    gingko::GenealogyNode * n9 = new gingko::GenealogyNode();
    n6->link(n3);
    n7->link(n3);
    n8->link(n4);
    n9->link(n5);
    gingko::GenealogyNode * na = new gingko::GenealogyNode();
    gingko::GenealogyNode * nb = new gingko::GenealogyNode();
    gingko::GenealogyNode * nc = new gingko::GenealogyNode();
    gingko::GenealogyNode * nd = new gingko::GenealogyNode();
    gingko::GenealogyNode * ne = new gingko::GenealogyNode();
    gingko::GenealogyNode * nf = new gingko::GenealogyNode();    
    na->link(n6);
    nb->link(n6);
    nc->link(n7);
    nd->link(n7);
    ne->link(n8);
    nf->link(n9);
    
    gingko::Tree tree;
    
    std::string label_a("a");
    tree.add_node(na, &label_a);    
    std::string label_b("b");
    tree.add_node(nb, &label_b);
    std::string label_c("c");
    tree.add_node(nc, &label_c);
    std::string label_d("d");
    tree.add_node(nd, &label_d);
    std::string label_e("e");
    tree.add_node(ne, &label_e);
    std::string label_f("f");
    tree.add_node(nf, &label_f);
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
    gingko::HaploidMarker h0;
    gingko::HaploidMarker h1;
    gingko::HaploidMarker h2;   
    h1.inherit(h0);
    h2.inherit(h0);
    gingko::HaploidMarker h3;
    gingko::HaploidMarker h4;
    gingko::HaploidMarker h5;
    h3.inherit(h1);
    h4.inherit(h1);
    h5.inherit(h2);
    gingko::HaploidMarker h6;
    gingko::HaploidMarker h7;
    gingko::HaploidMarker h8;
    gingko::HaploidMarker h9;
    h6.inherit(h3);
    h7.inherit(h3);
    h8.inherit(h4);
    h9.inherit(h5);
    gingko::HaploidMarker ha;
    gingko::HaploidMarker hb;
    gingko::HaploidMarker hc;
    gingko::HaploidMarker hd;
    gingko::HaploidMarker he;
    gingko::HaploidMarker hf;    
    ha.inherit(h6);
    hb.inherit(h6);
    hc.inherit(h7);
    hd.inherit(h7);
    he.inherit(h8);
    hf.inherit(h9);
    
    gingko::Tree tree;
    
    std::string label_a("a");
    tree.add_node(ha.node(), &label_a);    
    std::string label_b("b");
    tree.add_node(hb.node(), &label_b);
    std::string label_c("c");
    tree.add_node(hc.node(), &label_c);
    std::string label_d("d");
    tree.add_node(hd.node(), &label_d);
    std::string label_e("e");
    tree.add_node(he.node(), &label_e);
    std::string label_f("f");
    tree.add_node(hf.node(), &label_f);
    std::ostringstream s;
    tree.write_newick_tree(s);
    return s.str();
}

int main(int, char *) {

    std::string expected1 = "((((a:1, b:1):1, (c:1, d:1):1):1, e:3):1, f:4):1";
    std::string expected2 = "((((a:1, b:1):1, (c:1, d:1):1):1, e:3):1, f:4):9999";

    std::string tree1 = build_tree1();
    if (tree1 != expected1) {
        //std::cout << "FAIL" << std::endl;
        std::cerr << "\nExpecting:" << std::endl;
        std::cerr << expected1 << std::endl;
        std::cerr << "Observed:" << std::endl;
        std::cerr << tree1 << std::endl;
        //exit(1);
    }
    
    std::string tree2 = build_tree2();
    if (tree2 != expected2) {
        //std::cout << "FAIL" << std::endl;
        std::cerr << "\nExpecting:" << std::endl;
        std::cerr << expected2 << std::endl;
        std::cerr << "Observed:" << std::endl;
        std::cerr << tree2 << std::endl;
        //exit(1);
    }    
   
    std::string tree3 = build_tree3();
    if (tree3 != expected2) {
        //std::cout << "FAIL" << std::endl;
        std::cerr << "\nExpecting:" << std::endl;
        std::cerr << expected2 << std::endl;
        std::cerr << "Observed:" << std::endl;
        std::cerr << tree3 << std::endl;
        //exit(1);
    }    
    std::cout << "SUCCESS" << std::endl;
}
