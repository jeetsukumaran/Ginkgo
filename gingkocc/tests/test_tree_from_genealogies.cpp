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

#include "../biosys.h"
#include "../tree.h"
#include <iostream>
#include <string>
#include <sstream>


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
    tree.process_node(&na, &label_a);    
    std::string label_b("b");
    tree.process_node(&nb, &label_b);
    std::string label_c("c");
    tree.process_node(&nc, &label_c);
    std::string label_d("d");
    tree.process_node(&nd, &label_d);
    std::string label_e("e");
    tree.process_node(&ne, &label_e);
    std::string label_f("f");
    tree.process_node(&nf, &label_f);
    std::ostringstream s;
    tree.write_newick_tree(s);
    return s.str();
}

int main(int, char *) {

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
    std::string expected = "((((a:1,b:1):1,(c:1,d:1):1):1,e:3):1,f:4):1;";

    std::string tree1 = build_tree1();
    if (expected.compare(tree1) != 0) {
        std::cout << "FAIL" << std::endl;
        std::cerr << "Expecting:" << std::endl;
        std::cerr << expected << std::endl;
        std::cerr << "Observed:" << std::endl;
        std::cerr << tree1 << std::endl;
    }
    
    std::cout << "SUCCESS" << std::endl;
}
