#include "../include/pattern.h"

string encode_el(Pattern& p, VID_MAP& mapping) {
    assert(p.adj.size() > 0);

    //vertex reordering: original vid to new id mapping 
    ////BFS code
    ////For each edge (src id, src label, dest id, dest label, edge label),
    ////Find the lexicographical minimum "BFS" code; using a priority queue for BFS tree string
    ////(code, mapping tables)
    typedef std::tuple<int, int, int> TYPE_EDGE_CODE;
    typedef std::vector< TYPE_EDGE_CODE > TYPE_CODE;

    auto compare_code = [](const TYPE_CODE& lhs, const TYPE_CODE& rhs) {
        int size=min(lhs.size(), rhs.size());
        for(int i=0; i<size; i++) {
            if(lhs[i] != rhs[i]) {
                return lhs[i] < rhs[i];
            }
        }
        return lhs.size() > rhs.size();
    };
    typedef std::pair<int, bool> TYPE_MAPPING_TABLE_ELEMENT;
    typedef std::vector< TYPE_MAPPING_TABLE_ELEMENT > TYPE_MAPPING_TABLE;

    typedef std::pair<TYPE_CODE, TYPE_MAPPING_TABLE> TYPE_CODE_QUEUE_ELEMENT;
    typedef multimap<TYPE_CODE, TYPE_MAPPING_TABLE, decltype(compare_code)> TYPE_CODE_QUEUE;
    typedef TYPE_CODE_QUEUE::iterator TYPE_CODE_QUEUE_ITER;
    int pattern_size=p.getNumEdges();

    TYPE_CODE_QUEUE codes(compare_code);

    //Choose starting vertex
    {
        for (int i = 0; i < p.getNumVertices(); i++) {
            TYPE_MAPPING_TABLE renumber_mapping(p.getNumVertices(), TYPE_MAPPING_TABLE_ELEMENT(-1, false));
            renumber_mapping[i]=TYPE_MAPPING_TABLE_ELEMENT(0, false);

            TYPE_CODE code(1, TYPE_EDGE_CODE(-1,-1,-1));
            codes.insert( TYPE_CODE_QUEUE_ELEMENT(code, renumber_mapping) );
        }
    }
    //maximum length of codes so far: d_max

    auto is_permutable = []( vector< pair<int, int> >& list) {
        auto iter=list.end();
        while(iter != list.begin()) {
            int key=(iter-1)->first;
            auto cur_iter= std::find_if( list.begin(), iter, [key](const std::pair<int, int>& element){ return element.first == key;});
            if (std::next_permutation(cur_iter, iter))
                return true;
            iter=cur_iter;
        }
        return false;
    };

    auto min_mapped_vertex = [](TYPE_MAPPING_TABLE& renumber_mapping) {
        int min_val = MAX_INT, min_val_idx=-1;
        for(TYPE_MAPPING_TABLE::iterator iter=renumber_mapping.begin(); iter != renumber_mapping.end(); iter++) {
            if(iter->first != -1 && min_val > iter->first && !iter->second) {
                min_val = iter->first;
                min_val_idx = iter - renumber_mapping.begin();
            }

        }
        return min_val_idx;
    };

    size_t d_max=0;
    TYPE_CODE min_lex_code(1, TYPE_EDGE_CODE(MAX_INT, MAX_INT, MAX_INT) );
    TYPE_MAPPING_TABLE min_lex_code_mapping;

    while (!codes.empty()) {
        //Pop elements (BFS depth d)
        ////Find the minimum key value
        TYPE_CODE cur_key=codes.begin()->first;

        ////Skip if d < d_max (exists longer but lexicographically smaller pattern)
        if (cur_key.size() < d_max) {
            //Remove and skip
            codes.erase(codes.lower_bound(cur_key), codes.upper_bound(cur_key));
            continue;
        }
        //Update d_max
        d_max = cur_key.size();
        TYPE_CODE_QUEUE new_codes(compare_code);

        for (TYPE_CODE_QUEUE_ITER element = codes.lower_bound(cur_key); element != codes.upper_bound(cur_key); element++) {
            //Find the smallest vertex v_src not visited in the BFS tree
            int v_src=min_mapped_vertex(element->second);
            int starting_renumber_id=(*std::max_element(element->second.begin(), element->second.end(),  [](const TYPE_MAPPING_TABLE_ELEMENT& a, const TYPE_MAPPING_TABLE_ELEMENT& b) {return a.first < b.first;} )).first + 1;
            //If mapping is done!
            if(v_src == -1) {
                if(starting_renumber_id != p.getNumVertices()) {
                    for (int u=0; u<p.getNumVertices(); u++) {
                        if(element->second[ u ].first == -1) {
                            TYPE_MAPPING_TABLE cur_renumber_mapping(element->second);
                            cur_renumber_mapping[ u ].first=starting_renumber_id;
                            cur_renumber_mapping[ u ].second=false;
                            new_codes.insert(TYPE_CODE_QUEUE_ELEMENT(element->first, cur_renumber_mapping));
                        }
                    }
                }
                else if(min_lex_code > element->first) {
                    min_lex_code=element->first;
                    min_lex_code_mapping=element->second;
                }
                continue;
            }
            if(element->first.size() == pattern_size + 1) {
                if(min_lex_code > element->first) {
                    min_lex_code=element->first;
                    min_lex_code_mapping=element->second;
                }
                break;
            }

            int v_src_renumbered=element->second[v_src].first;
            element->second[v_src].second=true;

            //Calculate the updated code
            //Sort all edges by (edge label) order and then by (vertex id) order - already done
            TYPE_CODE cur_code(element->first);
            TYPE_CODE added_code;
            vector<pair<int, int>> cur_adj_permutable;

            int cur_renumber_id = starting_renumber_id;
            for(int u=0; u<p.adj[v_src].size(); u++) {
                //Is visited;
                if (element->second[ p.adj[v_src][u].second ].first != -1) {
                    int u_renumbered=element->second[ p.adj[v_src][u].second ].first; 
                    added_code.push_back( TYPE_EDGE_CODE( v_src_renumbered, u_renumbered, p.adj[v_src][u].first ) );
                }
                else {
                    added_code.push_back( TYPE_EDGE_CODE( v_src_renumbered, cur_renumber_id++, p.adj[v_src][u].first ) );
                    //Capture a subset of the adjacency list to permutate
                    cur_adj_permutable.push_back(p.adj[v_src][u]);
                }
            }
            if(added_code.size() > 0) {
                //added_code need to be sorted
                sort(added_code.begin(), added_code.end());
                cur_code.insert(cur_code.end(), added_code.begin(), added_code.end());
            } 
            do {
                //For tie cases (vertices with the same order); calculate permutation
                TYPE_MAPPING_TABLE cur_renumber_mapping(element->second);
                for (int u=0; u<cur_adj_permutable.size(); u++) {
                    assert(cur_renumber_mapping[ cur_adj_permutable[u].second ].first == -1);
                    cur_renumber_mapping[cur_adj_permutable[u].second].first=starting_renumber_id+u;
                    cur_renumber_mapping[cur_adj_permutable[u].second].second=false;
                }
                new_codes.insert(TYPE_CODE_QUEUE_ELEMENT(cur_code, cur_renumber_mapping));
            } while(is_permutable(cur_adj_permutable));
        }
        //Remove processed codes
        codes.erase(codes.lower_bound(cur_key), codes.upper_bound(cur_key));

        //Update new codes
        codes.insert(new_codes.begin(), new_codes.end());
    }

    mapping.resize(min_lex_code_mapping.size());
    for (int i = 0; i < min_lex_code_mapping.size(); i++)
        mapping[i] = min_lex_code_mapping[i].first;

    string code("");
    //(src id, src label, dest id, dest label, edge label)(src id, ...) in BFS order
    for (int i = 1; i < min_lex_code.size(); i++) {
        code += string("(") + to_string(get<0>(min_lex_code[i])) + string(",") + to_string(get<1>(min_lex_code[i])) + string(",") + to_string(get<2>(min_lex_code[i])) + string(")");
    }

    return code;
}

string encode_vl(Pattern& p, VID_MAP& mapping) {
    typedef std::pair<int, int> TYPE_EDGE_CODE;
    typedef int TYPE_MAPPING_TABLE_ELEMENT;
    typedef std::vector< TYPE_EDGE_CODE > TYPE_CODE;

    auto compare_code = [](const TYPE_CODE& lhs, const TYPE_CODE& rhs) {
        int size=min(lhs.size(), rhs.size());
        for(int i=0; i<size; i++) {
            if(lhs[i] != rhs[i]) {
                return lhs[i] < rhs[i];
            }
        }
        return lhs.size() > rhs.size();
    };
    typedef std::vector< TYPE_MAPPING_TABLE_ELEMENT > TYPE_MAPPING_TABLE;
    typedef std::pair<TYPE_CODE, TYPE_MAPPING_TABLE> TYPE_CODE_QUEUE_ELEMENT;
    typedef multimap<TYPE_CODE, TYPE_MAPPING_TABLE, decltype(compare_code)> TYPE_CODE_QUEUE;
    typedef TYPE_CODE_QUEUE::iterator TYPE_CODE_QUEUE_ITER;
    int pattern_size=p.getNumEdges();

    TYPE_CODE_QUEUE codes(compare_code);

    string code="";

    mapping.clear();
    mapping.resize(p.getNumVertices());

    map<int, vector<int>> order2vid;
    for(int i=0; i<p.getNumVertices(); i++) {
        order2vid[ p.vorder[i] ].push_back(i);
    }

    code+=string("(");
    vector<vector<int>> permutable_vertices;
    for(auto& o2v: order2vid) {
        for(int i=0; i<o2v.second.size(); i++)
            code+=to_string(o2v.first)+string(",");
        permutable_vertices.push_back(o2v.second);
    }
    code.pop_back();
    code+=string(")");

    while(1) {
        //Get a mapping and code
        TYPE_MAPPING_TABLE cur_mapping(p.getNumVertices());
        TYPE_CODE cur_code;
        int renumber_id=0;
        for(auto& uset: permutable_vertices) {
            for(auto& u: uset) {
                cur_mapping[u]=renumber_id++;
            }
        }
        for(auto& uset: permutable_vertices) {
            for(auto& u: uset) {
                vector<int> nbrs;
                for(auto& v: p.adj[u]) {
                    nbrs.push_back( cur_mapping[v.second] );
                }
                sort(nbrs.begin(), nbrs.end());
                for(auto& v: nbrs) {
                    cur_code.push_back( TYPE_EDGE_CODE( cur_mapping[u], v ) );
                }
            }
        }
        //Push it to queue 
        codes.insert( TYPE_CODE_QUEUE_ELEMENT(cur_code, cur_mapping) );

        //Get next permutation
        int cur_level=permutable_vertices.size()-1;
        while( !next_permutation( permutable_vertices[cur_level].begin(), permutable_vertices[cur_level].end()) ) {
            cur_level--;
            if(cur_level == -1) break;
        }
        if(cur_level == -1) break;
        for(int l=cur_level+1; l<permutable_vertices.size(); l++)
            sort( permutable_vertices[l].begin(), permutable_vertices[l].end() );
    }

    assert(codes.size() > 0);
    TYPE_CODE min_code=codes.begin()->first;
    for(int i=0; i<min_code.size(); i++) {
        code += string("(") + to_string( min_code[i].first ) + string(",") + to_string( min_code[i].second ) + string(")");
    }

    mapping=codes.begin()->second;

    return code;
}

string encode_evl(Pattern& p, VID_MAP& mapping) {
    assert(p.adj.size() > 0);

    //vertex reordering: original vid to new id mapping 
    ////BFS code
    ////For each edge (src id, src label, dest id, dest label, edge label),
    ////Find lexicographical minimum "BFS" code; using a priority queue for BFS tree string
    ////(code, mapping tables)
    typedef std::tuple<int, int, int, int, int> TYPE_EDGE_CODE;
    typedef std::vector< TYPE_EDGE_CODE > TYPE_CODE;

    auto compare_code = [](const TYPE_CODE& lhs, const TYPE_CODE& rhs) {
        int size=min(lhs.size(), rhs.size());
        for(int i=0; i<size; i++) {
            if(lhs[i] != rhs[i]) {
                return lhs[i] < rhs[i];
            }
        }
        return lhs.size() > rhs.size();
    };
    typedef std::pair<int, bool> TYPE_MAPPING_TABLE_ELEMENT;
    typedef std::vector< TYPE_MAPPING_TABLE_ELEMENT > TYPE_MAPPING_TABLE;

    typedef std::pair<TYPE_CODE, TYPE_MAPPING_TABLE> TYPE_CODE_QUEUE_ELEMENT;
    typedef multimap<TYPE_CODE, TYPE_MAPPING_TABLE, decltype(compare_code)> TYPE_CODE_QUEUE;
    typedef TYPE_CODE_QUEUE::iterator TYPE_CODE_QUEUE_ITER;
    int pattern_size=p.getNumEdges();

    TYPE_CODE_QUEUE codes(compare_code);

    //Choose starting vertex
    {
        for (int i = 0; i < p.getNumVertices(); i++) {
            TYPE_MAPPING_TABLE renumber_mapping(p.getNumVertices(), TYPE_MAPPING_TABLE_ELEMENT(-1, false));
            renumber_mapping[i]=TYPE_MAPPING_TABLE_ELEMENT(0, false);

            TYPE_CODE code(1, TYPE_EDGE_CODE(-1,-1,-1,-1,-1));
            codes.insert( TYPE_CODE_QUEUE_ELEMENT(code, renumber_mapping) );
        }
    }
    //maximum length of codes so far: d_max

    auto is_permutable = []( vector< tuple<int, int, int> >& list) {
        auto iter=list.end();
        while(iter != list.begin()) {
            int key1=get<0>(*(iter-1));
            int key2=get<1>(*(iter-1));
            auto cur_iter= std::find_if( list.begin(), iter, [key1, key2](const std::tuple<int, int, int>& element){ return get<0>(element) == key1 && get<1>(element) == key2;});
            if (std::next_permutation(cur_iter, iter))
                return true;
            iter=cur_iter;
        }
        return false;
    };

    auto min_mapped_vertex = [](TYPE_MAPPING_TABLE& renumber_mapping) {
        int min_val = MAX_INT, min_val_idx=-1;
        for(TYPE_MAPPING_TABLE::iterator iter=renumber_mapping.begin(); iter != renumber_mapping.end(); iter++) {
            if(iter->first != -1 && min_val > iter->first && !iter->second) {
                min_val = iter->first;
                min_val_idx = iter - renumber_mapping.begin();
            }

        }
        return min_val_idx;
    };

    size_t d_max=0;
    TYPE_CODE min_lex_code(1, TYPE_EDGE_CODE(MAX_INT, MAX_INT, MAX_INT, MAX_INT, MAX_INT) );
    TYPE_MAPPING_TABLE min_lex_code_mapping;

    while (!codes.empty()) {
        //Pop elements (BFS depth d)
        ////Find minimum key value
        TYPE_CODE cur_key=codes.begin()->first;

        ////Skip if d < d_max (exists longer but lexicographically smaller pattern)
        if (cur_key.size() < d_max) {
            //Remove and skip
            codes.erase(codes.lower_bound(cur_key), codes.upper_bound(cur_key));
            continue;
        }
        //Update d_max
        d_max = cur_key.size();

        TYPE_CODE_QUEUE new_codes(compare_code);

        for (TYPE_CODE_QUEUE_ITER element = codes.lower_bound(cur_key); element != codes.upper_bound(cur_key); element++) {
            //Find the smallest vertex v_src not visited in the BFS tree
            int v_src=min_mapped_vertex(element->second);

            int starting_renumber_id=(*std::max_element(element->second.begin(), element->second.end(),  [](const TYPE_MAPPING_TABLE_ELEMENT& a, const TYPE_MAPPING_TABLE_ELEMENT& b) {return a.first < b.first;} )).first + 1;

            //If mapping is done!
            if(v_src == -1) {
                if(starting_renumber_id != p.getNumVertices()) {
                    for (int u=0; u<p.getNumVertices(); u++) {
                        if(element->second[ u ].first == -1) {
                            TYPE_MAPPING_TABLE cur_renumber_mapping(element->second);
                            cur_renumber_mapping[ u ].first=starting_renumber_id;
                            cur_renumber_mapping[ u ].second=false;
                            new_codes.insert(TYPE_CODE_QUEUE_ELEMENT(element->first, cur_renumber_mapping));
                        }
                    }
                }
                else if(min_lex_code > element->first) {
                    min_lex_code=element->first;
                    min_lex_code_mapping=element->second;
                }
                continue;
            }
            if(element->first.size() == pattern_size + 1) {
                if(min_lex_code > element->first) {
                    min_lex_code=element->first;
                    min_lex_code_mapping=element->second;
                }
                break;
            }

            int v_src_renumbered=element->second[v_src].first;
            element->second[v_src].second=true;

            //Calculate updated code
            //Sort all edges by (edge label) order and then by (vertex id) order - already done
            TYPE_CODE cur_code(element->first);
            TYPE_CODE added_code;
            vector<tuple<int, int, int>> cur_adj_permutable;

            int cur_renumber_id = starting_renumber_id;
            for(int u=0; u<p.adj[v_src].size(); u++) {
                //Is visited;
                if (element->second[ p.adj[v_src][u].second ].first != -1) {
                    int u_renumbered=element->second[ p.adj[v_src][u].second ].first; 
                    added_code.push_back( TYPE_EDGE_CODE( v_src_renumbered, p.vorder[v_src], u_renumbered, p.vorder[u], p.adj[v_src][u].first ) );
                }
                else {
                    //Capture subset of the adjacency list to permutate
                    cur_adj_permutable.push_back( tuple<int, int, int>(p.vorder[u], p.adj[v_src][u].first, p.adj[v_src][u].second ) );
                    //added_code.push_back( TYPE_EDGE_CODE( v_src_renumbered, cur_renumber_id++, p.adj[v_src][u].first ) );
                    //cur_adj_permutable.push_back(p.adj[v_src][u]);
                }
            }

            if(added_code.size() > 0) {
                //added_code need to be sorted
                sort(added_code.begin(), added_code.end());
                cur_code.insert(cur_code.end(), added_code.begin(), added_code.end());
            }
            sort(cur_adj_permutable.begin(), cur_adj_permutable.end());
            for(int u=0; u<cur_adj_permutable.size(); u++) {
                int vl = get<0>(cur_adj_permutable[u]);
                int el = get<1>(cur_adj_permutable[u]);
                added_code.push_back( TYPE_EDGE_CODE( v_src_renumbered, p.vorder[v_src], cur_renumber_id++, vl, el) );
            }

            do {
                //For tie cases (vertices with the same order); calculate permutation
                TYPE_MAPPING_TABLE cur_renumber_mapping(element->second);

                for (int u=0; u<cur_adj_permutable.size(); u++) {
                    auto uid=get<2>(cur_adj_permutable[u]);
                    assert(cur_renumber_mapping[ uid ].first == -1);
                    cur_renumber_mapping[uid].first=starting_renumber_id+uid;
                    cur_renumber_mapping[uid].second=false;
                }
                new_codes.insert(TYPE_CODE_QUEUE_ELEMENT(cur_code, cur_renumber_mapping));
            } while(is_permutable(cur_adj_permutable));
        }
        //Remove processed codes
        codes.erase(codes.lower_bound(cur_key), codes.upper_bound(cur_key));

        //Update new codes
        codes.insert(new_codes.begin(), new_codes.end());
    }

    mapping.resize(min_lex_code_mapping.size());
    for (int i = 0; i < min_lex_code_mapping.size(); i++)
        mapping[i] = min_lex_code_mapping[i].first;

    string code("");
    //(src id, src label, dest id, dest label, edge label)(src id, ...) in BFS order
    for (int i = 1; i < min_lex_code.size(); i++) {
        code += string("(") + to_string(get<0>(min_lex_code[i])) + string(",") + to_string(get<1>(min_lex_code[i])) + string(",") + to_string(get<2>(min_lex_code[i])) + string(")");
    }

    return code;
}
