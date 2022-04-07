#include "../include/data_graph.h"
#include "../include/intersection.h"

#include <algorithm>
#include <fstream>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <cassert>
#include <experimental/filesystem>
#include <iostream>
#include <queue>

extern bool calculate_domains_context;

#ifdef PRINT_SEARCH
size_t CD_totalSearch = 0;
size_t CD_totalSearchLen = 0;

size_t totalSearch = 0;
size_t totalSearchLen = 0;
#endif

void DataGraph::ReadText(const char* fn) {
    std::cout << "DataGraph::ReadText " << fn << "\n";
    raw_.max_vl_ = raw_.max_el_ = -1;
    ifstream input(fn);
    string line;

    uint64_t num = 0;
    while (getline(input, line)) {
        auto tok = parse(line, " ");
        if (tok[0][0] == 'v') {
            vector<int> labels;
            for (int i = 2; i < tok.size(); i++) {
                int vl = stoi(tok[i]);
                if (vl < 0) continue;
                labels.push_back(vl);
                raw_.max_vl_ = std::max(raw_.max_vl_, vl);
            }
            raw_.vlabels_.push_back(labels);
        }
        if (tok[0][0] == 'e') {
            int src = stoi(tok[1]);
            int dst = stoi(tok[2]);
            for (int i = 3; i < tok.size(); i++) {
                int el = stoi(tok[i]);
                raw_.out_edges_.emplace_back(src, dst, el);
                raw_.in_edges_.emplace_back(dst, src, el);
                raw_.max_el_ = std::max(raw_.max_el_, el);
            }
        }
        ++num;
        if (num % 10000000 == 0) std::cout << "Parse " << num << " lines...\n";
    }
    raw_.el_cnt_.resize(raw_.max_el_ + 1);
    raw_.vl_cnt_.resize(raw_.max_vl_ + 1);
    raw_.el_rel_.resize(raw_.max_el_ + 1);
#ifdef INTERSECTION
    raw_.in_el_rel_.resize(raw_.max_el_ + 1);
#endif
    raw_.vl_rel_.resize(raw_.max_vl_ + 1);
    std::cout << "~DataGraph::ReadText " << fn << "\n";
}

void DataGraph::MakeBinary() {
    std::cout << "DataGraph::MakeBinary\n";
    int vnum = raw_.vlabels_.size();
    auto& out = raw_.out_edges_;
    auto& in = raw_.in_edges_;

    sort(out.begin(), out.end());
    sort(in.begin(), in.end());
    out.erase(unique(out.begin(), out.end()), out.end());
    in.erase(unique(in.begin(), in.end()), in.end());

    {
        auto& offsets = raw_.offset_;
        auto& elabels = raw_.label_;
        auto& adj_offsets = raw_.adj_offset_;
        auto& adj = raw_.adj_;

        int prev_u = -1, prev_el = -1;
        for (auto& e : out) {
            if (e.src != prev_u) {
                for (int u = prev_u + 1; u <= e.src; u++) {
                    offsets.push_back(elabels.size());
                }
                elabels.push_back(e.el);
                adj_offsets.push_back(adj.size()); 
                adj.push_back(e.dst);

                prev_u = e.src;
                prev_el = e.el;
            } else if (e.el != prev_el) {
                elabels.push_back(e.el);
                adj_offsets.push_back(adj.size()); 
                adj.push_back(e.dst);

                prev_el = e.el;
            } else {
                adj.push_back(e.dst);
            }
            raw_.el_cnt_[e.el]++;
        }

        for (int u = prev_u + 1; u <= vnum; u++)
            offsets.push_back(elabels.size());
        elabels.push_back(-1);
        adj_offsets.push_back(adj.size());
    }

    {
        auto& offsets = raw_.in_offset_;
        auto& elabels = raw_.in_label_;
        auto& adj_offsets = raw_.in_adj_offset_;
        auto& adj = raw_.in_adj_;

        int prev_u = -1, prev_el = -1;
        for (auto& e : in) {
            if (e.src != prev_u) {
                for (int u = prev_u + 1; u <= e.src; u++) {
                    offsets.push_back(elabels.size());
                }
                elabels.push_back(e.el);
                adj_offsets.push_back(adj.size()); 
                adj.push_back(e.dst);

                prev_u = e.src;
                prev_el = e.el;
            } else if (e.el != prev_el) {
                elabels.push_back(e.el);
                adj_offsets.push_back(adj.size()); 
                adj.push_back(e.dst);

                prev_el = e.el;
            } else {
                adj.push_back(e.dst);
            }
        }

        for (int u = prev_u + 1; u <= vnum; u++)
            offsets.push_back(elabels.size());
        elabels.push_back(-1);
        adj_offsets.push_back(adj.size());
    }

#if defined (INTERSECTION) && !defined(PRECOMPUTE_FR)
    for (int u = 0; u < vnum; u++) {
        for (auto& vl : raw_.vlabels_[u]) {
            raw_.vl_cnt_[vl]++;
            raw_.vl_rel_[vl].push_back(u);
        }
    }
#else
    int num_vl = 0;
    for (int u = 0; u < vnum; u++) {
        raw_.vl_offset_.push_back(num_vl);
        sort(raw_.vlabels_[u].begin(), raw_.vlabels_[u].end());
        num_vl += raw_.vlabels_[u].size();
        raw_.vl_.insert(raw_.vl_.end(), raw_.vlabels_[u].begin(), raw_.vlabels_[u].end()); 
        for (auto& vl : raw_.vlabels_[u]) {
            raw_.vl_cnt_[vl]++;
            raw_.vl_rel_[vl].push_back(u);
        }
    }
    raw_.vl_offset_.push_back(num_vl);
#endif

    {
        for (auto& e : out) {
#ifdef INTERSECTION
            raw_.el_rel_[e.el].push_back(e.src);
            raw_.in_el_rel_[e.el].push_back(e.dst);
#else
            raw_.el_rel_[e.el].push_back(make_pair(e.src, e.dst));
#endif
        }
#ifdef INTERSECTION
        for (int el = 0; el <= raw_.max_el_; el++) {
            auto& out = raw_.el_rel_[el];
            sort(out.begin(), out.end());
            out.erase(unique(out.begin(), out.end()), out.end());
            auto& in = raw_.in_el_rel_[el];
            sort(in.begin(), in.end());
            in.erase(unique(in.begin(), in.end()), in.end());
        }
#endif
    }
    std::cout << "~DataGraph::MakeBinary\n";
}

size_t DataGraph::BinarySize() {
    size_t ret = 0;

    size_t vn = raw_.vlabels_.size();

    ret += sizeof(int) * (vn + 1); 
    ret += sizeof(int);
    ret += sizeof(int) * raw_.label_.size() * 2;
    ret += sizeof(int);
    ret += sizeof(int) * raw_.adj_.size();

    ret += sizeof(int) * (vn + 1); 
    ret += sizeof(int);
    ret += sizeof(int) * raw_.in_label_.size() * 2;
    ret += sizeof(int);
    ret += sizeof(int) * raw_.in_adj_.size();

#ifdef INTERSECTION
    ret += sizeof(int) * (raw_.max_el_ + 2);
    for (int el = 0; el <= raw_.max_el_; el++)
        ret += sizeof(int) * raw_.el_rel_[el].size();
    ret += sizeof(int) * (raw_.max_el_ + 2);
    for (int el = 0; el <= raw_.max_el_; el++)
        ret += sizeof(int) * raw_.in_el_rel_[el].size();
#else
    ret += sizeof(int) * (raw_.max_el_ + 2);
    ret += sizeof(pair<int, int>) * raw_.out_edges_.size(); 
#endif

    ret += sizeof(int) * (raw_.max_vl_ + 2); 
    for (int vl = 0; vl <= raw_.max_vl_; vl++)
        ret += sizeof(int) * raw_.vl_rel_[vl].size();

#if !defined (INTERSECTION) || defined(PRECOMPUTE_FR)
    ret += sizeof(int) * (vn + 1); 
    ret += sizeof(int);
    ret += sizeof(int) * raw_.vl_.size();
#endif

    return ret;
}

bool DataGraph::HasBinary(const char* filename) {
#ifdef INTERSECTION
    string metadata = string(filename) + ".inter.meta";
#else
    string metadata = string(filename) + ".graph.meta";
#endif
    return std::experimental::filesystem::exists(metadata.c_str());
}

void DataGraph::WriteBinary(const char* filename) {
#ifdef INTERSECTION
    string fname = string(filename) + ".inter";
#else
    string fname = string(filename) + ".graph";
#endif
    std::cout << "DataGraph::WriteBinary" << fname << "\n";
    string metadata = fname + ".meta";
    FILE* fp = fopen(metadata.c_str(), "w");
    size_t encode_size = BinarySize();

    int vn = raw_.vlabels_.size();
    int en = raw_.out_edges_.size(); 

    fprintf(fp, "%d %d %d %d %zu\n", vn, en, raw_.max_vl_+1, raw_.max_el_+1, encode_size);
    for (int vl = 0; vl <= raw_.max_vl_; vl++)
        fprintf(fp, "%d ", raw_.vl_cnt_[vl]); 
    fprintf(fp, "\n");
    for (int el = 0; el <= raw_.max_el_; el++)
        fprintf(fp, "%d ", raw_.el_cnt_[el]); 
    fprintf(fp, "\n");
    fclose(fp);

    char* buffer = new char[encode_size];
    char* orig = buffer;
    int size[1] = {0};

    {
        assert(raw_.offset_.size() == vn + 1);
        assert(raw_.offset_.back() + 1 == raw_.label_.size());
        size_t offset_size = sizeof(int) * (vn + 1); 
        memcpy(buffer, raw_.offset_.data(), offset_size); 
        buffer += offset_size;

        size[0] = raw_.label_.size();
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);

        size_t label_size = sizeof(int) * raw_.label_.size(); 
        assert(raw_.label_.size() == raw_.adj_offset_.size());
        memcpy(buffer, raw_.label_.data(), label_size);
        buffer += label_size;
        memcpy(buffer, raw_.adj_offset_.data(), label_size);
        buffer += label_size;

        size[0] = raw_.adj_.size();
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);

        size_t adj_size = sizeof(int) * raw_.adj_.size();
        assert(raw_.adj_offset_.back() == raw_.adj_.size());
        memcpy(buffer, raw_.adj_.data(), adj_size);
        buffer += adj_size;
    }

    {
        assert(raw_.in_offset_.size() == vn + 1);
        assert(raw_.in_offset_.back() + 1 == raw_.in_label_.size());
        size_t offset_size = sizeof(int) * (vn + 1); 
        memcpy(buffer, raw_.in_offset_.data(), offset_size); 
        buffer += offset_size;

        size[0] = raw_.in_label_.size();
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);

        size_t label_size = sizeof(int) * raw_.in_label_.size(); 
        assert(raw_.in_label_.size() == raw_.in_adj_offset_.size());
        memcpy(buffer, raw_.in_label_.data(), label_size);
        buffer += label_size;
        memcpy(buffer, raw_.in_adj_offset_.data(), label_size);
        buffer += label_size;

        size[0] = raw_.in_adj_.size();
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);

        size_t adj_size = sizeof(int) * raw_.in_adj_.size();
        assert(raw_.in_adj_offset_.back() == raw_.in_adj_.size());
        memcpy(buffer, raw_.in_adj_.data(), adj_size);
        buffer += adj_size;
    }

    {
        assert(raw_.el_rel_.size() == raw_.max_el_ + 1);
        size[0] = 0;  
        for (int el = 0; el <= raw_.max_el_; el++) {
            memcpy(buffer, size, sizeof(int)); 
            size[0] += raw_.el_rel_[el].size();
            buffer += sizeof(int);
        }
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);
        for (int el = 0; el <= raw_.max_el_; el++) {
#ifdef INTERSECTION
            memcpy(buffer, raw_.el_rel_[el].data(), sizeof(int) * raw_.el_rel_[el].size()); 
            buffer += sizeof(int) * raw_.el_rel_[el].size();
#else
            memcpy(buffer, raw_.el_rel_[el].data(), sizeof(pair<int, int>) * raw_.el_rel_[el].size()); 
            buffer += sizeof(pair<int, int>) * raw_.el_rel_[el].size();
#endif
        }
    }

#ifdef INTERSECTION
    {
        assert(raw_.in_el_rel_.size() == raw_.max_el_ + 1);
        size[0] = 0;  
        for (int el = 0; el <= raw_.max_el_; el++) {
            memcpy(buffer, size, sizeof(int)); 
            size[0] += raw_.in_el_rel_[el].size();
            buffer += sizeof(int);
        }
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);
        for (int el = 0; el <= raw_.max_el_; el++) {
            memcpy(buffer, raw_.in_el_rel_[el].data(), sizeof(int) * raw_.in_el_rel_[el].size()); 
            buffer += sizeof(int) * raw_.in_el_rel_[el].size();
        }
    }
#endif

    {
        assert(raw_.vl_rel_.size() == raw_.max_vl_ + 1);
        size[0] = 0;
        for (int vl = 0; vl <= raw_.max_vl_; vl++) {
            memcpy(buffer, size, sizeof(int)); 
            size[0] += raw_.vl_rel_[vl].size();
            buffer += sizeof(int);
        }
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);
        for (int vl = 0; vl <= raw_.max_vl_; vl++) {
            memcpy(buffer, raw_.vl_rel_[vl].data(), sizeof(int) * raw_.vl_rel_[vl].size()); 
            buffer += sizeof(int) * raw_.vl_rel_[vl].size();
        }
    }

#if !defined (INTERSECTION) || defined(PRECOMPUTE_FR)
    {
        assert(raw_.vl_offset_.size() == vn + 1);
        size_t offset_size = sizeof(int) * raw_.vl_offset_.size();
        memcpy(buffer, raw_.vl_offset_.data(), offset_size); 
        buffer += offset_size;

        size[0] = raw_.vl_.size();
        memcpy(buffer, size, sizeof(int));
        buffer += sizeof(int);

        size_t label_size = sizeof(int) * raw_.vl_.size();
        memcpy(buffer, raw_.vl_.data(), label_size);
        buffer += label_size;
    }
#endif

    assert((buffer - orig) == encode_size);
    FILE* f = fopen(fname.c_str(), "w");
    fwrite(orig, 1, encode_size, f);
    cout << "wrote " << (buffer - orig) << " bytes to file " << fname << endl;
    fclose(f);
    std::cout << "~DataGraph::WriteBinary" << fname << "\n";
}

void DataGraph::ReadBinary(const char* filename) {
#ifdef INTERSECTION
    string fname = string(filename) + ".inter";
#else
    string fname = string(filename) + ".graph";
#endif
    std::cout << "DataGraph::ReadBinary" << fname << "\n";
    string metadata = fname + ".meta";
    FILE* fp = fopen(metadata.c_str(), "r");
    size_t encode_size;
    fscanf(fp, "%d %d %d %d %zu", &vnum_, &enum_, &vl_num_, &el_num_, &encode_size);

    vl_cnt_.resize(vl_num_);
    el_cnt_.resize(el_num_);
    for (int vl = 0; vl < vl_num_; vl++)
        fscanf(fp, "%d ", &vl_cnt_[vl]); 
    for (int el = 0; el < el_num_; el++)
        fscanf(fp, "%d ", &el_cnt_[el]); 
    fclose(fp);

    FILE* f = fopen(fname.c_str(), "r");
    char* buffer = new char[encode_size];
    fread(buffer, 1, encode_size, f);
    fclose(f);
    char* orig = buffer;

    {
        offset_ = (const int*) buffer;
        assert(offset_[0] == 0);
        buffer += sizeof(int) * (vnum_ + 1);

        int* size = (int*) buffer;
        assert(offset_[vnum_] + 1 == size[0]);
        buffer += sizeof(int);

        label_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];
        adj_offset_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];

        size = (int*) buffer;
        buffer += sizeof(int);

        adj_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];
    }

    {
        in_offset_ = (const int*) buffer;
        assert(in_offset_[0] == 0);
        buffer += sizeof(int) * (vnum_ + 1);

        int* size = (int*) buffer;
        assert(in_offset_[vnum_] + 1 == size[0]);
        buffer += sizeof(int);

        in_label_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];
        in_adj_offset_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];

        size = (int*) buffer;
        buffer += sizeof(int);

        in_adj_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];
    }

    {
        el_rel_offset_ = (const int*) buffer; 
        assert(el_rel_offset_[0] == 0);
        buffer += sizeof(int) * (el_num_ + 1);

#ifdef INTERSECTION
        el_rel_ = (const int*) buffer;
        buffer += sizeof(int) * el_rel_offset_[el_num_];
#else
        el_rel_ = (const pair<int, int>*) buffer;
        buffer += sizeof(pair<int, int>) * el_rel_offset_[el_num_];
#endif
    }

#ifdef INTERSECTION
    {
        in_el_rel_offset_ = (const int*) buffer; 
        assert(in_el_rel_offset_[0] == 0);
        buffer += sizeof(int) * (el_num_ + 1);

        in_el_rel_ = (const int*) buffer;
        buffer += sizeof(int) * in_el_rel_offset_[el_num_];
    }
#endif

    {
        vl_rel_offset_ = (const int*) buffer; 
        assert(vl_rel_offset_[0] == 0);
        buffer += sizeof(int) * (vl_num_ + 1);

        vl_rel_ = (const int*) buffer;
        buffer += sizeof(int) * vl_rel_offset_[vl_num_];
    }

#if !defined (INTERSECTION) || defined(PRECOMPUTE_FR)
    {
        vl_offset_ = (const int*) buffer;
        assert(vl_offset_[0] == 0);
        buffer += sizeof(int) * (vnum_ + 1);

        int* size = (int*) buffer;
        buffer += sizeof(int);

        vl_ = (const int*) buffer;
        buffer += sizeof(int) * size[0];
    }
#endif

    assert((buffer - orig) == encode_size);
    std::cout << "~DataGraph::ReadBinary" << fname << "\n";
}

void DataGraph::ClearRawData() {
    raw_.vl_cnt_.clear();
    raw_.el_cnt_.clear();
    raw_.vlabels_.clear();
    raw_.out_edges_.clear();
    raw_.in_edges_.clear();
    raw_.offset_.clear();
    raw_.label_.clear();
    raw_.adj_offset_.clear();
    raw_.adj_.clear();
    raw_.in_offset_.clear();
    raw_.in_label_.clear();
    raw_.in_adj_offset_.clear();
    raw_.in_adj_.clear();
    raw_.vl_offset_.clear();
    raw_.vl_.clear();
    raw_.vl_rel_.clear();
    raw_.el_rel_.clear();
#ifdef INTERSECTION
    raw_.in_el_rel_.clear();
#endif
}

int DataGraph::GetNumVertices() {
    return vnum_;
}

int DataGraph::GetNumVertices(int vl) {
    return vl == -1 ? vnum_ : vl_cnt_[vl];
}

int DataGraph::GetNumEdges() {
    return enum_;
}

int DataGraph::GetNumEdges(int el) {
    return el_cnt_[el];
}

int DataGraph::GetNumVLabels(int v) {
    if (v == -1)
        return vl_num_; 
    else {
        range r = GetVLabels(v);
        return r.end - r.begin;
    }
}

int DataGraph::GetNumELabels(int v, bool dir) {
    if (v == -1)
        return el_num_; 
    else {
        range r = GetELabels(v, dir);
        return r.end - r.begin;
    }
}

range DataGraph::GetVLabels(int v) {
    range r;
    r.begin = vl_ + vl_offset_[v];
    r.end   = vl_ + vl_offset_[v+1];
    return r;
}

bool DataGraph::HasVLabel(int v, int vl) {
    int begin = vl_offset_[v];
    int end   = vl_offset_[v+1];

    int res = search(vl_, begin, end, vl);
    return res != -1;
}

range DataGraph::GetELabels(int v, bool dir = true) {
    const int* offset = dir ? offset_ : in_offset_;
    const int* label  = dir ? label_  : in_label_;

    range r;
    r.begin = label + offset[v]; 
    r.end   = label + offset[v+1]; 
    return r;
}

bool DataGraph::HasELabel(int v, int el, bool dir = true) {
    range r = GetAdj(v, el, dir);
    return r.end != r.begin;
}

int DataGraph::GetELabelIndex(int v, int el, bool dir = true) {
    const int* offset = dir ? offset_ : in_offset_;
    const int* label  = dir ? label_ : in_label_; 
    const int* adj_o  = dir ? adj_offset_ : in_adj_offset_; 
    const int* adj    = dir ? adj_   : in_adj_;

    int begin = offset[v];
    int end   = offset[v+1];
    int res = search(label, begin, end, el);
    return res == -1 ? -1 : res - begin;
}

range DataGraph::GetAdj(int v, int el, bool dir = true) {
    const int* offset = dir ? offset_ : in_offset_;
    const int* label  = dir ? label_ : in_label_; 
    const int* adj_o  = dir ? adj_offset_ : in_adj_offset_; 
    const int* adj    = dir ? adj_   : in_adj_;

    range r;
    r.begin = r.end = NULL;

    int begin = offset[v];
    int end   = offset[v+1];
#ifdef PRINT_SEARCH
    if (calculate_domains_context) {
        CD_totalSearchLen += (end - begin);
        CD_totalSearch++;
    }
    else {
        totalSearchLen += (end - begin);
        totalSearch++;
    }
#endif
    int res = search(label, begin, end, el);
    if (res == -1) return r;
    r.begin = adj + adj_o[res];
    r.end   = adj + adj_o[res+1];
    return r;
}

int DataGraph::GetAdjSize(int v, int el, bool dir = true) {
    range r = GetAdj(v, el, dir);
    return r.end - r.begin;
}

bool DataGraph::HasEdge(int u, int v, int el, bool dir = true) {
    int u_adj_size = GetAdjSize(u, el, dir);
    int v_adj_size = GetAdjSize(v, el, !dir);

    range r = u_adj_size < v_adj_size ? GetAdj(u, el, dir) : GetAdj(v, el, !dir); 
    int target = u_adj_size < v_adj_size ? v : u;

    int s = 0; 
    int e = r.end - r.begin;

    return search(r.begin, s, e, target) != -1;
}

#ifndef INTERSECTION
vector<int> DataGraph::GetRandomEdge(int el) {
    vector<int> ret;
    int begin = el_rel_offset_[el]; 
    int end   = el_rel_offset_[el+1]; 

    if (begin == end)
        return ret;
    int r = rand();
    r %= (end - begin);
    r += begin;
    ret.resize(2);
    ret[0] = el_rel_[r].first;
    ret[1] = el_rel_[r].second;
    return ret;
}

vector<int> DataGraph::GetEdge(int el, int i) {
    vector<int> ret;
    int begin = el_rel_offset_[el]; 
    int end   = el_rel_offset_[el+1]; 

    if (begin == end)
        return ret;
    assert(i < end - begin);
    int r = i + begin;
    ret.resize(2);
    ret[0] = el_rel_[r].first;
    ret[1] = el_rel_[r].second;
    return ret;
}

double DataGraph::GetJVD(int el, int col) {
    pair<int, int> key(el, col);
    if (jvd_map_.find(key) != jvd_map_.end())
        return jvd_map_[key];

    double jvd = 0;
    int begin = el_rel_offset_[el]; 
    int end   = el_rel_offset_[el+1]; 

    if (begin != end) {
        set<int> V;
        for (int r = begin; r < end; r++) {
            V.insert(col == 0 ? el_rel_[r].first : el_rel_[r].second); 
        }
        jvd = (double)V.size()/(end - begin);
        assert(jvd <= 1);
    }
    jvd_map_[key] = jvd;
    return jvd;
}
#endif

vector<int> DataGraph::GetRandomEdge(int v, int el, bool dir) {
    vector<int> ret;
    range r = GetAdj(v, el, dir);
    if (r.begin == r.end)
        return ret;
    int rv = rand();
    rv %= (r.end - r.begin);
    int other = r.begin[rv];
    ret.resize(2);
    ret[0] = dir ? v : other;
    ret[1] = dir ? other : v;
    return ret;
}

vector<int> DataGraph::GetRandomVertex(int vl) {
    vector<int> ret;

    int begin = vl_rel_offset_[vl]; 
    int end   = vl_rel_offset_[vl+1]; 

    if (begin == end)
        return ret;
    int r = rand();
    r %= (end - begin);
    r += begin; 
    ret.push_back(vl_rel_[r]);
    return ret;
}

vector<int> DataGraph::GetVertex(int vl, int i) {
    vector<int> ret;

    int begin = vl_rel_offset_[vl]; 
    int end   = vl_rel_offset_[vl+1]; 

    if (begin == end)
        return ret;
    assert(i < end - begin);
    int r = i + begin; 
    ret.push_back(vl_rel_[r]);
    return ret;
}

range DataGraph::GetVertices(int vl) {
    //assert(vl != -1);
    range r;
    r.begin = vl_rel_ + vl_rel_offset_[vl]; 
    r.end   = vl_rel_ + vl_rel_offset_[vl+1];
    return r;
}

#ifdef INTERSECTION
range DataGraph::GetVertices(int el, bool dir) {
    const int* offset = dir ? el_rel_offset_ : in_el_rel_offset_; 
    const int* rel    = dir ? el_rel_        : in_el_rel_; 

    range r;
    r.begin = rel + offset[el];
    r.end   = rel + offset[el+1];
    return r;
}

int DataGraph::GetNumVertices(int el, bool dir) {
    const int* offset = dir ? el_rel_offset_ : in_el_rel_offset_; 
    return offset[el+1] - offset[el];
}

void DataGraph::GroupELabels(vector<int>& label2group, vector<vector<int>>& group2label) {
    vector<int> new_el_cnt(group2label.size(), 0);
    for (int i = 0; i < group2label.size(); i++)
        for (int l : group2label[i])
            new_el_cnt[i] += el_cnt_[l];
    for (int i = 0; i < group2label.size(); i++)
        el_cnt_[i] = new_el_cnt[i];
    el_num_ = group2label.size();

    {
        int* new_el_rel_offset = new int[group2label.size()+1];
        int* new_el_rel = new int[el_rel_offset_[label2group.size()]];
        int cur_offset = 0;

        for (int i = 0; i < group2label.size(); i++) {
            new_el_rel_offset[i] = cur_offset; 

            vector<range> ranges;
            for (auto l : group2label[i]) {
                range r = GetVertices(l, true);
                if (r.end != r.begin)
                    ranges.push_back(r);
            }

            typedef pair<int, int> elem;
            priority_queue<elem, vector<elem>, greater<elem>> q;
            for (int j = 0; j < ranges.size(); j++) {
                q.emplace(ranges[j].begin[0], j);
                ranges[j].begin++;
            }
            int cur_id = -1;
            while (!q.empty()) {
                pair<int, int> t = q.top();
                q.pop();
                if (t.first > cur_id) {
                    cur_id = t.first;
                    new_el_rel[cur_offset++] = cur_id;
                }
                if (ranges[t.second].end != ranges[t.second].begin) {
                    q.emplace(ranges[t.second].begin[0], t.second); 
                    ranges[t.second].begin++;
                }
            }
        }
        new_el_rel_offset[group2label.size()] = cur_offset;
        el_rel_offset_ = new_el_rel_offset;
        el_rel_ = new_el_rel;
    }
    {
        int* new_in_el_rel_offset = new int[group2label.size()+1];
        int* new_in_el_rel = new int[in_el_rel_offset_[label2group.size()]];
        int cur_offset = 0;

        for (int i = 0; i < group2label.size(); i++) {
            new_in_el_rel_offset[i] = cur_offset; 

            vector<range> ranges;
            for (auto l : group2label[i]) {
                range r = GetVertices(l, false);
                if (r.end != r.begin)
                    ranges.push_back(r);
            }

            typedef pair<int, int> elem;
            priority_queue<elem, vector<elem>, greater<elem>> q;
            for (int j = 0; j < ranges.size(); j++) {
                q.emplace(ranges[j].begin[0], j);
                ranges[j].begin++;
            }
            int cur_id = -1;
            while (!q.empty()) {
                pair<int, int> t = q.top();
                q.pop();
                if (t.first > cur_id) {
                    cur_id = t.first;
                    new_in_el_rel[cur_offset++] = cur_id;
                }
                if (ranges[t.second].end != ranges[t.second].begin) {
                    q.emplace(ranges[t.second].begin[0], t.second); 
                    ranges[t.second].begin++;
                }
            }
        }
        new_in_el_rel_offset[group2label.size()] = cur_offset;

        in_el_rel_offset_ = new_in_el_rel_offset;
        in_el_rel_ = new_in_el_rel;
    }
    vector<vector<int>> vids(group2label.size());
    {
        int* new_offset = const_cast<int*>(offset_); 
        int* new_label  = const_cast<int*>(label_);
        int* new_adj_o  = const_cast<int*>(adj_offset_);
        int* new_adj    = const_cast<int*>(adj_);

        int cur_offset = 0, cur_adj_offset = 0;
        for (int v = 0; v < vnum_; v++) {
            int begin = offset_[v];
            int end   = offset_[v+1];

            for (int o = begin; o < end; o++) {
                int el = label_[o];
                assert(el >= 0 && el < label2group.size());
                int group = label2group[el];
                //pruned label
                if (group == -1)
                    continue;
                for (int adj_o = adj_offset_[o]; adj_o < adj_offset_[o+1]; adj_o++) {
                    vids[group].push_back(adj_[adj_o]);
                }
            }

            new_offset[v] = cur_offset;

            for (int i = 0; i < group2label.size(); i++) {
                if (vids[i].size() > 0) {
                    new_label[cur_offset] = i;
                    new_adj_o[cur_offset] = cur_adj_offset;

                    cur_offset++;

                    std::sort(vids[i].begin(), vids[i].end());
                    vids[i].erase(std::unique(vids[i].begin(), vids[i].end()), vids[i].end());

                    memcpy(new_adj + cur_adj_offset, vids[i].data(), sizeof(int) * vids[i].size()); 
                    cur_adj_offset += vids[i].size();

                    vids[i].clear();
                }
            }
        }

        new_offset[vnum_] = cur_offset;
        new_label[cur_offset] = -1; 
        new_adj_o[cur_offset] = cur_adj_offset; 
    }
    {
        int* new_in_offset = const_cast<int*>(in_offset_); 
        int* new_in_label  = const_cast<int*>(in_label_);
        int* new_in_adj_o  = const_cast<int*>(in_adj_offset_);
        int* new_in_adj    = const_cast<int*>(in_adj_);

        int cur_offset = 0, cur_adj_offset = 0;
        for (int v = 0; v < vnum_; v++) {
            int begin = in_offset_[v];
            int end   = in_offset_[v+1];

            for (int o = begin; o < end; o++) {
                int el = in_label_[o];
                int group = label2group[el];
                //pruned label
                if (group == -1)
                    continue;
                for (int in_adj_o = in_adj_offset_[o]; in_adj_o < in_adj_offset_[o+1]; in_adj_o++)
                    vids[group].push_back(in_adj_[in_adj_o]);
            }

            new_in_offset[v] = cur_offset;

            for (int i = 0; i < group2label.size(); i++) {
                if (vids[i].size() > 0) {
                    new_in_label[cur_offset] = i;
                    new_in_adj_o[cur_offset] = cur_adj_offset;

                    cur_offset++;

                    std::sort(vids[i].begin(), vids[i].end());
                    vids[i].erase(std::unique(vids[i].begin(), vids[i].end()), vids[i].end());

                    memcpy(new_in_adj + cur_adj_offset, vids[i].data(), sizeof(int) * vids[i].size()); 
                    cur_adj_offset += vids[i].size();

                    vids[i].clear();
                }
            }
        }

        new_in_offset[vnum_] = cur_offset;
        new_in_label[cur_offset] = -1; 
        new_in_adj_o[cur_offset] = cur_adj_offset; 
    }
}

void DataGraph::GroupVLabels(vector<int>& label2group, vector<vector<int>>& group2label) {
    vl_num_ = group2label.size();

    {
        int* new_vl_rel_offset = new int[group2label.size()+1];
        int* new_vl_rel = new int[vl_rel_offset_[label2group.size()]];
        int cur_offset = 0;

        for (int i = 0; i < group2label.size(); i++) {
            new_vl_rel_offset[i] = cur_offset; 

            vector<range> ranges;
            for (auto l : group2label[i]) {
                range r = GetVertices(l);
                if (r.end != r.begin)
                    ranges.push_back(r);
            }

            int union_size = do_union(ranges, new_vl_rel + cur_offset);
            vl_cnt_[i] = union_size;
            cur_offset += union_size;
        }
        new_vl_rel_offset[group2label.size()] = cur_offset;
        vl_rel_offset_ = new_vl_rel_offset;
        vl_rel_ = new_vl_rel;
    }
    for (int i = 0; i < group2label.size(); i++) {
        range r = GetVertices(i);
        assert(GetNumVertices(i) == r.end - r.begin); 
    }
}


#endif
