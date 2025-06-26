/*
X-Set Heuristic Program for minimizing XOR's in executing a binary
matrix-vector product for erasure coding.

Version 1.0 

James S. Plank
Professor
EECS Department
University of Tennessee
Knoxville, TN 37996
plank@cs.utk.edu
http://web.eecs.utk.edu/~jplank


Copyright (c) 2010, James S. Plank
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in
   the documentation and/or other materials provided with the
   distribution.

 - Neither the name of the University of Tennessee nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <list>
#include <set>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
using namespace std;

#define MW (1111)
#define MW_SS (1112)
#define MW_SMALLEST_SUM (1113)
#define MW_SQ (1114)
#define UBER_XSET (1115)
#define MW_MATCHING (1116)
#define SUBEX (1117)

int Technique = 0;
int Thresh = 0;
int NStart = 0;
int Give_Up_When = 0;

#define BMS_SET(C) elts[(C)/8] |= (1 << ((C)%8))
#define BMS_CLEAR(C) elts[(C)/8] &= (0xff ^ (1 << (C)%8))
#define BMS_ISSET(C) ((elts[(C)/8] & (1 << ((C)%8))) ? 1 : 0)

#define BM_SET(BM, C) BM.elts[(C)/8] |= (1 << ((C)%8))
#define BM_CLEAR(BM, C) BM.elts[(C)/8] &= (0xff ^ (1 << (C)%8))
#define BM_ISSET(BM, C) ((BM.elts[(C)/8] & (1 << ((C)%8))) ? 1 : 0)

#define BMP_SET(BM, C) BM->elts[(C)/8] |= (1 << ((C)%8))
#define BMP_CLEAR(BM, C) BM->elts[(C)/8] &= (0xff ^ (1 << (C)%8))
#define BMP_ISSET(BM, C) ((BM->elts[(C)/8] & (1 << ((C)%8))) ? 1 : 0)

#define VIT(i, v) for (i = 0; i < v.size(); i++) 
#define IT(it, ds) for (it = ds.begin(); it != ds.end(); it++)
#define ITP(it, ds) for (it = ds->begin(); it != ds->end(); it++)
#define FUP(i, n) for (i = 0; i < n; i++)

#define O1(v) cout << v << endl
#define O2(v1, v2) cout << v1 << " " << v2 << endl
#define O3(v1, v2, v3) cout << v1 << " " << v2 << " " << v3 << endl
#define OVEC(v) { int iii; VIT(iii, v) cout << v[iii] << " " ; cout << endl; }

typedef map<int, string> ISmap;
typedef map<int, int> IImap;
typedef map<string, int> SImap;
typedef map<string, string> SSmap;

typedef map <string, int> SIMap;
typedef SIMap::iterator SIMit;

typedef vector <int> IVec;
typedef vector <char> CVec;
typedef vector <string> SVec;
typedef vector <IVec> IIVec;

class Bitset {
  public:
    string Key;
    void Print(int cols, int newline);
    void SetKey();
    Bitset(int c);
    Bitset(class Bitset *b1, class Bitset *b2);   /* Create the bitset from the XOR of two bitsets */
    Bitset(string &s, int c);
    void XOR_With(class Bitset *b);
    int Diff(class Bitset *b);
    int Ones();
    IVec Diff_Inds(class Bitset *b);
    int Is_Equal_To(class Bitset *b);
    CVec elts;
};

typedef vector <Bitset> BVec;
typedef vector <Bitset *> BPVec;

Bitset *ZERO;
typedef pair <int, int> IPair;

// Structure to track how each element was computed in X-Sets
struct XSetComputationStep {
    int result_index;        // Index in sched array
    int operand1_index;      // First operand index (-1 if it's a data element)
    int operand2_index;      // Second operand index (-1 if single operand or data element)
    bool is_target;          // True if this is a target element
    int target_index;        // Index in original target list (-1 if not a target)
};

typedef map <int, Bitset *> IBMap;
typedef IBMap::iterator IBMit;

typedef map <string, Bitset *> SBMap;
typedef SBMap::iterator SBMit;

typedef map <string, IPair> SPMap;
typedef SPMap::iterator SPMit;

/* typedef list <IBMap> IBMList;
typedef IBMList::iterator IBMLit; */

/* typedef list <SPMap> SPMList;
  typedef SPMList::iterator SPMLit; */

class XSet {
  public: 
    string id;
    IBMap indices;
    SPMap sums;    /* or "edges" */
};
  
typedef map <string, XSet> XMap;
typedef XMap::iterator XMit;

class MXS {
  public: 
    Bitset *t;
    int ind;
    int sz;
    
    XMap XSets;
    SIMap MX_Global_Edge_Weights;
    SIMap MX_Local_Edge_Weights;
    void Print();
};

typedef list <MXS *> MList;
typedef MList::iterator MLit;
typedef map <string, MXS *> SMMap;
typedef SMMap::iterator SMMit;

typedef set <string> SSet;

typedef map <int, MXS *> IMXMap;
typedef IMXMap::iterator IMXMit;
typedef map <string, IMXMap> SIMXMap;
typedef SIMXMap::iterator SIMXMit;

string ITOS(int i) { char buf[20]; string k; sprintf(buf, "%d", i); k = buf; return k; }
typedef map <string, string> SSMap;
typedef SSMap::iterator SSMit;

typedef set <int> ISet;
typedef ISet::iterator ISit;

class Schedule {
  public:
    BPVec targets;
    BPVec sched;
    SIMap schedtree;
    SMMap  sets;
    vector <IVec> Sched_By_Bits;     /* Sched_By_Bits[i] contains a vector of indices of schedule elements that have i bits. */
                                     /* I'm ignoring Sched_By_Bits[1].  That's just the data */
    vector <SSet> Edges_By_Weight;   /* element i contains a set of b->Key's that have that weight */
    SIMXMap Sched_Edge_Weights;      /* This is a map of all sums.  Key = b->key.  Val = map of targets in which sum appears */
    int cols;
    IVec Enumerator;                 /* This is for enumerating all starting points */
    Bitset *EN_Cannot_Use;              /* These are the starting points */
    vector <XSetComputationStep> computation_history;  // Track how each element was computed
    void Initialize_Enumerator(int last_product);
    int Enumerator_Done();
    void Enumerator_Next();
    void Read();
    void Copy(class Schedule *S);
    void Create_Initial_Sets();
    void Add_Element(Bitset *b);
    void Add_Element_With_History(Bitset *b, int op1_idx, int op2_idx);
    void Add_Element_From_GetBest(Bitset *b);
    void Add_Element_To_MX_Set(MXS *mx, XMit xptr, int index, int global_weight);
    void Add_Global_Edge(MXS *mx, string &s);
    void Subtract_Element_From_MX_Set(MXS *mx, XMit xptr, int index, int global_weight);
    void Subtract_Global_Edge(MXS *mx, string &s);
    void Clear_Sets_And_Potentials(MXS *mx);
    void Add_New_MX(MXS *mx, IVec &inds);                    /* Adds a new MX set.  If the set's size is == mx->sz, 
                                                                it adds the sums to global_edge_w */
    void Remove_Global_Edge_Weights(MXS *mx, XMit xptr); /* If xptr is no longer of size mx->sz, remove all the sums from 
                                                                global_edge_weights */
    void Add_Global_Edge_Weights(MXS *mx, XMit xptr);    /* If xptr is now of size mx->sz, add all the sums into global_eddge_weights */
    Bitset *Find_Two();
    vector <Bitset *> GetBest(int maxelts, class Schedule *me);
    void Print();
    void PrintMathematical();
    void Check();
    string Key;
    void  SetKey();
    ~Schedule();
};

/* BEGIN MATCHING CODE */

class Vertex {
  public:
    string name;
    list <class Edge *> edges;
    class Blossom *blossom;
    class Blossom *in_blossom;
    class Blossom *tmp;
    int exposed;
    int marked;
    int d_to_root;
    Vertex *root;
    class Edge *edge_to_root;
    class Edge *matching_edge;
    int in_forest;
    list <Vertex *>::iterator vertices_iterator;
    multimap <double, Vertex *>::iterator unmarked_iterator;

    Vertex(string n);
};

class Edge {
   public:
     string name;
     Vertex *v1;
     Vertex *v2;
     Edge *backedge;
     int marked;
     class Blossom *blossom;
     list <Edge *>::iterator v1_iterator;
     list <Edge *>::iterator matching_iterator;
     list <Edge *>::iterator edge_list_iterator;
};

class Blossom {
   public:
     list <Vertex *> vertices;   /* The order is in the order of the cycle */
     list <Edge *> cycle;
};

typedef list <Vertex *>::iterator LVIT;

class Graph {
  public:
    map <string, Vertex *> v_by_name;
    list <Vertex *> vertices;
    list <Edge *> edges;
    list <Edge *> matching;
    map <string, Vertex *> forest;

    void Add_Edge(string v1, string v2);
    int Find_Augmenting_Path();
    void Add_Edge_To_Matching(Edge *e);  
    Edge *Create_Edge_And_Backedge(Vertex *v1, Vertex *v2);
    void Lift(Vertex *vb);
    void Find_Matching();
};

void Graph::Add_Edge(string vv1, string vv2)
{
  Vertex *v1, *v2;
  map <string, Vertex *>::iterator vit;

  vit = v_by_name.find(vv1);
  if (vit == v_by_name.end()) {
    v1 = new Vertex(vv1);
    vertices.push_front(v1);
    v1->vertices_iterator = vertices.begin();
    v_by_name[vv1] = v1;
  } else {
    v1 = vit->second;
  }

  vit = v_by_name.find(vv2);
  if (vit == v_by_name.end()) {
    v2 = new Vertex(vv2);
    vertices.push_front(v2);
    v2->vertices_iterator = vertices.begin();
    v_by_name[vv2] = v2;
  } else {
    v2 = vit->second;
  }
  Create_Edge_And_Backedge(v1, v2);
}

Edge *Create_Edge(Vertex *v1, Vertex *v2)
{
  Edge *e;

  e = new Edge;
  e->v1 = v1;
  e->v2 = v2;
  e->blossom = NULL;
  e->name = "[" + v1->name + ":" + v2->name + "]";
  v1->edges.push_front(e);
  e->v1_iterator = v1->edges.begin();
  return e;
}

Edge *Graph::Create_Edge_And_Backedge(Vertex *v1, Vertex *v2)
{
  Edge *e1, *e2;

  e1 = Create_Edge(v1, v2);
  e2 = Create_Edge(v2, v1);
  e1->backedge = e2;
  e2->backedge = e1;
  edges.push_front(e1);
  e1->edge_list_iterator = edges.begin();
  edges.push_front(e2);
  e2->edge_list_iterator = edges.begin();
  return e1;
}

Vertex::Vertex(string n)
{
  name = n;
  exposed = 1;
  blossom = NULL;
  in_blossom = NULL;
  tmp = NULL;
}

Vertex *Find_Head_node(Vertex *v1, Vertex*v2)
{
  while(v1->d_to_root > v2->d_to_root) v1 = v1->edge_to_root->v2;
  while(v2->d_to_root > v1->d_to_root) v2 = v2->edge_to_root->v2;
  while (v1 != v2) {
    v1 = v1->edge_to_root->v2;
    v2 = v2->edge_to_root->v2;
  }
  return v1;
}

void Graph::Add_Edge_To_Matching(Edge *e)
{
  Edge *be;

  // cout << "Add_Edge_To_Matching(" << e->name << ")\n";
  matching.push_front(e);
  e->matching_iterator = matching.begin();
  be = e->backedge;
  matching.push_front(be);
  be->matching_iterator = matching.begin();
  e->v1->matching_edge = e;
  be->v1->matching_edge = be;

  if (e->v1->root == e->v1) {
    e->v1->exposed = 0;
  } else {
    matching.erase(e->v1->edge_to_root->matching_iterator);
    matching.erase(e->v1->edge_to_root->backedge->matching_iterator);
    Add_Edge_To_Matching(e->v1->edge_to_root->v2->edge_to_root->backedge);
  }

  if (e->v1->root != e->v2->root) { /* This code only gets executed on the initial call, not the recursive ones.  */
    if (e->v2->root == e->v2) {
      e->v2->exposed = 0;
    } else {
      matching.erase(e->v2->edge_to_root->matching_iterator);
      matching.erase(e->v2->edge_to_root->backedge->matching_iterator);
      Add_Edge_To_Matching(e->v2->edge_to_root->v2->edge_to_root->backedge);
    }
  }

  /* When we return, the augmenting path is done, so don't worry about 
     roots or edges to roots right now.  Fix that in a later iteration. */
}

int Graph::Find_Augmenting_Path()
{
  int retval;
  Vertex *v1, *v2, *vhead, *vb, *v;
  Blossom *b;
  Edge *e1, *e, *be1;
  list <Edge *> tmplist;
  list <Vertex *>::iterator vit;
  list <Edge *>::iterator mit;
  list <Edge *>::iterator eit;
  list <Edge *>::iterator v1eit;
  multimap <double, Vertex *> unmarked_v;
  multimap <double, Vertex *>::iterator uvit;

  // cout << "Finding Augmenting Path\n";
  // Print_Matching();

  forest.clear();

  for (vit = vertices.begin(); vit != vertices.end(); vit++) {
    v1 = *vit;
    v1->tmp = NULL;
    v1->in_blossom = NULL;
    if (v1->exposed) {
      v1->in_forest = 1;
      v1->d_to_root = 0;
      v1->root = v1;
      v1->edge_to_root = NULL;
      v1->marked = 0;
      v1->unmarked_iterator = unmarked_v.insert(make_pair(static_cast<double>(std::rand()) / (RAND_MAX + 1.0), v1));
      forest[v1->name] = v1;
    } else {
      v1->in_forest = 0;
    }
  }


  for (mit = edges.begin(); mit != edges.end(); mit++) {
    e1 = *mit;
    e1->marked = 0;
    e1->blossom = NULL;
  }

  for (mit = matching.begin(); mit != matching.end(); mit++) {
    e1 = *mit;
    e1->marked = 1;
  }

  /* 1. Find an unmarked vertex v in F s.t. distance v->root is even */

  while (1) {
    if (unmarked_v.empty()) return 0;
    uvit = unmarked_v.begin();
    v1 = uvit->second;
 
    // cout << "Considering unmarked vertex " << v1->name << " d_to_root: " << v1->d_to_root << endl;
    if (v1->d_to_root%2 == 0) {

    /* Traverse v1's edges to find unmarked edges */
    
      for (v1eit = v1->edges.begin(); v1eit != v1->edges.end(); v1eit++) {
        e1 = *v1eit;
        if (e1->marked == 0) {
          // cout << "Considering unmarked edge " << e1->name << endl;

          /* Case 3: Has an important by-product -- this will add a vertex to the unmarked vertices pool.  */
          if (!e1->v2->in_forest) {
            // cout << "Case 3 -- Adding " << e1->v2->name << "(" << e1->v2->exposed << ") & " << e1->v2->matching_edge->v2->name << " to " << e1->v1->name << "'s tree.\n";
   

            e1->v2->root = v1->root;
            e1->v2->d_to_root = v1->d_to_root + 1;
            e1->v2->edge_to_root = e1->backedge;
            e1->v2->in_forest = 1;
            
            v2 = e1->v2->matching_edge->v2;
            v2->root = v1->root;
            v2->d_to_root = v1->d_to_root + 2;
            v2->edge_to_root = v2->matching_edge;
            v2->in_forest = 1;
            v2->unmarked_iterator = unmarked_v.insert(make_pair(static_cast<double>(std::rand()) / (RAND_MAX + 1.0), v2));

            e1->marked = 1;
            e1->backedge->marked = 1;
      
          /* Case 4 */
          } else if (e1->v2->d_to_root % 2 == 1) {
            // cout << "Case 4 -- Mark edge " << e1->name << endl;
            e1->marked = 1;
            e1->backedge->marked = 1;

          /* Case 5: We have found an augmenting path from v1->root to v2->root.
             We need to remove edges along that path that are currently matching, and 
             remove the others.  This will be done with a recursive procedure called
             add_edge_to_matching. */
          } else if (v1->root != e1->v2->root) { 
            // cout << "Case 5 -- Adding Edge " << e1->name << " to Matching\n";
            Add_Edge_To_Matching(e1);
            return 1;
          
          } else {
            // cout << "Case 6 -- Discovered a blossom\n";
            v2 = e1->v2;
            vb = new Vertex("*-"+v1->name);
            vertices.push_front(vb);
            vb->vertices_iterator = vertices.begin();
            b = new Blossom;
            vb->blossom = b;
            vhead = Find_Head_node(v1, v2); 

            /* Create the cycle on b's edge list */
            
            b->cycle.push_back(e1->backedge);
            for (v = v1; v != vhead; v = v->edge_to_root->v2) {
              b->cycle.push_back(v->edge_to_root);
            } 
            for (v = v2; v != vhead; v = v->edge_to_root->v2) {
              tmplist.push_front(v->edge_to_root->backedge);
            } 
            for (eit = tmplist.begin(); eit != tmplist.end(); eit++) {
              b->cycle.push_back(*eit);
            }

            /* Create the vertex list in the same order as the cycle */
            for (eit = b->cycle.begin(); eit != b->cycle.end(); eit++) {
              e = *eit;
              b->vertices.push_back(e->v1);
            }
            tmplist.clear();

            /* Get rid of all matching edges for the blossom */

            for (vit = b->vertices.begin(); vit != b->vertices.end(); vit++) {
              v = *vit;
              v->in_blossom = b;
              if (!v->exposed) {
                matching.erase(v->matching_edge->matching_iterator);
                v->exposed = 1;
              }
            } 

            /* Now, traverse the vertices, and delete each edge from g.edges.
               If an edge goes from inside the blossom to outside the blossom,
               then if you should, create the new edges to/from vb.  Manage
               matching edges properly. */

            for (vit = b->vertices.begin(); vit != b->vertices.end(); vit++) {
              v = *vit;
              vertices.erase(v->vertices_iterator);
              for (eit = v->edges.begin(); eit != v->edges.end(); eit++) {
                e = *eit;
                edges.erase(e->edge_list_iterator);
                e->blossom = b;
                if (e->v2->in_blossom != b) {
                  if (e->v2->tmp != b) {
                    be1 = Create_Edge_And_Backedge(vb, e->v2);
                    e->v2->tmp = b;
                  }
                  if (!e->v2->exposed && e->v2->matching_edge == e->backedge) {
                    be1 = *(e->v2->edges.begin());
                    be1 = be1->backedge;
                    matching.erase(e->backedge->matching_iterator);
                    matching.push_front(be1);
                    be1->matching_iterator = matching.begin();
                    matching.push_front(be1->backedge);
                    be1->backedge->matching_iterator = matching.begin();
                    vb->matching_edge = be1;
                    e->v2->matching_edge = be1->backedge;
                    vb->exposed = 0;
                  }
                  e->v2->edges.erase(e->backedge->v1_iterator);
                  edges.erase(e->backedge->edge_list_iterator);
                }
              }
            } 
            retval = Find_Augmenting_Path();
            Lift(vb);
            return retval;
          }
        }
      }
    }

    unmarked_v.erase(uvit);
    v1->marked = 1;
  }
}

void Graph::Lift(Vertex *vb)
{
  list <Vertex *>::iterator vit;
  list <Edge *>::iterator eit;
  Blossom *b;
  Vertex *connecting_vertex;
  Vertex *v;
  Edge *e;
  int do_match;

  // cout << "Lifting " << vb->name << endl;
  // Print_Graph();

  connecting_vertex = NULL;

  b = vb->blossom;
  if (b == NULL) return;

  /* Restore all vertices and edges to their rightful places. */

  for (vit = b->vertices.begin(); vit != b->vertices.end(); vit++) {
    v = *vit;
    vertices.push_front(v);
    v->vertices_iterator = vertices.begin();
  }

  for (vit = b->vertices.begin(); vit != b->vertices.end(); vit++) {
    v = *vit;
    for (eit = v->edges.begin(); eit != v->edges.end(); eit++) {
      e = *eit;
      edges.push_front(e);
      e->edge_list_iterator = edges.begin();
      e->blossom = NULL;
      if (e->v2->in_blossom != b) {
        e->v2->edges.push_front(e->backedge);
        e->backedge->v1_iterator = e->v2->edges.begin();
        edges.push_front(e->backedge);
        e->backedge->edge_list_iterator = edges.begin();
        if (!vb->exposed && e->v2 == vb->matching_edge->v2) {
          matching.erase(vb->matching_edge->matching_iterator);
          matching.erase(e->v2->matching_edge->matching_iterator);
          vb->exposed = 1;
          matching.push_front(e);
          e->matching_iterator = matching.begin();
          e->v1->matching_edge = e;
          e->v1->exposed = 0;
          connecting_vertex = e->v1;
          matching.push_front(e->backedge);
          e->backedge->matching_iterator = matching.begin();
          e->v2->matching_edge = e->backedge;
        } 
      }
    }
  }

  for (vit = b->vertices.begin(); vit != b->vertices.end(); vit++) {
    v = *vit;
    v->in_blossom = NULL;
    v->tmp = NULL;
  }

  vertices.erase(vb->vertices_iterator);
  for (eit = vb->edges.begin(); eit != vb->edges.end(); eit++) {
    e = *eit;
    e->v2->edges.erase(e->backedge->v1_iterator);
    edges.erase(e->edge_list_iterator);
    edges.erase(e->backedge->edge_list_iterator);
  }

  if (connecting_vertex == NULL) {
    eit = b->cycle.begin();
    connecting_vertex = (*eit)->v1;
  } else {
    for (eit = b->cycle.begin(); (*eit)->v1 != connecting_vertex; eit++) ;
  }

  e = *eit;
  do_match = 0;
  while (e->v2 != connecting_vertex) {
    if (do_match) {
      matching.push_front(e);
      e->matching_iterator = matching.begin();
      e->v1->exposed = 0;
      e->v1->matching_edge = e;

      matching.push_front(e->backedge);
      e->backedge->matching_iterator = matching.begin();
      e->v2->exposed = 0;
      e->v2->matching_edge = e->backedge;
    }
    do_match = 1-do_match;
    eit++;
    if (eit == b->cycle.end()) eit = b->cycle.begin();
    e = *eit;
  }

  delete b;

  for (eit = vb->edges.begin(); eit != vb->edges.end(); eit++) {
    e = *eit;
    delete e->backedge;
    delete e;
  }
  delete vb;
  
}

void Graph::Find_Matching()
{
  while (Find_Augmenting_Path()) {
  }
}

/* END MATCHING CODE */

void Schedule::SetKey()
{
  int i;
  SIMap sm;
  SIMit smit;

  for (i = cols; i < sched.size(); i++) {
    sm[sched[i]->Key]++;
  }
  Key = "";
  IT(smit, sm) Key += smit->first;
}


Schedule::~Schedule()
{
  int i;
  SMMit sit;
  MXS *mx;

  IT(sit, sets) {
    mx = sit->second;
    delete mx;
  }
  VIT(i, sched) delete sched[i];
  VIT(i, targets) delete targets[i];
}

void Schedule::Copy(Schedule *S)
{
  int i;
  Bitset *empty;
  SMMit sit;
  MXS *mx, *newmx;
  XMit xit;
  IBMit mit;

  cols = S->cols;
  empty = new Bitset(cols);
  VIT(i, S->targets) targets.push_back(new Bitset(empty, S->targets[i]));

  VIT(i, S->sched) {     /* Sched will be their own pointers */
    sched.push_back(new Bitset(empty, S->sched[i]));
    schedtree[sched[i]->Key] = i;
  }
  Sched_Edge_Weights = S->Sched_Edge_Weights;
  Edges_By_Weight = S->Edges_By_Weight;
  Sched_By_Bits = S->Sched_By_Bits;

  IT(sit, S->sets) {
    mx = sit->second;
    newmx = new MXS;
    newmx->ind = mx->ind;
    newmx->t = targets[mx->ind];
    newmx->sz = mx->sz;
    newmx->XSets = mx->XSets;
    IT(xit, newmx->XSets) {
      IT(mit, xit->second.indices) {
        mit->second = sched[mit->first];
      }
    }
    newmx->MX_Global_Edge_Weights = mx->MX_Global_Edge_Weights;
    newmx->MX_Local_Edge_Weights = mx->MX_Local_Edge_Weights;
    sets[newmx->t->Key] = newmx;
  }
  delete empty;
}


void MXS::Print()
{
  int i;
  IBMit mit;
  SPMit psit;
  XMit xptr;
  string istring;

  istring = "Target: " + t->Key + " mx->sz: " + ITOS(sz) + " ";
  IT(xptr, XSets) {
    cout << istring;
    i = istring.size();
    istring.clear();
    istring.resize(i, ' ');
    cout << xptr->first << "(";
    IT(mit, xptr->second.indices) {
      if (mit != xptr->second.indices.begin()) cout << ",";
      cout << mit->second->Key;
    }
    cout << "){";
    IT(psit, xptr->second.sums) {
      if (psit != xptr->second.sums.begin()) cout << ",";
      cout << psit->first;
    }
    cout << "}\n";
  }
}

void Schedule::Check()
{
  SMMit sit;
  MXS *mx;
  XMit xptr;
  int i, ones;
  IBMit indit;
  SIMap potcounter, gpotcounter;
  SIMit pit;
  SPMit sumit;

  // Print();
  IT(sit, sets) {
    /* First check -- is there an XSet with the right size? */
    mx = sit->second;
    i = 0;
    IT(xptr, mx->XSets) {
      if (mx->sz == xptr->second.indices.size()) i++;
    }
    if (i == 0) {
      printf("Problem with MXS %s.  No xset of size %d\n", mx->t->Key.c_str(), mx->sz);
      mx->Print();
      exit(1);
    }

    /* Second check -- Are the keys correct? */

    IT(xptr, mx->XSets) {
      if (xptr->first != xptr->second.id) {
        printf("Problem MXS %s.  xptr->first %s doesn't equal xptr->second.id\n",
           mx->t->Key.c_str(), xptr->first.c_str(), xptr->second.id.c_str());
        mx->Print();
      }
      ones = 0;
      VIT(i, xptr->second.id) if (xptr->second.id[i] == '1') ones++;
      if (xptr->second.indices.size() != ones) {
        printf("Problem MXS %s.  XSet id %s doesn't match number of indices ",
           mx->t->Key.c_str(), xptr->second.id.c_str());
           cout << xptr->second.indices.size() << endl;
        mx->Print(); exit(1);
      }
      IT(indit, xptr->second.indices) {
        if (xptr->second.id[indit->first] != '1') {
          printf("Problem MXS %s.  XSet id %s bit not set at index %d\n",
           mx->t->Key.c_str(), xptr->second.id.c_str(), indit->first);
          mx->Print(); exit(1);
        }
      }
    }

    /* Third check -- Are the sum keys list correct? */

    IT(xptr, mx->XSets) {
      IT(sumit, xptr->second.sums) {
        VIT(i, sumit->first) if (sumit->first[i] != '0' && sumit->first[i] != '1') {
          printf("Problem MXS %s - %s: Bad sum key: %s\n", mx->t->Key.c_str(), xptr->first.c_str(), sumit->first.c_str());
          mx->Print(); exit(1);
        }
      }
    }
    /* Fourth check -- Check potentials against the mx edge weight lists */

    potcounter.clear();
    gpotcounter.clear();
    IT(xptr, mx->XSets) {
      IT(sumit, xptr->second.sums) {
        potcounter[sumit->first]++;
        if (xptr->second.indices.size() == mx->sz) gpotcounter[sumit->first]++;
      }
    }
    IT(pit, potcounter) {
      if (pit->second != mx->MX_Local_Edge_Weights[pit->first]) {
        printf("Problem MXS %s.  MX_Local_Edge_Weights is off\n", mx->t->Key.c_str());
        mx->Print(); 
        printf("Sum, Calculated, Stored\n");
        IT(pit, potcounter) O3(pit->first, pit->second, mx->MX_Local_Edge_Weights[pit->first]);
        exit(1);
      }
    }
    IT(pit, gpotcounter) {
      if (pit->second != mx->MX_Global_Edge_Weights[pit->first]) {
        printf("Problem MXS %s.  MX_Global_Edge_Weights is off\n", mx->t->Key.c_str());
        mx->Print(); exit(1);
      }
    }

  }
}

void Schedule::Print()
{
  int i;
  SIMXMit pit;
  SMMit mit;

  if (!sets.empty()) cout << "------------------------------\n";
  VIT(i, sched) sched[i]->Print(cols, 1);
  if (!sets.empty()) printf("MXSets:\n");
  IT(mit, sets) mit->second->Print();
  if (!Sched_Edge_Weights.empty()) printf("Potentials:\n");
  IT(pit, Sched_Edge_Weights) O2(pit->first, pit->second.size());
  for (i = 1; i < Edges_By_Weight.size(); i++) {
    printf("Edges_By_Weight[%d].size() = %lu\n", i, (unsigned long) Edges_By_Weight[i].size());
  }
  cout.flush();
}

void Schedule::PrintMathematical()
{
  int i, j;
  
  printf("Schedule:\n");
  
  // Print initial data elements (s0 = x0, s1 = x1, etc.)
  for (i = 0; i < cols; i++) {
    printf("s%d = x%d\n", i, i);
  }
  
  // Print computed intermediate elements
  for (i = 0; i < computation_history.size(); i++) {
    XSetComputationStep& step = computation_history[i];
    
    if (step.operand1_index == -1) {
      // This shouldn't happen for computed elements
      continue;
    } else if (step.operand2_index == -1) {
      // Single operand operation (shouldn't occur in XOR, but just in case)
      printf("s%d = s%d\n", step.result_index, step.operand1_index);
    } else {
      // Two operand XOR operation
      printf("s%d = s%d + s%d\n", step.result_index, step.operand1_index, step.operand2_index);
    }
  }
  
  // Print target assignments (bi = si)
  for (i = 0; i < computation_history.size(); i++) {
    XSetComputationStep& step = computation_history[i];
    if (step.is_target && step.target_index >= 0) {
      printf("b%d = s%d\n", step.target_index, step.result_index);
    }
  }
  
  printf("\nTotal XOR operations: %lu\n", computation_history.size());
}

void Schedule::Remove_Global_Edge_Weights(MXS *mx, XMit xptr)
{
  SPMit ppit;
  string s;

  IT(ppit, xptr->second.sums) { s = ppit->first; Subtract_Global_Edge(mx, s); }
}

void Schedule::Add_Global_Edge_Weights(MXS *mx, XMit xptr)
{
  SPMit ppit;
  string s;

  IT(ppit, xptr->second.sums) { s = ppit->first; Add_Global_Edge(mx, s); }
}

void Schedule::Add_Element_To_MX_Set(MXS *mx, XMit xptr, int index, int global_weight)
{
  Bitset *b;
  IBMit spit;
  IPair ip;

  IT(spit, xptr->second.indices) {
    b = new Bitset(spit->second, sched[index]); 
    if (schedtree.find(b->Key) == schedtree.end()) {
      if (xptr->second.sums.find(b->Key) != xptr->second.sums.end()) {
        /* If this happens, you have issues, because the sum of four elements equals zero. */
      } else {
        ip.first = index;
        ip.second = spit->first;
        xptr->second.sums.insert(make_pair(b->Key, ip));
        mx->MX_Local_Edge_Weights[b->Key]++;
        if (global_weight) Add_Global_Edge(mx, b->Key);
      }
    }
    delete b;
  }
  xptr->second.indices.insert(make_pair(index, sched[index]));
}

void Schedule::Add_Global_Edge(MXS *mx, string &s)
{
  int ew;

  mx->MX_Global_Edge_Weights[s]++;
  if (mx->MX_Global_Edge_Weights[s] == 1) {
    ew = Sched_Edge_Weights[s].size();
    if (ew != 0) Edges_By_Weight[ew].erase(Edges_By_Weight[ew].find(s));
    Sched_Edge_Weights[s].insert(make_pair(mx->ind, mx));
    ew++;
    if (Edges_By_Weight.size() <= ew) Edges_By_Weight.resize(ew+1);
    Edges_By_Weight[ew].insert(s);
  }
}

/* It's ok to call this on a set that doesn't have the element, although I don't. */

void Schedule::Subtract_Element_From_MX_Set(MXS *mx, XMit xptr, int index, int global_weight)
{
  Bitset *b;
  IBMit spit;
  SPMit pit;
  SIMit mpit;

  spit = xptr->second.indices.find(index);
  if (spit == xptr->second.indices.end()) return;               
  xptr->second.indices.erase(spit);

  IT(spit, xptr->second.indices) {
    b = new Bitset(spit->second, sched[index]); 
    pit = xptr->second.sums.find(b->Key);
    if (pit != xptr->second.sums.end()) {
      xptr->second.sums.erase(pit);                  /* Erase from sum's edge set */
      mpit = mx->MX_Local_Edge_Weights.find(b->Key);   /* Decrement from mx->MX_Global_Edge_Weights */
      mpit->second--;
      if (mpit->second == 0) { 
        mx->MX_Local_Edge_Weights.erase(mpit);
      }
      if (global_weight) Subtract_Global_Edge(mx, b->Key);
    }
    delete b;
  }
}

void Schedule::Subtract_Global_Edge(MXS *mx, string &s)
{
  SIMit mpit;
  SIMXMit eit;
  int ew;

  mpit = mx->MX_Global_Edge_Weights.find(s);   /* Decrement from mx->MX_Global_Edge_Weights */
  mpit->second--;
  if (mpit->second == 0) {             /* If it's no longer in mx->MX_Global_Edge_Weights, then modifify schedule.Sched_Edge_Weights */
    mx->MX_Global_Edge_Weights.erase(mpit);
    eit = Sched_Edge_Weights.find(s);
    ew = eit->second.size();
    Edges_By_Weight[ew].erase(Edges_By_Weight[ew].find(s));
    eit->second.erase(eit->second.find(mx->ind));
    ew--;
    if (eit->second.size() == 0) {
      Sched_Edge_Weights.erase(eit);
    } else {
      Edges_By_Weight[ew].insert(s);
    }
  }
}

void Schedule::Add_New_MX(MXS *mx, IVec &inds)
{
  XSet xs;
  XMit xptr;
  int i;
  int maxid;

  /* If the new set reduces the size, then you need to go through all the old
     sets of min size and remove their edge weights from the global edge weights. */

  if (inds.size() < mx->sz) {
    printf("Crazy -- adding additional MX set that's smaller!\n");
    IT(xptr, mx->XSets) {
      if (xptr->second.indices.size() == mx->sz) Remove_Global_Edge_Weights(mx, xptr);
    }
    mx->sz = inds.size();
  }

  /* Construct the new id */

  maxid = -1;
  VIT(i, inds) if (inds[i] > maxid) maxid = inds[i];
  xs.id.resize(maxid+1, '0');
  VIT(i, inds) xs.id[inds[i]] = '1';
  
  if (mx->XSets.find(xs.id) == mx->XSets.end()) {
    mx->XSets[xs.id] = xs;
    xptr = mx->XSets.find(xs.id);
    VIT(i, inds) Add_Element_To_MX_Set(mx, xptr, inds[i], (inds.size() == mx->sz));
  }
}

void Bitset::Print(int cols, int newline)
{
  int i;

  FUP(i, cols) cout << (BMS_ISSET(i)) ? '1' : '0';
  if (newline) cout << endl;
}
  
void Bitset::SetKey()
{
  int i, j;
  
  if (Key.size() == 0) FUP(i, elts.size()*8) Key.push_back((BMS_ISSET(i)) ? '1' : '0');
}
  
Bitset::Bitset(int c)
{
  int sz;

  sz = c/8;
  if (c%8 != 0) sz++;
  elts.resize(sz, 0);
}

Bitset::Bitset(string &s, int c)
{
  int sz, i;

  sz = c/8;
  if (c%8 != 0) sz++;
  elts.resize(sz, 0);
  FUP(i, c) if (s[i] == '1') BMS_SET(i);
  SetKey();
}

Bitset::Bitset(Bitset *b1, Bitset *b2)
{
  int i;

  elts.resize(b1->elts.size());
  VIT(i, b1->elts) elts[i] = (b1->elts[i] ^ b2->elts[i]);
  SetKey();
}

void Bitset::XOR_With(Bitset *b)
{
  int i;
  string s;

  VIT(i, elts) elts[i] ^= b->elts[i];
  FUP(i, elts.size()*8) s.push_back((BMS_ISSET(i)) ? '1' : '0');
  Key = s;
}

int Bitset::Diff(Bitset *b)
{
  int i, j, d, x;

  d = 0;
  VIT(i, elts) {
    x = elts[i] ^ b->elts[i];
    FUP(j, 8) if (x & (1 << j)) d++;
  }
  return d;
}

int Bitset::Ones()
{
  int i, j, d;

  d = 0;
  VIT(i, elts) FUP(j, 8) if (elts[i] & (1 << j)) d++;
  return d;
}

IVec Bitset::Diff_Inds(Bitset *b)
{
  int i, j, d, x;
  IVec rv;

  d = 0;
  VIT(i, elts) {
    x = elts[i] ^ b->elts[i];
    FUP(j, 8) if (x & (1 << j)) rv.push_back(i*8+j);
  }
  return rv;
}

int Bitset::Is_Equal_To(Bitset *b)
{
  int i;

  VIT(i, elts) if (elts[i] != b->elts[i]) return 0;
  return 1;
}

void Schedule::Read()
{
  int tc;
  string s;
  int i;
  Bitset *b;

  cols = -1;
  while (getline(cin, s)) {
    tc = 0;
    VIT(i, s) {
      if (s[i] == '0' || s[i] == '1') {
        tc++;
      } else if (!isspace(s[i])) {
        cerr << "Bad Line -- only accept 0's, 1's and spaces: " << s << endl;
        exit(1);
      }
    }
    if (tc > 0) {
      if (cols == -1) {
        cols = tc;
      } else if (tc != cols) {
         cerr << "All rows must have the same number of columns\n";
         exit(1);
      }
      b = new Bitset(cols);
      tc = 0;
      VIT(i, s) {
        if (s[i] == '0' || s[i] == '1') {
          if (s[i] == '1') BMP_SET(b, tc);
          tc++;
        }
      }
      b->SetKey();
      targets.push_back(b);
    }
  }

  FUP(i, cols) {
    b = new Bitset(cols);
    BMP_SET(b, i);
    b->SetKey();
    sched.push_back(b);
    schedtree[b->Key] = sched.size()-1;
  }

  Sched_By_Bits.resize(cols+1);
}

Bitset *Schedule::Find_Two()
{
  XMit xptr;
  SMMit sit;
  IBMit bit;
  Bitset *b1, *b2;
  MXS *mx;
  int idx1, idx2;

  IT(sit, sets) {
    mx = sit->second;
    if (mx->sz == 2) {
      IT(xptr, mx->XSets) if (xptr->second.indices.size() == 2) {
        bit = xptr->second.indices.begin();
        idx1 = bit->first;
        b1 = bit->second;
        bit++;
        idx2 = bit->first;
        b2 = bit->second;
        Bitset *result = new Bitset(b1, b2);
        
        // Track the computation history
        XSetComputationStep step;
        step.result_index = sched.size();
        step.operand1_index = idx1;
        step.operand2_index = idx2;
        step.is_target = false;
        step.target_index = -1;
        
        // Check if this completes a target
        for (int k = 0; k < targets.size(); k++) {
          if (result->Is_Equal_To(targets[k])) {
            step.is_target = true;
            step.target_index = k;
            break;
          }
        }
        
        computation_history.push_back(step);
        return result;
      }
      printf("Internal error -- mx->sz = 2, but no sets have a size of 2\n");
      exit(1);
    }
  }
  return NULL;
}
  
typedef multimap <int, string> ISMmap;
typedef ISMmap::iterator ISMmit;

vector <Bitset *>Schedule::GetBest(int maxelts, Schedule *me)
{
  SSet::iterator eit;
  SMMit sit;
  IMXMit ewit;
  SIMXMit sewit;
  int i, top, done, minxsz, sum, newtop, newn, maxn, maxtop;
  vector <Bitset *> retval;
  string s, key;
  MXS *mx;
  Bitset *b;
  Schedule *news;
  SVec pots;
  Graph g;
  XMit xptr;
  SPMit spmit;
  SSMap edges;
  SSMit edit;
  Edge *e;

  top = Edges_By_Weight.size()-1;

  /* MW is simple -- just pick maxelts with the highest edge weight */

  if (Technique == MW) {
    IT(eit, Edges_By_Weight[top]) {
      s = *eit;
      retval.push_back(new Bitset(s, cols));
      if (retval.size() == maxelts) return retval;
    }
    return retval;
  }

  if (Technique == MW_MATCHING || Technique == SUBEX) {
    IT(eit, Edges_By_Weight[top]) {
      s = *eit;
      if (top == 1) {
        retval.push_back(new Bitset(s, cols));
        return retval;
      }
      sewit = Sched_Edge_Weights.find(s);
      IT(ewit, sewit->second) {
        mx = ewit->second;
        IT(xptr, mx->XSets) {
          spmit = xptr->second.sums.find(s);
          if (spmit != xptr->second.sums.end()) {
            key = ITOS(spmit->second.first) + " " + ITOS(spmit->second.second);
            if (edges.find(key) == edges.end()) {
              g.Add_Edge(ITOS(spmit->second.first), ITOS(spmit->second.second));
              edges[key] = s;
              key = ITOS(spmit->second.second) + " " + ITOS(spmit->second.first);
              edges[key] = s;
            }
          }
        }
      }
    }
    g.Find_Matching();
    while (!g.matching.empty()) {
      e = *(g.matching.begin());
      key = e->v1->name + " " + e->v2->name;
      edit = edges.find(key);
      if (edit != edges.end()) {
        s = edit->second;
        edges.erase(edit);
        key = e->v2->name + " " + e->v1->name;
        edit = edges.find(key);
        edges.erase(edit);
        retval.push_back(new Bitset(s, cols));
        if (retval.size() == maxelts) return retval;
      }
      g.matching.erase(g.matching.begin());
    }
    return retval;
  }

  if (Technique == UBER_XSET) {
    minxsz = cols+1;
    IT(sit, sets) if (sit->second->sz < minxsz) minxsz = sit->second->sz;
    for (i = top; i > 0; i--) {
      IT(eit, Edges_By_Weight[i]) {
        sewit = Sched_Edge_Weights.find(*eit);
        IT(ewit, sewit->second) {
          mx = ewit->second;
          if (mx->sz == minxsz) {
            s = *eit;
            retval.push_back(new Bitset(s, cols));
            if (retval.size() == maxelts) return retval;
          }
        }
      }
      if (retval.size() > 0) return retval;
    }
    fprintf(stderr, "Error in UBER_XSET - fell out of for loop\n");
    exit(1);
  }

  /* Traverse the best edge weight and find edges whose MX sets are closest to being done.  
     Put maxelts of them on the list */
 
  if (Technique == MW_SS) {
    minxsz = cols+1;
    IT(eit, Edges_By_Weight[top]) {
      IT(ewit, Sched_Edge_Weights[*eit]) {
        mx = ewit->second;
        if (mx->sz < minxsz) minxsz = mx->sz;
      }
    }
    IT(eit, Edges_By_Weight[top]) {
      done = 0;
      IT(ewit, Sched_Edge_Weights[*eit]) {
        if (!done) {
          mx = ewit->second;
          if (mx->sz == minxsz) {
            s = *eit;
            retval.push_back(new Bitset(s, cols));
            if (retval.size() == maxelts) return retval;
            done = 1;
          }
        }
      }
    }
    return retval;
  }

  /* This one is heavyweight -- for each potential new value, go ahead and process it and
     count up the sums of the MXS's.  Choose the best. */

  
  if (Technique == MW_SMALLEST_SUM || Technique == MW_SQ) {
    minxsz = (cols+1)*sets.size();
    maxtop = -1;
    if (top == 1) {
      eit = Edges_By_Weight[top].begin();
      s = *eit;
      retval.push_back(new Bitset(s, cols));
      return retval;
    }
    IT(eit, Edges_By_Weight[top]) {
      s = *eit;
      b = new Bitset(s, cols);
      news = new Schedule;
      news->Copy(me);
      news->Add_Element(b);
      for (b = news->Find_Two(); b != NULL; b = news->Find_Two()) news->Add_Element(b);
      if (Technique == MW_SMALLEST_SUM) {
        sum = 0;
        IT(sit, sets) sum += sit->second->sz;
        if (sum < minxsz) {
          minxsz = sum;
          pots.clear();
        }
        if (sum == minxsz) pots.push_back(s);
      } else {
        newtop = news->Edges_By_Weight.size()-1;
        newn = news->Edges_By_Weight[newtop].size();
        if (newtop > maxtop || (newtop == maxtop && newn > maxn)) {
          maxtop = newtop;
          maxn = newn;
          pots.clear();
        }
        if (newtop == maxtop && newn == maxn) pots.push_back(s);
      }
      delete news;
    }
    VIT(i, pots) {
      retval.push_back(new Bitset(pots[i], cols));
      if (retval.size() == maxelts) return retval;
    }
    return retval;
  }

  return retval;
}

void Schedule::Clear_Sets_And_Potentials(MXS *mx)
{
  XMit xptr;
  int global;

  while (!mx->XSets.empty()) {
    xptr = mx->XSets.begin();
    global = (xptr->second.indices.size() == mx->sz);
    while (!xptr->second.indices.empty()) Subtract_Element_From_MX_Set(mx, xptr, xptr->second.indices.begin()->first, global);
    mx->XSets.erase(xptr);
  }

  if (!mx->MX_Global_Edge_Weights.empty() || !mx->MX_Local_Edge_Weights.empty()) {
    fprintf(stderr, "Internal error: Clear_Sets_And_Potentials called and mx->MX_Global_Edge_Weights is not empty\n");
    exit(1);
  }
}

void Schedule::Add_Element_With_History(Bitset *b, int op1_idx, int op2_idx)
{
  // Track computation history
  XSetComputationStep step;
  step.result_index = sched.size();
  step.operand1_index = op1_idx;
  step.operand2_index = op2_idx;
  step.is_target = false;
  step.target_index = -1;
  
  // Check if this completes a target
  int target_idx = -1;
  for (int k = 0; k < targets.size(); k++) {
    if (b->Is_Equal_To(targets[k])) {
      target_idx = k;
      step.is_target = true;
      step.target_index = k;
      break;
    }
  }
  
  computation_history.push_back(step);
  Add_Element(b);
}

void Schedule::Add_Element_From_GetBest(Bitset *b)
{
  // For elements from GetBest, try to find the operands
  int op1_idx = -1, op2_idx = -1;
  string s = b->Key;
  
  SIMXMit sewit = Sched_Edge_Weights.find(s);
  if (sewit != Sched_Edge_Weights.end()) {
    IMXMit ewit;
    IT(ewit, sewit->second) {
      MXS *mx = ewit->second;
      XMit xptr;
      IT(xptr, mx->XSets) {
        SPMit spmit = xptr->second.sums.find(s);
        if (spmit != xptr->second.sums.end()) {
          op1_idx = spmit->second.first;
          op2_idx = spmit->second.second;
          break;
        }
      }
      if (op1_idx != -1) break;
    }
  }
  
  // Track computation history
  XSetComputationStep step;
  step.result_index = sched.size();
  step.operand1_index = op1_idx;
  step.operand2_index = op2_idx;
  step.is_target = false;
  step.target_index = -1;
  
  // Check if this completes a target
  for (int k = 0; k < targets.size(); k++) {
    if (b->Is_Equal_To(targets[k])) {
      step.is_target = true;
      step.target_index = k;
      break;
    }
  }
  
  computation_history.push_back(step);
  Add_Element(b);
}

void Schedule::Add_Element(Bitset *b)
{
  Bitset *b2, *b3;
  SMMit sit;
  SPMit psit;
  IBMit ssit;
  MXS *mx;
  XMit xptr;
  SVec xptrs_to_modify;
  int i, id, j, newsz, gonna_be_smaller, global, ones, tones, to_test, maxid;
  IPair ip;
  string s, oldkey, newkey;
  IVec inds;
  
  // printf("Before adding %s:\n", b->Key.c_str()); Print(); Check(); 

  id = sched.size();
  sched.push_back(b);
  schedtree[b->Key] = id;
  
  /* Do a first pass to kill of any targets that are done */

  if (sets.find(b->Key) != sets.end()) {
    sit = sets.find(b->Key);
    mx = sit->second;
    Clear_Sets_And_Potentials(mx);
    delete mx;
    sets.erase(sit);
  }

  /* First see if the new product will improve any sets */
 
  IT(sit, sets) {
    mx = sit->second;
    xptrs_to_modify.clear();
    if (mx->MX_Local_Edge_Weights.find(b->Key) != mx->MX_Local_Edge_Weights.end()) {  

      if (mx->MX_Global_Edge_Weights.find(b->Key) != mx->MX_Global_Edge_Weights.end()) {
        newsz = mx->sz-1;
      } else {
        newsz = mx->sz;
      }

      /* Traverse the sets.  First, figure out if the global edges need to be modified, for good or bad.
         Then, modify the XOR set if need be. */

      IT(xptr, mx->XSets) {
        gonna_be_smaller = (xptr->second.sums.find(b->Key) != xptr->second.sums.end());
        if ((!gonna_be_smaller) && xptr->second.indices.size() == mx->sz && newsz < mx->sz) {
          Remove_Global_Edge_Weights(mx, xptr);
        } else if (gonna_be_smaller && xptr->second.indices.size() == mx->sz +1 && newsz == mx->sz) {
          Add_Global_Edge_Weights(mx, xptr);
        } 

        if (xptr->second.sums.find(b->Key) != xptr->second.sums.end()) {
          ip = xptr->second.sums.find(b->Key)->second;
          global = (xptr->second.indices.size()-1 == newsz);
          Subtract_Element_From_MX_Set(mx, xptr, ip.first, global);
          Subtract_Element_From_MX_Set(mx, xptr, ip.second, global);
          Add_Element_To_MX_Set(mx, xptr, id, global);
          xptrs_to_modify.push_back(xptr->second.id);
        }
      }

      mx->sz = newsz;
    } 

    /* If you modified any X Sets, you need to update their keys in XSets.  You may have discovered
     a duplicate.  In that case, you need to nuke the XSet, which is pretty heavyweight -- you need
     to go ahead and subtract everything from the set so that all the Edge/Sum pointers are valid. */

    VIT(i, xptrs_to_modify) {
      xptr = mx->XSets.find(xptrs_to_modify[i]);
      oldkey = xptr->second.id;
      newkey.clear();
      maxid = -1;
      IT(ssit, xptr->second.indices) if (ssit->first > maxid) maxid = ssit->first;
      newkey.resize(maxid+1, '0');
      IT(ssit, xptr->second.indices) newkey[ssit->first] = '1';
      if (mx->XSets.find(newkey) == mx->XSets.end()) {
        xptr->second.id = newkey;
        mx->XSets[newkey] = xptr->second;
      } else {
        global = (xptr->second.indices.size() == mx->sz);
        while (!xptr->second.indices.empty()) {
          Subtract_Element_From_MX_Set(mx, xptr, xptr->second.indices.begin()->first, global);
        }
      }
      mx->XSets.erase(xptr);
    }
  
  }

  /* Clean up the mx sets.  In other words, if there are any XSets in the MX sets that are too big, delete them. */

  // printf("Cleaning up MX sets\n"); Check();

  IT(sit, sets) {
    mx = sit->second;
    xptrs_to_modify.clear();
    IT(xptr, mx->XSets) {
      i = xptr->second.indices.size();
      if (i > mx->sz && i > mx->sz + Thresh) {  /* This extra test is so you can have negative Thresh's */
        xptrs_to_modify.push_back(xptr->first);
      }
    }
    VIT(i, xptrs_to_modify) {
      xptr = mx->XSets.find(xptrs_to_modify[i]);
      while (!xptr->second.indices.empty()) {
        Subtract_Element_From_MX_Set(mx, xptr, xptr->second.indices.begin()->first, 0);
      }
      mx->XSets.erase(xptr);
    }
  }
  // printf("And here\n"); Check();
  
  /* Next, see if the new node can be part of any MX sets that are equal to or smaller than the current ones.
     This is kind of redundant, because if I just did it above, then this may generate the same set.  Ha ha. */

  Initialize_Enumerator(id);
  
  while (!Enumerator_Done()) {

    /* Indexes are in Enumerator and EN_Cannot_Use.  Sentinel value is cols-1 */
    IT(sit, sets) {
      mx = sit->second;
  
      inds.clear();
      VIT(i, Enumerator) if (Enumerator[i] >= cols) inds.push_back(Enumerator[i]);
      b2 = new Bitset(mx->t, sched[inds[0]]);
      for (i = 1; i < inds.size(); i++) b2->XOR_With(sched[inds[i]]);
  
      /* Now, cruise through Sched_By_Bits and see if any currently scheduled items can reduce the size */
  
      ones = b2->Ones();
      to_test = ones*2 - 2;
      if (to_test > cols) to_test = cols;
      j = 0;
      while (to_test > 1) {
        if (Sched_By_Bits[to_test].size() >= j) {
          to_test--;
          j = 0;
        } else {
          i = Sched_By_Bits[to_test][j];
          j++;
          tones = b2->Diff(sched[i]);
          if ((!(BMP_ISSET(EN_Cannot_Use, i))) && tones < ones-1) {
            inds.push_back(i);
            ones = tones;
            b2->XOR_With(sched[i]);
            if (ones*2 - 2 < to_test) {
              to_test = ones*2-2;
              j = 0;
            }
          }
        }
      }

    /* When we get to the end, add the new MX set if need be. */

      ones += inds.size();
      if ((ones <= mx->sz + Thresh)) {
        FUP(i, cols) if (BMP_ISSET(b2, i)) inds.push_back(i);
        Add_New_MX(mx, inds);
      }
      delete b2;
//    mx->Print();
    }
    Enumerator_Next();
  }
  
  /* Put the new element into "Sched_By_Bits" */

  Sched_By_Bits[b->Ones()].push_back(id);
  
  /* Get rid of 0 entries at the end of Edges_By_Weight */

  i = Edges_By_Weight.size()-1;
  while (i > 0 && Edges_By_Weight[i].size() == 0) { Edges_By_Weight.pop_back(); i--; }

  // printf("After adding %s:\n", b->Key.c_str()); Print(); Check(); 
  // Check();

}

void Schedule::Initialize_Enumerator(int last_product)
{
  Enumerator.clear();
  EN_Cannot_Use = NULL;
  if (NStart == 0 || Edges_By_Weight.size()-1 <= Give_Up_When) return;
  EN_Cannot_Use = new Bitset(sched.size());
  BMP_SET(EN_Cannot_Use, last_product);
  Enumerator.resize(NStart, cols-1);
  Enumerator[0] = last_product;
}


int Schedule::Enumerator_Done()
{
  return (EN_Cannot_Use == NULL);
}

void Schedule::Enumerator_Next()
{
  int i;
  int li;

  if (NStart == 0) return;
  li = -1;
  for (i = 1; i < Enumerator.size(); i++) {
    if (Enumerator[i]+1 < Enumerator[i-1]) li = i;
  }
  if (li == -1) {
    delete EN_Cannot_Use;
    EN_Cannot_Use = NULL;
    return;
  }
  for (i = li; i < Enumerator.size() && Enumerator[i] >= cols; i++) BMP_CLEAR(EN_Cannot_Use, Enumerator[i]);
  Enumerator[li]++;
  BMP_SET(EN_Cannot_Use, Enumerator[li]);
  for (i = li+1; i < Enumerator.size() && Enumerator[i] >= cols; i++) Enumerator[i] = cols-1;
}

void Schedule::Create_Initial_Sets()
{
  int i, i1;
  MXS *mx;
  IVec inds;

  VIT(i, targets) {
    if (schedtree.find(targets[i]->Key) == schedtree.end()) {
      mx = new MXS;
      mx->t = targets[i];
      mx->ind = i;

      inds.clear();
      for (i1 = 0; i1 < cols; i1++) if (BMP_ISSET(mx->t, i1)) inds.push_back(i1);
      mx->sz = inds.size();
      Add_New_MX(mx, inds);
          
      sets[mx->t->Key] = mx;
    }
  }
}

void usage() 
{
    cerr << "usage: X-Sets thresh nstart give-up-when selection_technique\n"; 
    cerr << "       selection_technique = MW|MW_SS|MW_SMALLEST_SUM|MW_SQ|UBER_XSET|MW_MATCHING|SUBEX\n"; 
    exit(1); 
}

int main(int argc, char **argv)
{
  string s;
  Bitset *b;
  int i, j;
  Schedule *S;
  vector <Bitset *> to_try;

  if (argc != 5) usage();

  Thresh = atoi(argv[1]);
  NStart = atoi(argv[2]);
  Give_Up_When = atoi(argv[3]);

  S = new Schedule;

  std::srand(std::time(0));

  if (strcmp(argv[4], "MW") == 0) {
    Technique = MW;
  } else if (strcmp(argv[4], "MW_SS") == 0) {
    Technique = MW_SS;
  } else if (strcmp(argv[4], "UBER_XSET") == 0) {
    Technique = UBER_XSET;
  } else if (strcmp(argv[4], "MW_SMALLEST_SUM") == 0) {
    Technique = MW_SMALLEST_SUM;
  } else if (strcmp(argv[4], "MW_SQ") == 0) {
    Technique = MW_SQ;
  } else if (strcmp(argv[4], "MW_MATCHING") == 0) {
    Technique = MW_MATCHING;
  } else if (strcmp(argv[4], "SUBEX") == 0) {
    Technique = SUBEX;
  } else usage();

  S->Read();
  ZERO = new Bitset(S->cols);

  S->Create_Initial_Sets();

  if (Technique == SUBEX) {
    while (S->sets.size() > 0) {
      to_try = S->GetBest(S->targets.size(), S);
      VIT(i, to_try) {
        b = to_try[i];
        printf("Targets Left %2lu.  New Element ", S->sets.size()); b->Print(S->cols, 0);
        printf(" Edge Weight: %lu (Matching %d)\n", S->Sched_Edge_Weights[b->Key].size(), i);
        fflush(stdout);
        S->Add_Element(to_try[i]);
      }
    }
    S->Print();
    printf("%lu\n", S->sched.size());
    exit(0);
  }

  for (b = S->Find_Two(); b != NULL; b = S->Find_Two())  {
    printf("Targets Left %2lu.  New Element ", S->sets.size()); b->Print(S->cols, 0);
    printf(" Edge Weight: %lu (TH 4.1)\n", S->Sched_Edge_Weights[b->Key].size());
    fflush(stdout);
    S->Add_Element(b);
  }

  while (S->sets.size() > 0) {
    to_try = S->GetBest(1, S);  /* This will still work if you call it with something larger, just to see how many sets you get */
    b = to_try[0];
    printf("Targets Left %2lu.  New Element ", S->sets.size()); b->Print(S->cols, 0);
    printf(" Edge Weight: %lu\n", S->Sched_Edge_Weights[b->Key].size());
    fflush(stdout);

    VIT(i, to_try) if (i > 0) delete to_try[i]; /* This is in case you did call it with something larger */

    S->Add_Element_From_GetBest(b);
    for (b = S->Find_Two(); b != NULL; b = S->Find_Two()) {
      printf("Targets Left %2lu.  New Element ", S->sets.size()); b->Print(S->cols, 0);
      printf(" Edge Weight: %lu (TH 4.1)\n", S->Sched_Edge_Weights[b->Key].size());
      S->Add_Element(b);
      fflush(stdout);
    }
  }

  S->Print();
  printf("\n");
  S->PrintMathematical();
  // printf("%lu\n", S->sched.size());
  exit(0);

}