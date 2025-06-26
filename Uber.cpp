/*
Uber-CSHR Heuristic Program for minimizing XOR's in executing a binary
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
#include <cstring>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <list>
#include <set>
#include <iostream>
#include <cstdio>
#include <cstdlib>
using namespace std;

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

string ITOS(int i) { char buf[20]; string k; sprintf(buf, "%d", i); k = buf; return k; }

class Target {
  public:
    Bitset *b;
    int Min;
    IVec From;
    list <class Target *>::iterator ptr;
};

// New structure to track how each element was computed
struct ComputationStep {
    int result_index;        // Index in Sched array
    int operand1_index;      // First operand index (-1 if it's a data element)
    int operand2_index;      // Second operand index (-1 if single operand or data element)
    bool is_target;          // True if this is a target element
    int target_index;        // Index in original target list (-1 if not a target)
};

class Schedule {
  public:
    int Cols;
    int L;
    int Tori;
    list <Target *> Targets_Left;
    vector <Bitset *> Targets_Done;
    vector <Bitset *> Sched;
    IVec TD_Indices;     /* I hate doing this, but it makes Try_To_Improve easier */
    IVec Sched_Indices;   
    vector <ComputationStep> computation_history;  // Track how each element was computed
    vector <Bitset *> original_targets;            // Store original targets in order
    void Read();
    void Print();
    void PrintMathematical();
    void Do_Uber();
    void Try_To_Improve();
};

void Schedule::Print()
{
  int i, j;
  list <Target *>::iterator it;
  Target *t;

  if (!Targets_Left.empty()) {
    printf("-----------------\n");
    printf("Schedule so far:\n\n");
  }
  VIT(i, Sched) Sched[i]->Print(Cols, 1);
  if (!Targets_Left.empty()) {
    printf("\nTargets Left\n\n");
    IT(it, Targets_Left) {
      t = *it;
      t->b->Print(Cols, 0);
      printf(" %3d", t->Min);
      VIT(j, t->From) {
        cout << " " << t->From[j] << " "; Sched[t->From[j]]->Print(Cols, 0);  
      }
      printf("\n");
    }
    printf("\n");
  }
}

void Schedule::PrintMathematical()
{
  int i, j;
  
  printf("Schedule:\n");
  
  // Print initial data elements (s0 = x0, s1 = x1, etc.)
  for (i = 0; i < Cols; i++) {
    printf("s%d = x%d\n", i, i);
  }
  
  // Print computed intermediate elements
  for (i = 0; i < computation_history.size(); i++) {
    ComputationStep& step = computation_history[i];
    
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
    ComputationStep& step = computation_history[i];
    if (step.is_target && step.target_index >= 0) {
      printf("b%d = s%d\n", step.target_index, step.result_index);
    }
  }
  
  printf("\nTotal XOR operations: %lu\n", computation_history.size());
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
  int i, j;
  Bitset *b;
  Target *t;
  list <Target *>::iterator it;

  cout.flush();
  Cols = -1;
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
      if (Cols == -1) {
        Cols = tc;
      } else if (tc != Cols) {
         cerr << "All rows must have the same number of columns\n";
         exit(1);
      }
      b = new Bitset(Cols);
      tc = 0;
      VIT(i, s) {
        if (s[i] == '0' || s[i] == '1') {
          if (s[i] == '1') BMP_SET(b, tc);
          tc++;
        }
      }
      b->SetKey();
      t = new Target;
      t->b = b;
      original_targets.push_back(b);  // Store original targets in order
      Targets_Left.push_back(t);
    }
  }

  FUP(i, Cols) {
    b = new Bitset(Cols);
    BMP_SET(b, i);
    b->SetKey();
    Sched.push_back(b);
  }
 
  IT(it, Targets_Left) { 
    t = *it;
    t->ptr = it;
  }

  it = Targets_Left.begin();
  while (it != Targets_Left.end()) {
    t = *it;
    if (t->b->Ones() <= 1) {
      it++;
      Targets_Left.erase(t->ptr);
      delete t;
    } else {
      t->Min = t->b->Ones()-1;
      FUP(j, Cols) {
        if (BMP_ISSET(t->b, j)) {
          t->From.push_back(j);
          j = Cols;
        }
      }
      it++;
    }
  }
}

void Schedule::Try_To_Improve()
{
  int i, j, l, done;
  int lindex, d;
  Bitset *last, *b;
  IVec p;
  list <Target *>::iterator it;
  Target *t;
  IVec *Indices;

  if (Tori == 'T') {
    Indices = &TD_Indices;
  } else {
    Indices = &Sched_Indices;
  }

  lindex = Indices->size()-1;
  last = Sched[(*Indices)[lindex]];
  IT(it, Targets_Left) {
    t = *it;
    d = last->Diff(t->b);
    if (d < t->Min) {
      t->Min = d;
      t->From.clear();
      t->From.push_back((*Indices)[lindex]);
    }
  }

  for (l = 1; l < L && lindex >= l; l++) {
    
    p.resize(l);
    for (i = 0; i < l; i++) p[i] = i;
    b = new Bitset(last, Sched[(*Indices)[p[0]]]);
    for (i = 1; i < l; i++) b->XOR_With(Sched[(*Indices)[p[i]]]);
    do {
      IT(it, Targets_Left) {
        t = *it;
        d = b->Diff(t->b);
        if (d + l < t->Min) {
          t->Min = d + l;
          t->From.clear();
          VIT(j, p) t->From.push_back((*Indices)[p[j]]);
          t->From.push_back((*Indices)[lindex]);
        }
      }
      for (i = 0; i < l && p[l-i-1] == lindex-i-1; i++) ;
      if (i == l) {
        done = 1;
      } else {
        done = 0;
        j = l-i-1;
        b->XOR_With(Sched[(*Indices)[p[j]]]);
        p[j]++;
        b->XOR_With(Sched[(*Indices)[p[j]]]);
        for (j++; j < l; j++) {
          b->XOR_With(Sched[(*Indices)[p[j]]]);
          p[j] = p[j-1]+1;
          b->XOR_With(Sched[(*Indices)[p[j]]]);
        }
      }
    } while (!done);
    p.clear();
    delete b;
  }
}

void Schedule::Do_Uber()
{
  int minones;
  int i, j, fbs;
  Bitset *b, *b2;
  IVec tmp;
  Target *t;
  list <Target *>::iterator it;

  while (Targets_Left.size() > 0) {
    minones = Cols+1;
    IT(it, Targets_Left) {
      if ((*it)->Min < minones) { t = *it; minones = t->Min; }
    }

    if (t->From.size() == 1) { 
      b = Sched[t->From[0]];
      FUP(i, Cols) {
        if ((BMP_ISSET(b, i) && !BMP_ISSET(t->b, i)) ||
            (!BMP_ISSET(b, i) && BMP_ISSET(t->b, i))) {
          j = i;
          i = Cols;
        }
      }
      b2 = new Bitset(Cols);
      BMP_SET(b2, j);
      b2->XOR_With(b);
      Sched_Indices.push_back(Sched.size());
      
      // Track computation: xj + s[From[0]]
      ComputationStep step;
      step.result_index = Sched.size();
      step.operand1_index = j;  // Data element
      step.operand2_index = t->From[0];  // Existing schedule element
      step.is_target = false;
      step.target_index = -1;
      computation_history.push_back(step);
      
      Sched.push_back(b2);
      t->From[0] = Sched.size()-1;
    } else {
      b = new Bitset(Sched[t->From[0]], Sched[t->From[1]]);
      Sched_Indices.push_back(Sched.size());
      
      // Track computation: s[From[0]] + s[From[1]]
      ComputationStep step;
      step.result_index = Sched.size();
      step.operand1_index = t->From[0];
      step.operand2_index = t->From[1];
      step.is_target = false;
      step.target_index = -1;
      computation_history.push_back(step);
      
      Sched.push_back(b);
      t->From[0] = Sched.size()-1;
      fbs = t->From.size();
      t->From[1] = t->From[fbs-1];
      t->From.resize(fbs-1);
    }
    t->Min--;
    if (t->Min == 0) {
      // Find which original target this corresponds to
      int target_idx = -1;
      for (int k = 0; k < original_targets.size(); k++) {
        if (t->b->Is_Equal_To(original_targets[k])) {
          target_idx = k;
          break;
        }
      }
      
      // Mark the last computation as producing a target
      if (!computation_history.empty()) {
        computation_history.back().is_target = true;
        computation_history.back().target_index = target_idx;
      }
      
      Targets_Done.push_back(t->b);
      Targets_Left.erase(t->ptr);
      TD_Indices.push_back(Sched.size()-1);
      delete t;
    }
    
    if (t->Min == 0 || Tori == 'I') Try_To_Improve();
    
    // Print();
  }
}

void usage() 
{
    cerr << "usage: Uber I|T L\n"; 
    exit(1); 
}

int main(int argc, char **argv)
{
  int i, j;
  Schedule *S;

  if (argc != 3) usage();

  S = new Schedule;
  if (strcmp(argv[1], "T") != 0 && strcmp(argv[1], "I") != 0) usage();
  S->Tori = argv[1][0];
  S->L = atoi(argv[2]);

  S->Read();

  S->Do_Uber();
  S->Print();
  printf("\n");
  S->PrintMathematical();
  // printf("%lu\n", S->Sched.size());
  exit(0);
}