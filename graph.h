/**********************/
// Date: July 10, 2018
// Name: Ruben Navarro
// Degrees of Separation
// Header file graph.h for Graph Search and Survey
/**********************/

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <entry.h>
#include <vector.h>
#include <list.h>
#include <primes.h>
#include <genalg.h>
#include <deque.h>

namespace fsu
{
  template < typename N >
  class ALUGraph
  {
  public:
    typedef N                                   Vertex;
    typedef typename fsu::List<Vertex>          SetType;
    typedef typename SetType::ConstIterator     AdjIterator;
     
    void   SetVrtxSize     (N n);
    size_t VrtxSize        () const;
    void   AddEdge         (Vertex from, Vertex to);
    bool   HasEdge         (Vertex from, Vertex to) const;
    size_t EdgeSize        () const;              // Theta (|V| + |E|)
    size_t OutDegree       (Vertex v) const;
    size_t InDegree        (Vertex v) const;
    void Shuffle           (); 

    void Clear             ();
    void Dump              (std::ostream& os);
    
    AdjIterator Begin      (Vertex x) const;  
    AdjIterator End        (Vertex x) const;

    ALUGraph               ();
    explicit ALUGraph      (N n);

  protected:
    fsu::Vector < SetType > al_;
  };

  template < typename N >
  class ALDGraph : public ALUGraph <N>
  {
  public:
    typedef N                                   Vertex;
    typedef typename ALUGraph<N>::SetType       SetType;
    typedef typename ALUGraph<N>::AdjIterator   AdjIterator;

    // void SetVrtSize          (N n);
    // size_t VrtxSize          () const;
    void AddEdge              (Vertex from, Vertex to);
    // bool HasEdge (Vertex from, Vertex to) const;
    size_t EdgeSize () const;         // Theta (|V| + |E|)
    // size_t OutDegree (Vertex v) const;
    size_t InDegree (Vertex v) const; // Theta (|V| + |E|)

    // void Clear ();
    // void Dump (std::ostream& os); // Theta (|V| + |E|)

    // AdjIterator Begin (Vertex x) const;
    // AdjIterator End   (Vertex x) const;

    ALDGraph          ();
    explicit ALDGraph (N n);

    // creates d as the reverse directed graph of *this 
    void Reverse (ALDGraph& d) const;
  };

  template < typename N >
  ALUGraph<N>::ALUGraph () 
  {}

  template < typename N >
  ALUGraph<N>::ALUGraph(N n) : al_(n) // set vector with n object
  { }

  template < typename N >
  typename ALUGraph<N>::AdjIterator ALUGraph<N>::Begin(Vertex x) const
  {
    return al_[x].Begin();  // return location where iterator begins  
  }

  template < typename N > 
  typename ALUGraph<N>::AdjIterator ALUGraph<N>::End(Vertex x) const
  {
    return al_[x].End();   // returns location where iterator ends    
  }

  template < typename N >
  void ALUGraph<N>::SetVrtxSize (N n)
  {
    al_.SetSize(n);       // Set Vertex Size by sending n into vector
  }

  template < typename N >
  size_t ALUGraph<N>::VrtxSize () const
  {
    return  al_.Size();   // returns size of vertex
  }

  template < typename N >
  void ALUGraph<N>::AddEdge(Vertex from, Vertex to)
  {
    al_[from].Insert(to);
    al_[to].Insert(from);
  }

  template < typename N >
  bool ALUGraph<N>::HasEdge (Vertex from, Vertex to) const
  {
    AdjIterator i = al_[from].Includes(to);
    if (i == End(from))
      return 0;
    return 1;
  }

  template < typename N >
  size_t ALUGraph<N>::EdgeSize () const // Theta (|V| + |E|)
  {
    size_t esize = 0;
    for (Vertex v = 0; v < al_.Size(); ++v)
      esize += al_[v].Size();
    return esize/2;
  }

  template < typename N >
  size_t ALUGraph<N>::OutDegree (Vertex v) const
  {
    return al_[v].Size();  // retuns size of row/list in vertex 
  }

  template < typename N >
  size_t ALUGraph<N>::InDegree (Vertex v) const
  {
    size_t indegree = 0;
    AdjIterator j;
    for (Vertex x = 0; x < VrtxSize(); ++x)
    {
      for (j = this->Begin(x); j != this->End(x); ++j)
      {
        if (v == *j) ++indegree;
      }
    }
    return indegree; 
  }

  template < typename N >
  void ALUGraph<N>::Clear ()
  {
    SetType list;
    list.Clear();
    al_.Clear();
  }
  
  template < typename N >
  void ALUGraph<N>::Dump (std::ostream& os) // Theta (|V| + |E|)
  { 
    AdjIterator j;
    for (Vertex v = 0; v < VrtxSize(); ++v)
    {
      os << '[' << v << "]->";
      j = this->Begin(v);
      if (j != this->End(v))
      {
        os << *j;
        ++j;
      }
      for (; j != this->End(v); ++j)
      {
        os << ',' << *j;
      }
        os << '\n';
    }
  } 

  template < typename N >
  void ALUGraph<N>::Shuffle()
  {
    for (Vertex v = 0; v < VrtxSize(); ++v)
      al_[v].Shuffle();
  }

  // ALGGraph
  template < typename N > 
  ALDGraph<N>::ALDGraph() : ALUGraph<N>::ALUGraph() // initialize to use ALUGraph 
  {}

  template < typename N >
  ALDGraph<N>::ALDGraph(N n) : ALUGraph<N>::al_(n) // initialize vector to object n 
  {} 

  template < typename N >
  void ALDGraph<N>::AddEdge (Vertex from, Vertex to)
  {
    ALUGraph<N>::al_[from].Insert(to);
  }

  template < typename N >
  size_t ALDGraph<N>::EdgeSize () const // Theta (|V| + |E|)
  {
    size_t esize = 0;
    for (Vertex v = 0; v < ALUGraph<N>::al_.Size(); ++v)
    esize += ALUGraph<N>::al_[v].Size();
    return esize;
  }

  template < typename N >
  size_t ALDGraph<N>::InDegree (Vertex v) const // Theta (|V| + |E|)
  {
    size_t indegree = 0;
    AdjIterator j;
    for (Vertex x = 0; x < ALUGraph<N>::VrtxSize(); ++x)
    {
      for (j = this->Begin(x); j != this->End(x); ++j)
      {
        if (v == *j) ++indegree;
      }
    }
    return indegree;
  }
  template < typename N >
  void ALDGraph<N>::Reverse (ALDGraph& d) const
  {
    AdjIterator i;   
    d.SetVrtxSize(ALUGraph<N>::VrtxSize());              // sets vertex size 
    for (Vertex x = 0; x < ALUGraph<N>::al_.Size(); ++x) // traverse through size of vertex
    {
      for (i = this->Begin(x); i!=this->End(x);++i)      // iterate through list 
      {
        d.AddEdge(*i,x);                                 // add edge to d
      }
    }
  }
} // end of namespace
#endif
 
