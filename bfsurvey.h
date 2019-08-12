/**********************/
// Date: July 10, 2018
// Name: Ruben Navarro
// Degrees of Separation
/**********************/

#ifndef _BFSURVEY_H
#define _BFSURVEY_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <deque.h>
#include <primes.h>
#include <genalg.h>
#include <list.h>
#include <vector.h>
#include <entry.h>

namespace fsu
{

  template < class G >
  class BFSurvey
  {

 public:
    typedef G                              Graph;
    typedef typename Graph::Vertex         Vertex;
    typedef typename Graph::AdjIterator    AdjIterator;
    
    BFSurvey    ( const Graph& g );
    BFSurvey    ( const Graph& g , Vertex start );
    void Search ();
    void Search ( Vertex v );
    void Reset  ();
    void Reset  ( Vertex start );
                                     
    const fsu::Vector< Vertex >& Distance () const { return distance_; }
    const fsu::Vector< Vertex >& DTime    () const { return dtime_; }
    const fsu::Vector< Vertex >& Parent   () const { return parent_; }
    const fsu::Vector< char >&   Color    () const { return color_; }

    size_t VrtxSize () const { return g_.VrtxSize(); }
    size_t EdgeSize () const { return g_.EdgeSize(); }
    
    size_t InfiniteTime     () const { return forever_; }
    size_t InfiniteDistance () const {return  infinity_; }
    Vertex NullVertex       () const { return null_; }   

    //trace support
    bool traceQue;
    void ShowQueSetup (std::ostream& os) const;
    void ShowQue      (std::ostream& os) const;

  private:
    const Graph&          g_;
    Vertex                start_;       // default is vertex 0
    size_t                time_;        // global sequencing clock
    size_t                infinity_;     // unreachable distance = 1+|E|    
    size_t                forever_;     // unreachable time = |V|
    Vertex                null_;        // undefined vertex = |V|

  public:
    fsu::Vector < Vertex >  distance_; // distance from search origin
    fsu::Vector < Vertex >  dtime_;    // discovery time
    fsu::Vector < Vertex >  parent_;   // for BFS tree
    fsu::Vector < char >    color_;    // using chars 'w'=white, 'g'=grey,

  private:
    fsu::Deque < Vertex > conQ_;        // control queue

};

  template < class G >
  BFSurvey<G>::BFSurvey (const Graph& g)
    : g_(g), start_(0), time_(0),
      infinity_(1+g_.EdgeSize()), forever_((Vertex)g_.VrtxSize()), null_((Vertex)g_.VrtxSize()),
      distance_ (g_.VrtxSize(), infinity_),
      dtime_    (g_.VrtxSize(), forever_),
      parent_   (g_.VrtxSize(), null_),
      color_    (g_.VrtxSize(), 'w'), // 'w' = white
      conQ_() 
  {}

  template < class G >
  BFSurvey<G>::BFSurvey (const Graph& g, Vertex start)
    : g_(g), start_(start), time_(0),
      infinity_ (1+g_.EdgeSize()), forever_(g_.VrtxSize()), null_((Vertex)g_.VrtxSize()),
      distance_ (g_.VrtxSize(), infinity_),
      dtime_    (g_.VrtxSize(), forever_),
      parent_   (g_.VrtxSize(), null_),
      color_    (g_.VrtxSize(), 'w'), // 'w' = white
      conQ_()
  {}
  
  template < class G >
  void BFSurvey<G>::Reset()
  {
    time_ = 0;
    conQ_.Clear();
    if (color_.Size() != g_.VrtxSize()) // g_ has changed vertex size
    {
      infinity_ = 1+g_.EdgeSize(); // unreachable distance
      forever_ = g_.VrtxSize(); // unreachable time
      null_ = (Vertex)g_.VrtxSize(); // undefined parent
      distance_.SetSize (g_.VrtxSize(), infinity_);
      dtime_.SetSize    (g_.VrtxSize(), forever_);
      parent_.SetSize   (g_.VrtxSize(), null_);
      color_.SetSize    (g_.VrtxSize(), 'w'); // w = white
    }
    else
    {
      for (Vertex x = 0; x < g_.VrtxSize(); ++x)
      {
        distance_[x] = infinity_; //unreachable distance
        dtime_[x] = forever_; //unreachable time
        parent_[x] = null_; //undefined parent
        color_[x] = 'w'; // 'w' = white
      }
    }
  }

  template < class G >
  void BFSurvey<G>::Reset ( Vertex start )
  {
    start_ = start;
    Reset();
  }        

  template < class G >
  void BFSurvey<G>::Search( Vertex v )
  {
    distance_[v] = 0;
    dtime_[v] = time_++;
    conQ_.PushBack(v);
    //  if (traceQue) { ShowQue(std::cout); }
    color_[v] = 'g'; // 'g' = grey
    Vertex front;
    AdjIterator i;
    while (!conQ_.Empty())
    {
      front = conQ_.Front();
      // add all unvisited neighbors of front to queue
      for (i = g_.Begin(front); i != g_.End(front); ++i)
      {
        if ('w' == color_[*i]) // 'w' = white = unvisited
        {
          distance_[*i] = distance_[front] + 1;
          dtime_[*i] = time_++;
          parent_[*i] = front;
          color_[*i] = 'g'; // 'g' = grey
          conQ_.PushBack(*i);
          //          if (traceQue) { ShowQue(std::cout); }
        }
      }
      // remove front of queue
      conQ_.PopFront();
      //  if (traceQue) { ShowQue(std::cout); }
      color_[front] = 'b'; // 'b' = black
    }
  }

  template < class G >
  void BFSurvey<G>::Search()
  {
    Reset();
    if (traceQue) { ShowQueSetup(std::cout);  ShowQue(std::cout); }
    for (Vertex v = start_; v < g_.VrtxSize(); ++v)
    {
      if (color_[v] == 'w') // 'w' = white = unvisited
        Search(v);
    }
    for (Vertex v = 0; v < start_; ++v)
    {
      if (color_[v] == 'w') // 'w' = white = unvisited
        Search(v);
    }
  }

  template < class G >
  void BFSurvey<G>::ShowQueSetup (std::ostream& os) const
  {
    os << "\n  conQueue\n"
    << "  <-------\n";
  }

  template < class G >
  void BFSurvey<G>::ShowQue (std::ostream& os) const
  {
    os << "  ";
    if (conQ_.Empty())
      os << "NULL";
    else
      conQ_.Display(os, ' ');
    os << '\n';
  } 
}
#endif
