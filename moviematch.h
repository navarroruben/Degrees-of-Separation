/**********************/
// Date: July 10, 2018
// Name: Ruben Navarro
// Degrees of Separation
/**********************/

#ifndef _MOVIEMATCH_H
#define _MOVIEMATCH_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <xstring.h>
#include <graph.h>
#include <list.h>
#include <vector.h>
#include <graph_util.h>
#include <hashtbl.h>
#include <gheap.h>
#include <hashclasses.h>
#include <gbsearch.h>
#include <survey_util.h>
#include <bfsurvey.h>
#include <string>
 
// types used

typedef uint32_t                                Vertex;
typedef fsu::String                             Name;
typedef fsu::ALUGraph   <Vertex>                Graph;
typedef fsu::BFSurvey   <Graph>                 BFS;
typedef hashclass::KISS <Name>                  Hash;
typedef fsu::HashTable  <Name,Vertex,Hash>      AA;        // associative array 
typedef fsu::Vector     <Name>                  Vector;
typedef fsu::Vector     <char>                  Color;

class MovieMatch
{
  public:
       MovieMatch    ();
  bool Load          (const char* filename);
  bool Init          (const char* actor);
  void Shuffle       ();
  long MovieDistance (const char* actor); 
  void ShowPath      (std::ostream& os) const;
  void ShowStar      (Name name, std::ostream& os) const;
  void Hint          (Name name, std::ostream& os) const;
  void Dump          (std::ostream& os) const;

  private:
  static void Line (std::istream& is, Vector& movie);
  // symbol graph private members
    Graph  g_;
    Vector name_;
    Color color_;    // Vector to hold color coding     
    AA     vrtx_;

  // BFS survey private members
    Name baseActor_;
    BFS bfs_;

  // Store path private member
     fsu::List<Vertex> path_;
}; // end of MovieMatch     

//template < typename T >

// Predicate class
class LessThan
{
public:
  bool operator () (const fsu::String&s1, const fsu::String& s2) const
  {
    size_t i = 0;                               // counter for while loop 
    size_t counter = 0;                         // counter to hold string size

    if (s1.Size() > s2.Size())                  // check if string 1 > string 2 
      counter = s2.Size();                      // assign length of string2 to counter
    else
      counter = s1.Size();                      // assign length of string1 to counter 
   
    while (i < counter)
    {
      if (tolower(s1[i]) > tolower(s2[i]))      // ignore case and compare strings
        return false;
      else if (tolower(s1[i]) < tolower(s2[i]))
        return true;
      else
        ++i;                         
    }
    return true;
  }
};

MovieMatch::MovieMatch() : g_(), name_(), color_('\0'), vrtx_(), baseActor_(0), bfs_(g_), path_() 
{}

void MovieMatch::Line (std::istream& is, fsu::Vector<Name>& movie)
{
  movie.Clear();                                  
  fsu::String line;
  char delim = '/';
  line.GetLine(is);
  char* name_buffer = new char [1+line.Size()];
  size_t pos = 0, next = 0;
  while (pos < line.Size())
  {
    next = line.Position(delim,pos);
    for (size_t i = pos; i < next; ++i)
    {
      name_buffer[i - pos]= line[i];
    }
    name_buffer[next-pos] = '\0';
    movie.PushBack(name_buffer);
    pos = ++next; // skip delimiter
  }
delete [] name_buffer;
}

bool MovieMatch::Load (const char* filename)
{
  Name s = "";                           // String object
  Vertex v = 0;                          // vertex for data
  std::ifstream ifs;                     // infile stream object
  ifs.open(filename);                    // open file 
  if (!ifs)                              // exit if file does not open 
    return false;
  
  std::cout << "\n Loading database " << filename << " (first read) "; 
  while (!ifs.eof())
  {
    s.GetNext(ifs, '/');                // Get next word using until delimiter reached 
    if (!s.Size())                      // break out of look if size returns 0 
      break;

    if (name_.Size() > 2 * vrtx_.NumBuckets()) // Check size vs number of buckets  
      vrtx_.Rehash();                          // Rehash for decreased load time
       
    if (!vrtx_.Retrieve(s, v))                 // check if key not in database
    {
      name_.PushBack(s);                       // Push string into vector 
      vrtx_.Insert(s,name_.Size() - 1);        // insert key and data into AA
    }
  }    
  ifs.close();                                 // close file 
  
  Vector s2;                                   // second vector for color coding
  size_t m = 0;                                // counter to hold movie colors
  size_t a = 0;                                // counter to hold actor colors
  color_.SetSize(name_.Size());                // set size of vector to size of name vector
  g_.SetVrtxSize(vrtx_.Size());                // set graph size to vector size 
  ifs.open(filename);                          // open file name 
  
  if (!ifs)
    return false;

  std::cout << "...(second read)";              // inform user of second read 

  while (!ifs.eof())
  {
    Line(ifs,s2);                               // send file object and second vector into line function
  
  if (!s2.Size())
    break;

  for (size_t i = 1; i < s2.Size(); ++i)       // loop through vector 2 size 
  {  
    g_.AddEdge(vrtx_[s2[0]], vrtx_[s2[i]]);    // add an edge from first element of s2 to graph from i  
    color_[vrtx_[s2[0]]] = 'm';                // color code vector for movie  
    color_[vrtx_[s2[i]]] = 'a';                // color code vector for actor
  }
  }

  for (size_t j = 0; j < color_.Size(); ++j)    // loop through color vector
  {
    if (color_[j] == 'm')                        
      ++m;                                       // increment movie count if location is m
    else if (color_[j] == 'a')                   
      ++a;                                       // increment actor count if location is a
  }
  std::cout << " ... done." << std::endl;
  std::cout << " " << m << " movies and " << a << " actors read from " << filename << std::endl;  
  ifs.close();                                   // close file
  vrtx_.Rehash();                                // call rehash for decreased load times
   return true;
}

bool MovieMatch::Init (const char* actor)
{
  baseActor_ = actor;                            // assign passed in actor to base actor
  bfs_.Reset(vrtx_[baseActor_]);                                  
  Vertex v = 0; 
    
  if (!vrtx_.Retrieve(baseActor_, v) || color_[v] == 'm')    // check if baseactor is not in aa or vertex is not a movie    
    return false;
  else
  { 
    bfs_.Search(vrtx_[baseActor_]);
  }                                                          // bfsearch for baseactor 
  return true; 
}  

void MovieMatch::Shuffle ()
{
  g_.Shuffle();                                             // call graph shuffle
  bfs_.Reset();                                             // reset bfs 
  bfs_.Search(vrtx_[baseActor_]);                           // bfsearch for baseactor
}  

long MovieMatch::MovieDistance (const char* actor)
{
  Vertex v = 0;
  path_.Clear();
  if (vrtx_.Retrieve(actor, v))                              // check if actor in db
    {
      if (bfs_.Color()[v] == 'b')                            // determine if actor reachable 
       {
         if (bfs_.Distance()[v] % 2 != 0)                     
           return -1;                                        // return -1 if not actor
         else
         {
           if (actor == baseActor_)
             return 0;

           path_.Insert(vrtx_[actor]);                       // insert actor vertex to path/list
           while (bfs_.Parent()[v] != vrtx_[baseActor_])      
           {
             path_.Insert(bfs_.Parent()[v]);                 // Inset parentinto path/list
             v = bfs_.Parent()[v];                           // assign parent vertex as current vertex
           }
           path_.Insert(vrtx_[baseActor_]);                  // insert baseactor vertex to path/list 
           return (path_.Size()/2);                          // return half path size
         }
       }
       else
       {
         return -2;                                          // return -2 name not reachable from baseacotor
       }
    } 
  return -3;                                                 // return -3 when name not in DB 
}

void MovieMatch::ShowPath (std::ostream& os) const
{
  fsu::List<Vertex>::ConstIterator i;                // const iterator for path traversal
  os << '\n';
  size_t p = 0;                                      // p counter for on off printing
  for (i = path_.Begin(); i != path_.End(); ++i)
  {
    if (p % 2 == 0)                                  // display name
      os << " " << name_[*i] << '\n';
    else
      os << "   | " << name_[*i] << '\n';            // displace name with seperator
    ++p;                                             // increment p
  }  
  os << '\n';
}

void MovieMatch::ShowStar (Name name, std::ostream& os) const
{
  Vertex v = 0;
  typename Graph::AdjIterator i;
  Vector starvect_;                                    // vector for star data  
  size_t x = 0; 
  if (vrtx_.Retrieve(name, v))                         // retrieve vertex associated with name
    {
      os << '\n';
      os << " " << name << '\n';                      
      for(i = g_.Begin(v); i != g_.End(v); ++i)      
      {
        starvect_.PushBack(name_[*i]);                // Push data into star vector
      }
  typename Vector::Iterator j = starvect_.Begin();
  typename Vector::Iterator k = starvect_.End();

  g_heap_sort(j,k);                                   // call to heap sort 

      for ( i = g_.Begin(v); i != g_.End(v); ++i)
      {     
        os << "   | " << starvect_[x++] << '\n';      // display star vector 
      }
      os << '\n';
    }
}

void MovieMatch::Hint (Name name, std::ostream& os) const
{
  
  Vector hintvect;                                     // create hint vector
  LessThan p;                                          // predicate class object
  hintvect.SetSize(name_.Size());                      // set size of hint vect to name vector size
  hintvect = name_;                                    // copy name vector to hint vector
  typename Vector::Iterator i = hintvect.Begin();      // set i iterator to beginning of hint vector
  typename Vector::Iterator j = hintvect.End();        // set j iterator to beginning of hint vector
  g_heap_sort(i, j, p);                                // call heap sort
  typename Vector::Iterator ubi;                       // upper bounds iterator
  typename Vector::Iterator lbi;                       // lower bounds iterator
  
  if (name.Length() > 6)                               // check if length of name is larger than 6
  {
    name.SetSize(6);                                   // truncate name
    ubi = g_upper_bound(i, j, name, p);                // call upper bounds 
    name = name + "zz";                                // append string with "zz"
    lbi = g_lower_bound(i, j, name, p);                // call lower bounds
  }
  else
  {
    ubi = g_upper_bound(i, j, name, p);                // if name smaller than 6 characters
    name = name + "zz";                                // call upper bounds and append "zz" 
    lbi = g_lower_bound(i, j, name, p);                // call lower bounds  
  }

  typename Vector::Iterator uppGuard = hintvect.Begin(); // upper bounds guard iterator 
  typename Vector::Iterator lowGuard = hintvect.End();   // lower bounds guard iterator

  for (size_t i = 0; i < 3; ++i)                         // suggestion count
  {
    if (ubi != uppGuard)                                 // check is upperbound is not beginning of vect
    {
      --ubi;                                             // decrement upper bounds 
    }
    else
    {/* do nothing */}
    
    if (lbi != lowGuard)                                  // check is lowerbound not ending of vect
    {
      ++lbi;                                              // increment lowerbounds
    }
    else
    {/*do nothing */}
  }
    
  for (; ubi != lbi; ++ubi)                               // traverse through upper to lower bounds  
    os << *ubi << std::endl;                              // display suggestions
}
void MovieMatch::Dump (std::ostream& os) const
{
   ShowAL(g_,os);
   WriteData(bfs_,os);
   vrtx_.Dump(os);
   for (size_t i = 0; i < name_.Size(); ++i)
   {
    os << "name_[" << i << "] = " << name_[i] << '\t';
    os << "vrtx_[" << name_[i] << "] = " << vrtx_[name_[i]] << '\n';
  }
   vrtx_.Analysis(std::cout);
}
#endif
