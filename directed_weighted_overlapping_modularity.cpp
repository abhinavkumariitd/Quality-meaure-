#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <string>
#include <iterator>
#include <stdio.h>
#include <list>
using namespace std;

class Graph
    {
    public:
        set<int> V; //only vertices with positive degree are stored
         map<int,set<int>>in_neighbors, out_neighbors;
        map<int, float> indegree, outdegree, strength;
        map<int, map<int, float> > weight;
        void read_edgelist(string&, bool, bool);
        inline int order(){ return V.size(); }
        int ecount();
        float get_weight(int, int);
        void print_graph();
        friend float overlapping_weighted_modularity(Graph&, map<int, set<int> >&, int, map<int, set<int> >& , int &);
    };

void Graph::read_edgelist(string& edgefile, bool weighted, bool directed)
{
    ifstream fin;
    fin.open(edgefile);
    if(!fin.is_open())
    {
        cout<<"The file containing the edgelist could not be opened."<<endl;
        exit(1);
    }
    string line;
    while ( getline(fin, line ) )
    {
        istringstream is( line );
        int u, v;
        is>>u;
        is>>v;
        if(u == v)
            continue;
        float w;
        if(weighted == true)
        {
            if(is.str() == " " || is.eof())
            {   cout<<endl<<"The edge list has missing weights."<<endl;
                exit(1);
            }
            is>>w;
        }
        else
            w = 1;
        weight[u][v] = w;
       if (directed==true)
         {
           if(outdegree.find(u) == outdegree.end())
                outdegree[u] = w;
            else
                 outdegree[u] = outdegree[u]+w;
            if(indegree.find(v) == indegree.end())
                indegree[v] = w;
            else
                 indegree[v] = indegree[v]+w;
          }
        V.insert(u);
        V.insert(v);
        in_neighbors[v].insert(u);
       out_neighbors[u].insert(v);
    }
}
int Graph::ecount()
{
    map<int, map<int, float> >::iterator mi;
    int degree_sum = 0;
    for(mi = weight.begin(); mi != weight.end(); ++mi)
        degree_sum += mi->second.size();
    return degree_sum/2;
}
float Graph::get_weight(int u, int v)
{
    if(weight[u].find(v) != weight[u].end())
        return weight[u][v];
    else
        return 0;
}

void read_communities(string& comfile, map<int, set<int> >&, map<int, set<int> >&);
void print_communities(map<int, set<int> >&);
float overlapping_weighted_modularity(Graph&, map<int, set<int> >&, map<int, set<int> >&, int & );
void print_set(set<int>&);
void print_vector(vector<int>& );
void usage();

int main(int argc, char* argv[])
{
    string network_file, community_file;
    bool weighted = false, directed=true, suppress_info = false;
    if(argc < 3 || argc > 5)
        usage();
    else
    {
        network_file = string(argv[1]);
        community_file = string(argv[2]);
    }
    int i=3;
    while(i <= argc-1)
    {
        string arg = string(argv[i]);
        if(arg == "-w")
        {
            weighted = true;
            i++;
        }
        else
            if(arg == "-si")
            {
                suppress_info = true;
                i++;
            }
            else
                usage();
    }
    Graph g;
    map<int, set<int> > coms;
    g.read_edgelist(network_file, weighted, directed);
   // g.print_graph();
    map<int, set<int> > memberships;
    read_communities(community_file, coms, memberships);
    //print_communities(coms);
    int ghost_communities = 0;
    //print_communities(C2);
    float Q_wo = overlapping_weighted_modularity(g, coms, memberships, ghost_communities);
    if(!suppress_info)
    {
        cout<<"-------------------------------------------------------------"<<endl;
        cout.setf(ios::left, ios::adjustfield);
        cout<<setw(25)<<"actual communities"<<"= "<<setw(20)<<coms.size()<<endl;
        cout<<setw(25)<<"ghost communities"<<"= "<<setw(20)<<ghost_communities<<endl;
        cout<<setw(25)<<"Q_wo"<<"= "<<setw(20)<<Q_wo<<endl;
        cout<<"--------------------------------------------------------------"<<endl;
    }
    else
        cout<<" "<<Q_wo<<" ";
    return 0;
}

void read_communities(string& comfile, map<int, set<int> >& coms, map<int, set<int> >& memberships)
{
    std::ifstream fin(comfile);
    string line;
    if(!fin.is_open())
    {   cout<<"Community file could not be opened."<<endl;
        exit(1);
    }
    int i = 1;
    while ( std::getline(fin, line ) )
    {
        istringstream is( line );
        int v;
        while(is>>v)
        {
            coms[i].insert(v);
            memberships[v].insert(i);
        }
        i++;
   }
}

void print_communities(map<int, set<int> >& C)
{
    cout<<"community"<<endl;
    map<int, set<int> >::iterator sitr;
    for(sitr = C.begin(); sitr != C.end(); ++sitr)
    {
        copy((*sitr).second.begin(), (*sitr).second.end(), ostream_iterator<int>(cout, " "));
        cout<<endl;
    }
}

float overlapping_weighted_modularity(Graph& g, map<int, set<int> >& coms, map<int, set<int> >& memberships, int &ghost_communities)
{
    map<int, set<int> >::iterator ci;
    set<int>::iterator si;
    float indeg, deg, ind_frac,rho,indeg_c, outdeg_c, indeg_inside_c,outdeg_inside_c;
    float Q_wo = 0;
    for(ci = coms.begin(); ci != coms.end(); )
    {
        bool flag = false;
        for(auto p = ci->second.begin(); p != ci->second.end(); ++p)
            if(g.V.find(*p) == g.V.end())
            {
                ci = coms.erase(ci);
                flag = true;
                ghost_communities++;
                break;
            }
        if(flag == false)
            ++ci;
    }
   // print_communities(coms);
    int p=0;
   for(ci = coms.begin(); ci != coms.end(); ++ci)
    {
       if(ci->second.size()>1)
         {
           int total_memberships = 0;
           indeg_inside_c=0,outdeg_inside_c=0,indeg_c=0,outdeg_c=0;

          for(si = ci->second.begin(); si != ci->second.end(); ++si)
           {
               // cout<<"node = "<<*si<<endl;
               total_memberships += memberships[*si].size();

              for(auto ki=g.in_neighbors[*si].begin(); ki != g.in_neighbors[*si].end(); ++ki)
               {
                 if(ci->second.find(*ki) != ci->second.end())
                   indeg_inside_c += g.weight[*ki][*si];
               }

              //  cout<<"indeg_inside_c=="<<indeg_inside_c<<endl;

             for(auto ki=g.out_neighbors[*si].begin();ki!=g.out_neighbors[*si].end();++ki)
               {
                 if(ci->second.find(*ki) != ci->second.end())
                    outdeg_inside_c += g.weight[*si][*ki];
                }
                  //cout<<"outdeg_inside_c=="<<outdeg_inside_c<<endl;
                  indeg_c+=g.indegree[*si];
                  outdeg_c+=g.outdegree[*si];
                //cout<<"indeg_c= "<<indeg_c<<"  outdeg_c = "<<outdeg_c<<endl;

            }

         //cout<<"total memberships = "<<total_memberships<<endl;
         //cout<<"com size = "<<coms.size()<<endl;
         //cout<<"node in  community = "<<ci->second.size()<<endl;
         if(indeg_c!=0 &&outdeg_c!=0)
          Q_wo += (((ci->second.size()*coms.size() - total_memberships)*((outdeg_c*indeg_inside_c)+(indeg_c*outdeg_inside_c)))/(ci->second.size()*coms.size()*2*indeg_c*outdeg_c));
          else
            Q_wo+=0;
         // cout<<"Q_wo ==="<<Q_wo<<endl;
          //cout<< "\n\n\n";
         }


    }
    return Q_wo/coms.size();
}
void usage()
{
    cout<<"Please, follow the syntax as given below:"<<endl;
    cout<<"-----------------------------------------------------------"<<endl;
    cout<<"./Q_dwo network_file community_file -w -si"<<endl;
    cout<<"-----------------------------------------------------------"<<endl;
    cout<<"Here, network_file is the file holding the directed network in edgelist form,"<<endl;
    cout<<"and community_file is the file holding the communities, each written on a separate line."<<endl;
    cout<<"Use flag w if your network is weighted."<<endl;
    //cout<<"Use flag d if your network is directed."<<endl;
    cout<<"Use flag si if you do not want other than required information."<<endl;
    exit(1);
}


    void Graph::print_graph()
{
    cout<<"No. of vertices: "<<order()<<endl;
    cout<<"No. of edges: "<<ecount()<<endl;
    cout<<"\n"<<"vertex"<<setw(8)<<"degree"<<endl;
    cout<<"\n"<<"node"<<setw(10)<<"nbrs"<<setw(10)<<"weight"<<endl;
    map<int, map<int, float> >::iterator it;
    for(it=weight.begin(); it != weight.end(); ++it)
    {
        cout<<it->first<<setw(5)<<"--- "<<endl;
        map<int, float>::iterator jt;
        for(jt = it->second.begin(); jt != it->second.end(); ++jt)
            cout<<setw(10)<<jt->first<<setw(10)<<jt->second<<endl;
    }


cout<<"indegree sum"<<endl;
      for(auto it=indegree.begin();it!=indegree.end();++it)
      {
          cout<<it->first<<"--"<<it->second<<endl;
      }
cout<<"outdegree sum"<<endl;
      for(auto it=outdegree.begin();it!=outdegree.end();++it)
      {
          cout<<it->first<<"--"<<it->second<<endl;
      }

      cout<<"strength"<<endl;
      for(auto it=strength.begin();it!=strength.end();++it)
      {
          cout<<it->first<<"--"<<it->second<<endl;
      }


}


void print_vector(vector<int>& v)
{
    cout<<endl;
    copy(v.begin(),v.end(),ostream_iterator<int>(cout, " "));
    cout<<endl;
}

void print_set(set<int> &s)
{
    cout<<endl;
    copy(s.begin(),s.end(),ostream_iterator<int>(cout, " "));
    cout<<endl;
}
