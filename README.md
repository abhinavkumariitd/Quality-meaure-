It contains quality measures for evaluating the strengths of overlapping or disjoint communities in directed weighted/unweighted networks To use the code, first compile it on a linux/unix terminal as follows:

g++ -std=c++11 directed_weighted_overlapping_modularity.cpp -o Qdwo

Then an executable code Qdwo would be produced. This can be used as follows:

./Qdwo network-file community-file -w 

If the network is unweighted, then the code can be called as follows:

./Qdwo network-file community-file 

Here, network-file is the file containing the directed network in edgelist format(in two/three columns), where the first two columns indicate the list of edges, and the third one indicates the corresponding weights. And, community-file is the file containing the community structure of the network, with each line containg a community. Note that both the files must be in the same directory where the executable code Qdwo lies. The flag -w is needed only if the network is weighted, that is, if the network-file contains three columns.


 
