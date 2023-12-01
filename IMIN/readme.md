   ## Example
        At first, there are three txt files in the 'dataset/test' folder: (i) attribute.txt; (ii) graph.txt; (iii) rumorSet_1.txt.
    
        The sample network graph.txt, in this case, contains only 4 nodes and 4 edges and is formated as follows:
            
            1 2 0.3 (an edge from node 1 to node 2, with a probability of 0.3)
            
            2 3 0.5
            
            1 3 0.4
            
            1 4 0.2
        
        The attribute.txt contains the number of edges and nodes in graph.txt, i.e.,
    
            n=4 (number of nodes)
            m=4 (number of edges)
        
        The rumorSet_1.txt contains the set of seed nodes, i.e.,
    
            1
    
        Next we convert graph.txt to binary file:
            Execute:
                ./el2bin graph.txt graph_ic
        
        Then we can start running our algorithm:
            Execute:
                ./IMIN.out -dataset dataset/test -k 10 -rumorNum 1 -algo SandIMIN -epsilon 0.2 -gamma 0.1 -beta 0.1
                ./IMIN.out -dataset dataset/test -k 10 -rumorNum 1 -algo SandIMIN- -epsilon 0.2 -gamma 0.1 -beta 0.1
            Arguments:
                -dataset:
                    path to the dataset directory
                -epsilon, -gamma, -beta:
                    the parameter 
                -k:
                    number of blockers
                -rumorNum:
                    number of seed nodes of misinformation
                -algo:
                    algorithm (SandIMIN/SandIMIN-)





