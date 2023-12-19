# SandIMIN/SandIMIN-

   This project implements the SandIMIN and SandIMIN- algorithms for the following paper:
   <ul>
    <li>Jinghao Wang, Yanping Wu, Xiaoyang Wang, Ying Zhang, Lu Qin, Wenjie Zhang, Xuemin Lin, "Efficient Influence Minimization via Node Blocking", 2024.</li>
    </ul>

## Dataset
   

All the datasets in this paper can be downloaded in <strong>SNAP</strong>.

## Compile 

<code>g++ -O3 -o IMIN Sandwich.cpp sfmt/SFMT.c</code>

## Preliminary
<ul>
<li>There are three txt files in the 'dataset/test' folder: (i) attribute.txt; (ii) graph.txt; (iii) rumorSet_1.txt.</li>
<li>
The sample network graph.txt, in this case, contains only 4 nodes and 4 edges and is formated as follows:</li>
        
        1 2 0.3 (an edge from node 1 to node 2, with a probability of 0.3)
            
        2 3 0.5
            
        1 3 0.4
            
        1 4 0.2

<li>
The attribute.txt contains the number of edges and nodes in graph.txt, i.e.,</li>

    n=4 (number of nodes)

    m=4 (number of edges)

<li>
The rumorSet_1.txt contains the set of seed nodes, i.e.,</li>

    1
</ul>



## Run 

<ul>
    <li>Format the graph</li>
    <code>./el2bin graph.txt graph_ic</code>
    <li>Execute the Program</li>
    <code>./IMIN.out -dataset dataset/test -k 10 -rumorNum 1 -algo SandIMIN -epsilon 0.2 -gamma 0.1 -beta 0.1</code>
    <li>Arguments</li>
    <ul>
    <li>dataset: path to the dataset directory</li>
    <li>epsilon, gamma, beta: the parameter</li>
    <li>rumorNum: number of seed nodes of misinformation</li>
    <li>algo: algorithm (SandIMIN/SandIMIN-)</li>
    </ul>
</ul>






