// CPCS 324 Algorithms & Data Structures 2
// Graph data structure starter - Transitive Closure Package



// -----------------------------------------------------------------------
// simple graph object with linked-list edge implementation and minimal fields
// extra vertex and edge member fields and methods to be added later as needed
//


var _v = [], _e = [];   


// -----------------------------------------------------------------------

function main_graph()   
{
    
	// create an undirected graph 

    var g = new Graph();
    

    // set input graph properties (label, directed etc.)

    g.label = 'Exercise 9.2: 1b (Levitin, 3rd edition)';


    // use global input arrays _v and _e to initialize its internal data structures

    g.readGraph(_v, _e);

    g.makeAdjMatrix(); 


    // use printGraph() method to check graph

    g.printGraph();


    // perform depth-first search and output stored result

    g.connectedComp = g.topoSearch(g.dfs);

    document.write("dfs_push: ", g.dfs_push, "<br><br>");


    // report connectivity status if available

    document.write(g.componentInfo());


    // perform breadth-first search and output stored result

    g.connectedComp = g.topoSearch(g.bfs);

    document.write("bfs_order: ", g.bfs_order, "<br>");


    
	// output transitive closure matrix using Warshall-Floyd
    
    document.write('<br>Transitive closure<br>');
    
    g.warshallFloyd();
    
    if(g.adjMatrix[0][0].tc === undefined)
    {
        document.write("Transitive closure is computed for directed graph only<br>");
    }
    else
    {
        for(var i=0; i<g.nv; i++){
            for(var j=0; j<g.nv; j++){
                document.write(g.adjMatrix[i][j].tc,',');
            }
            document.write('<br>');
        }
    }


    // output distance matrix using Warshall-Floyd
	
	document.write('<br>Distance matrix <br>');
    
    g.warshallFloyd();
    	
    if(g.adjMatrix[0][0].dist === undefined)
    {
        document.write("Distance matrix is computed for weighted graph only<br>");
    }
    else
    {
        for(var i=0; i<g.nv; i++){
            for(var j=0; j<g.nv; j++){
                document.write(g.adjMatrix[i][j].dist,',');
            }
            document.write('<br>');
        }
    }

    
    //print minimal spanning tree

    document.write('<br>MST by Prim2 (linear PQ)<br>');
    
    g.prim();
    if(g.weighted){
        print_tree(g.VT);
    }
    else
    {
        document.write('minimal spanning tree works for weighted graph only<br>');
    }

    //Print shortest path by Dijkstra from vertex 0. Shortest path problems works for weighted graph only

    document.write('<br>Shortest paths by Dijkstra from vertex 0<br>');
    g.Dijkstra(0);
    if(g.weighted){
        print_tree(g.VT);
    }
    else
    {
        document.write('Shortest path problems works for weighted graph only<br>');
    }

    
    //Print distance matrix from Dijkstra. Distance matrix is computed for weighted graph only
    document.write('<br>Distance matrix from Dijkstra<br>');
    
    var disMDijkstra = [];
    
    if(g.weighted){
        for(i=0; i<g.nv; i++){
            //The shortes path from vertex i to all other vertices
            g.Dijkstra(i);             
            disMDijkstra[i] = [];

            //fill the dist of row i
            for(j=0; j<g.nv; j++){
                disMDijkstra[i][g.VT[j].tree] = g.VT[j].distance;
            }

            //print row i
            document.write(disMDijkstra[i],'<br>');
        }
    }
    else
    {
        document.write('Distance matrix is computed for weighted graph only<br>');
    }
}


// -----------------------------------------------------------------------

function Vertex(v)
{
	// base base property fields
	this.label = v.label;  

		
	// --------------------
	// more student fields next
	this.visit = false; 
	this.adjacent = new List(); 
		
	// --------------------
	// more student methods next
	this.adjacentById = adjacentByIdImpl;
    this.incidentEdges = incidentEdgesImpl; 
    this.insertAdjacent = insertAdjacentImpl;
    this.vertexInfo = vertexInfoImpl; 
}

// -----------------------------------------------------------------------

function Edge(vert_i,weight)
{
	// base property fields	
	this.target_v = vert_i;  

	
	// --------------------
	// more student fields next
	if (!(weight === undefined)) 
    {
        this.weight = weight;
    }
    else
    {
        this.weight = null;
    }
	
}


// -----------------------------------------------------------------------

function Graph()
{	
	// base property fields

	this.vert = [];
	this.nv = 0;  
    this.ne = 0; 
    this.digraph = false; 
    this.dfs_push = []; 
    this.bfs_order = []; 
    this.label = "";

	
	// base member methods
	
	this.listVerts = listVertsImpl; 
    this.readGraph = better_input;
    this.add_edge = addEdgeImpl3;
    this.printGraph = printGraphImpl; 
    this.dfs = dfsImpl;
    this.bfs = bfsImpl; 
	
	// more student fields next

	this.weighted = false; 
    this.connectedComp = 0; 
    this.adjMatrix = [];
    this.VT = []; 
  	
	
	// --------------------
	// more student methods next 
    this.makeGraph = makeGraphImpl;
    this.makeAdjMatrix = makeAdjMatrixImpl3; 
    this.isConnected = isConnectedImpl; 
    this.componentInfo = componentInfoImpl;
    this.topoSearch = topoSearchImpl;
    
    /**
       Return the minimum spanning tree (MST)
      @method
    */
    this.prim = primImpl2;
    
    /**
       Return the single source shortest path (SSSP) of a specific vertex
      @method
    */
    this.Dijkstra = Dijkstra_; 
	

	// transitive closure package (requirements in line comments, to be removed and replaced by JSDOCs) 
	/**
       Return true if two vertices have a path in digraph otherwise false
      @method
    */

    this.hasPathImpl = hasPathImpl;                 
    
    /**
      Return distance of the shortest path between two vertices in weighted graph
      @method
    */

    this.shortestPath = shortestPathImpl;                
    
    /**
      Return true if directed graph is acyclic otherwise false
      @method
    */

    this.isDAG = isDAGImpl;                    
    
      /**
      fill the tc (transitive closure) field of the adjacent matrix if the graph is directed.
      fill the dist (distance matrix) field of the adjacent matrix if the graph is weighted.
      @method
    */

    this.warshallFloyd = warshallFloydImpl;            

    /**
      Return TC matrix for digraph based of DFS
      @method
    */

	this.dfsTC = dfsTCImpl;                     
		

}


// -----------------------------------------------------------------------
// functions used by methods of Graph and ancillary objects

// -----------------------------------------------------------------------
// begin student code section
// -----------------------------------------------------------------------

// transitive closure package 

/**
 * function thah helps to print the output of prim and Dijkstra algorithms
   @author Ansam Ali
   @param {object[]} VT an array of 
*/

function print_tree(VT){
    for(i=0; i<VT.length-1; i++){
        document.write((VT[i].distance===undefined? '':VT[i].distance) ,'(',VT[i].parent,',',VT[i].tree,'),');
    }
    document.write((VT[VT.length-1].distance===undefined? '':VT[VT.length-1].distance) ,'(',VT[VT.length-1].parent,',',VT[VT.length-1].tree,').<br>');
}
/**
 * @descriptionThis is Dijkstra Algorithm, This algorithm is applicable to undirected and directed graphs
 *                  with nonnegative weights only.
 * @async
 * @method 
 * @param {integer} s start vertex... email - abahkali0005@stu.kau.edu.sa
 * @author Ansam Ali 
 * 
 */

function Dijkstra_(s){

    //chech if the graph is weighted, otherwise the shortest path will not be computed. 
    if(this.weighted){
        var pq = new PQueue();
        var penultimateV = [];
        var dist = [];

        //initialize queue
        for(var i=0; i<this.nv; i++){
            dist[i] = Infinity;        //dist[i] is the shortest path length to vertex i. initially infinity(unreachable) 
            penultimateV[i] = "-";           //penultimateV[i] is the next-to-last vertex of vertex i.
            pq.insert(i, dist[i]);

            //mark all vertices as unvisited. that is, not in the tree
            this.vert[i].visit = false;
        }

        //intalize the tree by the source
        dist[s] = 0;
        pq.insert(s, dist[s]);
        

        //get the next vertex to be added to the tree
        for(var i=0; i<this.nv; i++){
            u = pq.deleteMin();
            this.vert[u].visit = true;

            this.VT[i] = {
                tree: u,                    //tree represents the vertex.
                parent: penultimateV[u],     //parent represents the next-to-last vertex in path.
                distance: dist[u]      //distance represents the total distance from the source to the current vertex.
            } 
        

            //to update the tree
            //get adjacents of the newly added vertex
            var adjacents = this.vert[u].incidentEdges();
            for(var j=0; j<adjacents.length; j++){
                var id = adjacents[j].adjVert_i;
                var weight = adjacents[j].edgeWeight;

                //if there's a shorter path to a vertex, consider that path (update the queue, dist, and penultimateV arrays) 
                if(dist[u] + weight < dist[id] && !this.vert[id].visit){
                    dist[id] = dist[u] + weight;
                    penultimateV[id] = u;
                    pq.insert(id,dist[id]);
                }

            }
        }
    }
}



function primImpl2(){

  
 //chech if the graph is weighted, otherwise there is no minimal spanning tree. 
 
 if(this.weighted){

    var pq = new PQueue();  // create qeueue 

    var penultimateV = []; // array contains penulimate

    var weight = []; //  array contains weight of edges

    this.VT=[]; //create an array of json objects (it's minimum spanning tree)

    var u;



    // mark all vertices unvisited
    // intitalize penulimate array by '-' that mean's no penulimate until now.
    // intitalize weight array by infinity.

    for (var j = 0; j < this.nv; j++)

    {

        this.vert[j].visit = false;

        penultimateV[j] = "-";

        weight[j] = Infinity;

    }



    // edges array contains all the edges (v,u)
    var edges =this.vert[0].incidentEdges();

    for(var i=0;i<edges.length;i++){

        //insert dest. vertex and weight of edge in the queue.
        pq.insert(edges[i].adjVert_i, edges[i].edgeWeight); 

        // make a src. vertex penultimateV of dest. vertex.
        penultimateV[edges[i].adjVert_i]=0; 

        //add weight of edge to weight array.
        weight[edges[i].adjVert_i] = edges[i].edgeWeight;

    }



    //intitalize Vt by the first vertex and mark it as true.
    this.vert[0].visit = true;

    this.VT[0] = {

        //vertex
        tree: 0,

        //penultimateV of this vertex
        parent: penultimateV[0]

    } 



    //find a way to the rest n-1 vertices

    for(var i=1; i<this.nv; i++){

        //delete minimum-weight edge from queue
        u =pq.deleteMin();

        //mark a vertex as visit
        this.vert[u].visit = true;

       
        //add vertex that delete from qeueue into vt
        this.VT[i] = {

            tree: u,

            parent: penultimateV[u]

        }  



       // edges array contains all the edges (v,u)
        edges=this.vert[u].incidentEdges();

        for(var j=0;j<edges.length;j++){


         //if vertex doesn't mark as visit which mean's don't add into VT & weight of edge less than weight of edge in qeueue 
         //then insert or update  dest. vertex and its weight ,
        //update penultimateV and update weight of edge in weight array. 
            if(!this.vert[edges[j].adjVert_i].visit && edges[j].edgeWeight < weight[edges[j].adjVert_i]){

                pq.insert(edges[j].adjVert_i,edges[j].edgeWeight);

                penultimateV[edges[j].adjVert_i]=u;

                weight[edges[j].adjVert_i] = edges[j].edgeWeight;

            }
            

        }

    }

 }

}

/**
 * This is an imlementation for Prim's algorithm which finds the minimal spanning tree(MSP)
 * It finds the connected acyclic subgraph (tree) that contains all vertices with smallest total weight.
   @author Ansam Ali
   @deprecated Obsolete first implementation. Use {@link primImpl2} instead.
   @implements Graph#prim
   @returns {object} array of edges of the minimal spanning tree 
*/

function primImpl(){

    // mark all vertices unvisited
    for (var j = 0; j < this.nv; j++)
    {
        this.vert[j].visit = false;
    }
    
    //intitalize Vt by the first vertex and mark it as true, Et by empty array.
    Vt = [this.vert[0]];
    this.vert[0].visit = true;
    Et = [];
    var src;                    //id of the source of the edge

    //find a way to the rest n-1 vertices
    for(i=1; i<this.nv; i++){

        var min_weight = Infinity;      //initalize the minimum weight with Infinity
        var min_weight_edge = null;     //no edge found till now

        
        //check fringe vertices Vt vertices
        for(j=0; j<Vt.length; j++){
            edges = Vt[j].incidentEdges();
            for(k=0; k<edges.length; k++){
                //update the minimum edge information if found
                if(edges[k].edgeWeight < min_weight && !this.vert[edges[k].adjVert_i].visit){   
                    min_weight = edges[k].edgeWeight;   
                    min_weight_edge = edges[k];
                    src = Vt[j];                
                }
            }
        }
        
        //add the target vertex to Vt and mark it as true
        Vt[i] = this.vert[min_weight_edge.adjVert_i];
        this.vert[min_weight_edge.adjVert_i].visit = true;

        Et[i-1] = {
            source: src,
            dest: this.vert[min_weight_edge.adjVert_i],
            weights: min_weight_edge.edgeWeight
        };
    }
    return Et;
}
/**
 * This func is to check if there is a path between 2 given vertices in a directed graph.
   @author Ansam Ali
   @implements Graph#hasPathImpl
   @param {number} u_i Source vertex id
   @param {number} v_i Target vertex id
   @returns {boolean}  Is there a path from vertex with id u_i to vertex with id v_i in a directed graph. error message if undirected.
*/

function hasPathImpl(u_i, v_i)
{
    //if the grapgh is undirected, tc will be undefiend ==> return error message.
    //if the grapgh is directed, and tc !=1 ==> no path ==> return false.
    //otherwise return true

    if(this.adjMatrix[u_i][v_i].tc === undefined){
        return "<br>This method doesn't work for undirected graph<br>";
    }
	return (this.adjMatrix[u_i][v_i].tc == 1 ? true : false);
}



function shortestPathImpl(u_i,v_i)
{
    //if the dist is undefiend ==> unweightd graph ==> can't compute the shortest path.
    var d = this.adjMatrix[u_i][v_i].dist;
	return (d === undefined? "<br>This method doesn't work for unweighted graph<br>":d);
}


function isDAGImpl()
{
    if(this.adjMatrix[0][0].tc === undefined){
        return "<br>This method doesn't work for undirected graph<br>";
    }

	for (var i = 0; i < this.nv; i++) 
		{
		if (this.hasPathImpl(i, i))     //If there's a path from a vertix to itself, then there's a cycle. Otherwise not.
		{
			return false;
		}
	}
	return true;
}



function warshallFloydImpl() 
{
this.makeAdjMatrix();

floyd = [];
warshall = [];

for (var i = 0; i < this.nv; i++) 
{
    //make a copy of the adjacent matrix array
	floyd[i] = this.adjMatrix[i].slice();
	warshall[i] = this.adjMatrix[i].slice();


	for (var j = 0; j < this.nv; j++)
	{
        //if the value on the adjacent matrix at i,j is zero, there is no path from i to j
        //so the distance in infinity (unreachable).
        //except if i=j. There is a path of length 0 from any vertex to itself. 
		if (this.adjMatrix[i][j] == 0 && (i != j)) 
		{
			floyd[i][j] = Infinity;
		}
	}
}

    for (var k = 0; k < this.nv; k++) 
    {
        for (var i = 0; i < this.nv; i++) 
        {
            for (var j = 0; j < this.nv; j++) 
            {
                //If there's a path i->k and k->j  ===> there's a path i->j
                warshall[i][j] = (warshall[i][j] || (warshall[i][k] && warshall[k][j])) ? 1 : 0;
                
                //The shortest path i->j is the shortest between (i->k) + (k->j) and path (i->j)  
                floyd[i][j] = Math.min(floyd[i][j], floyd[i][k] + floyd[k][j]);

                this.adjMatrix[i][j] = {
                    tc: this.digraph? warshall[i][j]:undefined, //transitive closure is computed for undirected graph
                    dist: this.weighted? floyd[i][j]:undefined  //shortest distance matric is computed for weighted graph
                }
            }
        }

    }
       
} 


function dfsTCImpl()
{
    if(this.digraph){
        TC = [];

        for (var i = 0; i < this.nv; i++)
        {
            // mark all vertices unvisited
            for (var j = 0; j < this.nv; j++)
            {
                this.vert[j].visit = false;
            }

            //create + init the corresponding row

            TC[i] = [];

            for (var j = 0; j < this.nv; j++)
            {
                TC[i][j] = 0;
            }


            //dfs the current vertex v to visit all vertexs u such that there is a path (v,u) for all u
            this.dfs_push = []; //reset the dfs push order array
            this.dfs(i);  
            
            for (var j = 1; j < this.dfs_push.length ; j++) //start from 1 to exclude the vertex itself, we'll deal with cycle later
            {
                TC[i][this.dfs_push[j]]=1;	//set all vetices reachable by vertix i to 1
            }
        }

        for (var i = 0; i < this.nv; i++)

            {   //solution of the cycle problem. if a-->b and b-->a then a-->a, and we have all relations from a vertex to another vertex 
                for (var j = 0; j < this.nv; j++)   
                {
                    if ((TC[i][j] == 1) && (TC[j][i] == 1))
                    {
                        TC[i][i] = 1;
                    }
                }
            }

        return TC;
    }
    
    return null;
}
// -----------------------------------------------------------------------
// published docs section (ref. assignment page)
// use starter6-based P1M1 code as-is (fixes/improvements OK)
// no JSDOC comments in this section (docs already published)
// -----------------------------------------------------------------------

function addEdgeImpl(u_i, v_i)
{
    // fetch edge vertices using their id, where u: source vertex, v: target vertex
    var u = this.vert[u_i];
    var v = this.vert[v_i];


    // insert (u,v), i.e., insert v (by id) in adjacency list of u

    u.adjacent.insert(v_i);


    // insert (v,u) if undirected graph (repeat above but reverse vertex order)

    if (!this.digraph)
    {
        v.adjacent.insert(u_i);
    }
}


function addEdgeImpl2(u_i, v_i, weight)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.vert[u_i];
    var v = this.vert[v_i];


    // insert (u,v), i.e., insert v in adjacency list of u
    // (first create edge object using v_i as target, then pass object)

    var v_edge = new Edge(v_i);

    if (!(weight === undefined))
    {

        v_edge.weight = weight;
    }

    u.adjacent.insert(v_edge);


    // insert (v,u) if undirected graph (repeat above but reverse vertex order)

    if (!this.digraph)
    {

        var u_edge = new Edge(u_i);


        if (!(weight === undefined))
        {

            u_edge.weight = weight;

        }

        v.adjacent.insert(u_edge);
    }

}


function addEdgeImpl3(u_i, v_i, weight)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.vert[u_i];
    var v = this.vert[v_i];


    // insert (u,v), i.e., insert v in adjacency list of u
    // (first create edge object using v_i as target, then pass object)

    u.insertAdjacent(v_i, weight);


    // insert (v,u) if undirected graph (repeat above but reverse vertex order)

    if (!this.digraph)
    {

        v.insertAdjacent(u_i, weight);
    }
}


function printGraphImpl()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "WEIGHTED, " : "", this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ", this.ne, " EDGES:</p>");


    // report connectivity status if available

    document.write(this.componentInfo());


    // list vertices
    this.listVerts();

    document.write('<br>');
}


function listVertsImpl()
{
    var i, v; // local vars
    for (i = 0; i < this.nv; i++)
    {
        v = this.vert[i];
        document.write("VERTEX: ", i, v.vertexInfo());
    }
}


function better_input(v, e)
{
    // set vertex and edge count fields
    this.nv = v.length;
    this.ne = e.length;

    // input vertices into internal vertex array
    for (var i = 0; i < this.nv; i++)
    {
        this.vert[i] = new Vertex(v[i]);
    }

    //check if the graph is weighted or not
    this.weighted = e[0].w === undefined ? false : true;

    // input vertex pairs from edge list input array
    // remember to pass vertex ids to add_edge()
    for (var i = 0; i < this.ne; i++)
    {
        this.add_edge(e[i].u, e[i].v, e[i].w);
    }

    // double edge count if graph undirected
    if (!this.digraph)
    {
        this.ne = e.length * 2;
    }
}


function better_output()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "WEIGHTED, " : "", this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ", this.ne, " EDGES:</p>");


	// report connectivity status if available
	
	switch(this.connectedComp){
		case 0: document.write('<br> no connectivity info <br><br>');
		case 1: document.write('<br> CONNECTED <br><br>');
		default: document.write('<br> DISCONNECTED: '.concat(this.connectedComp,'<br>'));
	}

    // list vertices
    this.listVerts();
}


function topoSearchImpl(fun)
{
    connectedComp = 0;

    // mark all vertices unvisited
    for (var i = 0; i < this.nv; i++)
    {
        this.vert[i].visit = false;
    }

    // traverse unvisited connected components
    for (var i = 0; i < this.nv; i++)
    {
        if (!this.vert[i].visit)
        {
            fun == this.dfs ? this.dfs(i) : this.bfs(i);

            connectedComp++;
        }
    }

    return connectedComp;

}


function dfsImpl(v_i)
{
    // get landing vert by id then process
    var v = this.vert[v_i];
    v.visit = true;
    len = this.dfs_push.length;
    this.dfs_push[len] = v_i;

    // recursively traverse unvisited adjacent vertices
    var w = v.adjacentById();
    for (var j = 0; j < w.length; j++)
    {
        if (!this.vert[w[j]].visit)
        {
            this.dfs(w[j]);
        }
    }
}


function bfsImpl(v_i)
{
    // get vertex v by its id
    var v = this.vert[v_i];

    // process v 
    v.visit = true;

    // initialize queue with v
    var q = new Queue();
    q.enqueue(v_i);


    // while queue not empty
    while (!q.isEmpty())
    {

        // dequeue and process a vertex, u
        var u_i = q.dequeue();
        var u = this.vert[u_i];
        this.bfs_order[this.bfs_order.length] = u_i; //fill the bfs_order array when dequeu the vertex


        // queue all unvisited vertices adjacent to u
        var w = u.adjacentById();
        for (var i = 0; i < w.length; i++)
        {
            if (!this.vert[w[i]].visit)
            {
                this.vert[w[i]].visit = true;
                q.enqueue(w[i]);
            }
        }
    }
}


function makeAdjMatrixImpl()
{
    // initially create row elements and zero the adjacency matrix
    for (var i = 0; i < this.nv; i++)
    {

        this.adjMatrix[i] = [];

        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // for each vertex, set 1 for each adjacent  

        var w = this.vert[i].adjacentById();

        for (var j = 0; j < w.length; j++)
        {
            this.adjMatrix[i][w[j]] = 1;
        }
    }
}


function makeAdjMatrixImpl2()
{
    // initially create row elements and zero the adjacency matrix
    for (var i = 0; i < this.nv; i++)
    {

        this.adjMatrix[i] = [];

        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // for each vertex, set 1 for each adjacent if unweighted, its weight if graph is weighted
        var adj = this.vert[i].adjacent.traverse();
        for (var j = 0; j < adj.length; j++)
        {
            var edge = adj[j];
            this.adjMatrix[i][edge.target_v] = this.weighted ? edge.weight : 1;
        }
    }
}


function makeAdjMatrixImpl3()
{
    for (var i = 0; i < this.nv; i++)
    {

        this.adjMatrix[i] = [];

        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // for each vertex, set 1 for each adjacent if unweighted, its weight if graph is weighted
        var enodes = this.vert[i].incidentEdges();
        for (var j = 0; j < enodes.length; j++)
        {
            var edge_node = enodes[j];
            this.adjMatrix[i][edge_node.adjVert_i] = this.weighted ? edge_node.edgeWeight : 1;
        }
    }
}


function isConnectedImpl()
{
    connectedComp = this.topoSearch(this.dfs); //Do topological search (DFS, or BFS) for the graph to see the number of connected components if not set
    return connectedComp == 1; //The graph is connected if it has only one connected component
}


function componentInfoImpl()
{
	
	switch(this.connectedComp){
		case 0: return ('no connectivity info <br><br>');
		case 1: return ('CONNECTED <br><br>');
		default: return ('DISCONNECTED: '.concat(this.connectedComp,'<br>'));
	}

}


function adjacentByIdImpl()
{
    var traversal = this.adjacent.traverse();

    var traversal_array = [];

    for (var i = 0; i < traversal.length; i++)
    {
        traversal_array[i] = traversal[i].target_v;
    }

    return traversal_array;
}


function incidentEdgesImpl()
{
    var edges = this.adjacent.traverse(); //get array of incident edges objects

    var edges_info = [];

    for (var i = 0; i < edges.length; i++)
    { //create an array of json objects which contains edge information.
        edges_info[i] = {
            adjVert_i: edges[i].target_v,
            edgeLabel: "",
            edgeWeight: edges[i].weight
        };
    }

    return edges_info;
}


function vertexInfoImpl()
{
    var info = "".concat(" {", this.label, "} - VISIT: ", this.visit, " - ADJACENCY: ", this.adjacentById(), "<br>");

    return info;

}


function insertAdjacentImpl(v_i, weight)
{
    this.adjacent.insert(new Edge(v_i, weight));
}


//not implemented now
function makeGraphImpl()
{
}

// -----------------------------------------------------------------------
// Implementing functions of priority queue 


function isEmptyImpl()
{
    // check if PQ is empty or not
	return (this.pq.isEmpty());
	
}


function deleteMinImpl()
{
       // pointer ponits to the first node in the PQ 
       var pointer = this.pq.first;
       // the node of item with minimum value of key is the first node of PQ (intial value)
       var minimum = this.pq.first;
       // initialize the previous node of node of item with minimum value with null value
       var previous = null ;
       // search for the node of item with the minimum value of key 
       while(pointer != null){
           // check if the key of the current item is smaller than the key of the minimum item
         if((pointer.next != null) && (pointer.next.item.prior < minimum.item.prior)){
           // update the value of node of item with the minimum value of key
           minimum = pointer.next ;
           // update the value of node that has previous item of the minimum item 
           previous = pointer ;
         }
         // to move to the next node in PQ
         pointer = pointer.next;
       }
       // check if node of item with minimum value of key is the first one in PQ
       if(minimum == this.pq.first){ 
          // return the data item of item with the minimum value of key and delete it
          return this.pq.deleteFirst().item;
       }else{ // if node of item with minimum value of key is not the first one in PQ 
          // previous node points to the node next to the node with minimum item 
          previous.next = minimum.next;
          // return the data item of item with the minimum value of key
          return minimum.item.item;
       }
}


function insertImpl(item, key)
{
    // creat new PQNode to store the item and key 
    var item = new PQNode(item,key);
    // pointer ponits to the first item in the PQ
    var pointer = this.pq.first;
    // boolean variable : true in insert case and false in update case
    var toInsert = true;
    // check if the PQ is empty so it must be the insert case 
    if (this.isEmpty()) {
        // insert the item in the first of the PQ
        this.pq.insert(item);
        // change toInsert variable to false to not insert again the first position
        toInsert = false;

    } else { // if PQ is  not empty so it might be either the insert case or update case
        // search for the received item if it exists then it must be the update case 
        while (pointer != null) { 
            // if the data item of the received item equal to some data item for an item in PQ
            if (item.item == pointer.item.item) { 
                // update the priority value 
                pointer.item.prior = key;
                //change toInsert variable to false to not insert the updated item
                toInsert = false;
                // break the loop
                break; 
            }
            // to go to the next item in the PQ
            pointer = pointer.next;
        }

    }
    // insert case
    if (toInsert) {  
        // insert item in PQ
        this.pq.insert(item);
    }
    
}