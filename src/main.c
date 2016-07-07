/*********************************************************************
Team: AK-Group - Deparment of Information Systems, University of Engineering and Technology
ACM SigMod Programming Contest 2016
@Department of Information Systems - University of Engineering and Technology - VNU Hanoi, Vietnam
*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <getopt.h>
//#include <limits.h>

//#define EXEC_FILE
#define UNSIGNED_SIZE 	4
#define MAX_QUERY		2000000

//build the graph
void build_graph(unsigned num_nodes, unsigned num_edges);

//add an edge to the graph
static inline void add_edge(unsigned src, unsigned dst);

//delete an edge from the graph
static inline void del_edge(unsigned src, unsigned dst);

//compute the shortest distance from src to dst
static inline int shortest_distance(unsigned src, unsigned dst);

//execute all consecutive queries
static inline void 	exec_queries();

//graph data
unsigned  max_edges = 80000000, max_nodes=20000000, max_2nodes;//100M nodes
//50M nodes => needs 2x4x4x50M = 1,6G for in/out queues with CPU 4T; 50M for maps
//50M nodes => needs 2x4x24x50M = 9,6G for in/out queues with CPU 24T; 300M for maps
//100M nodes => needs 2x4x24x100M = 2x9.6G for in/out queues with CPU 24T; 600M for maps

unsigned num_edges;				//num of edges
unsigned num_nodes, num_2nodes;				//num of nodes
//infos for alls incoming nodes of each vertex
unsigned *incoming_edges;		//list of incoming nodes for all vertex
unsigned *incoming_index;			//num of incoming nodes for each vertex
// o => incoming_start, 1 => num_incoming
//unsigned *incoming_starts;		//location of incoming nodes in the list

//infos for alls outoging nodes of each vertex
unsigned *outgoing_edges;		//list of outgoing nodes for all vertex
unsigned *outgoing_index;			//num of outgoing nodes for each vertex
// 0 => outgoing_starts , 1 => num_outgoing
//unsigned *outgoing_starts;		//location of outgoing nodes in the list

//for handling query 
unsigned *gmap; 				//global buffer for map; used for marking all travelled vertexs
unsigned *gin_queue; 			//incoming queue
unsigned *gout_queue;			//outgoing queue

unsigned long query_params[MAX_QUERY*2];// 0 -> query index, 1 -> src; 2 -> dst
unsigned num_query;
int 	 distances[MAX_QUERY];

unsigned num_threads;// number of CPU cores
unsigned map_isize;// size of map in 4bytes

char actions[MAX_QUERY];
unsigned long vertices[MAX_QUERY];
unsigned query_index[MAX_QUERY];
unsigned num_actions = 0;

//set bit at index in the array addr
static inline void set_bit(unsigned index, unsigned *addr)
{
	unsigned p = index >> 5, x = index & 31;
	addr[p] |= 1 << x;
}

//test bit at index in the array addr
static inline int test_bit(unsigned index, const unsigned *addr)
{
	unsigned p = index >> 5, x = index & 31;
	return (addr[p] >> x) & 1;
}

//clear map at positions specified by queue array 
//__attribute__((vector))
static inline void reset_map(unsigned *map, unsigned* queue, unsigned top)
{
	for (unsigned i = 0; i < top; i++)	map[queue[i]>>5] = 0;
	//map[queue[0:top]>>5] = 0;
}

//compare qsort
int node_cmp(const void * a, const void * b)
{
	return (*(unsigned*)a - *(unsigned*)b);
}

//build_graph graph from ini-file
//build_graph graph from ini-file
void build_graph(unsigned max_nodes, unsigned max_edges)
{
	unsigned *src, *dst;
	//allocate all buffer
	incoming_index = (unsigned *)aligned_alloc(16,UNSIGNED_SIZE*max_nodes*2);
	memset(incoming_index, 0, UNSIGNED_SIZE*max_nodes*2);
	
	outgoing_index = (unsigned *)aligned_alloc(16,UNSIGNED_SIZE*max_nodes*2);
	memset(outgoing_index, 0, UNSIGNED_SIZE*max_nodes*2);
	//edges 
	src = (unsigned *)aligned_alloc(16, UNSIGNED_SIZE*max_edges);
	dst = (unsigned *)aligned_alloc(16, UNSIGNED_SIZE*max_edges);

	char *line = NULL;
	unsigned n =0, max = 0, maxin=0, maxout=0;
	int res;
	size_t bufsize=0;
	
	//fetch the list of edges
	FILE *fp = stdin;
#ifdef EXEC_FILE	
	fp = fopen("init-file.txt", "r");	
	//fp = fopen("soc-pokec-relationships.txt", "r");  
#endif
	while (1) {
		res = getline(&line,&bufsize,fp);
		if (res == -1) break;
		if (line[0] == 'S') break;
	
		res = sscanf(line, "%u %u", &src[n], &dst[n]);
		if ( !res || res == EOF ) {
			continue;
		} 
		if (src[n]>max) max = src[n];if (dst[n]>max) max = dst[n];
		
		//Count all outgoing/incoming for each node
		++outgoing_index[(src[n] << 1) + 1];
		++incoming_index[(dst[n] << 1) + 1];//*2
		
		++n;
	}
	//fclose(fp);
	free(line);	
	num_edges = n; num_nodes = max+1; num_2nodes = num_nodes << 1;
	
	max = 0;
	incoming_index[0] = outgoing_index[0] = 0;
	
	//compute the position of incomings/outgoings from a node in array of edges
	for (unsigned i = 0; i< max_nodes-1; i++){	
		n = i<< 1;	
		incoming_index[n+2] = incoming_index[n] + incoming_index[n+1];
		// if ((i>BUCKET_INDEX) && ((i & BUCKET_MODULO) == 0)) {
		// 	incoming_index[n+2] += BUCKET_SIZE;
		// }
		
		outgoing_index[n+2] = outgoing_index[n] + outgoing_index[n+1];
		// if ((i>BUCKET_INDEX) && ((i & BUCKET_MODULO) == 0)){
		// 	outgoing_index[n+2] += BUCKET_SIZE ;
		// }	
		if (maxin < incoming_index[n+1]) maxin = incoming_index[n+1];
		if (maxout < outgoing_index[n+1]) maxout = outgoing_index[n+1];	
	}

	//allocate buffer for edges
	incoming_edges = (unsigned *)aligned_alloc(16, UNSIGNED_SIZE * (incoming_index[(max_nodes<<1)-2] + 4000000  ));//add more 100.000 nodes
	outgoing_edges = (unsigned *)aligned_alloc(16, UNSIGNED_SIZE * (outgoing_index[(max_nodes<<1)-2]+ 4000000));//add more 100.000 nodes	
	
	//assign incomings/outgoings for each node
	unsigned *ni = 	(unsigned *)aligned_alloc(16, UNSIGNED_SIZE * num_nodes);
	memset(ni, 0, UNSIGNED_SIZE * num_nodes);
	unsigned *no = 	(unsigned *)aligned_alloc(16, UNSIGNED_SIZE * num_nodes);
	memset(no, 0, UNSIGNED_SIZE * num_nodes);
	for (unsigned i =0; i< num_edges; i++){
		unsigned u = src[i]<<1, v = dst[i] << 1;
		incoming_edges[ incoming_index[v] + ni[dst[i]] ] = src[i];
		++ni[dst[i]];

		outgoing_edges[ outgoing_index[u] + no[src[i]] ] = dst[i];
		++no[src[i]];		
	}
	free(ni); free(no); free(src); free(dst);
	//sort incomings/outgoings for each node
	for (unsigned i = 0; i < num_nodes; i++) {
		n = i<< 1;
		qsort(incoming_edges + incoming_index[n], incoming_index[n+1] , UNSIGNED_SIZE, node_cmp);
		qsort(outgoing_edges + outgoing_index[n], outgoing_index[n+1] , UNSIGNED_SIZE, node_cmp);
	}

	//remove duplicate	
	for (unsigned i = 0; i < num_nodes; i++) {
		unsigned i2 = i<<1;
		n = outgoing_index[i2]; max = n + outgoing_index[i2+1];
		for (unsigned j = n+1; j < max; j++){
			if (outgoing_edges[j] != outgoing_edges[j-1]){
				++n;
				outgoing_edges[n] = outgoing_edges[j];
			}			
		}
		if (n < max-1 )	{
			outgoing_index[i2+1] = n - outgoing_index[i2] + 1;
		}
		
		n = incoming_index[i2]; max = n + incoming_index[i2+1];
		for (unsigned j = n+1; j < max; j++){
			if (incoming_edges[j] != incoming_edges[j-1]){
				++n;
				incoming_edges[n] = incoming_edges[j];
			}			
		}
		if (n < max-1 )	{
			incoming_index[i2+1] = n - incoming_index[i2]+1;
		}
	}

	//compute the size of map for marking travelled nodes
	//map size in 4bytes
	map_isize = max_nodes / 32 + (max_nodes % 32 ==0? 0: 1);
	unsigned long MAP_BUFFER = (unsigned long)num_threads *2*map_isize*UNSIGNED_SIZE ;
	gmap = (unsigned *)aligned_alloc(16, MAP_BUFFER);
	memset(gmap, 0, MAP_BUFFER);
	
	unsigned long queuesize = (unsigned long) UNSIGNED_SIZE *max_nodes *num_threads;
	gin_queue = (unsigned *)aligned_alloc(16,queuesize);
	gout_queue = (unsigned *)aligned_alloc(16,queuesize);	

	//statistics
	fprintf(stderr, "Graph n=%u e=%u MI=%u MO=%u MS=%u GMS=%lu QS=%lu\n", num_nodes, num_edges,maxin, maxout, map_isize, MAP_BUFFER, queuesize);
	//fprintf(stderr, "InOut1 %u ,In0 %u, In1 %u, Out0 %u, Out1 %u\n", num_inout1, num_in0, num_in1, num_out0, num_out1 );		
	//fprintf(stderr, "end %u %u %u\n", incoming_index[num_nodes<<1],incoming_index[(num_nodes<<1)+1], incoming_index[(num_nodes<<1) +2] );
	//exit(1);
}

//check if vertex is in the graph
//__attribute__((vector))
static inline int find_vertex(unsigned *nodes, unsigned num, unsigned v)
{
	for (unsigned i = 0; i < num; i++){
		if (nodes[i] == v)	{	
			return i;
		}		
	}
	return -1;	
}

//add a new edge
//__attribute__((vector))
static inline void add_edge(unsigned src, unsigned dst)
{	
	unsigned u = src << 1, v = dst << 1, u1=u+1, v1=v+1;
	//Verify if this edge already in the graph
	if (incoming_index[v1] < outgoing_index[u1]){
		if (find_vertex(incoming_edges + incoming_index[v], incoming_index[v1], src) > -1) 
			return;
	}
	else {
		if (find_vertex(outgoing_edges + outgoing_index[u], outgoing_index[u1], dst) > -1) 
			return;		
	}
	
	//exec all query in queue
	if (num_query > 0) exec_queries(); 

	//unsigned  max =  outgoing_index[u] + outgoing_index[u1];
	unsigned num_index = outgoing_index[u1];
	memcpy(outgoing_edges + outgoing_index[num_2nodes], 
				outgoing_edges + outgoing_index[u],
					(num_index<<2));
	outgoing_index[u] = outgoing_index[num_2nodes];
	outgoing_edges[outgoing_index[u] + num_index] = dst;
	outgoing_index[u1] += 1;
	outgoing_index[num_2nodes] = outgoing_index[u] + outgoing_index[u1];

	num_index = incoming_index[v1];
	memcpy(incoming_edges + incoming_index[num_2nodes], 
			incoming_edges + incoming_index[v],
				(num_index << 2));
	incoming_index[v] = incoming_index[num_2nodes];
	incoming_edges[incoming_index[v] + num_index] = src;
	incoming_index[v1] += 1;
	incoming_index[num_2nodes] = incoming_index[v] + incoming_index[v1]; 
}
//remove an edge in the graph
//__attribute__((vector))
static inline void del_edge(unsigned src, unsigned dst)
{	
	unsigned u = src << 1, v = dst << 1, u1=u+1, v1=v+1;
	if (outgoing_index[u1]< incoming_index[v1]){
		int n;
		unsigned *p =  outgoing_edges + outgoing_index[u];
		//verify if this edge is in the graph or not
		n = find_vertex(p, outgoing_index[u1], dst);
		if (n<0) return; //quit if this edge is not found

		//do the search if having queries
		if (num_query > 0) exec_queries();  

		--outgoing_index[u1];
		unsigned  size = (outgoing_index[u1] - n); 
		p = p+n;

		for (unsigned i=0; i<size; i++) 
			p[i] = p[i+1];
		
		p = incoming_edges + incoming_index[v];
		n = find_vertex(p, incoming_index[v1], src);
		--incoming_index[v1];
		size = (incoming_index[v+1] - n); 
		p = p+n;

		for (unsigned i=0; i<size; i++) 
			p[i] = p[i+1];

	} else {
		int n;
		unsigned *p = incoming_edges + incoming_index[v];
		n = find_vertex(p, incoming_index[v1], src);
		if (n<0) return; //quit if this edge is not found
		
		if (num_query > 0) exec_queries();  
		
		--incoming_index[v1];
		unsigned size = (incoming_index[v1] - n); 
		p = p+n;

		for (unsigned i=0; i<size; i++) 
			p[i] = p[i+1];

		p =  outgoing_edges + outgoing_index[u];	
		n = find_vertex(p, outgoing_index[u1], dst);	
		--outgoing_index[u1];
		size = (outgoing_index[u1] - n); 
		p = p+n;

		for (unsigned i=0; i<size; i++) 
			p[i] = p[i+1];
		//p[0:size] = p[1:size];
	}
}
//do a exec_queries
//__attribute__((vector))
static inline void exec_queries()
{	
	cilk_for(unsigned i=0; i< num_query; i++){	
		long uv = vertices[query_index[i]];		
		distances[i]  = shortest_distance((unsigned)(uv >> 32), (unsigned) uv);			
	}
	//write out results
	for (unsigned i=0; i< num_query; i++) fprintf(stdout, "%d\n", distances[i]);
	num_query = 0;	
}

//find shortest length from src to dst, bitmap is firstly used in the global buffer resetted to 0
// if all global buffer is consumed, it is then resetted by using last queues 
//__attribute__((vector))
static inline int shortest_distance(unsigned src, unsigned dst)
{	
	if (src == dst) return 0;	
	unsigned u = src << 1, v= dst << 1;//for indexing
	unsigned out_total = outgoing_index[u+1], in_total = incoming_index[v+1];
	if ((out_total == 0) || (in_total == 0) )  return -1;

	unsigned 	current_in=1, current_out=1, 
				last_in = 0, last_out = 0, 
				in_cost = 0, out_cost = 0,
				count = 1, 
				i,j, //for index
				e, q, // edge & queue
				top,  max,				
				tid = __cilkrts_get_worker_number();
	unsigned long queuesize =  (unsigned long) tid * max_nodes;
	unsigned *in_queue = gin_queue + queuesize;
	unsigned *out_queue = gout_queue + queuesize;	
	
	unsigned *in_map, *out_map;
	in_map = gmap + (map_isize << 1)*tid;
	out_map = in_map + map_isize;	
		
	//init in-queue	
	in_queue[0] = dst;
	set_bit(dst, in_map);

	//init out-queue
	out_queue[0] = src;	
	set_bit(src, out_map);
	
	while (count>0) {		
		if (out_total < in_total)
		{	
			out_total = 0; 
			out_cost += 1; 
			top = current_out;
			for (i=last_out; i<top; i++)
			{
				q = out_queue[i]<<1; max = outgoing_index[q] + outgoing_index[q+1];
				for (j=outgoing_index[q]; j < max; j++) 
				{	
					e = outgoing_edges[j];
					if (!test_bit(e, out_map) ) {						
						if (e == dst) 
						{
							reset_map(out_map, out_queue, current_out);
							reset_map(in_map, in_queue, current_in);
							return out_cost;
						}
						if (test_bit(e, in_map)) 
						{
							reset_map(out_map, out_queue, current_out);
							reset_map(in_map, in_queue, current_in);
							return in_cost + out_cost;
						}
						out_queue[current_out++] = e;				
						set_bit(e, out_map); 
						e = (e << 1) + 1;
						out_total +=  outgoing_index[e] ;                    
					}				
				}//end for outgoing_edges					
			}//end for out_queue
			last_out = top; 
			count = current_out - last_out;
			out_total += count;			
		}
		else 
		{
			in_total = 0; 
			in_cost += 1;
			top = current_in;
			for (i=last_in; i<top; i++)
			{
				q = in_queue[i] << 1;max = incoming_index[q]  + incoming_index[q+1];
				for (j=incoming_index[q]; j < max; j++) 
				{
					e = incoming_edges[j];
					if (!test_bit(e,in_map)) {						
						if (e == src) 
						{
							reset_map(out_map, out_queue, current_out);
							reset_map(in_map, in_queue, current_in);
							return in_cost;												
						}
						if (test_bit(e, out_map)) 
						{
							reset_map(out_map, out_queue, current_out);
							reset_map(in_map, in_queue, current_in);
							return in_cost + out_cost;
						} 
						in_queue[current_in++] = e;				
						set_bit(e, in_map);
                        e = (e << 1) + 1;	
						in_total += incoming_index[e] ;
					}				
						
				}// end for imcomings				
			}// end for queue
			last_in = top;			
			count = current_in - last_in;
			in_total += count;
		}
		
	}
	reset_map(out_map, out_queue, current_out);
	reset_map(in_map, in_queue, current_in);
	return -1;
}

//free all allocated buffers
void free_buffer()
{
	free(outgoing_index); free(outgoing_edges); 
	free(incoming_index); free(incoming_edges);
	free(gin_queue); free(gout_queue); 
	//free(actions); free(src); free(dst);
	free(gmap);
}

//usage info
void usage() {
	printf("Usage: akgroup [options]\n");
	printf("Program Options:\n");
	printf("  -n  --max_nodes <N>  Graph has maximum N nodes - Default N = 20.000.000\n");
	printf("  -e  --max_edges <E>  Graph has maximum E edges - Default E= 20.000.000\n");
	printf("  -?  --help           This is program to perform actions on a large-scale directed graph.\n");
}

static inline void run_batch() {
	num_query = 0;
	//unsigned n,m;
	long uv;
	for (unsigned i = 0; i < num_actions; i++) {
		//n = i<<1;
		if(actions[i] == 'A'){
			uv= vertices[i];
			add_edge( (unsigned) (uv >> 32), (unsigned) uv );	
		}else if(actions[i] == 'D'){
			uv= vertices[i];
			del_edge( (unsigned) (uv >> 32), (unsigned) uv);	
		}else if(actions[i] == 'Q'){
			//query_params[num_query] = vertices[i];
			query_index[num_query] = i;
			++num_query;
		}
	}
	num_actions = 0;
}

//main function
int main(int argc, char** argv) {
	//checking to run parallel 
	num_threads = __cilkrts_get_nworkers();
	//__cilkrts_set_param("nworkers",nworkers);
	
	int opt;
	static struct option options[] = {
		{"max_nodes", optional_argument, 0, 'n'},
		{"max_edges", optional_argument, 0, 'e'},
		{"help", optional_argument, 0, '?'},
		{0,0,0,0}
	};

	while ((opt = getopt_long(argc, argv, "n:e:?", options, NULL)) != EOF) {
		switch (opt) {
		case 'n':
				max_nodes = atoi(optarg);
				break;
		case 'e':
				max_edges = atoi(optarg);
				break;
		case '?':default:
			usage();
			return 1;
		}
	}
	// end parsing of commandline options

	if (optind > argc) {
		usage();
		return 1;
	}
	
	//ingesting paragraph
	build_graph(max_nodes, max_edges);
	
	//starting workload
	char *line = NULL;
	size_t linesize=0;
	int res;
	//unsigned n;
	unsigned long u, v;	
	//stdout buffer
	char out_buff[100000];memset(out_buff, 0, 100000);
	setvbuf(stdout, out_buff, _IOFBF, 100000);
	//actions = (char*)malloc(1000000);	src = (unsigned*)malloc(1000000);	dst = (unsigned*)malloc(1000000);
	num_actions = 0;
	num_query = 0;max_2nodes = max_nodes << 1;
	FILE *fp = stdin;
#ifdef EXEC_FILE	
	fp = fopen("workload-file.txt", "r");	
	//fp = fopen("workload-pokec.txt", "r");	
#endif
	fprintf(stdout,"R\n");	fflush(stdout);
	while (1) {
		res = getline(&line,&linesize,fp);
		if (res < 1) {
			run_batch();
			exec_queries();		
			fflush(stdout);	
			break;
		}
		if (line[0] == 'F'){
			run_batch();
			exec_queries();				
			fflush(stdout);
			continue;
		}	
		res = sscanf(line, "%c %lu %lu", &actions[num_actions], &u, &v);
		vertices[num_actions] = (u << 32) | v; 
		num_actions++;

	}	//end of while
	free_buffer(); 
	return 0;
}
