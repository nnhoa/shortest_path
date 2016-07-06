# shortest_path
Compute the shortest distance on the large-scale, directed, unweighted dynamic graph

1. Team Name:
	AK-Group - Deparment of Information Systems, University of Engineering and Technology, Hanoi, Vietnam.
	
2. Team members:
	- Pham Hai Dang, dangph@vnu.edu.vn, University of Engineering and Technology - Vietnam National University, Information Technology, BSc.

	- DU Phuong Hanh, hanhdp@vnu.edu.vn, Deparment of Information Systems, University of Engineering and Technology - Vietnam National University, Vietnam, PhD Student.

	- Vu Ba Duy, duyvb@vnu.edu.vn, University of Engineering and Technology - Vietnam National University, MSc.

3. Advisors
	- Assoc. Prof. Nguyen Ngoc Hoa, hoa.nguyen@vnu.edu.vn, Deparment of Information Systems, University of Engineering and Technology - Vietnam National University, Vietnam.
	
4. Brief Explanation
	- Graph is organized as arrays of firstly sorted incoming & outcoming nodes; 

	- Each 32 nodes has a reserved buffer for adding new nodes, this buffer is called BUCKET. By default, BUCKET size is 8, that means we can add more atleast 8 nodes in the interval of 32 vertices

	- Searching algorithm is bi-directional; 2 bitmap arrays are used for remarking travelled incoming/outgoing nodes.

	- Global incoming queue and outgoing queue are used for each call of searching shortest length. Then, each searching thread will use only proper in/out queues determined by an interval of globals in/out queues.

	- Each searching thread will also own the proper in/out map slots, computed from the global in/out maps. Once the searching is finished, these in/out slots will be cleared for the next search.

	- Cilk is used for performing the searching algo in parallel of queries 

5. Third Party Code:
- Intel® Cilk for parallel processing.

 6. Usage: akgroup [options]
Program Options:
  -n  --max_nodes <N>  Graph has maximum N nodes - Default N = 20.000.000
  -e  --max_edges <E>  Graph has maximum E edges - Default E= 20.000.000
  -?  --help           This is program to perform actions on a large-scale directed graph.

 Constraints
	- Currently, we use the default maximum number of vertices is 20M. If the graph has the vertex bigger than 20M, please run this program with the -n parameter for changing the max-nodes.
