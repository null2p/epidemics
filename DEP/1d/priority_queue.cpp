#include<algorithm>
#include<iostream>
#include<cmath>
#include<vector>
#include<string>

static inline std::string do_padding (unsigned index, unsigned mlength);
static inline void printer (std::vector<double> const & tree, unsigned index, unsigned mlength);
static inline void print_tree (std::vector<double> & tree);

struct Element{
	int reaction;
	double tau;
};

static inline bool tau_ascending(Element const& lhs, Element const& rhs){
	return lhs.tau < rhs.tau;
};


class IdxPriorityQueue{
	public:
		Element *heap_element;
		int * heap_index;
		int size;

		IdxPriorityQueue(int size, double * tau_arr) : size(size)
		{
			heap_element = new Element[size];
			heap_index = new int[size]{};

			for(int i =0; i<size ; i++){
			      	heap_element[i].tau = tau_arr[i];
				heap_element[i].reaction = i;
			}

			std::sort(heap_element,heap_element+size,&tau_ascending);

			//reaction number 'A' is on 'i th' site in queue ( heap_index[A] = i );
			for( int i = 0; i < size ; i ++ ){
				heap_index[heap_element[i].reaction] = i;
			}


		}
		~IdxPriorityQueue()
		{
			delete[] heap_element;
			delete[] heap_index;
		}

		void print_tau(){
			/*
			int next_line = 0;
			int number_cout = 0;
			std::cout<<"\t\t\t\t";
			for(int i=0;i<size;i++){
			      	std::cout << heap_element[i].tau<<" ";
				number_cout ++;
				if(number_cout == std::pow(2,next_line)){
				
					number_cout = 0;
					next_line ++;
					std::cout<<std::endl;
					for(int j=0; j<4-next_line; j++) std::cout<<"\t";
				}
			}
				std::cout<<std::endl;
				*/
			std::vector<double> tau_v(size);
			for(int i =0; i< size ; i ++) tau_v[i] = heap_element[i].tau;
			std::make_heap(tau_v.begin(), tau_v.end(), std::greater<double>());
			std::cout<<"\n";
			print_tree(tau_v);
		}
		void print_reaction(){
			for(int i=0;i<size;i++)  std::cout << heap_index[i]<<" ";
			std::cout<<std::endl;
		}

		void swap(int new_i, int new_j)
		{
			int tmp_reaction = heap_element[new_i].reaction;
			double tmp_tau = heap_element[new_i].tau;

			heap_element[new_i].reaction = heap_element[new_j].reaction;
			heap_element[new_i].tau = heap_element[new_j].tau;
			
			heap_element[new_j].reaction = tmp_reaction;
			heap_element[new_j].tau = tmp_tau;

			heap_index[heap_element[new_i].reaction] = new_i;
			heap_index[heap_element[new_j].reaction] = new_j;
			
		}
		void update(int i, double new_tau)
		{
			heap_element[i].tau = new_tau;
			update_aux(i);
		}
		void update_aux(int i)
		{
			heap_up(i);
			heap_down(i);
		}
		void heap_down(int i)
		{
			double tau_min_child = 987654321;
			int idx_min_child = size-1;

			int first_child = 2*i+1;
			int second_child = 2*i+2;

			if(first_child >= size) return;
			else if(second_child >= size){
				idx_min_child = first_child;
			       	tau_min_child = heap_element[first_child].tau;
			}	
			else{
				if(heap_element[first_child].tau > heap_element[second_child].tau){
					idx_min_child = second_child;
					tau_min_child = heap_element[second_child].tau;
				}
				else{
					idx_min_child = first_child;
					tau_min_child = heap_element[first_child].tau;
				}
			}

			if(heap_element[i].tau > tau_min_child)
			{
				swap(i, idx_min_child);
				heap_down(idx_min_child);
				//update_aux(idx_min_child);
			}

		}
		void heap_up(int i)
		{
			if( i == 0 ) return;

			int parent = (i-1)/2;

			if(heap_element[i].tau < heap_element[parent].tau){
				swap(i,parent);
				heap_up(parent);
				//update_aux(parent);
			}
		}


};

std::string do_padding (unsigned index, unsigned mlength){
  std::string padding;
  if (int((index-1)/2) != 0){
    return (int((index-1)/2) % 2 == 0) ?
    (do_padding(int((index-1)/2),mlength) + std::string(mlength+4,' ') + " ")  :
    (do_padding(int((index-1)/2),mlength) + std::string(mlength+3,' ') + " |") ;
  }
  return padding;
}


void printer (std::vector<double> const & tree, unsigned index, unsigned mlength){
  auto last = tree.size() - 1 ;
  auto  left = 2 * index + 1 ;
  auto  right = 2 * index + 2 ;
  std::cout << " " << tree[index] << " ";
  if (left <= last){
    auto llength = std::to_string(tree[left]).size();
    std::cout << "---" << std::string(mlength - llength,'-');
    printer(tree,left,mlength);
    if (right <= last) {
      auto rlength = std::to_string(tree[right]).size();
      std::cout << "\n" << do_padding(right,mlength) << std::string(mlength+ 3,' ') << " | ";
      std::cout << "\n" << do_padding(right,mlength) << std::string(mlength+ 3,' ') << " └" <<
      std::string(mlength - rlength,'-');
      printer(tree,right,mlength);
    }
  }
}


void print_tree (std::vector<double> & tree){
  unsigned mlength = 0;
  for (double & element : tree){
    auto clength = std::to_string(element).size();
    if (clength > mlength) {
      mlength = std::to_string(element).size();
    }
  }
  std::cout <<  std::string(mlength- std::to_string(tree[0]).size(),' ');
  printer(tree,0,mlength);
}



