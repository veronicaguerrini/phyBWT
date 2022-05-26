#include "Tools.h"

#if DEBUG
	#ifndef CK
		#define CK 0
	#endif

	#ifndef CK_2
		#define CK_2 1
	#endif
#endif

#ifndef PRINT_TABLE
	#define PRINT_TABLE 1
#endif

#ifndef USE_LIST
	#define USE_LIST 1
#endif



dataTypelenSeq k_min=16; //minLCP
double tau=1.0/2.0;
uint max_it=7;

dataTypeNChar n=0; //eBWT size

dataTypelenSeq max_minLCP=0;

              //A B C D E F G H I J K L M N O P Q R S T U V W X Y
int ORD[25] = {0,0,1,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,3,0,0,0,0,0};
#define ord(c) (ORD[c-65])
//terminator character at the end of the reads
char TERM = '#';

dataTypeNSeq numData=0;

vector<string> NAME;
vector<pair<string,leaves>> to_build_newick;

int Close(dataTypeNChar ind_start, dataTypeNChar ind_fin, dataTypeNSeq *bufferCDA,dataTypedimAlpha *buffereBWT, std::bitset<SIZE_BITSET> &isIn, dataTypeNChar *n_clust){
	
	//ind_start: starting index of the eBWT positional cluster
	//ind_fin: final index (NOT included)
	//isIn stores colors in the cluster
	//B_Symb stores symbols in the cluster
	std::bitset<ALF> B_Symb;
	
	for(dataTypeNChar index=ind_start; index<ind_fin; index++){
		isIn[bufferCDA[index]]=1;
		if(buffereBWT[index]==TERM){
			B_Symb[ALF-1]=1;
		}
		else{
			B_Symb[ord(buffereBWT[index])]=1;
		}
	}
	
	//Check positional cluster maximal
	if(B_Symb.count()>1){
		#if CK 
			fprintf(stderr,"Close return 1 -->\nCDA:\t");				
			for(dataTypeNChar index=ind_start; index<ind_fin; index++)
				fprintf(stderr,"%d", (int)bufferCDA[index]);
			fprintf(stderr,"\neBWT:\t");
			for(dataTypeNChar index=ind_start; index<ind_fin; index++)
				fprintf(stderr,"%c",buffereBWT[index]);
			fprintf(stderr,"\n");
		#endif
		(*n_clust)++;
		return 1;
	}
	else
		return 0;
}

int UpdateTable(vector<dataTypeNChar> &partial_sum,dataTypelenSeq w, std::bitset<SIZE_BITSET> &isIn, leaves &inSet){
	//w is the weight to use for updating partial_sum
	#if CK
		fprintf(stderr,"\nUpdateTable --> isIn: ");				
		for(dataTypeNSeq s=0; s<inSet.size(); s++)
			fprintf(stderr,"%d,", (int)isIn[s]);
		fprintf(stderr,"\n-->c_weight=%d\n",w);
	#endif
	
	//B_vect stores which elements of inSet are in the cluster
	std::bitset<SIZE_BITSET> B_vect;
	
	//the cluster votes the candidate part if ALL the elements appear over the threshold tau
	bool vote=true;
	double abund=0.0;

	//scan inSet
	dataTypeNSeq k=0;
	while ( vote && k<inSet.size()){
		dataTypeNSeq card=inSet[k].count();
		dataTypeNSeq i=(isIn&inSet[k]).count();
		abund=(double)i/(double)card;
		#if CK
			fprintf(stderr,"\t--> inSet[%d]: ",k);				
			for(dataTypeNSeq s=0; s<numData; s++)
				fprintf(stderr,"%d,", (int)inSet[k][s]);
			fprintf(stderr,"\n-->abund=%lf\n",abund);
		#endif
		if(abund>=tau)
			B_vect[k]=1;
		else if(abund>0)
			vote=false;

		k++;
	}
	#if CK
		fprintf(stderr,"vote -->%d\n",vote);
	#endif
	
	//Increase table partial_sum
	if(vote && (B_vect.count()<inSet.size()) && (B_vect.count()>1)) {
		dataTypeNChar index=B_vect.to_ulong();	
		if(max_minLCP<w)	max_minLCP=w;
		partial_sum[index]+=w; //Increase vector partial_sum
		return 1;
	}
	else return 0;
	
}

dataTypeNChar restrictDS(dataTypelenSeq *out_lcp, dataTypeNSeq *out_cda,dataTypedimAlpha *out_ebwt, dataTypelenSeq *in_lcp,dataTypeNSeq *in_cda,dataTypedimAlpha *in_ebwt,vector<dataTypeNSeq> colors){
	
	dataTypeNChar numchar=0; //output variable
	
	dataTypelenSeq inherited_lcp=MAX_LCP;
	bool prev_removed=false;
	
	std::pair<std::vector<dataTypeNSeq>::iterator,std::vector<dataTypeNSeq>::iterator> bounds;
	
	//Scan in_DS
	for(dataTypeNChar indexbuffer=0; indexbuffer<n; indexbuffer++){
		bounds=std::equal_range (colors.begin(), colors.end(), in_cda[indexbuffer]);  
		//if in_cda[indexbuffer] IS present--> write out_cda, out_ebwt and out_lcp
		if(bounds.second-bounds.first==1){
			out_cda[numchar]=in_cda[indexbuffer];
			out_ebwt[numchar]=in_ebwt[indexbuffer];
			//Change the LCP entry
			if(prev_removed)	out_lcp[numchar]=(in_lcp[indexbuffer]<inherited_lcp)?in_lcp[indexbuffer]:inherited_lcp;
			else	out_lcp[numchar]=in_lcp[indexbuffer];
			numchar++;
			prev_removed=false;
			inherited_lcp=MAX_LCP;
		}	
		//if in_cda[indexbuffer] is NOT --> do NOT write, but SAVE inherited_lcp
		else{
			prev_removed=true;
			inherited_lcp=(in_lcp[indexbuffer]<inherited_lcp)?in_lcp[indexbuffer]:inherited_lcp;
		}
	}
	
	return numchar;
}

dataTypeNChar ClusterAnalysis(vector<dataTypeNChar> &partial_sum,leaves &in_set,dataTypelenSeq* bufferLCP, dataTypeNSeq *bufferCDA,dataTypedimAlpha *buffereBWT, dataTypeNChar *n_clust){
	
	dataTypeNChar nClusters=0;
	dataTypeNChar numcharLCP=n;
	
	dataTypelenSeq *b_lcp;
	dataTypeNSeq *b_da;
	dataTypedimAlpha *b_ebwt;
	//Ordered vector of CDA elements in in_set
	vector<dataTypeNSeq> colors;
	for(dataTypeNSeq i=0; i<in_set.size(); i++){
		for(dataTypeNSeq s=0; s<numData; s++)
			if(in_set[i][s]==1) colors.push_back(s);
	}
	sort(colors.begin(),colors.end());
	
	#if CK
		fprintf(stderr,"ClusterAnalysis --> in_set.size=%lu, colors=[",in_set.size());
		for(dataTypeNSeq i=0; i<colors.size(); i++)
			fprintf(stderr,"%d,", (int)colors[i]);
		fprintf(stderr,"]\n");
	#endif
	
	if(colors.size()<numData){
		//Create restricted LCP array,eBWT and CDA 
		b_lcp = (dataTypelenSeq*) malloc (sizeof(dataTypelenSeq)*n);
		b_da = (dataTypeNSeq*) malloc (sizeof(dataTypeNSeq)*n);
		b_ebwt = (dataTypedimAlpha*) malloc (sizeof(dataTypedimAlpha)*n);
		
		numcharLCP = restrictDS(b_lcp,b_da,b_ebwt,bufferLCP,bufferCDA,buffereBWT,colors);
	}
	else
	{
		b_lcp=bufferLCP;
		b_da=bufferCDA;
		b_ebwt=buffereBWT;
	}
	
	
	bool started=false;
	dataTypelenSeq pred = -1;
	bool grow = true;
	dataTypeNChar tail_len = 0;
	dataTypeNChar pos_init = 0;
	dataTypelenSeq c_weight=k_min;
		
	//Scan LCP to detect maximal eBWT positional clusters
	for(dataTypeNChar indexbuffer=1; indexbuffer<numcharLCP; indexbuffer++){		
		if(b_lcp[indexbuffer]>=k_min) {
			if (not started){
				started=true, grow=true;
				pos_init=indexbuffer-1;
				c_weight=b_lcp[indexbuffer];
				tail_len=0;
			}
			else //cluster is open
			{
				//Case1: increasing sequence of LCP values
				if (grow){if(b_lcp[indexbuffer]<pred) grow = false;}
				else{//Case2: non increasing LCP values 
					if (b_lcp[indexbuffer]<pred) tail_len=0; 
					else if (b_lcp[indexbuffer]==pred) tail_len++;
					else { //close cluster and start a new one
						std::bitset<SIZE_BITSET> isIn;
						if(Close(pos_init,indexbuffer-1-tail_len,b_da,b_ebwt,isIn,n_clust)==1){	//Close returns 1 if the eBWT positional cluster is maximal							
							c_weight=(c_weight>b_lcp[indexbuffer-tail_len-2])?b_lcp[indexbuffer-tail_len-2]:c_weight; //set c_weight
							//Update table partial_sum
							nClusters+=UpdateTable(partial_sum,c_weight,isIn,in_set);
						}
						//Start a new cluster
						started=true, grow=true;
						pos_init=indexbuffer-1-tail_len;
						c_weight=b_lcp[indexbuffer-tail_len];
						tail_len=0;
					}
				}
			}//end-else if(not started)		
			pred=b_lcp[indexbuffer];
		}
		else if (started){ //close cluster
			std::bitset<SIZE_BITSET> isIn;
			if(Close(pos_init,indexbuffer,b_da,b_ebwt,isIn,n_clust)==1){
				c_weight=(c_weight>b_lcp[indexbuffer-1])?b_lcp[indexbuffer-1]:c_weight; 
				nClusters+=UpdateTable(partial_sum,c_weight,isIn,in_set);
			}
			started=false;
		}//end-else if(b_lcp[indexbuffer]>=k_min)
			
	}//end-for
			
	//Close a possibly last cluster 
	if(started){
		std::bitset<SIZE_BITSET> isIn;
		if(Close(pos_init,numcharLCP,b_da,b_ebwt,isIn,n_clust)==1){
			c_weight=(c_weight>b_lcp[numcharLCP-1])?b_lcp[numcharLCP-1]:c_weight; 
			nClusters+=UpdateTable(partial_sum,c_weight,isIn,in_set);
		}
	}	
	
	if(colors.size()<numData){
		delete [] b_lcp;
		delete [] b_da;
		delete [] b_ebwt;
	}
	
	return nClusters;
}

void SortTable(vector<dataTypeNChar> &v_in,vector<pair<dataTypeNChar,dataTypeNChar>> &v_out){
	//Scan table
	for(dataTypeNChar i=0; i<v_in.size(); i++){
		if(v_in[i]>0)	v_out.push_back(make_pair(v_in[i],i));
	}
	v_in.clear();
	sort(v_out.begin(), v_out.end()); //By default sorts on basis of first element of pairs (i.e. score).
}

#if PRINT_TABLE
void PrintTable(vector<pair<dataTypeNChar,dataTypeNChar>> &v_in, leaves &set_ele,FILE *OutC){
	for(dataTypeNChar i=v_in.size(); i>0; i--){
		std::bitset<SIZE_BITSET> B(v_in[i-1].second);//B bitset for v[i-1].second	
		for(dataTypeNSeq s=0; s<set_ele.size(); s++){
			if(B[s]==1){
				for(dataTypeNSeq k=0; k<set_ele[s].size(); k++){
					if(set_ele[s][k]==1)
						fprintf(OutC,"%s,", &NAME[k][0]);
				}
			}
		}
		fprintf(OutC,"\t%lu\n",v_in[i-1].first);
	}
}
#endif

void BuildOutputPartition(vector<pair<dataTypeNChar,dataTypeNChar>> &table,leaves &input_set, vector<leaves> &out_vect){
	
	dataTypeNChar num_it=0;
	bool invalid_set=false;
	bool stop=false;
	
	#if USE_LIST
		//list2 is a list of compatible extensions of the output list
		leaves list2;
	#endif
	
	//Bitvector storing partitioned elements of input_set
	//input_set.size=#elements to be partitioned;
	std::bitset<SIZE_BITSET> B;
	
	dataTypeNChar diff_curr=0;
	dataTypeNChar diff_next=0;
	if(table.size()>1)
		diff_next=table[table.size()-1].first-table[table.size()-2].first;
	
	
	while((!stop) && (num_it<table.size()) && (B.count()<input_set.size()-1)){
		//max_ele is the maximum combination set among the remaining
		dataTypeNChar max_ele=table[table.size()-1-num_it].second;
        
        //max_ele is a valid part for the output partition if
        //1. its elements are not already in a part (i.e. B[s]==0)
        //2. the first element in list2 (if any) with non-empty intersection is a superset of max_ele
        std::bitset<SIZE_BITSET> bit_max_ele(max_ele);
        leaves new_ele;
        
		#if CK_2
			fprintf(stderr,"findMax --> num_it=%lu, maxele=%lu, bit_max_ele=",num_it,max_ele);
			for(dataTypeNSeq s=0; s<input_set.size(); s++)
				fprintf(stderr,"%d,",(int)bit_max_ele[s]);
		#endif
		
		
		if((B&bit_max_ele).count()==0)
			invalid_set=false;
		else
			invalid_set=true;
        
		#if CK_2
			fprintf(stderr," invalid_set=%d\n",invalid_set);
		#endif
		
		#if USE_LIST
		if(invalid_set){//Check if bit_max_ele is a compatible extension 
			if ((B&bit_max_ele)==B){
				list2.push_back(bit_max_ele);
				#if CK_2
					fprintf(stderr,"Put bit_max_ele in list2, list2.size=%d\n",list2.size());
				#endif
			}
		}
		#endif
		
		if(not invalid_set){
			//Candidate part to be inserted into out_vect
			#if USE_LIST
			//Check if it is compatible
			bool takeIt=true;
			dataTypeNSeq i=0; 
			while(takeIt && i<list2.size()){
				dataTypeNSeq num=(bit_max_ele&list2[i]).count();
                if(num>0){
                    if(num<bit_max_ele.count())
                        takeIt=false;
                    i=list2.size();
                }
				else
					i++;
			}
			#if CK_2
				fprintf(stderr,"takeIt=%d\n",takeIt);
			#endif
			if(takeIt)
			#endif
			{//the part can be inserted
				for(dataTypeNSeq s=0; s<input_set.size(); s++){
					if(bit_max_ele[s]==1)
						new_ele.push_back(input_set[s]);
				}
				out_vect.push_back(new_ele);
				//modify B accordingly
				B=(B|bit_max_ele); //OR
				#if CK_2
				fprintf(stderr,"B= ");
				for(dataTypeNSeq s=0; s<numData; s++)
					fprintf(stderr,"%d,",(int)B[s]);
				fprintf(stderr,"\n");
				#endif
			}
		}
		
		num_it++;
		//Check break condition
		diff_curr=diff_next;
		if(num_it<table.size()-1){
			diff_next=table[table.size()-1-num_it].first-table[table.size()-2-num_it].first;
		}
		
		if((num_it>max_it) && (diff_curr<=diff_next)){
			stop=true;
		}
	}//end-while
	
	//Add singleton of input_set to our output partition
	for(dataTypeNSeq s=0; s<input_set.size(); s++){
		if(B[s]==0){
			leaves new_ele_singleton;
			new_ele_singleton.push_back(input_set[s]);
			out_vect.push_back(new_ele_singleton);
		}
		
	}
	
	fprintf(stderr,"END BuildOutputPartition --> out_vect.size=%lu\n",out_vect.size());
}
#if PRINT_TABLE
void Partition(leaves in_set,vector<leaves> &out_set,dataTypelenSeq* bufferLCP, dataTypeNSeq *bufferCDA,dataTypedimAlpha *buffereBWT, FILE *f_out){
#else
void Partition(leaves in_set,vector<leaves> &out_set,dataTypelenSeq* bufferLCP, dataTypeNSeq *bufferCDA,dataTypedimAlpha *buffereBWT){
#endif
	//partial_sum is the vector storing scores for each combination set
	vector<dataTypeNChar> partial_sum;
	dataTypeNChar sizeTable=pow(2,in_set.size());
	for(dataTypeNChar i=0; i<sizeTable; i++) 
		partial_sum.push_back(0);
	
	//Detect positional clusters and write Table partial_sum
	dataTypeNChar n_clust=0;
	dataTypeNChar n_clust_vote=ClusterAnalysis(partial_sum,in_set,bufferLCP,bufferCDA,buffereBWT,&n_clust);
	fprintf(stderr,"Partition -> Number maximal positional clusters=%lu, of which voting=%lu\n",n_clust,n_clust_vote);
	
	vector<pair<dataTypeNChar,dataTypeNChar>> table;
	SortTable(partial_sum,table);
	
	#if PRINT_TABLE
		fprintf(f_out,"printTable\n");
		PrintTable(table,in_set,f_out);
	#endif
	
	//Take the first top combination sets (as far as they are compatible) and create a partition
	BuildOutputPartition(table,in_set,out_set);
	
	fprintf(stderr,"END Partition --> out_set.size=%lu\n",out_set.size());
	
}

void UpdateEdges(leaves &P){
	//Create node
	string conc="(";
	bitset<SIZE_BITSET> vect;
	leaves bit_node; 
	dataTypeNSeq i=0;
	while(i<P.size()){
		dataTypeNSeq j=0;
		bool stop=false;
		while(j<to_build_newick.size() && !stop){
			if((to_build_newick[j].second)[0]==P[i]){
				conc+=to_build_newick[j].first;
				stop=true;
				dataTypeNSeq k=1;
				while((to_build_newick[j].second).size()>k){
					bit_node.push_back((to_build_newick[j].second)[k]);
					k++;
				}
			}
			else
				j++;
		}
		if(stop)
			to_build_newick.erase(to_build_newick.begin()+j);
		else{
			conc+="*";
			bit_node.push_back(P[i]);
		}
		vect=(vect|P[i]);
		conc+=",";
		i++;
	}
	conc.pop_back();
	conc+=")";
	bit_node.insert(bit_node.begin(),vect);
	to_build_newick.push_back(make_pair(conc,bit_node));
}

int main(int argc, char **argv) {
  
	fprintf(stderr,"USE_LIST %d\n",USE_LIST);
	
	//time_t t_refine=0, t_total=0;
	time_t t_total=0;
    //clock_t c_refine=0, c_total=0;
    clock_t c_total=0;
	
	if( argc < 4) {
		fprintf(stderr,"Error usage %s fileFasta fileInfo fileOutput k_min tau t\n",argv[0]); 
		exit(1);
	}

	string fileFasta=argv[1];
	string fileInfo=argv[2];
	string output_new=argv[3];

	sscanf(argv[4], "%hhu", &k_min);
	sscanf(argv[5], "%lf", &tau);
	sscanf(argv[6], "%u", &max_it);
	
	fprintf(stdout,"k_min: %d\ntau: %lf\nt: %u\n",k_min,tau,max_it);
	
	fprintf(stderr,"sizeof(dataTypeNChar)=%lu bytes.\n",sizeof(dataTypeNChar));
	fprintf(stderr,"sizeof(dataTypelenSeq)=%lu bytes.\n",sizeof(dataTypelenSeq));
	fprintf(stderr,"sizeof(dataTypeNSeq)=%lu bytes.\n",sizeof(dataTypeNSeq));
	
	std::ifstream in_list;
	in_list.open(fileInfo.c_str(), std::ifstream::in);
	if (!in_list.is_open()){
		fprintf(stderr,"Error opening file %s\n",fileInfo.c_str());
		exit (EXIT_FAILURE);
	}
	
	while(getline(in_list,fileInfo,'\t')){
		NAME.push_back(fileInfo);
		getline(in_list,fileInfo,'\n');
		numData++;
	}
	
	fprintf(stdout,"numData: %u\n",numData);
	
    //Files .lcp, .ebwt and .cda 
	string fnLCP=fileFasta+".lcp";
	string fnCDA=fileFasta+".cda";
	string fnEBWT=fileFasta+".ebwt";
	
	FILE *InLCP= fopen(fnLCP.c_str(), "rb");
	FILE *InCDA= fopen(fnCDA.c_str(), "rb");
	FILE *InEBWT= fopen(fnEBWT.c_str(), "rb");
	
    if(InLCP==NULL){
		fprintf(stderr,"Error opening file %s\n",fnLCP.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    if(InCDA==NULL){
        fprintf(stderr,"Error opening file %s\n",fnCDA.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    if (InEBWT==NULL){
        fprintf(stderr,"Error opening file %s\n",fnEBWT.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
	}
	
	//Set the eBWT length
    fseek(InEBWT, 0, SEEK_END);
	n = ftell(InEBWT)/sizeof(dataTypedimAlpha);

	//To upload in main mamory LCP, DA and eBWT
	dataTypeNChar numcharLCP, numchar;
	dataTypelenSeq *bufferLCP;
	dataTypeNSeq *bufferCDA;
	dataTypedimAlpha *buffereBWT;
	
	bufferLCP = (dataTypelenSeq*) malloc (sizeof(dataTypelenSeq)*n);
	bufferCDA = (dataTypeNSeq*) malloc (sizeof(dataTypeNSeq)*n);
	buffereBWT = (dataTypedimAlpha*) malloc (sizeof(dataTypedimAlpha)*n);
	
	//Upload in main memory LCP/DA/eBWT file
	fseek(InLCP, 0, SEEK_SET);
	numcharLCP=fread(bufferLCP,sizeof(dataTypelenSeq),n,InLCP);
	fseek(InCDA, 0, SEEK_SET);
	numchar=fread(bufferCDA,sizeof(dataTypeNSeq),n,InCDA);
	assert(numcharLCP==numchar);
	fseek(InEBWT, 0, SEEK_SET);
	numchar=fread(buffereBWT,sizeof(dataTypedimAlpha),n,InEBWT);
	assert(numcharLCP==numchar);
	
	//Close files
	fclose(InLCP);
    fclose(InCDA);
    fclose(InEBWT);
	
	//Output file(s)
	FILE *OutNew= fopen(output_new.c_str(), "w");
    #if PRINT_TABLE
	std::stringstream ssout;
	ssout << "Table_k_" << (int)k_min << "_tau" << (int)(tau*10) << ".txt\0";
	
	string fnOUT=ssout.str();
	FILE *OutP= fopen(fnOUT.c_str(), "w");
    if ((OutP==NULL)){
		std::cerr << "Error opening " << fnOUT << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
    #endif
    
	time_start(&t_total, &c_total); //start time

	//We have a queue of subsets of colors that must be partitioned
	std::deque<leaves> job_list;
	
	//At the beginning the queue job_list has only one element --> the set of all the sequences as singletons
	leaves startSet;
	for(dataTypeNSeq i=0; i<numData; i++){
		bitset<SIZE_BITSET> o;
		o[i]=1;
		startSet.push_back(o);
		//Inizialize to_build_newick
		leaves oo;
		oo.push_back(o);
		to_build_newick.push_back(make_pair(NAME[i],oo));
	}
	job_list.push_back(startSet);
	
	dataTypeNSeq step=1;
	
	while(job_list.size()>0){
		//Pick the first element and remove it from the queue
		leaves P = job_list.front();
		job_list.pop_front();
		
		if(P.size()>1)//Call Partition
		{ 
			vector<leaves> outputPartition;
			if(P.size()>2){
                #if PRINT_TABLE
				Partition(P,outputPartition,bufferLCP,bufferCDA,buffereBWT,OutP);
                #else
                Partition(P,outputPartition,bufferLCP,bufferCDA,buffereBWT);
                #endif
				fprintf(stderr,"Partition step: %d -> outputPartition.size=%lu\n",step,outputPartition.size());
			}
			
			if(P.size()==2 || P.size()==outputPartition.size()){
				#if CK_2
					fprintf(stdout,"END step: %d -> (",step);
					dataTypeNSeq i=0;
					while(i<P.size()){
						fprintf(stdout,"(");
						for(dataTypeNSeq j=0; j<P[i].size(); j++){
							if(P[i][j]==1)	fprintf(stdout,"%s,",&NAME[j][0]);
						}
						fprintf(stdout,"),");
						i++;
					}
					fprintf(stdout,")\n");
				#endif
				//Update to_build_newick
				UpdateEdges(P);
			}
			else
			{
				leaves upperSubTree;
				for(dataTypeNSeq i=0; i<outputPartition.size(); i++){
					//Put each part in queue
					job_list.push_back(outputPartition[i]);
					
					//bitset part corresponds to the part in outputPartition[i]
					bitset<SIZE_BITSET> part;
					for(dataTypeNSeq j=0; j<outputPartition[i].size(); j++)
						part=(part|outputPartition[i][j]);
					
					upperSubTree.push_back(part);
				}
				
				//Put the whole partition in job_list
				job_list.push_back(upperSubTree);
			}
			
		}
		else{
			#if CK_2
				fprintf(stdout,"END step: %d -> Leaf: (",step);
				for(dataTypeNSeq i=0; i<P[0].size(); i++){
					if(P[0][i]==1)
					fprintf(stdout,"%s,",&NAME[i][0]);
				}
				fprintf(stdout,")\n");
			#endif
		}
		step++;
	}
	
	#if CK_2
		cout << "to_build_newick.size=" << to_build_newick.size() << endl;
	#endif
	
	while(to_build_newick.size()!=1){
        int last=to_build_newick.size()-1;
        assert((to_build_newick[last].second).size()==1);
        string to_sub=to_build_newick[last].first;
        bool sub=false;
        int i=last-1;
        while(i>=0 && !sub){
            if((to_build_newick[i].second[0]&to_build_newick[last].second[0])==to_build_newick[last].second[0]){
                uint num=1;
                while(num<(to_build_newick[i].second).size() && to_build_newick[last].second[0]!=to_build_newick[i].second[num])
                    num++;
                if(num<(to_build_newick[i].second).size()){
                    sub=true;
                    (to_build_newick[i].second).erase((to_build_newick[i].second).begin()+num);
                    std::size_t pos=0;
                    while(num>0 && pos!=std::string::npos){
                        pos = (to_build_newick[i].first).find("*",pos+1);
                        num--;
                    }
                    (to_build_newick[i].first).replace(pos,1,to_sub);
                }
            }
            i--;
        }
        assert(sub);
        to_build_newick.pop_back();
	}
	fprintf(OutNew,"%s;\n",&((to_build_newick[0].first)[0]));
    
    #if PRINT_TABLE
	//Close output file
	fclose(OutP);
    #endif
	fclose(OutNew);

    fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}

