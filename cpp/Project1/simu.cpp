#include "simu.h"

double proba(int nei_sum, double p, double q, int spin){
    double gamma = 1-2*q;
    short f = 1;
    if(nei_sum == 0)
        f = 0;
    else if( nei_sum < 0)
        f = -1;
    else
        f = 1;
    return 0.5*(1-p)*(1-gamma*spin*f)+p*0.5;
}
std::vector<short> make_ran_lattice(int N)
{
	std::vector<short> lat(N);
	for_each(lat.begin(), lat.end(), [lat](short& b) mutable {b = (randZero_One_Uni()) < 0.5 ? -1 : 1; });
	return lat;
}

void simulationAverages(int N, int num_nei, int mcs, int av_num, double pstart,double pmax, int pnum, std::vector<double> qs,std::string type, bool plot_in_time ){
    // FULLY CONNECTED
    // N - lattice sites
    // num_nei - neighbors number
    // mcs - Monte Carlo steps
    // av_num - number of averages
    // pstart - starting probability p
    // pmax - ending p
    // pnum - number of ps
  
//#pragma omp parallel for
    std::for_each(std::execution::par, std::begin(qs),std::end(qs),[&](auto&& q)
	//int qs_num = qs.size();
	//#pragma omp parallel for
    //for(int qi=0; qi < qs_num; qi++)
	{
    	//double q = qs[qi];
        ofstream file;
    	string fname = type +",q="+ to_string(q) + ",N="+to_string(N)+ ",mcs" + to_string(mcs) + ",avNum"+to_string(av_num) + ".dat";
        file.open(fname);
    	if(!file.is_open())
    	{
    		cout << "couldn't open a file" << endl;
    	}
        double p = pstart;
        double pstep = 1.0*(pmax - p)/pnum;
        if(p == pmax)
            pstep = 1;
        std::cout<< "-------------- starting simulation for q = " + std::to_string(q) + "--------------" << std::endl;
    	//#pragma omp parallel while reduction(+:av_abs_m) private(p)
        while(p <= pmax){
            double av_abs_m = 0;
        	double av_m = 0;
        	double av_m2 = 0;
        	double av_m4 = 0;
        	int m =0;
            for(int av = 0; av < av_num;av++){
                std::vector<short> nums = (make_ran_lattice(N));
                //m = std::accumulate(nums.begin(), nums.end(), 0);	
                for(int mc = 0; mc < mcs*N; mc++){
                    const int spin_num = randInt_Uni(N-1); // choose the spin to propose a flip
                    int spin = nums[spin_num];
                    int nei_sum = 0;
                    for(int j = 0; j<num_nei;j++){
                        int nei_choice = spin_num;
                        //print("spin_num = " + str(spin_num))
                        while(nei_choice == spin_num){
                        nei_choice = randInt_Uni(N-1); // choose the spin to propose a neighbor
                            //print(nei_choice)
                        }
                        nei_sum += nums[nei_choice];
                    }
                    //print(nei_sum)
                    double prob = proba(nei_sum,p,q,spin);
                    if(randZero_One_Uni() <= prob){
                        nums[spin_num] = -spin;
                        //m -= 2*spin;
                    }
                }
                m = std::accumulate(nums.begin(), nums.end(), 0);
            	av_abs_m += abs(1.0*m/N);
            	av_m += 1.0*m/N;
            	av_m2 += ((1.0*m/N)*(1.0*m/N));
            	av_m4 += (1.0*m/N)*(1.0*m/N)*(1.0*m/N)*(1.0*m/N);
            }
        	av_abs_m /= av_num;
            av_m /= av_num;
            av_m2 /= av_num;
            av_m4 /= av_num;
            file << p << "\t" << (av_abs_m) << "\t" << N*(av_m2-av_abs_m*av_abs_m) << "\t" << 1-(av_m4/(3*av_m2*av_m2)) <<  endl;
        	//std:: cout << "for q = " + std::to_string(q) + " p = " + std::to_string(p) << " ->  ->  -> m = " +to_string(av_abs_m/N/av_num)  << std::endl;
            p +=pstep;
        }
        //------------------------------------------------------------------------------------            
        std::cout<< "-------------- finished simulation for q = " + std::to_string(q) + "--------------" << std::endl;
        //print(ps)
        //print(av_abs_ms)
    	file.close();
    });
	
}
void simulationAveragesNonInfty(int Lx,int Ly, int mcs, int av_num, double pstart,double pmax, int pnum, std::vector<double> qs,std::string type,std::string savename){
    // FULLY CONNECTED
    // N - lattice sites
    // num_nei - neighbors number
    // mcs - Monte Carlo steps
    // av_num - number of averages
    // pstart - starting probability p
    // pmax - ending p
    // pnum - number of ps
    int N = Lx*Ly;
    general::lattice2D* lat;
	if(type == "square") lat = new general::square_lattice(Lx,Ly);
    else lat = new general::triangle_lattice(Lx,Ly);
 
	

	
//#pragma omp parallel for
    std::for_each(std::execution::par, std::begin(qs),std::end(qs),[&](auto&& q)
	//int qs_num = qs.size();
	//#pragma omp parallel for
    //for(int qi=0; qi < qs_num; qi++)
	{
    	//double q = qs[qi];
        ofstream file;
    	ofstream fileTime;
    	string fname = savename+",q="+ to_string(q) + ",N="+to_string(N)+ ",mcs" + to_string(mcs) + ",avNum,NONINFTY"+to_string(av_num) + ".dat";
        file.open(fname);
    	if(!file.is_open())
    	{
    		cout << "couldn't open a file" << endl;
    	}
    	file << "p" << "\t" << "(av_abs_m)" << "\t"<<"susceptibility"<<"\t" << "binder cumulant" << endl;
        double p = pstart;
        double pstep = 1.0*(pmax - p)/pnum;
        if(p == pmax)
            pstep = 1;
        std::cout<< "-------------- starting simulation for q = " + std::to_string(q) + "--------------" << std::endl;
    	//#pragma omp parallel while reduction(+:av_abs_m) private(p)
        while(p <= pmax){
            double av_abs_m = 0;
        	double av_m = 0;
        	double av_m2 = 0;
        	double av_m4 = 0;
        	int m =0;
            for(int av = 0; av < av_num;av++){
            	m = 0;
                std::vector<std::vector<short>> nums(Lx,std::vector<short>(Ly,0));
            	for (std::vector<short>& num : nums)
                {
	                num = (make_ran_lattice(Ly));
            	}
                for(int mc = 0; mc < mcs*N; mc++){
                    int spin_num_x = randInt_Uni(Lx-1); // choose the spin to propose a flip -x
                	int spin_num_y = randInt_Uni(Ly-1); // choose the spin to propose a flip -x
                    int spin = nums[spin_num_x][spin_num_y];
                    int nei_sum = 0;
                	int nei_num = lat->get_nn_number(spin_num_x,spin_num_y);
                    for(int j = 0; j < nei_num;j++){
                        auto nei_choice = lat->get_nn(spin_num_x,spin_num_y,j);
                    	int x = std::get<0>(nei_choice);
                    	int y = std::get<1>(nei_choice);
                        nei_sum += nums[x][y];
                    }
                    //std::cout<<(nei_sum);
                    double prob = proba(nei_sum,p,q,spin);
                    if(randZero_One_Uni() <= prob){
                        nums[spin_num_x][spin_num_y] = -spin;
                        //m -= 2*spin;
                    }
                }
            	for (auto& num : nums)
                {
            		for (short col : num)
                    {
            			m += col;
            		}

            	}	
                av_abs_m += abs(1.0*m/N);
            	av_m += 1.0*m/N;
            	av_m2 += ((1.0*m/N)*(1.0*m/N));
            	av_m4 += (1.0*m/N)*(1.0*m/N)*(1.0*m/N)*(1.0*m/N);
            }
        	av_abs_m /= av_num;
            av_m /= av_num;
            av_m2 /= av_num;
            av_m4 /= av_num;
            file << p << "\t" << (av_abs_m) << "\t" << N*(av_m2-av_abs_m*av_abs_m) << "\t" << 1-(av_m4/(3*av_m2*av_m2)) <<  endl;
        	//std:: cout << "for q = " + std::to_string(q) + " p = " + std::to_string(p) << " ->  ->  -> m = " +to_string(av_abs_m/N/av_num)  << std::endl;
            p +=pstep;
        }
        //------------------------------------------------------------------------------------            
        std::cout<< "-------------- finished simulation for q = " + std::to_string(q) + "--------------" << std::endl;
        //print(ps)
        //print(av_abs_ms)
    	file.close();
    });
	//delete[] lat;
}


//simulationAverages(50, 4, 500, 50, 0,1, 100, qs,savename = "simulation.png")
            
    



double randZero_One_Uni()
{
    return std::uniform_real_distribution<double>(0, 1)(eng);
}
int randInt_Uni(int L)
{
    return std::uniform_int_distribution<int>(0, L)(eng);
}


/// <summary>
/// Makes Euclidean modulo
/// </summary>
/// <param name="a">first value</param>
/// <param name="b">second value</param>
/// <returns>a%b</returns>
int myModuloEuclidean(int a, int b)
{
	int m = a % b;
	if (m < 0) {
		// m += (b < 0) ? -b : b; // avoid this form: it is UB when b == INT_MIN
		m = (b < 0) ? m - b : m + b;
	}
	return m;
}
