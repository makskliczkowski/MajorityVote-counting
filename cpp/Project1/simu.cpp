#include "simu.h"
#include "xoshiro_pp.h"

int sign(int x) {
    return (x > 0) - (x < 0);
}
inline double proba(int nei_sum, double p, double q, int spin){
    return 0.5*((1-p)*(1-(1-2*q)*spin*sign(nei_sum))+p);
}
std::vector<short> make_ran_lattice(int N,XoshiroCpp::Xoshiro256PlusPlus& eng)
{
	std::vector<short> lat(N);
	for_each(lat.begin(), lat.end(), [lat,eng](short& b) mutable {b = (randZero_One_Uni(eng)) < 0.5 ? -1 : 1; });
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
    //std::random_device rd{};
	//std::mt19937 eng{rd()};
	const std::uint64_t seed = static_cast<uint64_t>(std::time(nullptr));
	XoshiroCpp::Xoshiro256PlusPlus eng(Random_SeedInit(seed));
//#pragma omp parallel for
	const int qs_num = qs.size();
    for(int qi=0; qi < qs_num; qi++)
	{
    	const double q = qs[qi];
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
                std::vector<short> nums = (make_ran_lattice(N,eng));
                //m = std::accumulate(nums.begin(), nums.end(), 0);	
                for(int mc = 0; mc < mcs*N; mc++){
                    const int spin_num = randInt_Uni(N-1,eng); // choose the spin to propose a flip
                    int spin = nums[spin_num];
                    int nei_sum = 0;
                    for(int j = 0; j<num_nei;j++){
                        int nei_choice = spin_num;
                        //print("spin_num = " + str(spin_num))
                        while(nei_choice == spin_num){
                        nei_choice = randInt_Uni(N-1,eng); // choose the spin to propose a neighbor
                            //print(nei_choice)
                        }
                        nei_sum += nums[nei_choice];
                    }
                    //print(nei_sum)
                    double prob = proba(nei_sum,p,q,spin);
                    if(randZero_One_Uni(eng) <= prob){
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
    }
	
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
    const int N = Lx*Ly;
    general::lattice2D* lat;
	if(type == "square") lat = new general::square_lattice(Lx,Ly);
    else lat = new general::triangle_lattice(Lx,Ly);

	
	const std::uint64_t seed = static_cast<uint64_t>(std::time(nullptr));
	XoshiroCpp::Xoshiro256PlusPlus eng(Random_SeedInit(seed));

	const int qs_num = qs.size();
    for(int qi=0; qi < qs_num; qi++)
	{
    	const double q = qs[qi];
    	int done = 0;
        ofstream file;
    	ofstream fileTime;
    	string fname = savename+",q="+ to_string(q) + ",N="+to_string(N)+ ",mcs" + to_string(mcs) + ",avNum,NONINFTY"+to_string(av_num) + ".dat";
        file.open(fname);
    	if(!file.is_open())
    	{
    		cout << "couldn't open a file" << endl;
    	}
    	file << "p" << "\t" << "(av_abs_m)" << "\t"<<"susceptibility"<<"\t" << "binder cumulant" << endl;
        //double p = pstart;
        double pstep = 1.0*(pmax - pstart)/pnum;
        if(pstart == pmax)
            pstep = 1;
        std::cout<< "-------------- starting simulation for q = " + std::to_string(q) + " and N = " + std::to_string(N) + "--------------" << std::endl;
    	std::vector<std::vector<std::tuple<double,double>>> averages(3);
    	#pragma omp parallel for num_threads(5) shared(done)
        for(int pi = 0;pi< pnum; pi++){
        	double p = pstart + pi*pstep;
            double av_abs_m = 0;// 0
        	//double av_m = 0; 
        	double av_m2 = 0; // 1
        	double av_m4 = 0; // 2
        	int m = 0;
        	std::vector<std::vector<short>> nums(Lx,std::vector<short>(Ly));
        	for (std::vector<short>& num : nums)
            {
            	for(auto & elem : num)
            	{
            		elem = randZero_One_Uni( eng) >= 0.5? 1:-1;
            	}
                //num = (make_ran_lattice(Ly));
            }
        	const int nei_num = lat->get_nn_number(0,0);
        	const int mcsteps = mcs*N;
        	for(int mc = 0; mc < mcsteps; mc++){
                const int spin_num = randInt_Uni(N-1,eng); // choose the spin to propose a flip -x
                //const int spin_num_y = randInt_Uni(Ly-1,eng); // choose the spin to propose a flip -x
        		const int spin_num_x = lat->get_coordinates(spin_num,0);
        		const int spin_num_y = lat->get_coordinates(spin_num,1);
                const short spin = nums[spin_num_x][spin_num_y];
                int nei_sum = 0;
                for(int j = 0; j < nei_num;j++){
                    auto nei_choice = lat->get_nn(spin_num_x,spin_num_y,j);
                    nei_sum += nums[std::get<0>(nei_choice)][std::get<1>(nei_choice)];
                }
                //std::cout<<(nei_sum);
                if(XoshiroCpp::DoubleFromBits(eng()) <= proba(nei_sum,p,q,spin))
                    nums[spin_num_x][spin_num_y] = - spin;
                    //m -= spin * 2;
            }
        	const int corr_time = 100*N;
			for(int av = 0; av < av_num;av++){
				for(int mc = 0; mc < corr_time; mc++){
					const int spin_num = randInt_Uni(N-1,eng); // choose the spin to propose a flip -x
        			const int spin_num_x = lat->get_coordinates(spin_num,0);
        			const int spin_num_y = lat->get_coordinates(spin_num,1);
	                const short spin = nums[spin_num_x][spin_num_y];
	                int nei_sum = 0;
	                for(int j = 0; j < nei_num;j++){
	                    auto nei_choice = lat->get_nn(spin_num_x,spin_num_y,j);
	                    nei_sum += nums[std::get<0>(nei_choice)][std::get<1>(nei_choice)];
	                }
	                //std::cout<<(nei_sum);
	                if(XoshiroCpp::DoubleFromBits(eng()) <= proba(nei_sum,p,q,spin))
	                    nums[spin_num_x][spin_num_y] = -spin;
	                    //m -= spin * 2;
                }
				m=0;
				for (std::vector<short>& num : nums)
	            {
            		for(auto & elem : num)
            		{
            			m+=static_cast<int>(elem);
            		}
	                //num = (make_ran_lattice(Ly));
				}
            	//double mx=1.0*m;
                av_abs_m += 1.0*abs(m)/N;
            	//averages[0][pi] +=abs(m);
            	//averages[1][pi] +=m*m;
            	//averages[2][pi] +=m*m*m*m;
            	//av_m += m/xN;
            	av_m2 += 1.0*m*m/N;
            	av_m4 += 1.0*m*m*m*m/(N*N);	
            }
        	//averages[0][pi] /= av_num;
            //averages[1][pi] /=av_num;
            //averages[2][pi] /=av_num;
        	av_abs_m = av_abs_m/(1.0*av_num);
            //av_m /= av_num;
            av_m2 = av_m2/(1.0*N * av_num);
            av_m4 = av_m4/(1.0*N *N * av_num);
            averages[0].emplace_back(p,av_abs_m);
        	averages[1].emplace_back(p,av_m2);
        	averages[2].emplace_back(p,av_m4);
        	
        	//std:: cout << "for q = " + std::to_string(q) + " p = " + std::to_string(p) << " ->  ->  -> m = " +to_string(av_abs_m/N/av_num)  << std::endl;
            p +=pstep;
        	done++;
        	if(done%10 == 0)
        	{
        		cout << "Finished " << (1.0 * done / pnum) *100 << "%\n";
        	}
        }

        for(int pi = 0; pi<pnum; pi++)
        {
        	double p = std::get<0>(averages[0][pi]);//pstart + pi*pstep;
        	double abs_m = std::get<1>(averages[0][pi]);
        	double m2 = std::get<1>(averages[1][pi]);
        	double m4 = std::get<1>(averages[2][pi]);
	        file << p << "\t" << abs_m << "\t" << N*(m2-abs_m*abs_m) << "\t" << 1.0-(m4/(3.0*m2*m2)) <<  endl;
        }
    	

        //------------------------------------------------------------------------------------            
        std::cout<< "-------------- finished simulation for q = " + std::to_string(q) + " and N = " + std::to_string(N) + "--------------" << std::endl;
        //print(ps)
        //print(av_abs_ms)
    	file.close();
    }
	//delete[] lat;
}


//simulationAverages(50, 4, 500, 50, 0,1, 100, qs,savename = "simulation.png")
            
    



inline double randZero_One_Uni(XoshiroCpp::Xoshiro256PlusPlus& eng)
{
    return std::uniform_real_distribution<double>(0, 1)( eng);
}
inline int randInt_Uni(int L,XoshiroCpp::Xoshiro256PlusPlus& eng)
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
uint64_t Random_SeedInit(uint64_t n)
{
    std::vector<uint64_t> s(16,0);
	for (int i = 0; i < 16; i++)
	{
		n ^= n >> 12;   // a
		n ^= n << 25;   // b
		n ^= n >> 27;   // c

		// 2685821657736338717 = 72821711 * 36882155347, from Pierre L'Ecuyer's paper
		s[i] = n * 2685821657736338717LL;
	}
    return std::accumulate(s.begin(),s.end(), 0.0);
}