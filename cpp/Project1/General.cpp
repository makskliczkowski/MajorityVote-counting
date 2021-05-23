#include "General.h"
using namespace std;


/// <summary>
/// Function to split string into substrings by delimiter
/// </summary>
/// <param name="s">String to be split</param>
/// <param name="delimiter">Delimiter to be used</param>
/// <returns>Vector of substrings</returns>
vector<string> general::split_str(std::string s, std::string delimiter) {
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	string token;
	vector<string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
		token = s.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	}

	res.push_back(s.substr(pos_start));
	return res;
}



/* SQUARE LATTICE */

general::square_lattice::square_lattice()
{
	this->Lx = 0;
	this->Ly = 0;
}

general::square_lattice::square_lattice(int Lx, int Ly)
{
	this->Lx = Lx;
	this->Ly = Ly;
	this->Ns = Lx * Ly;
	this->calculate_nn_pbc();
	//this->calculate_nnn();
	this->type = general::lattice2D::lattice_types::square;

}

general::square_lattice::~square_lattice()
{
}

general::square_lattice::square_lattice(const square_lattice& A)
{
	this->Ns = A.Ns;
	this->type = A.type;
	this->nearest_neighbours = A.nearest_neighbours;
	this->next_nearest_neighbours = A.next_nearest_neighbours;
	this->type = general::lattice2D::lattice_types::square;
	/* Square lattice params */
	this->Lx = A.Lx;
	this->Ly = A.Ly;
}

general::square_lattice::square_lattice(square_lattice&& A) noexcept
{
	this->Ns = A.Ns;
	this->type = A.type;
	this->nearest_neighbours = std::move(A.nearest_neighbours);
	this->next_nearest_neighbours = std::move(A.next_nearest_neighbours);
	this->type = general::lattice2D::lattice_types::square;
	/* Square lattice params */
	this->Lx = A.Lx;
	this->Ly = A.Ly;
}

int general::square_lattice::get_Lx()
{
	return this->Lx;
}

int general::square_lattice::get_Ly()
{
	return this->Ly;
}

/// <summary>
/// Nearest with periodic boundary conditions in 2D
/// </summary>
void general::square_lattice::calculate_nn_pbc()
{
		/* Two dimensions */
	this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int,int>>> (Ly,std::vector<std::tuple<int, int>>(4)));
	for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			/* LEFT AND TOP ARE SET TWO FIRST TWO */
			this->nearest_neighbours[x][y][3] = std::make_tuple(myModuloEuclidean(x+1,Lx), y); // right
			this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x-1,Lx), y); //left
			this->nearest_neighbours[x][y][2] = std::make_tuple(x, myModuloEuclidean(y+1,Ly)); // bottom
			this->nearest_neighbours[x][y][1] = std::make_tuple(x, myModuloEuclidean(y-1,Ly)); // top
		}
	}
}
/// <summary>
/// Next nearest with periodic boundary conditions in 2D
/// </summary>
void general::square_lattice::calculate_nnn_pbc()
{
	/* Two dimensions */
	this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int,int>>> (Ly,std::vector<std::tuple<int, int>>(4)));
	for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			/* LEFT AND TOP ARE SET TWO FIRST TWO */
			this->nearest_neighbours[x][y][3] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
			this->nearest_neighbours[x][y][0] = std::make_tuple(x - 1 == -1 ? Lx-1 : x - 1, y); //left
			this->nearest_neighbours[x][y][2] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom
			this->nearest_neighbours[x][y][1] = std::make_tuple(x, y - 1 == -1 ? Ly-1 : y - 1); // top
		}
	}
}
/// <summary>
/// nearest with free boundary conditions in 2D
/// </summary>
void general::square_lattice::calculate_nn()
{
	/* Two dimensions */
	this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int,int>>> (Ly,std::vector<std::tuple<int, int>>(4)));
	/*for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			this->nearest_neighbours[x][y][3] = std::make_tuple(x + 1 == Lx ? -1 : x + 1, y); // right
			this->nearest_neighbours[x][y][0] = std::make_tuple(x - 1, y); //left
			this->nearest_neighbours[x][y][2] = std::make_tuple(x, y + 1 == Ly ? -1 : y + 1); // bottom
			this->nearest_neighbours[x][y][1] = std::make_tuple(x, y - 1); // top
		}
	}*/
	for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			/* LEFT AND TOP ARE SET TWO FIRST TWO */
			this->nearest_neighbours[x][y][3] = std::make_tuple(myModuloEuclidean(x+1,Lx), y); // right
			this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x-1,Lx), y); //left
			this->nearest_neighbours[x][y][2] = std::make_tuple(x, myModuloEuclidean(y+1,Ly)); // bottom
			this->nearest_neighbours[x][y][1] = std::make_tuple(x, myModuloEuclidean(y-1,Ly)); // top
		}
	}
}
/// <summary>
/// Next neares with free
/// not used 
/// </summary>
void general::square_lattice::calculate_nnn()
{

}



/* GENERAL LATTICE */

int general::lattice2D::get_Ns()
{
	return this->Ns;
}

general::lattice2D::lattice_types general::lattice2D::get_type()
{
	return this->type;
}

std::string general::lattice2D::getString_type()
{
	std::string lat;
	switch (this->get_type()) {
	case general::lattice2D::lattice_types::square:
		lat = "square";
		break;
	case general::lattice2D::lattice_types::triangle:
		lat = "triangle";
		break;
	case general::lattice2D::lattice_types::hexagonal:
		lat = "hexagonal";
	default:
		std::cout << "Don't know that option for lattice, excuse me\n";
		std::cout << "Try square,triangle,hexagonal -> got'em\n";
		exit(-1);
	}
	return lat;
}

std::tuple<int, int> general::lattice2D::get_nn(int x, int y, int nei_num)
{
	return this->nearest_neighbours[x][y][nei_num];
}

std::tuple<int, int> general::lattice2D::get_nnn(int x, int y, int nei_num)
{
	return std::tuple<int, int>();
}

int general::lattice2D::get_nn_number(int x, int y)
{
	return this->nearest_neighbours[x][y].size();
}

int general::lattice2D::get_nnn_number(int x, int y)
{
	return this->next_nearest_neighbours[x][y].size();
}

/* TRAINGULAR LATTICE */

general::triangle_lattice::triangle_lattice()
{
	this->Ly = 0;
	this->Lx = 0;
	this->type = general::lattice2D::lattice_types::triangle;
}

general::triangle_lattice::triangle_lattice(int Lx, int Ly)
{
	this->Lx = Lx;
	this->Ly = Ly;
	this->Ns = Lx * Ly;
	this->calculate_nn_pbc();
	this->calculate_nnn();
	this->type = general::lattice2D::lattice_types::triangle;
}

general::triangle_lattice::~triangle_lattice()
{

}

general::triangle_lattice::triangle_lattice(const triangle_lattice& A)
{
	this->Ns = A.Ns;
	this->type = A.type;
	this->nearest_neighbours = A.nearest_neighbours;
	this->next_nearest_neighbours = A.next_nearest_neighbours;
	this->Lx = A.Lx;
	this->Ly = A.Ly;
	this->type = general::lattice2D::lattice_types::triangle;
}

general::triangle_lattice::triangle_lattice(triangle_lattice&& A) noexcept
{
	this->Ns = A.Ns;
	this->type = A.type;
	this->nearest_neighbours = std::move(A.nearest_neighbours);
	this->next_nearest_neighbours = std::move(A.next_nearest_neighbours);
	this->Lx = A.Lx;
	this->Ly = A.Ly;
	this->type = general::lattice2D::lattice_types::triangle;
}

int general::triangle_lattice::get_Lx()
{
	return this->Lx;
}

int general::triangle_lattice::get_Ly()
{
	return this->Ly;
}

void general::triangle_lattice::calculate_nn_pbc()
{
	this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int, int>>>(Ly, std::vector<std::tuple<int, int>>(6)));
	for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			/* LEFT AND TOP ARE SET TWO FIRST TWO */
			/* ZIGZAG */
			if (y % 2 == 0) {
				this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(x, myModuloEuclidean(y - 1,Ly)); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y + 1 == Ly ? 0 : y + 1); // bottom_right
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, myModuloEuclidean(y - 1,Ly)); // top_right 
			}
			else {
				this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(myModuloEuclidean(x - 1,Lx), myModuloEuclidean(y - 1,Ly)); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y + 1 == Ly ? 0 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom_rigt
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // top_right
			}

		}
	}
	
}

void general::triangle_lattice::calculate_nnn_pbc()
{
	this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int, int>>>(Ly, std::vector<std::tuple<int, int>>(6)));
	for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			/* LEFT AND TOP ARE SET TWO FIRST TWO */
			/* ZIGZAG */
			if (y % 2 == 0) {
				this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(x, myModuloEuclidean(y - 1,Ly)); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y + 1 == Ly ? 0 : y + 1); // bottom_right
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, myModuloEuclidean(y - 1,Ly)); // top_right 
			}
			else {
				this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(myModuloEuclidean(x - 1,Lx), myModuloEuclidean(y - 1,Ly)); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y + 1 == Ly ? 0 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom_rigt
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // top_right
			}

		}
	}
}

void general::triangle_lattice::calculate_nn()
{
	/* Two dimensions */
	this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int, int>>>(Ly, std::vector<std::tuple<int, int>>(6)));
	/*for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {

			if (y % 2 == 0) {
				this->nearest_neighbours[x][y][0] = std::make_tuple(x - 1, y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(x, y - 1); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(x, y + 1 == Ly ? -1 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x + 1 == Lx ? -1 : x + 1, y + 1 == Ly ? -1 : y + 1); // bottom_right
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? -1 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x + 1 == Lx ? -1 : x + 1, y - 1); // top_right 
			}
			else {
				this->nearest_neighbours[x][y][0] = std::make_tuple(x - 1, y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(x - 1, y - 1); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(x-1, y + 1 == Ly ? -1 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x, y + 1 == Ly ? -1 : y + 1); // bottom_rigt
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? -1 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x, y + 1 == Ly ? -1 : y + 1); // top_right
			}

		}
	}*/
	//this->nearest_neighbours = neighborsStructure(Lx, std::vector<std::vector<std::tuple<int, int>>>(Ly, std::vector<std::tuple<int, int>>(6)));
	for (int x = 0; x < this->Lx; x++) {
		for (int y = 0; y < this->Ly; y++) {
			/* LEFT AND TOP ARE SET TWO FIRST TWO */
			/* ZIGZAG */
			if (y % 2 == 0) {
				this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(x, myModuloEuclidean(y - 1,Ly)); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y + 1 == Ly ? 0 : y + 1); // bottom_right
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, myModuloEuclidean(y - 1,Ly)); // top_right 
			}
			else {
				this->nearest_neighbours[x][y][0] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y); // left
				this->nearest_neighbours[x][y][1] = std::make_tuple(myModuloEuclidean(x - 1,Lx), myModuloEuclidean(y - 1,Ly)); // top_left
				this->nearest_neighbours[x][y][2] = std::make_tuple(myModuloEuclidean(x - 1,Lx), y + 1 == Ly ? 0 : y + 1); // bottom_left
				this->nearest_neighbours[x][y][3] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // bottom_rigt
				this->nearest_neighbours[x][y][4] = std::make_tuple(x + 1 == Lx ? 0 : x + 1, y); // right
				this->nearest_neighbours[x][y][5] = std::make_tuple(x, y + 1 == Ly ? 0 : y + 1); // top_right
			}

		}
	}

}

void general::triangle_lattice::calculate_nnn()
{
}
