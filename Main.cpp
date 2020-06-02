/*
	ABOUT PROGRAM
The our model consist water that forced into a porous medium containig oil. The idea is to use the water 
to recover as much least resistance. We consider a lattice of size 2L*L, with the water(the invader) initialy
occuping left edge. The resistance to the invader is given by uniformly distributed random nember between 0 and 1, 
which assign to each site in the lattice, and held fixed throughout the invasion. Sites that are nearest neighbors
of the invader site are perimeter sites. At each step time the perimeter with lowest random number is occupied by 
the invader and the oil(defender) is diplaced.
The invading cluster growth until a path forms wich connects the left and right edges of lattice.

** We use periodic condition in TOP and Bottom of edges for minimze boundary effects.
** All quantities are measured only over central L*L region of lattice.
** site(i,j), store random number. If the site at (i,j) is occuiped, then site(i,j) increased by 1, if is perimeter,
increased by 2. Then it inserted into proper orederd position in the perimeter list perX, perY.
The perimeter list are ordered with the site with the largest random number at the begining of the list.


*/

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>

using namespace std;



// define structures
struct Location
{
    int x, y;
    Location() : x(0), y(0)
    {}
    Location(int x, int y) :x(x), y(y)
    {}

    Location& operator = (const Location& loc)
    {
        x = loc.x;
        y = loc.y;
        return *this;
    }
};
// End define Structers


// Fileds
const int maxSizeY = 1000, maxSizeX = 2*maxSizeY;
int Ly, Lx;
const int maxPerimterSize = maxSizeX * maxSizeY;
Location perimeter[maxPerimterSize];
double site[maxSizeX][maxSizeY];
int n_per;		// the number of perimtere sites in list.
int P[maxSizeX/2+1];  //  P(r)dr that represent the probabilty of site that ha random number between r and r+dr.
const double dr = 0.05;
int L_r = (int)(1/dr)+1; // the lenght of P array.
int nr[maxSizeX/2+1]; // the array that  how much sites that  their random value is r .
int n_occuiped ; // total number of occuiped site.
int n;	// number sites in middle  half of lattice.
double P_Final[maxSizeX/2+1];

const int nx[4] = {1, 0, -1, 0}, ny[4] = {0, 1, 0, -1};
		
//


// Prototype Methods
inline int Rand( int nMin, int nMax );
inline double Rand();

void Assign();
void Insert(int x, int y);
int IndexBinarySearch(int x, int y);
int IndexLinearSearch(int x, int y);
void Invade();
void Average();
void ResetAll();
void Work1();
void Work2();
//



// calculate P(r)dr that represent the probabilty of site that ha random number between r and r+dr.
void Work1()
{
	clock_t t1 = clock();
	
	const int iteration = 2000;
	Ly = 200;
	Lx = 2 * Ly;

	cout << "0%" ;
	for(int it = 0 ; it < iteration; it++)
	{
		ResetAll();
		Assign();
		Invade();
		Average();

		for(int i_r = 0; i_r < L_r; i_r++)
			if(nr[i_r] >0)
				P_Final[i_r] += (double)P[i_r] / nr[i_r] ;

		if((it+1) % 20 == 0)
			cout << "\b\b\b" <<100*(it+1)/iteration << "%" ;
	}
	for(int i_r = 0; i_r < L_r; i_r++)
		P_Final[i_r] /= iteration;

	t1 = clock() - t1;
	cout << "\n\n\a ** Elapsed time = " << t1/CLOCKS_PER_SEC  << endl;


	ofstream outP("P(r).txt");
	outP << "r" << "\t" << "P(r)" << endl;
	for(int i_r = 0; i_r < L_r; i_r++)
		outP << dr * i_r << "\t" << P_Final[i_r] << endl;
	outP.close();
};


// calculate M(L) : the numbr of occuiped site in middle half of latice.   M(L) = L^d
void Work2()
{
	clock_t t1 = clock();

	const int iterations = 2000;
	const int lMin = 20, lMax = 200, dl = 20;
	int i_min, i_max, x, y;

	int M;

	ofstream outM("M(L).txt");
	outM << "L" << "\t" << "M(L)" << endl;

	cout <<"0%";
	for(int L = lMin; L <= lMax; L += dl)
	{
		M = 0;
		Ly = L;
		Lx = 2 * L;
		for(int i = 0; i < iterations; i++)
		{
			ResetAll();
			Assign();
			Invade();
			i_min = L/2;
			i_max = Lx - L/2 - 1;
			
			for ( x = i_min; x <= i_max; x++)
				for(y = 0; y < Ly; y++)
					if(site[x][y] >= 1 && site[x][y] < 2)
						n_occuiped++;

			M += n_occuiped;
		}
		
		outM << L << "\t" << (double)M / iterations << endl;

		cout <<"\b\b\b" << 100*(L-lMin + dl)/(lMax-lMin+dl) << "%";
	}
	outM.close();
	t1 = clock() - t1;
	cout << "\n\n\aElapsed time = " << t1/CLOCKS_PER_SEC  << endl;
};

void ResetAll()
{
	n_occuiped = 0;
	n_per = 0;
	n = 0;

	int i;
	for( i =0 ; i < L_r; i++)
	{
		P[i] = 0;
		nr[i] = 0;
	}
};


int main()
{
	srand(unsigned(time(NULL)));

	
	Work2();

	

	//cout << "** Total number of occuiped sites in middle half of latice(" <<  Lx << "*" << Ly << ") = " << n_occuiped << endl;


	return 0;
}


void Assign()
{
	int i, j;
	for( j = 0; j < Ly; j++ )
		site[0][j] = 1;			// occuping firs columns

	for( j = 0; j < Ly; j++ )
		for( i = 1; i < Lx; i++)		// assign random number
			site[i][j] = Rand();

	// sites in second column are inital perimters sites. 
	i = 1;
	n_per = 0;
	for( j = 0; j < Ly; j++)
	{
		site[i][j] += 2;
		n_per++;
		Insert(i, j); // sort perimeter sites. order perimiter list.
		
	}
}



void Insert(int x, int y)
{
	int i_insert = IndexBinarySearch(x, y);  // index of insert 
	
	// move site with smaler value to next high array index
	for(int i = n_per - 1; i > i_insert; i--)
		perimeter[i] = perimeter[i-1];

	perimeter[i_insert].x = x;
	perimeter[i_insert].y = y;

}


int IndexBinarySearch(int x, int y)
{
	int iFirst = 0, iEnd = n_per - 2;
	if(iEnd  < 0  || site[x][y] >= site[perimeter[iFirst].x][perimeter[iFirst].y])
		return 0;
	else if(site[x][y] <= site[perimeter[iEnd].x][perimeter[iEnd].y])
		return iEnd + 1;


	int iMid = (iFirst + iEnd) / 2;
	int xMid, yMid;

	do
	{
		if( iEnd - iFirst <= 1 )
			return iFirst + 1;
		
		xMid = perimeter[iMid].x;
		yMid = perimeter[iMid].y;
		if(site[x][y] >= site[xMid][yMid])
			iEnd = iMid;
		else
			iFirst = iMid;
		iMid = (iFirst + iEnd) / 2;
	}while(true);

	return 0;

}

int IndexLinearSearch(int x, int y)
{
	if(n_per == 1)
		return 0;

	int tempX, tempY, i;
	for( i = 0; i < n_per - 2; i++ )
	{
		tempX = perimeter[i].x;
		tempY = perimeter[i].y;
		if(site[x][y] >= site[tempX][tempY])
			return i;
	}
	return i;
}

void Invade()
{
	int x, y, xNew, yNew, i;

	do
	{
		x = perimeter[n_per-1].x;
		y = perimeter[n_per-1].y;
		n_per--;

		site[x][y]--;

		for( i=0; i<4; i++)
		{
			xNew = x + nx[i];
			yNew = y + ny[i];

			if(yNew >= Ly)	   yNew = 0;
			else if(yNew < 0 ) yNew = Ly - 1;

			if(site[xNew][yNew] < 1)
			{
				site[xNew][yNew] += 2;
				n_per++;
				Insert(xNew, yNew);
			}
		}
	}while(x < Lx-1);			// until clutser reach the right edge
}


void Average()
{
	int i_min = Ly/2;
	int i_max = Lx - Ly/2 -1;
	n = (i_max - i_min + 1)* Ly;  
	int x, y, i_r;

	for ( x = i_min; x <= i_max; x++)
	{
		for(y = 0; y < Ly; y++)
		{
			i_r = (int)(L_r * fmod(site[x][y], 1.0));
			nr[i_r]++;
			if(site[x][y] >= 1 && site[x][y] < 2)		
			{
				n_occuiped++;
				P[i_r]++;
			}
			
		}
	}

	
}


inline int Rand( int nMin, int nMax )
{
    int nRange = nMax - nMin;
    int nNum = rand() % nRange;
    return( nNum + nMin );
}
inline double Rand()
{
	return (double)rand() / RAND_MAX;

}
