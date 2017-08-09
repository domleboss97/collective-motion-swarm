/*
This program takes console input of properties of the 2d Vicsek model and generates
a CSV file of particle positions, direction vectors, and the global order parameter.
The first three lines of the output file are the population, the size of the environment
(in terms of the particle's radius of alignment) and eta, the noise parameter.
The csv file prints a line for each particle, and inserts a blank line between iterations.
*/

#include <iostream>
#include <iomanip>
#include <cmath>    //cos and sin
#include <vector>   //used to track particles within a cell
#include <algorithm>//min and max
#include <iterator> //vector pointer operations
#include <ctime>    //time-based seed
#include <random>   //uniform random distribution
#include <fstream>

#define PI 3.14159265

using namespace std;

//Contains data about each particle in simulation
class Swimmer{
    public:
        float angle;
        double x_pos, y_pos, x_sum, y_sum;
        int angle_inc;

        Swimmer()
        {
            angle = 0;
            x_pos = 0;
            y_pos = 0;
            x_sum = 0;
            y_sum = 0;
        }

        void set_angle(float a)
        {
            angle = a;
        }

        void set_x_position(double x)
        {
            x_pos = x;
        }

        void set_y_position(double y)
        {
            y_pos = y;
        }
};

//Structure that counts how many particles are in a 1x1 grid
struct Square{
    vector<int> neighbors;

    Square()
    {
        vector<int> neighbors(0);
    }
};

//Function which determines whether or not a particle is close enough to affect distance in vicsek model
double coord_distance(double x_1, double y_1, double x_2, double y_2)
{
    return sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2));
}

int t, b, j, k, l, m, n, p; //For-loop counters

double old_x_pos, old_y_pos, x_movement, y_movement, order_check_new, order_check_old = 0;
bool stop = false;

default_random_engine generator1(time(0)%91);
default_random_engine generator2(time(0)%83);


int main()
{
    int x_cell, y_cell, new_x_cell, new_y_cell; //For tracking cell
    double noise;
    double order_coeff, order_sum_x = 0, order_sum_y = 0, time_track = 0;
    vector<int>::iterator it;
    vector<int> neighbor_list; //Used to keep track of the particles within a close enough range to potentially affect direction

    uniform_real_distribution<double> angles (0,360.0);//Uniform real distribution for determining random velocity angle

    ofstream fout("output.csv", ios::out);
    ofstream fout3("order.csv", ios::out);

    int popsize, fieldsize, timesteps;
    double eta;

    cout << "Population size?" << endl;
    cin >> popsize;

    cout << "Field size?" << endl;
    cin >> fieldsize;

    cout << "Maximum iterations?" << endl;
    cin >> timesteps;

    cout << "Noise coefficient?" << endl;
    cin >> eta;

    fout << popsize << endl << fieldsize << endl << eta << endl;

    fout << fixed;
    fout << setprecision(2);

    uniform_real_distribution<double> positions (0.0, fieldsize); //Uniform real distribution for determining random particle position

    vector<Swimmer> parts(popsize);
    vector< vector<Square> > cells(fieldsize, vector<Square>(fieldsize));

    time_track = clock();

    //Initialize properties of the swimmer objects
    for (b = 0; b < popsize; b++)
    {
        float a = angles(generator1);
        double x = positions(generator1);
        double y = positions(generator2);
        parts[b].set_angle(a);
        parts[b].set_x_position(x);
        parts[b].set_y_position(y);

        x_cell = x;
        y_cell = y;

        cells[x_cell][y_cell].neighbors.push_back(b);
    }

    //Greater outer loop with as many iterations as time steps.
    for (t = 1; stop == false; t++)
    {
        order_sum_x = 0;
        order_sum_y = 0;

        //Add angle to order calculation, then reset angle calculations at each new timestep
        for (b = 0; b < popsize; b++)
        {
            //cout << "angle " << i << " =: " << parts[i].angle << endl;
            order_sum_x += cos(parts[b].angle * PI / 180.0);
            order_sum_y += sin(parts[b].angle * PI / 180.0);

            parts[b].x_sum = 0;
            parts[b].y_sum = 0;
            parts[b].angle_inc = 0;
            parts[b].angle = fmod(parts[b].angle, 360.0);
        }

        order_coeff = sqrt(pow(order_sum_x/popsize, 2) + pow(order_sum_y/popsize, 2));

        fout3 << order_coeff << endl;

        //Loop repeated for each particle in which neighboring-cell particles are located
        for (b = 0; b < popsize; b++)
        {
            m = parts[b].x_pos - 1; //floored result from x position
            n = parts[b].y_pos - 1; //floored result from y position

            if (m < 0)              //loop back m/n if too large/small
                m += fieldsize;
            else if (m > fieldsize - 1)
                m -= fieldsize;

            if (n < 0)
                n += fieldsize;
            else if (n > fieldsize - 1)
                n -= fieldsize;

            for (j = 0; j <= 2; j++)//increment X dimension
            {
                for (k = 0; k <= 2; k++)//increment Y dimension
                {
                    x_cell = (m + j) % fieldsize;
                    y_cell = (n + k) % fieldsize;

                    for (p = 0; p < cells[x_cell][y_cell].neighbors.size(); p++)//Checks size of array of particles within neighboring cell
                    {
                        neighbor_list.push_back(cells[x_cell][y_cell].neighbors[p]); //Adds indeces of particles from neighboring cells to neighbor_list, which tracks
                    }
                }
            }

            for (p = 0; p < neighbor_list.size(); p++)//Cycle through each element of neighbor_list
            {
                l = neighbor_list[p];//set l as the index of the particle held in neighbor_list
                //If statement checks if particle is within 1 unit of another particle, and accounts for wraparound effect of board
                if (l != b  &&((coord_distance(    parts[l].x_pos,     parts[l].y_pos,         parts[b].x_pos,     parts[b].y_pos) < 1)
                            || (coord_distance(min(parts[l].x_pos,     parts[b].x_pos)+10,     parts[l].y_pos, max(parts[l].x_pos,     parts[b].x_pos), parts[b].y_pos) < 1)
                            || (coord_distance(    parts[l].x_pos, min(parts[l].y_pos,         parts[b].y_pos)+10, parts[b].x_pos, max(parts[l].y_pos, parts[b].y_pos)) < 1)
                            || (coord_distance(min(parts[l].x_pos,     parts[b].x_pos)+10, min(parts[l].y_pos,     parts[b].y_pos)+10, max(parts[l].x_pos, parts[b].x_pos), max(parts[l].y_pos, parts[b].y_pos)) < 1)))

                {
                    parts[b].x_sum += cos(parts[l].angle * PI / 180.0);
                    parts[b].y_sum += sin(parts[l].angle * PI / 180.0);
                    parts[b].angle_inc++;
                }
            }

            neighbor_list.clear();//Clears neighbor_list and resets capacity to 0
        }

        for (b = 0; b < popsize; b++)
        {
            old_x_pos = parts[b].x_pos; //Stores previous position in placeholder variable
            old_y_pos = parts[b].y_pos;

            noise = (angles(generator1) - 180.0) * eta; //noise term is a function of present velocity direction and uniform real distr. of angles

            if (parts[b].angle_inc > 0)
            {
                parts[b].angle = (atan2(parts[b].y_sum, parts[b].x_sum) * 180.0 / PI) + noise; //When all contributing angles are summed and divided by the # of contributing angles, average angle achieved
            }
            else
            {
                parts[b].angle += noise;
            }

            x_movement = cos((parts[b].angle) * PI / 180.0);
            y_movement = sin((parts[b].angle) * PI / 180.0);


            parts[b].x_pos += x_movement; //move in x and y direction based on all factors contributing to velocity angle
            parts[b].y_pos += y_movement;

            if (parts[b].x_pos > fieldsize)
            {
                parts[b].x_pos -= fieldsize; //account for particle landing outside of field
            }
            else if (parts[b].x_pos < 0)
            {
                parts[b].x_pos += fieldsize;
            }

            if (parts[b].y_pos > fieldsize)
            {
                parts[b].y_pos -= fieldsize;
            }
            else if (parts[b].y_pos < 0)
            {
                parts[b].y_pos += fieldsize;
            }

            if      ((static_cast<int>(old_x_pos) != static_cast<int>(parts[b].x_pos)) || (static_cast<int>(old_y_pos) != static_cast<int>(parts[b].y_pos))) //check to see if particle has moved to new cell
            {

                x_cell = static_cast<int>(old_x_pos);
                y_cell = static_cast<int>(old_y_pos);

                new_x_cell = static_cast<int>(parts[b].x_pos);
                new_y_cell = static_cast<int>(parts[b].y_pos);

                it = find(cells[x_cell][y_cell].neighbors.begin(), cells[x_cell][y_cell].neighbors.end(), b); //Search for particle index in previous cell and clears it
                cells[x_cell][y_cell].neighbors.erase(it);

                cells[new_x_cell][new_y_cell].neighbors.push_back(b); //Places particle index in new cell storage vector
            }

                fout << parts[b].x_pos << "," << parts[b].y_pos << "," << x_movement << "," << y_movement << "," << order_coeff << endl;

        }
        fout  << endl;
        if (t % 500 == 0)
        {

            order_check_new = order_coeff;


            cout << "New: " << setw(10) << order_check_new << ", Old: " << setw(10) << order_check_old << endl;

            if (abs(order_check_new - order_check_old) < 0.02)
            {
                stop = true;
            }
            else if (t >= timesteps)
            {
                stop = true;
            }

            order_check_old = order_coeff;
        }

    }

    time_track = clock() - time_track;
    cout << "time = " << time_track/CLOCKS_PER_SEC << endl;

    return 0;
}
